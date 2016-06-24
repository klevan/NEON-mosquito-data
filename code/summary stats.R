rm(list=ls())
library(dplyr)
library(graphics)
options(stringsAsFactors = F)

# NOTE: Run 'getting-the-data.R' first before running this script
# Change 'myPathToData' to your file location
myPathToData <- "~/GitHub/NEON-mosquito-data/data" 
samp <- read.csv(file = paste(myPathToData, "mos_samplingeffort_pub.csv", 
                              sep = "/"), 
                 header = TRUE, encoding = "UTF-8")

id <- read.csv(file = paste(myPathToData, "mos_identification_pub.csv", 
                            sep = "/"), 
               header = TRUE, encoding = "UTF-8")

arc <- read.csv(file = paste(myPathToData, "mos_archival_pub.csv", 
                             sep = "/"), 
                header = TRUE, encoding = "UTF-8")

path <- read.csv(file = paste(myPathToData, "mos_pathogenresults_pub.csv", 
                              sep = "/"), 
                 header = TRUE, encoding = "UTF-8")

# Sampling summary table
sampling <- summarize(group_by(samp, eventID), 
                      siteID = unique(siteID),
                      totalTrapHours = sum(trapHours), 
                      numTrapsDeployed = length(trapHours), 
                      startCollectDate = min(setDateTime), 
                      endCollectDate = max(collectDateTime))

sampling <- left_join(sampling,
                      summarize(group_by(id[grepl("species", id$taxonRank), ], 
                                         eventID),
                                estimatedAbundance = sum(estimatedAbundance),
                                richness = length(unique(scientificName))),
                      by = "eventID")

class(sampling) <- "data.frame"

sampling$startCollectDate <- as.Date(sampling$startCollectDate)
sampling$adjustedAbundance <- sampling$estimatedAbundance/sampling$totalTrapHours

# Species table
# EventID is unique to each bout of sampling
sppList <- sort(unique(id$scientificName[grepl("species", id$taxonRank)]))
spp <- unique.data.frame(sampling[, c("siteID", "eventID", "startCollectDate")])
spp <- cbind(spp, matrix(data = rep(0, length(sppList)*dim(spp)[1]), 
                             nrow = dim(spp)[1], ncol = length(sppList), 
                             dimnames = list(c(1:dim(spp)[1]), sppList)))
for (i in spp$eventID){
  rows <- which(id$eventID==i)
  spp$estimatedAbundance[spp$eventID==i] <- sum(id$estimatedAbundance[rows])
  spp$richness[spp$eventID==i] <- length(unique(id$scientificName[id$eventID==i & 
                                                                  grepl("species", 
                                                                        id$taxonRank)]))
  for (j in which(colnames(spp)%in%unique(id$scientificName))){
    spp[which(spp$eventID==i), j] <- round(sum(id$estimatedAbundance[id$eventID==i & 
                                                               grepl(colnames(spp)[j], 
                                                                     id$scientificName)])/
                                             spp$estimatedAbundance[spp$eventID==i], 2)
  }
}


# Pathogen status
pathogenStatus <- summarize(group_by(path[path$testNumber==1, ], testingID),
                            siteID = unique(siteID),
                            eventID = unique(eventID),
                            scientificName = unique(scientificName),
                            numMosquitoesPresent = sum(id$estimatedAbundance
                                                       [id$testingID == testingID]),
                            numMosquitoesTested = sum(poolSize),
                            relativeAbundance = round(numMosquitoesPresent/
                              sampling$estimatedAbundance[sampling$eventID == eventID], 
                              3))

pathogenStatus <- left_join(pathogenStatus, 
                            summarize(group_by(path[path$finalResult=="Y" &
                                                    path$testResult=="Positive", ], 
                                               testingID),
                                      maxPositiveMosquitoes = sum(poolSize),
                                      minPositiveMosquitoes = length(poolSize),
                                      testPathogenName = unique(testPathogenName)),
                            by = "testingID")

class(pathogenStatus) <- "data.frame"

pathogenStatus$maxPropPositive <- ((pathogenStatus$maxPositiveMosquitoes/
                                     pathogenStatus$numMosquitoesTested)*
                                     pathogenStatus$relativeAbundance)
pathogenStatus$minPropPositive <- ((pathogenStatus$minPositiveMosquitoes/
                                     pathogenStatus$numMosquitoesTested)*
                                     pathogenStatus$relativeAbundance)

sampling <- left_join(sampling, 
                      pathogenStatus[!is.na(pathogenStatus$maxPropPositive), 
                                     c(3, 11, 12)],
                      by = "eventID")

sampling[is.na(sampling$maxPropPositive), 10:11] <- 0

# Plotting Abundance and Richness
par(mfrow = c(2, 2))
for (i in unique(sampling$siteID)[c(3, 1, 2, 4)]){
  par(mar = c(5, 5, 2, 5))
  with(sampling[sampling$siteID==i, ], 
       plot(startCollectDate, estimatedAbundance, 
            type="l", col="darksalmon",
            ylab="Abundance", main = i, lwd = 2,
            xlim = c(min(sampling$startCollectDate),
                     max(sampling$startCollectDate)),
            xlab = "Date"))
  par(new = T)
  with(sampling[sampling$siteID==i, ], 
       plot(startCollectDate, richness,
            pch = 21, axes = F, xlab = NA, cex = 1.5,
            xlim = c(min(sampling$startCollectDate),
                     max(sampling$startCollectDate)),
            ylab = NA, cex.lab = 1, bg = "turquoise"))
  axis(side = 4)
  mtext(side = 4, line = 3, 'Richness')
  if(i == "STER"){
    legend("topleft",
           legend = c("Abundance", "Richness"),
           lty = c(1,0), pch=c(NA, 16), lwd = c(1.5, 0), 
           col = c("darksalmon", "turquoise"))
  }
}

# Species Composition
tenMosquitoes <- sort(sapply(spp[, 4:45], max), 
                      decreasing = T)[1:10]

for (i in unique(spp$siteID)[c(3,1,2,4)]){
  par(mar = c(5, 5, 2, 5))
  plot(spp[spp$siteID==i, "startCollectDate"], 
       spp[spp$siteID==i, 4], 
       type="n", ylim = c(0, 1), ylab="Relative abundance", 
       main = i, xlab = "Date")
  t = 1
  for (j in sort(names(tenMosquitoes))){
    par(new = T)
    col <- colSums(spp[, 4:45])
    plot(spp[spp$siteID==i, "startCollectDate"], 
         spp[spp$siteID==i, j], 
         type="l", col = rainbow(10)[t], 
         ylim = c(0, 1),
         xlim = c(min(sampling$startCollectDate),
                  max(sampling$startCollectDate)),
         ylab = NA, main = i, lwd = 3,
         xlab = NA, axes = F)
    t = t + 1
  }
  if (i == "STER"){
    legend(x = "topleft", legend = sort(names(tenMosquitoes)),
           col = rainbow(10), pch = rep(NA, 10),
           lwd = 2, lty = 1)
  }
}



# Relationship between Abundance & Richness
par(mfrow = c(1, 1))
plot(x = sampling$estimatedAbundance,
     y = sampling$richness,
     xlab = "Estimated mosquito abundance",
     ylab = "Richness", pch = 21, bg = "forestgreen",
     cex = 1.5)

# Pathogen status within the mosquito community
# Central Plains Experimental Range 
with(sampling[sampling$siteID=="CPER", ], 
     plot(x = startCollectDate,
          y = estimatedAbundance,
          xlab = "Date", main = "Central Plains Experimental Range (CPER)",
          ylab = "Estimated mosquito abundance", 
          type = "l", lwd = 2,
          col = "black"))
par(new = T)
with(sampling[sampling$siteID=="CPER", ], 
     plot(startCollectDate, maxPropPositive,
          pch = 21, axes = F, xlab = NA, cex = 1.5,
          ylim = c(0, .1),
          ylab = NA, cex.lab = 1, bg = "darksalmon"))
par(new = T)
with(sampling[sampling$siteID=="CPER", ], 
     plot(startCollectDate, minPropPositive,
          pch = 21, axes = F, xlab = NA, cex = 1.5,
          ylim = c(0, .1),
          ylab = NA, cex.lab = 1, bg = "forestgreen"))
abline(h = .001, lty = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Proportion Positive for West Nile')
legend("topleft",
       legend = c("Abundance (all mosquitoes)",  
                  "Expected infection rate",
                  "Measured max infection rate",
                  "Measured min infection rate"),
       lty = c(1, 2, 0, 0), pch=c(NA, NA, 16, 16), lwd = c(1.5, 1, 0, 0), 
       col = c("black", "black", "darksalmon", "forestgreen"))
box()

# North Sterling, CO 
with(sampling[sampling$siteID=="STER", ], 
     plot(x = startCollectDate,
          y = estimatedAbundance,
          xlab = "Date", main = "North Sterling, CO (STER)",
          ylab = "Estimated mosquito abundance", 
          type = "l", lwd = 2,
          col = "black"))
par(new = T)
with(sampling[sampling$siteID=="STER", ], 
     plot(startCollectDate, maxPropPositive,
          pch = 21, axes = F, xlab = NA, cex = 1.5,
          ylim = c(0, .5),
          ylab = NA, cex.lab = 1, bg = "darksalmon"))
par(new = T)
with(sampling[sampling$siteID=="STER", ], 
     plot(startCollectDate, minPropPositive,
          pch = 21, axes = F, xlab = NA, cex = 1.5,
          ylim = c(0, .5),
          ylab = NA, cex.lab = 1, bg = "forestgreen"))
abline(h = .001, lty = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Proportion Positive for West Nile')
legend("topleft",
       legend = c("Abundance (all mosquitoes)",  
                  "Expected infection rate",
                  "Measured max infection rate",
                  "Measured min infection rate"),
       lty = c(1, 2, 0, 0), pch=c(NA, NA, 16, 16), lwd = c(1.5, 1, 0, 0), 
       col = c("black", "black", "darksalmon", "forestgreen"))
box()


# Ordway-Swisher Biological Station
with(sampling[sampling$siteID=="OSBS", ], 
     plot(x = startCollectDate,
          y = estimatedAbundance,
          xlab = "Date", main = "Ordway-Swisher Biological Station (OSBS)",
          ylab = "Estimated mosquito abundance", 
          type = "l", lwd = 2,
          col = "black"))
par(new = T)
with(sampling[sampling$siteID=="OSBS", ], 
     plot(startCollectDate, maxPropPositive,
          pch = 21, axes = F, xlab = NA, cex = 1.5,
          ylim = c(0, .5),
          ylab = NA, cex.lab = 1, bg = "darksalmon"))
par(new = T)
with(sampling[sampling$siteID=="OSBS", ], 
     plot(startCollectDate, minPropPositive,
          pch = 21, axes = F, xlab = NA, cex = 1.5,
          ylim = c(0, .5),
          ylab = NA, cex.lab = 1, bg = "forestgreen"))
abline(h = .003, lty = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Proportion Positive for Eastern Equine')
legend("topleft",
       legend = c("Abundance (all mosquitoes)",  
                  "Expected infection rate",
                  "Measured max infection rate",
                  "Measured min infection rate"),
       lty = c(1, 2, 0, 0), pch=c(NA, NA, 16, 16), lwd = c(1.5, 1, 0, 0), 
       col = c("black", "black", "darksalmon", "forestgreen"))
box()
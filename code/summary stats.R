rm(list=ls())
library(dplyr)
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

# Relationship between Abundance & Richness
par(mfrow = c(1, 1))
plot(x = sampling$estimatedAbundance,
     y = sampling$richness,
     xlab = "Estimated mosquito abundance",
     ylab = "Richness", pch = 21, bg = "forestgreen",
     cex = 1.5)

# Plotting Abundance and Richness
for (i in unique(sampling$siteID)[c(3, 1, 2, 4)]){
  par(mfrow = c(2, 2), mar = c(5, 5, 2, 5))
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

# Pathogen status within the mosquito community
par(mfrow = c(1, 1))
plot(x = sampling$estimatedAbundance,
     y = sampling$richness,
     xlab = "Estimated mosquito abundance",
     ylab = "Richness", pch = 21, bg = "forestgreen",
     cex = 1.5)
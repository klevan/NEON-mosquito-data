rm(list=ls())
multipleCombine <- function(input, ply = llply){
  require(plyr)
  require(dplyr)
  ply(input, function(x){
    t <- read.table(x, header = TRUE, sep = ",",
                    stringsAsFactors = FALSE,
                    encoding = "UTF-8")
    t1 <- rbind(t) 
    return(t1)
  }
  )
}

# Getting the data
myPathToData <- "~/GitHub/NEON-mosquito-data/data" # Rename to your own directory location
allFiles <- list.files(path = myPathToData, full.names = T)
# Unzip the files
for (i in 1:length(allFiles[grepl("zip", allFiles)])){
  unzip(zipfile = allFiles[grepl("zip", allFiles)][i], list = FALSE, 
        junkpaths = TRUE, exdir = myPathToData, overwrite = TRUE)
}
# Combine files from multiple sites
allFiles <- list.files(path = myPathToData, full.names = T)
samp <- multipleCombine(allFiles[grepl("sampling", allFiles)], ply = ldply)
id <- multipleCombine(allFiles[grepl("identification", allFiles)], ply = ldply)
arc <- multipleCombine(allFiles[grepl("archival", allFiles)], ply = ldply)
path <- multipleCombine(allFiles[grepl("pathogen", allFiles)], ply = ldply)

# Assign taxonRank
id$taxonRank <- "species"
id$taxonRank[grepl("sp\\.", id$scientificName)] <- "genus"
id$taxonRank[grepl("(\\s[a-zA-Z]+\\s)", id$scientificName)] <- "subspecies"
id$taxonRank[grepl("Culicidae", id$scientificName)] <- "family"

# Find plotID
samp$plotID[samp$sampleID!=""] <- substr(samp$sampleID[samp$sampleID!=""], 1, 8)
id$plotID[id$sampleID!=""] <- substr(id$sampleID[id$sampleID!=""], 1, 8)
location <- unique.data.frame(samp[nchar(samp$plotID)>0, 
                                   c("plotID", "decimalLatitude", 
                                     "decimalLongitude")])
for (i in unique(location$plotID)){
  lat <- location$decimalLatitude[which(i==location$plotID)]
  lon <- location$decimalLongitude[which(i==location$plotID)]
  samp$plotID[samp$decimalLatitude==lat & samp$decimalLongitude==lon] <- i
  id$plotID[id$decimalLatitude==lat & id$decimalLongitude==lon] <- i
}
rm(location)

# Determine sampleID
id$sampleID <- substr(id$subsampleID, 1, 22)

# Remove extra files
allFiles <- list.files(path = myPathToData, full.names = T)
allFiles <- allFiles[!grepl("zip",allFiles) & !grepl("txt",allFiles) & !grepl("variables",allFiles)]
file.remove(allFiles) # tidying the number of files

# Write files that are ready for analysis
write.csv(samp, file = paste(myPathToData, "mos_samplingeffort_pub.csv", sep = "/"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(id, file = paste(myPathToData, "mos_identification_pub.csv", sep = "/"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(arc, file = paste(myPathToData, "mos_archival_pub.csv", sep = "/"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(path, file = paste(myPathToData, "mos_pathogenresults_pub.csv", sep = "/"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")
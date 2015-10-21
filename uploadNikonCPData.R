library(synapseClient)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(rGithubClient)

synapseLogin()

repo <- getRepo("MEP-LINCS/MEP_LINCS_Pilot", ref="branch", refName="accessSynapse")
thisScript <- getPermlink(repo, "uploadNikonCPData.R")

synapseRawDataDir <- "syn4997435"


# Take row of data frame and remove file name
# Convert to a list to use as Synapse annotations
toAnnotationList <- function(x) {
  as.list(x %>% select(-filename))
}

# Take row of data frame with filename and annots
# Upload to Synapse and set annotations
uploadToSynapse <- function(x, parentId) {
  annots <- toAnnotationList(x)
  obj <- File(x$filename, parentId=parentId)
  synSetAnnotations(obj) <- annots
  obj <- synStore(obj, forceVersion=FALSE, 
                  activityName="Upload", executed=thisScript)
  obj
}

for(cellLine in c("PC3", "MCF7", "YAPC")){
  for (ss in c("SS1", "SS2", "SS3")){
    dataDir <- paste("/Users/dane/Documents/MEP-LINCS",cellLine,ss,"RawData/v1", sep = "/")
    # Take file names and turn into basic annotation set
    # Replace this with a better way to get basic annotations from 
    # a standardized source
    dataFiles <- data.frame(filename=list.files(path=dataDir, full.names = TRUE), stringsAsFactors = FALSE) %>%
      mutate(level=0,
             CellLine=cellLine,
             StainingSet=ss,
             Filename=as.character(filename),
             basename=str_replace(filename, ".*/", "")) %>% 
      mutate(basename=str_replace(basename, "\\.csv", "")) %>% 
      separate(basename, c("Barcode", "Well", "Location"))
    
    res <- dlply(dataFiles[,], .(filename), uploadToSynapse, parentId=synapseRawDataDir)
    
  }
}




#title: "MEP-LINCs Preprocessing"
#author: "Mark Dane"
# 9/30/2015

##Introduction

#   The MEP-LINCs dataset contains imaging data from a Nikon automated microscope that is analyzed with a CellProfiler pipeline.
# 
# Part of this preprocessing of the dataset will be deprecated when the merging of the data and metadata happens within the CellProfiler part of the pipeline. For now, the metadata about the ECM proteins is read from the GAL file and the metadata about the wells (cell line, stains and ligands) is read from Excel spreadsheets.

source("MEPLINCSFunctions.R")
library("limma")#read GAL file and strsplit2
library("MEMA")#merge, annotate and normalize functions
library("data.table")#fast file reads, data merges and subsetting
library("parallel")#use multiple cores for faster processing

#Select a staining set
ss <- "SS3"
#Select a CellLine
cellLine <- "PC3"
#select analysis version
analysisVersion <- "v1"

#Rules-based classifier thresholds for perimeter cells
neighborsThresh <- 0.4 #Gates sparse cells on a spot
wedgeAngs <- 20 #Size in degrees of spot wedges used in perimeter gating
outerThresh <- 0.5 #Defines out cells used in perimeter gating
neighborhoodNucleiRadii <- 7 #Defines the neighborhood annulus

#Filter out debris based on nuclear area
nuclearAreaThresh <- 50
nuclearAreaHiThresh <- 4000

#Only process a curated set of the data
curatedOnly <- TRUE
curatedCols <- "ImageNumber|ObjectNumber|_Area$|_Eccentricity|_Perimeter|_MedianIntensity_|_IntegratedIntensity_|_Center_|_PA_"

#Flag to control files updates
writeFiles <- TRUE

#Process a subset of the data to speed development
limitBarcodes <- 8

#Do not normalized to Spot level
normToSpot <- TRUE

#QA flags are used to enable analyses that require minimum cell and
#replicate counts

#Set a threshold for the lowSpotCellCount flag
lowSpotCellCountThreshold <- 5

#Set a threshold for the lowRegionCellCount flag
lowRegionCellCountThreshold <- .4

#Set a threshold for the loess well level QA Scores
lthresh <- 0.6

#Set a threshold for lowWellQA flag
lowWellQAThreshold <- .7

#Set a threshold for the lowSpotReplicates flag
lowReplicateCount <- 3

##Summary
# This script prepares cell-level data and metadata for the MEP LINCs Analysis Pipeline. 
# 
# In the code, the variable ss determines which staining set (SS1, SS2 or SS3) to merge and the variable cellLine determines the cell line (PC3,MCF7, etc). All .txt data files in the "./RawData" folder will be merged with the well (xlsx) and log (XML) data from the "./Metadata" folder.
# 
# The merging assumes that the actual, physical B row wells (B01-B04) have been printed upside-down. That is, rotated 180 degrees resulting in the spot 1, 1 being in the lower right corner instead of the upper left corner. The metadata is matched to the actual printed orientation.

# Read and clean spotmetadata

#Read in the spot metadata from the gal file
smd <- readSpotMetadata(paste0("./",cellLine,"/",ss,"/Metadata/20150515_LI8X001_v1.2.gal"))
#Relabel the column Name to ECMpAnnotID
setnames(smd, "Name", "ECMpAnnotID")

#Add the print order and deposition number to the metadata
ldf <- readLogData(paste0("./",cellLine,"/",ss,"/Metadata/20150512-112336.xml"))
spotMetadata <- merge(smd,ldf, all=TRUE)
setkey(spotMetadata,Spot)
#Make a rotated version of the spot metadata to match the print orientation
spotMetadata180 <- rotateMetadata(spotMetadata)
ARowMetadata <- data.table(spotMetadata,Well=rep(c("A01", "A02","A03","A04"),each=nrow(spotMetadata)))
BRowMetadata <- data.table(spotMetadata180,Well=rep(c("B01", "B02","B03","B04"),each=nrow(spotMetadata180)))

# The well metadata describes the cell line, ligands and staining endpoints that are all added on a per well basis. There is one mutlisheet .xlsx file for each plate. Each filename is the plate's barcode.


# The raw data from all wells in all plates in the dataset are read in and merged with their spot and well metadata. The number of nuclei at each spot are counted and a loess model of the spot cell count is added. Then all intensity values are normalized through dividing them by the median intensity value of the control well in the same plate.
# 
# Next, the data is filtered to remove objects with a nuclear area less than nuclearAreaThresh pixels or more than nuclearAreaHiThresh pixels.

#merge_normalize_QA, echo=FALSE}
#The next steps are to bring in the well metadata, the print order and the CP data

cellDataFiles <- dir(paste0("./",cellLine,"/", ss,"/RawData/",analysisVersion),full.names = TRUE)
splits <- strsplit2(strsplit2(cellDataFiles,split = "_")[,1],"/")

if(limitBarcodes) {
  barcodes <- unique(splits[,ncol(splits)])[1:limitBarcodes] 
} else barcodes <- unique(splits[,ncol(splits)])

expDTList <- lapply(barcodes, function(barcode){
  #browser()
  plateDataFiles <- grep(barcode,cellDataFiles,value = TRUE)
  wells <- unique(strsplit2(split = "_",plateDataFiles)[,2])
  wellDataList <- lapply(wells,function(well){
    #browser()
    wellDataFiles <- grep(well,plateDataFiles,value = TRUE)
    imageDataFile <- grep("Image",wellDataFiles,value=TRUE,
                          ignore.case = TRUE)
    nucleiDataFile <- grep("Nuclei",wellDataFiles,value=TRUE,
                           ignore.case = TRUE)
    if (ss %in% c("SS1","SS3")){
      cellsDataFile <- grep("Cell",wellDataFiles,value=TRUE,
                            ignore.case = TRUE)
      cytoplasmDataFile <- grep("Cytoplasm",wellDataFiles,value=TRUE,
                                ignore.case = TRUE)
    }
    #Read in csv data
    image <- convertColumnNames(fread(imageDataFile))
    setkey(image,CP_ImageNumber)
    nuclei <- convertColumnNames(fread(nucleiDataFile))
    if (curatedOnly) nuclei <- nuclei[,grep(curatedCols,colnames(nuclei)), with=FALSE]
    setkey(nuclei,CP_ImageNumber,CP_ObjectNumber)
    if (ss %in% c("SS1","SS3")){
      cells <- convertColumnNames(fread(cellsDataFile))
      if (curatedOnly) cells <- cells[,grep(curatedCols,colnames(cells)), with=FALSE]
      setkey(cells,CP_ImageNumber,CP_ObjectNumber)
      cytoplasm <- convertColumnNames(fread(cytoplasmDataFile))
      if (curatedOnly) cytoplasm <- cytoplasm[,grep(curatedCols,colnames(cytoplasm)), with=FALSE]
      setkey(cytoplasm,CP_ImageNumber,CP_ObjectNumber)
    }
    
    #Add the data location as a prefix in the column names
    setnames(nuclei,paste0("Nuclei_",colnames(nuclei)))
    if (ss %in% c("SS1","SS3")){
      setnames(cells,paste0("Cells_",colnames(cells)))
      setnames(cytoplasm,paste0("Cytoplasm_",colnames(cytoplasm)))
    }
    
    #Merge the cells, cytoplasm and nuclei data
    if (ss %in% c("SS1","SS3")){
      setkey(cells,Cells_CP_ImageNumber,Cells_CP_ObjectNumber)
      setkey(cytoplasm,Cytoplasm_CP_ImageNumber,Cytoplasm_CP_ObjectNumber)
      setkey(nuclei,Nuclei_CP_ImageNumber,Nuclei_CP_ObjectNumber)
      
      DT <- cells[cytoplasm[nuclei]]
      setnames(DT,"Cells_CP_ImageNumber","Spot")
      setnames(DT,"Cells_CP_ObjectNumber","ObjectNumber")
    } else {
      DT <- nuclei
      setnames(DT,"Nuclei_CP_ImageNumber","Spot")
      setnames(DT,"Nuclei_CP_ObjectNumber","ObjectNumber")
    }
    
    #Add the well name as a parameter
    DT <- DT[,Well := well]
    
    #Merge the data with its metadata based on the row it's in
    m <- regexpr("[[:alpha:]]",well)
    row <- regmatches(well,m)
    setkey(DT,Spot)
    DT <- switch(row, A = merge(DT,spotMetadata,all=TRUE),
                 B = merge(DT,spotMetadata180,all=TRUE))
    
    return(DT)
  })
  
  #Create the cell data.table with spot metadata for the plate 
  pcDT <- rbindlist(wellDataList, fill = TRUE)
  #Read the well metadata from a multi-sheet Excel file
  wellMetadata <- data.table(readMetadata(paste0("./",cellLine,"/",
                                                 ss,"/Metadata/",barcode,".xlsx")), key="Well")
  
  #merge well metadata with the data and spot metadata
  pcDT <- merge(pcDT,wellMetadata,by = "Well")
  pcDT <- pcDT[,Barcode := barcode]
  #Count the cells at each spot
  pcDT<-pcDT[,Spot_PA_SpotCellCount := .N,by="Barcode,Well,Spot"]
  
  pcDT <- pcDT[pcDT$Nuclei_CP_AreaShape_Area > nuclearAreaThresh,]
  pcDT <- pcDT[pcDT$Nuclei_CP_AreaShape_Area < nuclearAreaHiThresh,]
  
  return(pcDT)
})

cDT <- rbindlist(expDTList, fill = TRUE)
#TODO delete unwanted columns here such as Euler Number
densityRadius <- sqrt(median(cDT$Nuclei_CP_AreaShape_Area)/pi)

#Create short display names, then replace where not unique
cDT$ECMp <- gsub("_.*","",cDT$ECMpAnnotID)
cDT$Ligand <- gsub("_.*","",cDT$LigandAnnotID)

#Use entire AnnotID for ligands with same uniprot IDs
cDT$Ligand[grepl("NRG1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "NRG1",annotIDs = cDT$LigandAnnotID[grepl("NRG1",cDT$Ligand)])
cDT$Ligand[grepl("TGFB1",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "TGFB1",annotIDs = cDT$LigandAnnotID[grepl("TGFB1",cDT$Ligand)])
cDT$Ligand[grepl("CXCL12",cDT$Ligand)] <- simplifyLigandAnnotID(ligand = "CXCL12",annotIDs = cDT$LigandAnnotID[grepl("CXCL12",cDT$Ligand)])

cDT$MEP <- paste(cDT$ECMp,cDT$Ligand,sep = "_")

# After merging the metadata with the cell-level data, several types of derived parameters are added. These include:
#   
#   The origin of coordinate system is placed at the median X and Y of each spot and the local cartesian and polar coordinates are added to the dataset.
# 
# The number of nuclei within three nuclear radii of `r densityRadius ` around each nuclei is counted and stored as a neighbor count parameter. The neighbor count value is thresholded to classify each cell as Sparse or not.The distance from the local origin is used to classify each cell as an OuterCell or not. The Sparse, OutCell and Wedge classifications are used to classify each cell as a Perimeter cell or not. 
# 
# For staining set 2, each cell is classified as EdU+ or EdU-. The threshold for EdU+ is based on kmeans threshold of the mean EdU intensity from the control well of each plate.
# 
# The intensity values are normalized at each spot so that spot-level variations can be analyzed.
# 
# The cell level raw data and metadata is saved as Level 1 data. Normalized values are added to the dataset and saved as Level 2 data.

#Add the local polar coordinates and Neighbor Count
cDT <- cDT[,Nuclei_PA_Centered_X :=  Nuclei_CP_AreaShape_Center_X-median(Nuclei_CP_AreaShape_Center_X)]
cDT <- cDT[,Nuclei_PA_Centered_Y :=  Nuclei_CP_AreaShape_Center_Y-median(Nuclei_CP_AreaShape_Center_Y)]
cDT <- cDT[, Nuclei_PA_AreaShape_Center_R := sqrt(Nuclei_PA_Centered_X^2 + Nuclei_PA_Centered_Y^2)]
cDT <- cDT[, Nuclei_PA_AreaShape_Center_Theta := calcTheta(Nuclei_PA_Centered_X, Nuclei_PA_Centered_Y)]
cDT <- cDT[,Nuclei_PA_AreaShape_Neighbors := cellNeighbors(.SD, radius = densityRadius*neighborhoodNucleiRadii), by = "Barcode,Well,Spot"]

#Rules for classifying perimeter cells
cDT <- cDT[,Spot_PA_Sparse := Nuclei_PA_AreaShape_Neighbors < neighborsThresh]

#Add a local wedge ID to each cell based on conversations with Michel Nederlof
cDT <- cDT[,Spot_PA_Wedge:=ceiling(Nuclei_PA_AreaShape_Center_Theta/wedgeAngs)]

#Define the perimeter cell if it exists in each wedge
#Classify cells as outer if they have a radial position greater than a thresh
cDT <- cDT[,Spot_PA_OuterCell := labelOuterCells(Nuclei_PA_AreaShape_Center_R, thresh=outerThresh),by="Barcode,Well,Spot"]

#Require the cell not be in a sparse region
denseOuterDT <- cDT[!cDT$Spot_PA_Sparse  & cDT$Spot_PA_OuterCell]
denseOuterDT <- denseOuterDT[,Spot_PA_Perimeter := findPerimeterCell(.SD) ,by="Barcode,Well,Spot,Spot_PA_Wedge"]
setkey(cDT,Barcode,Well,Spot,ObjectNumber)
setkey(denseOuterDT,Barcode,Well,Spot,ObjectNumber)
cDT <- denseOuterDT[,list(Barcode,Well,Spot,ObjectNumber,Spot_PA_Perimeter)][cDT]
cDT$Spot_PA_Perimeter[is.na(cDT$Spot_PA_Perimeter)] <- FALSE

#Set 2N and 4N DNA status
cDT <- cDT[,Nuclei_PA_Cycle_State := kmeansDNACluster(Nuclei_CP_Intensity_IntegratedIntensity_Dapi), by="Barcode,Well"]
#Manually reset clusters for poorly classified wells
#This is based on review of the clusters after a prior run
cDT$Nuclei_PA_Cycle_State[cDT$Barcode=="LI8X00403" & cDT$Well=="A03" & cDT$Nuclei_CP_Intensity_IntegratedIntensity_Dapi >150] <- 2

cDT <- cDT[,Nuclei_PA_Cycle_2NProportion := calc2NProportion(Nuclei_PA_Cycle_State),by="Barcode,Well,Spot"]
cDT$Nuclei_PA_Cycle_4NProportion <- 1-cDT$Nuclei_PA_Cycle_2NProportion

#Add spot level normalizations for selected intensities
if(normToSpot){
  intensityNamesAll <- grep("_CP_Intensity_Median",colnames(cDT), value=TRUE)
  intensityNames <- grep("Norm",intensityNamesAll,invert=TRUE,value=TRUE)
  for(intensityName in intensityNames){
    #Median normalize the median intensity at each spot
    cDT <- cDT[,paste0(intensityName,"_SpotNorm") := medianNorm(.SD,intensityName),by="Barcode,Well,Spot"]
  }
}

#Create staining set specific derived parameters
if(ss %in% c("SS1")){
  
} else if (ss == "SS2"){
  
  cDT <- cDT[,Nuclei_PA_Gated_EduPositive := kmeansCluster(.SD, value="Nuclei_CP_Intensity_MedianIntensity_Edu", ctrlLigand = "FBS"), by="Barcode"]
  #Calculate the EdU Positive Percent at each spot
  cDT <- cDT[,Nuclei_PA_Gated_EduPositiveProportion := sum(Nuclei_PA_Gated_EduPositive)/length(ObjectNumber),by="Barcode,Well,Spot"]
  
} else if (ss == "SS3"){
  #Calculate a lineage ratio of luminal/basal or KRT19/KRT5
  cDT <- cDT[,Cytoplasm_PA_Intensity_LineageRatio := Cytoplasm_CP_Intensity_MedianIntensity_KRT19/Cytoplasm_CP_Intensity_MedianIntensity_KRT5]
  
} else stop("Invalid ss parameter")

# Eliminate Variations in the Endpoint metadata
endpointNames <- grep("End",colnames(cDT), value=TRUE)
endpointWL <- regmatches(endpointNames,regexpr("[[:digit:]]{3}|DAPI",endpointNames))
setnames(cDT,endpointNames,paste0("Endpoint",endpointWL))

#Identify parameters that shouldn't be normalized
normParameters <- grep("Sparse|Wedge|OuterCell|Spot_PA_Perimeter|Nuclei_PA_Cycle_State",colnames(cDT),value=TRUE,invert=TRUE)

#Save the un-normalized parameters to merge in later
mdKeep <- cDT[,grep("Barcode|^Well$|^Spot$|ObjectNumber|Sparse|Wedge|OuterCell|Spot_PA_Perimeter|Nuclei_PA_Cycle_State",colnames(cDT),value=TRUE), with = FALSE]

#Normalized each feature by dividing by the median of its plate's FBS value well
#cDT <- normDataset(cDT[,normParameters,with=FALSE])
#Normalize each feature by subtracting the median of its plate's FBS value
# and divding by its plates MAD
cDT <- normRZSDataset(cDT[,normParameters, with = FALSE])
cDT <- merge(cDT,mdKeep)

#The cell-level data is median summarized to the spot level and coefficients of variations on the replicates are calculated. The spot level data and metadata are saved as Level 3 data.

#### Level3 ####

#Summarize cell data to medians of the spot parameters
parameterNames<-grep(pattern="(Children|_CP_|_PA_|Barcode|^Spot$|^Well$)",x=names(cDT),value=TRUE)

#Remove any spot-normalized and cell level parameters
parameterNames <- grep("SpotNorm|^Nuclei_PA_Gated_EduPositive$|^Nuclei_PA_Gated_EduPositive_RZSNorm$",parameterNames,value=TRUE,invert=TRUE)

#Remove any raw parameters
parameterNames <- grep("Barcode|^Spot$|^Well$|Norm|Nuclei_CP_Intensity_MedianIntensity_Dapi$|Cytoplasm_CP_Intensity_MedianIntensity_Actin$|Cytoplasm_CP_Intensity_MedianIntensity_CellMask$|Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker$|Nuclei_CP_Intensity_MedianIntensity_H3$|Nuclei_CP_Intensity_MedianIntensity_Fibrillarin$|Nuclei_CP_Intensity_MedianIntensity_Edu$|Cytoplasm_CP_Intensity_MedianIntensity_KRT5$|Cytoplasm_CP_Intensity_MedianIntensity_KRT19$|Spot_PA_SpotCellCount$", parameterNames, value = TRUE)

cDTParameters<-cDT[,parameterNames,with=FALSE]

slDT<-cDTParameters[,lapply(.SD,numericMedian),keyby="Barcode,Well,Spot"]
slDTse <- cDTParameters[,lapply(.SD,se),keyby="Barcode,Well,Spot"]

#Add _SE to the standard error column names
setnames(slDTse, grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$",colnames(slDTse), value = TRUE, invert = TRUE),"_SE"))

#Merge back in the spot and well metadata
#TODO: Convert the logic to not name the metadata
metadataNames <- grep("(Row|Column|PrintOrder|Block|^ID$|Array|CellLine|Ligand|Endpoint|ECMp|MEP|Barcode|^Well$|^Spot$)", x=colnames(cDT), value=TRUE)
setkey(cDT,Barcode, Well,Spot)
mDT <- cDT[,metadataNames,keyby="Barcode,Well,Spot", with=FALSE]
slDT <- mDT[slDT, mult="first"]
#Merge in the standard errr values
slDT <- slDTse[slDT]
#Add a count of replicates
slDT <- slDT[,Spot_PA_ReplicateCount := .N,by="LigandAnnotID,ECMpAnnotID"]

#Add the loess model of the SpotCellCount on a per well basis
slDT <- slDT[,Spot_PA_LoessSCC := loessModel(.SD, value="Spot_PA_SpotCellCount", span=.5), by="Barcode,Well"]

#Add well level QA Scores to spot level data
slDT <- slDT[,QAScore := calcQAScore(.SD,threshold=lthresh,maxNrSpot = max(cDT$ArrayRow)*max(cDT$ArrayColumn),value="Spot_PA_LoessSCC"),by="Barcode,Well"]

#Add QA scores to cell level data
setkey(cDT,Barcode, Well, Spot)
setkey(slDT, Barcode, Well, Spot)
cDT <- cDT[slDT[,list(Barcode, Well, Spot, QAScore, Spot_PA_LoessSCC)]]
#The spot level data is median summarized to the replicate level and is stored as Level 4 data and metadata.

#Level4Data

#Summarize spot level data to MEP level by taking the medians of the parameters
#TODO convert spot level QA flags to proportions
mepNames<-grep("Norm|LigandAnnotID|ECMpAnnotID|Barcode|ReplicateCount", x=names(slDT),value=TRUE)
#remove the _SE values
mepNames <- grep("_SE",mepNames, value = TRUE, invert = TRUE)

mepKeep<-slDT[,mepNames,with=FALSE]
mepDT<-mepKeep[,lapply(.SD,numericMedian),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
mepDTse <- mepKeep[,lapply(.SD,se),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
#Add _SE to the standard error column names
setnames(mepDTse, grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(mepDTse), value = TRUE, invert = TRUE), paste0(grep("Barcode|^Well$|^Spot$|Ligand|ECMp",colnames(mepDTse), value = TRUE, invert = TRUE),"_SE"))

#Merge back in the replicate metadata
mDT <- slDT[,list(Well,CellLine,Ligand,Endpoint488,Endpoint555,Endpoint647,EndpointDAPI,ECMp,MEP),keyby="LigandAnnotID,ECMpAnnotID,Barcode"]
mepDT <- mDT[mepDT, mult="first"]
mepDT <- mepDTse[mepDT]


#Write QA flags into appropriate data levels
#Low cell count spots
cDT$QA_LowSpotCellCount <- cDT$Spot_PA_SpotCellCount < lowSpotCellCountThreshold
slDT$QA_LowSpotCellCount <- slDT$Spot_PA_SpotCellCount < lowSpotCellCountThreshold

#Manually flag low quality DAPI wells
cDT$QA_LowDAPIQuality <- FALSE
cDT$QA_LowDAPIQuality[cDT$Barcode=="LI8X00420"&cDT$Well=="B01"] <- TRUE
slDT$QA_LowDAPIQuality <- FALSE
slDT$QA_LowDAPIQuality[slDT$Barcode=="LI8X00420"& slDT$Well=="B01"] <- TRUE

#Flag spots below automatically loess QA threshold
cDT$QA_LowRegionCellCount <- cDT$Spot_PA_LoessSCC < lowRegionCellCountThreshold
slDT$QA_LowRegionCellCount <- slDT$Spot_PA_LoessSCC < lowRegionCellCountThreshold

#Flag wells below automatically calculated QA threshold
slDT$QA_LowWellQA <- FALSE
slDT$QA_LowWellQA[slDT$QAScore < lowWellQAThreshold] <- TRUE
cDT$QA_LowWellQA <- FALSE
cDT$QA_LowWellQA[cDT$QAScore < lowWellQAThreshold] <- TRUE

#Level 4
mepDT$QA_LowReplicateCount <- mepDT$Spot_PA_ReplicateCount < lowReplicateCount

#WriteData

if(writeFiles){
  #Write out cDT without normalized values as level 1 dataset
  level1Names <- grep("Norm",colnames(cDT),value=TRUE,invert=TRUE)
  write.table(format(cDT[,level1Names, with=FALSE], digits=4), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(cDT$CellLine),"_",ss,"_","Level1.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
  
  normParmameterNames <- grep("Norm",colnames(cDT), value=TRUE)
  rawParameterNames <- gsub("_?[[:alnum:]]*?Norm$", "", normParmameterNames)
  metadataNormNames <- colnames(cDT)[!colnames(cDT) %in% rawParameterNames]
  #Paste back in the QA and selected raw data
  
  level2Names <- c(metadataNormNames,
                   grep("Nuclei_CP_Intensity_MedianIntensity_Dapi$|Cytoplasm_CP_Intensity_MedianIntensity_Actin$|Cytoplasm_CP_Intensity_MedianIntensity_CellMask$|Cytoplasm_CP_Intensity_MedianIntensity_MitoTracker$|Nuclei_CP_Intensity_MedianIntensity_H3$|Nuclei_CP_Intensity_MedianIntensity_Firbillarin$|Nuclei_CP_Intensity_MedianIntensity_Edu$|Cytoplasm_CP_Intensity_MedianIntensity_KRT5$|Cytoplasm_CP_Intensity_MedianIntensity_KRT19$|Spot_PA_SpotCellCount$", colnames(cDT), value = TRUE))
  
  #Write out cDT with normalized values as level 2 dataset
  write.table(format(cDT[,level2Names, with = FALSE], digits=4), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(cDT$CellLine),"_",ss,"_","Level2.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
  
  write.table(format(slDT, digits = 4), paste0("./",cellLine,"/", ss, "/AnnotatedData/", unique(slDT$CellLine),"_",ss,"_","Level3.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
  
  write.table(format(mepDT, digits = 4), paste0("./",cellLine,"/",ss, "/AnnotatedData/", unique(slDT$CellLine),"_",ss,"_","Level4.txt"), sep = "\t",row.names = FALSE, quote=FALSE)
  
  #Write the pipeline parameters to  tab-delimited file
  write.table(c(
    ss=ss,
    cellLine = cellLine,
    analysisVersion = analysisVersion,
    neighborhoodNucleiRadii = neighborhoodNucleiRadii,
    neighborsThresh = neighborsThresh,
    wedgeAngs = wedgeAngs,
    outerThresh = outerThresh,
    nuclearAreaThresh = nuclearAreaThresh,
    nuclearAreaHiThresh = nuclearAreaHiThresh,
    curatedOnly = curatedOnly,
    curatedCols = curatedCols,
    writeFiles = writeFiles,
    limitBarcodes = limitBarcodes,
    normToSpot = normToSpot,
    lowSpotCellCountThreshold = lowSpotCellCountThreshold,
    lowRegionCellCountThreshold = lowRegionCellCountThreshold,
    lowWellQAThreshold = lowWellQAThreshold,
    lowReplicateCount =lowReplicateCount,
    lthresh = lthresh
  ),
  paste0("./",cellLine,"/",ss, "/AnnotatedData/", cellLine,"_",ss,"_","PipelineParameters.txt"), sep = "\t",col.names = FALSE, quote=FALSE)
}
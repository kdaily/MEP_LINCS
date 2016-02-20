# Make the pre-release manifest table
library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)

library(synapseClient)
synapseLogin()

annotatedFolder <- "syn5004301"
rawFolder <- "syn4997435"
reportFolder <- "syn5007815"

query <- paste('select id,name,versionNumber,level,Barcode,CellLine,StainingSet,Location,Well',
               'from file where parentId=="%s"')

annotFiles <- synQuery(sprintf(query, annotatedFolder), blockSize = 400)$collectAll()
rawFiles <- synQuery(sprintf(query, rawFolder), blockSize = 400)$collectAll()
reportFiles <- synQuery(sprintf(query, reportFolder), blockSize = 400)$collectAll()

allFiles <- rbind(annotFiles, rawFiles, reportFiles) 
colnames(allFiles) <- gsub(".*\\.", "", colnames(allFiles))

allFiles <- allFiles %>%
  filter(CellLine %in% c("PC3", "YAPC"), StainingSet %in% c("SS2")) %>% 
  mutate(versionNumber=as.numeric(versionNumber)) %>% 
  arrange(level, CellLine, StainingSet) %>% 
  select(id,versionNumber,name,level,Barcode,CellLine,StainingSet,Location,Well)

tableName <- "Pre-release manifest"
tblCols <- as.tableColumns(allFiles)
schema <- TableSchema(name=tableName, columns=tblCols$tableColumns, 
                      parent="syn2862345")
tbl <- synStore(Table(tableSchema = schema, values=allFiles))

# Raw data
allFiles %>% 
  filter(level == 0) %>% 
  mutate(id=paste(id, versionNumber, sep=".")) %>%
  dcast(Barcode + CellLine + StainingSet + Location ~ Well, value.var="id", fill=NA) %>% 
  arrange(CellLine, Location, Barcode, StainingSet) %>% 
  knitr::kable()

# Processed data and analysis
allFiles %>% 
  filter(level != 0) %>% 
  mutate(id=paste(id, versionNumber, sep=".")) %>%
  dcast(CellLine + StainingSet ~ level, value.var="id", fill=NA) %>% 
  knitr::kable()

library("rmarkdown")
#cellLine
#ss is staining set
#s is signal
#m is method

for(cellLine in c("PC3")){
  for(s in c("EdU")){
    for(m in c("NaiveRandRUV"))
      render(paste0("Mep-LINCS_QANorm_",cellLine,"_",s,"_",m,".Rmd"),
             output_file = paste0("Mep-LINCS_QANorm_",cellLine,"_",s,"_",m,".html"),
             output_format = "html_document") 
  }
}
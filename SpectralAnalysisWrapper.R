library("rmarkdown")

analysisVersion <- "v1"
for(cellLine in c("PC3")){
  for(ss in c("SS3")){
    render("MEP-LINCS_Spectral_Analysis.Rmd", output_file = paste0("Mep-LINCS_Spectral_Analysis_",cellLine,"_",analysisVersion,"_",ss,".html"), output_format = "html_document") 
  }
}

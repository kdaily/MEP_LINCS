library("rmarkdown")

for(cellLine in c("HCC1143", "HCC1954")){
  for(ss in c("SS3")){
    render("MEP-LINCS_Analysis.Rmd", output_file = paste0("Mep-LINCS_Analysis_",cellLine,"_",ss,".html"), output_format = "html_document") 
  }
}

library("rmarkdown")

PC3df <- data.frame(cellLine=rep(c("PC3"), 4),
                    ss=c("SS1", "SS2","SS3","SS2noH3"),
                    analysisVersion="av1.4",
                    rawDataVersion=c("v2","v2.1","v2.1", "v1"),
                    stringsAsFactors=FALSE)

MCF7df <- data.frame(cellLine=rep(c("MCF7"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     analysisVersion="av1.4",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

YAPCdf <- data.frame(cellLine=rep(c("YAPC"), 3),
                     ss=c("SS1", "SS2","SS3"),
                     analysisVersion="av1.4",
                     rawDataVersion=c("v2","v2","v2"),
                     stringsAsFactors=FALSE)

MCF10Adf <- data.frame(cellLine=rep(c("MCF10A"), 2),
                       ss=c("SS1","SS3"),
                       analysisVersion="av1.4",
                       rawDataVersion=c("v2","v2"),
                       stringsAsFactors=FALSE)

ssDatasets <- rbind(PC3df,MCF7df,YAPCdf,MCF10Adf)

renderQAReports <- function(x){
  cellLine <- x[["cellLine"]]
  ss <- x[["ss"]]
  rawDataVersion <- x[["rawDataVersion"]]
  analysisVersion <- x[["analysisVersion"]]
  render("MEP-LINCS_QA.Rmd",
         output_file = paste0("./QAReports/MEP-LINCS_QA_",
                              x[["cellLine"]],"_",
                              x[["ss"]],"_",
                              x[["rawDataVersion"]],"_",
                              x[["analysisVersion"]],".html"),
         output_format = "html_document")
}

tmp <- apply(ssDatasets[c(12),], 1, renderQAReports)

downloadData <- function(Species, Sample_Type, geneList) {
  # # Bug Testing
  # Species <- "hsapiens"
  # Sample_Type <- "Normal_Tissue"
  downloadFolder <- system.file("data", package = "correlationAnalyzeR")
  downloadFolder <- file.path(downloadFolder, "Correlation_Data",
                              Species, Sample_Type)
  for (i in 1:length(geneList)) {
    gene <- geneList[i]
    geneFile <- paste0(gene, ".RData")
    geneFilePath <- file.path(downloadFolder, geneFile)
    if (! file.exists(geneFilePath)) {
      s3URL <- file.path("http://packagedata.s3.amazonaws.com",
                         "Correlation_Data", Species, Sample_Type, geneFile)
      download.file(url = "http://packagedata.s3.amazonaws.com/geneCorrelationDataFiles/ATM.RData",
                    destfile = geneFilePath)
    }
  }
}









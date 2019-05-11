#' Download Gene Correlation Data
#'
#' Downloads gene correlation data when requested from AWS S3 bucket
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#'     Either "All", "Normal_Tissues", or "Tumor_Tissues".
#' @param geneList Vector of genes for which data will be extracted.
#'
#' @return Downloads files to package data storage.
#'
#' @examples
#' downloadData(Species = "hsapiens", Sample_Type = "Normal_Tissues",
#'     geneList = c("ATM", "BRCA1"))
#'
#' @export
downloadData <- function(Species, Sample_Type, geneList) {
  # # Bug Testing
  # Species <- "hsapiens"
  # Sample_Type <- "Normal_Tissues"
  dataFolder <- system.file("data", package = "correlationAnalyzeR")
  downloadFolder <- file.path(dataFolder, "Correlation_Data",
                              Species[1], Sample_Type[1])
  if (! dir.exists(downloadFolder)) {
    dir.create(file.path(dataFolder, "Correlation_Data"))
    dir.create(file.path(dataFolder, "Correlation_Data",
                         Species))
    dir.create(file.path(dataFolder, "Correlation_Data",
              Species, Sample_Type))
  }
  for (i in 1:length(geneList)) {
    gene <- geneList[i]
    geneFile <- paste0(gene, ".RData")
    geneFilePath <- file.path(downloadFolder, geneFile)
    if (! file.exists(geneFilePath)) {
      s3URL <- file.path("http://packagedata.s3.amazonaws.com",
                         "Correlation_Data", Species, Sample_Type, geneFile)
      download.file(url = s3URL,
                    destfile = geneFilePath)
    }
  }
}









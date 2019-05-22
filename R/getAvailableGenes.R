#' Get Available Genes
#'
#' Finds available genes within correlation data
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @return Dataframe containing available genes and NCBI descriptions.
#'
#' @examples
#' getAvailableGenes("hsapiens")
#'
#' @export
getAvailableGenes <- function(Species = c("hsapiens", "mmusculus")) {

  # # # Bug testing
  # Species <- "hsapiens"

  # Specify information about the download location and species type
  downloadFolder <- system.file("data", package = "correlationAnalyzeR")
  if (Species[1] == "hsapiens") {
    gene <- "A1BG"
    geneFile <- "A1BG.RData"
  } else if (Species[1] == "mmusculus") {
    gene <- "A1bg"
    geneFile <- "A1bg.RData"
  } else {
    stop("\ncorrelationAnalyzeR currently supports only Human and Mouse data.
         Please select either 'hsapiens' or 'mmusculus' for Species parameter.
         \n")
  }
  Sample_Type = "Normal_Tissues" # This is default behavior
  # Download a sample file which contains all gene identifiers
  downloadData(Species = Species[1],
               Sample_Type = Sample_Type,
               geneList = gene)
  # Construct file for sample data
  file <- file.path(downloadFolder, "Correlation_Data",
                    Species[1], Sample_Type, geneFile)
  # load sample data
  load(file)
  geneNames <- names(vec)
  geneNamesDF <- as.data.frame(geneNames)
  colnames(geneNamesDF)[1] <- "geneName"
  if (Species[1] == "hsapiens") {
    data("humanGenes", package = "correlationAnalyzeR")
    colnames(humanGenes)[1] <- "geneName"
    geneNamesDF <- merge(x = humanGenes, y = geneNamesDF, by = "geneName",
                         all.y = T)
  } else if (Species[1] == "mmusculus") {
    data("mouseGenes", package = "correlationAnalyzeR")
    colnames(MouseGenes)[1] <- "geneName"
    geneNamesDF <- merge(x = MouseGenes, y = geneNamesDF, by = "geneName",
                         all.y = T)
  }
  return(geneNamesDF)
}

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
#' correlationAnalyzeR::getAvailableGenes("hsapiens")
#'
#' @export
getAvailableGenes <- function(Species = c("hsapiens", "mmusculus")) {


  # Specify information about the download location and species type
  if (Species[1] == "hsapiens") {
    gene <- "A1BG"
  } else if (Species[1] == "mmusculus") {
    gene <- "A1bg"
  } else {
    stop("\ncorrelationAnalyzeR currently supports only Human and Mouse data.
         Please select either 'hsapiens' or 'mmusculus' for Species parameter.
         \n")
  }
  Sample_Type = "normal" # This is default behavior
  Tissue = "brain"
  # Download a sample file which contains all gene identifiers
  geneNamesDF <- correlationAnalyzeR::getCorrelationData(Species = Species[1],
                                                         Sample_Type = Sample_Type,
                                                         Tissue = Tissue,
                                                         geneList = gene)
  geneNamesDF <- as.data.frame(rownames(geneNamesDF))
  colnames(geneNamesDF)[1] <- "geneName"
  if (Species[1] == "hsapiens") {
    humanGenes <- correlationAnalyzeR::humanGenes
    colnames(humanGenes)[1] <- "geneName"
    geneNamesDF <- merge(x = humanGenes, y = geneNamesDF, by = "geneName",
                         all.y = T)
  } else if (Species[1] == "mmusculus") {
    mouseGenes <- correlationAnalyzeR::mouseGenes
    colnames(mouseGenes)[1] <- "geneName"
    geneNamesDF <- merge(x = mouseGenes, y = geneNamesDF, by = "geneName",
                         all.y = T)
  }
  return(geneNamesDF)
}

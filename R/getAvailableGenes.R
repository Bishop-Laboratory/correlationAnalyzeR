#' Get Available Genes
#'
#' Finds available genes within correlation data
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied.
#' @return A vector of genes with associated correlation data
#'
#' @examples
#' correlationAnalyzeR::getAvailableGenes("hsapiens")
#'
#' @export
getAvailableGenes <- function(Species = c("hsapiens", "mmusculus"), pool = NULL) {

  # # Bug testing
  # Species <- "hsapiens"
  # pool <- NULL

  if (! is.null(pool)) {
    if (! pool$valid) {
      pool <- NULL
    }
  }

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
                                                         geneList = gene, pool = pool)
  avGenes <- rownames(geneNamesDF)
  return(avGenes)
}

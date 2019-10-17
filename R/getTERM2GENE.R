#' Obtains TERM2GENE object for corGSEA
#'
#' Wrapper for msgidb::msigdbr() function
#'
#' @param Species Species to obtain gene names for.
#' Either 'hsapiens' or 'mmusculus'
#'
#' @param GSEA_Type Whether GSEA should consider all msigdb annotations,
#' or just those in the most popular categories. Should be one of either
#' 'simple' or 'complex'.
#'
#' @return A tbl object with columns "gs_name" and "gene_symbol"
#'
#' @examples
#' correlationAnalyzeR::getTERM2GENE(Species = "hsapiens", GSEA_Type = "simple")
#'
#' @export
getTERM2GENE <- function(GSEA_Type = c("simple", "complex"),
                         Species = c("hsapiens", "mmusculus")) {

  # Species = "mmusculus"
  # GSEA_Type = "simple"

  if (Species == "hsapiens") {
    msigSpec <- "Homo sapiens"
  } else {
    msigSpec <- "Mus musculus"
  }

  if (! GSEA_Type %in% c("simple", "complex")) {
    stop("\nPlease enter either 'simple' or 'complex' for GSEA_Type\n")
  } else if (GSEA_Type[1] == "simple") {
    categories <- list("H" = c("all"),
                       "C2" = c("all"),
                       "C5" = c("BP"),
                       "C6" = c("all"))
    msigList <- list()
    for (i in 1:length(categories)) {
      category <- names(categories)[i]
      subcategory <- categories[[i]]
      if (category != "C5") {
        msigList[[i]] <- msigdbr::msigdbr(msigSpec, category = category)
      } else {
        msigList[[i]] <- msigdbr::msigdbr(msigSpec,
                                          category = category,
                                          subcategory = subcategory)
      }
      m_dfFinal <- data.table::rbindlist(msigList)
      TERM2GENE <- unique(m_dfFinal[,c(1, 8)])
    }

  } else {
    m_df <- msigdbr::msigdbr(msigSpec)
    TERM2GENE <- unique(m_df[,c(1, 8)])
  }
  return(TERM2GENE)
}



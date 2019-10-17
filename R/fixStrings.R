#' Fix strings (helper function)
#'
#' Convert vector of GSEA (or other) names to publication-ready titles
#'
#' @param StringVec A vector of titles (usually GSEA) to clean for visualizations
#' @return A vector of cleaned string titles in the same order as the input.
#' @export


fixStrings <- function(StringVec) {
  StringVec <- gsub(StringVec, pattern = "_", replacement = " ")
  StringVec <- tolower(StringVec)
  StringVec <- stringr::str_to_title(StringVec)
  StringVec <- gsub(StringVec, pattern = "Iii", replacement = "III", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ii", replacement = "II", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Of ", replacement = " of ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " To ", replacement = " to ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " By ", replacement = " by ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " For ", replacement = " for ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Via ", replacement = " via ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " With ", replacement = " with ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Vs ", replacement = " vs ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " And ", replacement = " and ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Rna", replacement = "RNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "rna$", replacement = "RNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "rna ", replacement = "RNA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Dna", replacement = "DNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "dna", replacement = "DNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Tca", replacement = "TCA", ignore.case = FALSE)
  return(StringVec)
}

#' Get available tissue types
#'
#' Finds tissue types with available correlation data for a given species
#'
#' @param Species Species to obtain tissue types for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @param useBlackList Should black-listed tissue/disease categories for this species
#' be removed from the returned list? Improves the quality of analysis by removing
#' categories with low sample numbers and high observed variance.
#'
#' @return Vector containing available tissue types.
#'
#' @examples
#' correlationAnalyzeR::getTissueTypes("hsapiens")
#'
#' @export
getTissueTypes <- function(Species = c("hsapiens", "mmusculus"),
                           useBlackList = FALSE) {

  Species <- Species[1]
  con <- DBI::dbConnect(RMySQL::MySQL(), user = "public-rds-user", port = 3306,
                        dbname="bishoplabdb",
                        password='public-user-password',
                        host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
  tabs <- DBI::dbListTables(con)
  tabs <- tabs[grep(tabs, pattern = paste0("correlations_", Species))]
  tabs <- gsub(tabs, pattern = paste0("correlations_", Species, "_"),
               replacement = "")
  tabs <- gsub(tabs, pattern = "repiratory",
               replacement = "respiratory")
  if (useBlackList) {
    blackListHuman <- correlationAnalyzeR::blackListHuman
    blackListMouse <- correlationAnalyzeR::blackListMouse
    if (Species == "hsapiens") {
      blackList <- blackListHuman
    } else {
      blackList <- blackListMouse
    }

    tabs <- tabs[grep(x = tabs, pattern = paste(blackList, collapse = "|"), invert = TRUE)]
  }
  tabs <- strsplit(tabs, split = "_")
  tissues <- sapply(tabs, "[[", 2)
  types <- sapply(tabs, "[[", 1)
  result <- paste0(tissues, " - ", types)
  result <- result[order(result)]
  DBI::dbDisconnect(con)
  return(result)
}

#' Get tissue types
#'
#' Finds tissue types with available correlation data
#'
#' @param Species Species to obtain tissue types for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @return Vector containing available tissue types.
#'
#' @examples
#' getTissueTypes("hsapiens")
#'
#' @export
getTissueTypes <- function(Species = c("hsapiens", "mmusculus")) {

  # # # # # Bug testing
  # Species <- "hsapiens"
  Species <- Species[1]
  con <- DBI::dbConnect(RMySQL::MySQL(), user = "public-rds-user", port = 3306,
                        dbname="bishoplabdb",
                        password='public-user-password',
                        host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
  tabs <- DBI::dbListTables(con)
  tabs <- tabs[grep(tabs, pattern = paste0("correlations_", Species))]
  tabs <- gsub(tabs, pattern = paste0("correlations_", Species, "_"),
               replacement = "")
  tabs <- strsplit(tabs, split = "_")
  tissues <- sapply(tabs, "[[", 2)
  types <- sapply(tabs, "[[", 1)
  result <- paste0(tissues, " - ", types)
  result <- result[order(result)]
  DBI::dbDisconnect(con)
  return(result)
}

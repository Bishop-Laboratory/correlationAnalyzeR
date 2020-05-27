#' Get available tissue types
#'
#' Finds tissue types with available correlation data for a given species
#'
#' @param useBlackList Should black-listed tissue/disease categories for this species
#' be removed from the returned list? Improves the quality of analysis by removing
#' categories with low sample numbers and high observed variance.
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied.
#' @return Vector containing available tissue types.
#'
#' @examples
#' correlationAnalyzeR::getTissueTypes()
#'
#' @export
getTissueTypes <- function(#Species = c("hsapiens", "mmusculus"),
                           useBlackList = FALSE, pool = NULL) {

  # # bug testing
  # Species = "hsapiens"
  # useBlackList = TRUE
  # pool = NULL

  Species = "hsapiens"

  if (! is.null(pool)) {
    if (! pool$valid) {
      pool <- NULL
    } else {
      conn <- pool::poolCheckout(pool)
      doPool <- TRUE
      on.exit(pool::poolReturn(conn))
    }
  }

  if (is.null(pool)) {
    doPool <- FALSE
    conn <- NULL
    retryCounter <- 1
    while(is.null(conn)) {
      conn <- try(silent = T, eval({
        DBI::dbConnect(
          drv = RMySQL::MySQL(),
          user = "public-rds-user@m2600az-db01p.mysql.database.azure.com", port = 3306,
          dbname="correlation_analyzer",
          password='public-user-password',
          host="m2600az-db01p.mysql.database.azure.com"
        )
      }))
      if ("try-error" %in% class(conn)) {
        if (retryCounter == 3) {
          stop("Unable to connect to database. Check internet connection and please contanct",
               " package maintainer if you believe this is an error.")
        }
        warning(paste0("Failed to establish connection to database ... retrying now ... ",
                       (4-retryCounter), " attempts left."),
                immediate. = T)
        conn <- NULL
        retryCounter <- retryCounter + 1
        Sys.sleep(1)
      }
    }

    on.exit(DBI::dbDisconnect(conn))
  }

  Species <- Species[1]
  tabs <- DBI::dbListTables(conn)
  tabs <- tabs[grep(tabs, pattern = paste0("correlations_", Species))]
  tabs <- gsub(tabs, pattern = paste0("correlations_", Species, "_"),
               replacement = "")
  if (useBlackList) {
    blackList <- correlationAnalyzeR::blackListHuman
    tabs <- tabs[grep(x = tabs, pattern = paste(blackList, collapse = "|"), invert = TRUE)]
  }
  tissues <- gsub(tabs, pattern = "^([a-z]+)_([a-z_]+)$", replacement = "\\2")
  types <- gsub(tabs, pattern = "^([a-z]+)_([a-z_]+)$", replacement = "\\1")
  result <- paste0(tissues, " - ", types)
  result <- result[order(result)]
  return(result)
}

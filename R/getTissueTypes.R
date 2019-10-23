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
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied.
#' @return Vector containing available tissue types.
#'
#' @examples
#' correlationAnalyzeR::getTissueTypes("hsapiens")
#'
#' @export
getTissueTypes <- function(Species = c("hsapiens", "mmusculus"),
                           useBlackList = FALSE, pool = NULL) {
  # Species = "hsapiens"
  # useBlackList = TRUE


  if (is.null(pool)) {
    retryCounter <- 1
    cat("\nEstablishing connection to database ... \n")
    while(is.null(pool)) {
      pool <- try(silent = T, eval({
        pool::dbPool(
          drv = RMySQL::MySQL(),
          user = "public-rds-user", port = 3306,
          dbname="bishoplabdb",
          password='public-user-password',
          host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com"
        )
      }))
      if ("try-error" %in% class(pool)) {
        if (retryCounter == 3) {
          stop("Unable to connect to database. Check internet connection and please contanct",
               " package maintainer if you believe this is an error.")
        }
        warning(paste0("Failed to establish connection to database ... retrying now ... ",
                       (4-retryCounter), " attempts left."),
                immediate. = T)
        pool <- NULL
        retryCounter <- retryCounter + 1
      }
    }

    on.exit(function() {
      pool::poolClose(pool)
    })
  }

  Species <- Species[1]
  tabs <- pool::dbListTables(pool)
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
  tissues <- vapply(tabs, "[[", FUN.VALUE = "character", 2)
  types <- vapply(tabs, "[[", FUN.VALUE = "character", 1)
  result <- paste0(tissues, " - ", types)
  result <- result[order(result)]
  return(result)
}

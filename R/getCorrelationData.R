#' Get Gene Correlation Data
#'
#' Obtain correlation data by querying MySQL database
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#' Either "all", "normal", or "cancer". Can be a single value for all genes,
#' or a vector corresponding to geneList.
#' @param Tissue Which tissue type should gene correlations be derived from?
#' Default = "all". Can be a single value for all genes,
#' or a vector corresponding to geneList.
#' Run getTissueTypes() to see available tissue list.
#' @param geneList Vector of genes for which data will be extracted.
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied.
#' @return A correlation data frame object
#'
#' @examples
#'correlationAnalyzeR::getCorrelationData(Species = "hsapiens",
#'                                     Sample_Type = "normal",
#'                                     Tissue = "blood",
#'                                     geneList = c("ATM", "BRCA1"))
#'
#' @export
getCorrelationData <- function(Species, Sample_Type,
                               Tissue, geneList, pool = NULL) {

  # Species = "hsapiens"
  # Sample_Type = "normal"
  # Tissue = c("all", "female0asdreproductive")
  # geneList = c("A1BG", "A1BG")

  # Make sure all tissue entries are appropriate
  goodConditions <- correlationAnalyzeR::getTissueTypes(Species = Species,
                                                        pool = pool)
  goodTissues <- unique(gsub(goodConditions,
                             pattern = "(.*) - (.*)",
                             replacement = "\\1"))
  goodSamples <- unique(gsub(goodConditions,
                             pattern = "(.*) - (.*)",
                             replacement = "\\2"))

  if (! all(Sample_Type %in% goodSamples)) {
    stop("Sample type must be either 'normal' or 'cancer'")
  } else if (! all(Tissue %in% goodTissues)) {
    stop("Tissue types must be selected from available options.",
         " Run correlationAnalyzeR::getTissueTypes() to see available tissue - sample groups.")
  }


  if (! is.null(pool)) {
    if (! pool$valid) {
      pool <- NULL
    }
  }
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


  if (Species == "hsapiens") {
    geneNames <- correlationAnalyzeR::hsapiens_corrSmall_geneNames
  } else {
    geneNames <- correlationAnalyzeR::mmusculus_corrSmall_geneNames
  }
  # Queries from multiple db at once
  if (length(Tissue) == 1) {
    Tissue <- rep(Tissue, length(geneList))
  } else if (length(Tissue) > 1) {
    if (length(Tissue) != length(geneList)) {
      warning("Number of valid genes not equal
           to length of supplied Tissue vector. Using only ", Tissue[1])
      Tissue <- rep(Tissue[1], length(geneList))
    }
  }
  if (length(Sample_Type) == 1) {
    Sample_Type <- rep(Sample_Type, length(geneList))
  } else if (length(Sample_Type) > 1) {
    if (length(Sample_Type) != length(geneList)) {
      warning("Number of valid genes not equal
           to length of supplied Sample_Type vector. Using only ", Sample_Type[1])
      Sample_Type <- rep(Sample_Type[1], length(geneList))
    }
  }

  if (length(unique(Tissue)) == 1 & length(unique(Sample_Type)) == 1) {
    if (Tissue == "respiratory") {
      TissueNow <- "repiratory"
    } else {
      TissueNow <- Tissue
    }
    sql <- paste0("SELECT * FROM correlations_",
                  Species, "_",
                  tolower(unique(Sample_Type)), "_",
                  tolower(unique(TissueNow)),
                  " WHERE row_names IN ('",
                  paste(geneList, collapse = "','"), "')")
    resdf <- try(silent = T, eval({
      pool::dbGetQuery(pool, sql)
    }))
    if ("try-error" %in% class(resdf)) {
      pool::poolClose(pool)
      pool <- NULL
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
      resdf <- try(silent = T, eval({
        pool::dbGetQuery(pool, sql)
      }))
      if ("try-error" %in% class(resdf)) {
        stop("Unable to connect to the database at the moment. If you",
             " believe this is an error, please contact the package maintainer.")
      }
    }
    resdf2 <- stringr::str_split_fixed(resdf$values, stringr::fixed(","), n = Inf)
    resdf2 <- apply(t(resdf2), 1:2, as.numeric)
    resdf2 <- as.data.frame(resdf2)
    colnames(resdf2) <- resdf$row_names
    rownames(resdf2) <- geneNames
    if (length(geneList) > 1) {
      resdf2 <- resdf2[,order(match(colnames(resdf2), geneList))]
    }
  } else {
    resDfList <- list()
    for ( i in 1:length(geneList) ) {
      geneName <- geneList[i]
      TissueNow <- Tissue[i]
      if (TissueNow == "respiratory") {
        TissueNow2 <- "repiratory"
      } else {
        TissueNow2 <- TissueNow
      }
      Sample_TypeNow <- Sample_Type[i]
      sql <- paste0("SELECT * FROM correlations_",
                    Species, "_",
                    tolower(Sample_TypeNow), "_", tolower(TissueNow2),
                    " WHERE row_names IN ('",
                    geneName, "')")
      resdf <- try(silent = T, eval({
        pool::dbGetQuery(pool, sql)
      }))
      if ("try-error" %in% class(resdf)) {
        pool::poolClose(pool)
        pool <- NULL
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
        resdf <- try(silent = T, eval({
          pool::dbGetQuery(pool, sql)
        }))
        if ("try-error" %in% class(resdf)) {
          stop("Unable to connect to the database at the moment. If you",
               " believe this is an error, please contact the package maintainer.")
        }
      }
      resdf2 <- stringr::str_split_fixed(resdf$values, stringr::fixed(","), n = Inf)
      resdf2 <- apply(t(resdf2), 1:2, as.numeric)
      resdf2 <- as.data.frame(resdf2)
      colnames(resdf2) <- resdf$row_names
      rownames(resdf2) <- geneNames
      resDfList[[i]] <- resdf2
    }
    resdf2 <- dplyr::bind_cols(resDfList)
    rownames(resdf2) <- geneNames
  }
  return(resdf2)
}


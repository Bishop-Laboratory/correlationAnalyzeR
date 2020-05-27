#' Get VST values for tissues and gene of interest
#'
#' Downloads VST values for tissues of interest
#'
#' @param genesOfInterest A length-two vector with genes to compare.
#'
#' @param Tissues Which tissue type should VST be collected for? See available options
#' with getTissueTypes().
#'
#' @param Sample_Type Type of RNA Seq samples to obtain VST for? See available options
#' with getTissueTypes().
#'
#' @param useBlackList Should black-listed tissue/disease categories for this species
#' be removed from the returned list? Improves the quality of analysis by removing
#' categories with low sample numbers and high observed variance.
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied.
#' @return List of VST matrices for each selected tissue-disease combination.
#'
#' @examples
#' res <- getTissueVST(genesOfInterest = c("BRCA1", "ATM"),
#'                     Tissues = c("brain", "respiratory"),
#'                     Sample_Type = "all",
#'                     useBlackList = TRUE)
#' @export
getTissueVST <- function(genesOfInterest,
                         # Species = c("hsapiens", "mmusculus"),
                         Tissues = "all",
                         Sample_Type = c("normal", "cancer"),
                         useBlackList = TRUE, pool = NULL) {

  # genesOfInterest = c("ATM")
  # Species = "hsapiens"
  # Tissues = c("male_reproductive")
  # Sample_Type = c("all")
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
    # cat("\nEstablishing connection to database ... \n")
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

  samples <- correlationAnalyzeR::sampleVSTOrderHuman
  possibleGenes <- correlationAnalyzeR::humanGenesVST

  # Get samples for each tissue group
  possibleTissues <- correlationAnalyzeR::getTissueTypes(#Species = Species,
                                                         useBlackList = useBlackList,
                                                         pool = pool)
  possibleTissues1 <- gsub(possibleTissues, pattern = " - .*", replacement = "")
  possibleTissues2 <- gsub(possibleTissues, pattern = ".* - ", replacement = "")
  possibleRetrieval <- paste0(possibleTissues2, "_", possibleTissues1)
  if (any(Tissues == "all")) {
    ofInterest <- possibleRetrieval
  } else {
    ofInterest <- possibleRetrieval[grep(x = possibleRetrieval,
                                         pattern = paste0("_", Tissues, collapse = "|"))]
    if (! length(ofInterest)) {
      stop("No valid tissue types returned. Please check that your tissue types are ",
           "correct by running getTissueTypes()")
    }
  }
  if (all(Sample_Type != "all")) {
    ofInterest <- ofInterest[grep(x = ofInterest,
                                  pattern = paste0(Sample_Type, collapse = "|"))]
  }
  genesOfInterest <- unique(genesOfInterest)
  genesOfInterestBad <- genesOfInterest[which(! genesOfInterest %in% possibleGenes)]
  genesOfInterestFinal <- genesOfInterest[which(genesOfInterest %in% possibleGenes)]
  if (! length(genesOfInterestFinal)) {
    stop(paste0(genesOfInterestBad, collapse = ", "), " not found in VST data.",
         " view data(humanVSTGenes) or data(mouseVSTGenes) to see available gene list")
  } else if (length(genesOfInterestBad)) {
    warning(paste0(genesOfInterestBad, collapse = ", "), " not found in VST data.",
            " view data(humanVSTGenes) or data(mouseVSTGenes) to see available gene list")
  }



  # Gather VST across samples for genes of interest
  sql <- paste0("SELECT * FROM VSD_",
                Species,
                " WHERE row_names IN ('",
                paste(genesOfInterestFinal, collapse = "','"), "')")

  resdf <- try(silent = T, eval({
    DBI::dbGetQuery(conn, sql)
  }))
  if ("try-error" %in% class(resdf)) {
    stop("ERROR IN DB CONNECTION: ", resdf)
  }
  # Parse VST frame
  resdf2 <- stringr::str_split_fixed(resdf$values, stringr::fixed(","), n = Inf)
  resdf2 <- apply(t(resdf2), 1:2, as.numeric)
  resdf2 <- as.data.frame(resdf2)
  colnames(resdf2) <- resdf$row_names
  resdf2 <- cbind(samples, resdf2)
  rownames(resdf2) <- NULL
  # Return VST frame for each specified group
  resFrameList <- list()
  newList <- unlist(correlationAnalyzeR::human_grouplist, recursive = F)
  groups <- gsub(names(newList), pattern = "\\.", replacement = " - ")
  newList <- newList[! groups %in% correlationAnalyzeR::blackListHuman]
  groupNow2 <- c()
  for (i in 1:length(newList)) {
    samplesNow <- newList[[i]]
    groupNow <- gsub(names(newList)[i], pattern = "\\.", replacement = " - ")
    name2 <- gsub(names(newList)[i], pattern = "(.+)\\.(.+)", replacement = "\\2_\\1")
    groupNow2 <- c(groupNow2, gsub(name2, pattern = "\\.| ", replacement = "_"))
    frameNow <- resdf2[which(resdf2$samples %in% samplesNow),, drop = FALSE]
    resFrameList[[i]] <- frameNow
    names(resFrameList)[i] <- groupNow
  }

  keepInd <- which(groupNow2 %in% ofInterest)
  resFrameList <- resFrameList[keepInd]

  return(resFrameList)
}


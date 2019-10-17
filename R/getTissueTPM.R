#' Get TPM values for tissues and gene of interest
#'
#' Downloads TPM values for tissues of interest
#'
#' @param genesOfInterest A length-two vector with genes to compare.
#'
#' @param Species Species to obtain tissue types for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @param Tissues Which tissue type should TPM be collected for? See available options
#' with getTissueTypes().
#'
#' @param Sample_Type Type of RNA Seq samples to obtain TPM for? See available options
#' with getTissueTypes().
#'
#' @param useBlackList Should black-listed tissue/disease categories for this species
#' be removed from the returned list? Improves the quality of analysis by removing
#' categories with low sample numbers and high observed variance.
#'
#' @return List of TPM matrices for each selected tissue-disease combination.
#'
#' @examples
#' res <- getTissueTPM(genesOfInterest = c("BRCA1", "ATM"),
#'                     Species = "hsapiens",
#'                     Tissues = c("brain", "respiratory"),
#'                     Sample_Type = "all",
#'                     useBlackList = TRUE)
#' @export
getTissueTPM <- function(genesOfInterest,
                         Species = c("hsapiens", "mmusculus"),
                         Tissues = "all",
                         Sample_Type = c("all", "normal", "cancer"),
                         useBlackList = TRUE) {

  # genesOfInterest = c("BRCA1", "ATM")
  # Species = "hsapiens"
  # Tissues = "all"
  # Sample_Type = "all"
  # useBlackList = TRUE

  if (Species[1] == "hsapiens") {
    samples <- correlationAnalyzeR::sampleTPMOrderHuman
    possibleGenes <- correlationAnalyzeR::humanGenesTPM
  } else {
    samples <- correlationAnalyzeR::sampleTPMOrderMouse
    possibleGenes <- correlationAnalyzeR::mouseGenesTPM
  }
  # Get samples for each tissue group
  possibleTissues <- correlationAnalyzeR::getTissueTypes(Species = Species,
                                                         useBlackList = useBlackList)
  possibleTissues1 <- gsub(possibleTissues, pattern = " - .*", replacement = "")
  possibleTissues2 <- gsub(possibleTissues, pattern = ".* - ", replacement = "")
  possibleRetrieval <- paste0(possibleTissues2, "_", possibleTissues1)
  if (Tissues == "all") {
    ofInterest <- possibleRetrieval
  } else {
    ofInterest <- possibleRetrieval[grep(x = possibleRetrieval,
                                         pattern = paste0(Tissues, collapse = "|"))]
    if (! length(ofInterest)) {
      stop("No valid tissue types returned. Please check that your tissue types are ",
           "correct by running getTissueTypes()")
    }
  }
  ofInterest <- gsub(x = ofInterest, pattern = "respiratory", replacement = "repiratory")
  if (Sample_Type != "all") {
    ofInterest <- ofInterest[grep(x = ofInterest,
                                  pattern = paste0(Sample_Type, collapse = "|"))]
  }
  genesOfInterest <- unique(genesOfInterest)
  genesOfInterestBad <- genesOfInterest[which(! genesOfInterest %in% possibleGenes)]
  genesOfInterestFinal <- genesOfInterest[which(genesOfInterest %in% possibleGenes)]
  if (! length(genesOfInterestFinal)) {
    stop(paste0(genesOfInterestBad, collapse = ", "), " not found in TPM data.",
         " view data(humanTPMGenes) or data(mouseTPMGenes) to see available gene list")
  } else if (length(genesOfInterestBad)) {
    warning(paste0(genesOfInterestBad, collapse = ", "), " not found in TPM data.",
            " view data(humanTPMGenes) or data(mouseTPMGenes) to see available gene list")
  }

  con <- DBI::dbConnect(RMySQL::MySQL(), user = "public-rds-user", port = 3306,
                        dbname="bishoplabdb",
                        password='public-user-password',
                        host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
  sql <- paste0("SELECT * FROM ",
                Species, "_sample_group_key ",
                " WHERE row_names IN ('",
                paste(ofInterest, collapse = "','"), "')")
  res <- DBI::dbSendQuery(conn = con, statement = sql)
  resdf <- DBI::dbFetch(res, n=-1)
  DBI::dbClearResult(res)
  resdf2 <- data.frame(row.names = resdf$row_names, samples = resdf$samples)
  resdfList <- stats::setNames(split(resdf2, seq(nrow(resdf2))), rownames(resdf2))
  newList <- lapply(resdfList, FUN = function(x){
    newX <- as.character(x[,1])
    newX <- unlist(strsplit(newX, split = ","))
  })


  # Gather TPM across samples for genes of interest
  sql <- paste0("SELECT * FROM TPM_",
                Species,
                " WHERE row_names IN ('",
                paste(genesOfInterestFinal, collapse = "','"), "')")
  res <- DBI::dbSendQuery(conn = con, statement = sql)
  resdf <- DBI::dbFetch(res, n=-1)
  DBI::dbClearResult(res)
  DBI::dbDisconnect(con)

  # Parse TPM frame
  resdf2 <- stringr::str_split_fixed(resdf$values, stringr::fixed(","), n = Inf)
  resdf2 <- apply(t(resdf2), 1:2, as.numeric)
  resdf2 <- as.data.frame(resdf2)
  colnames(resdf2) <- resdf$row_names
  resdf2$samples <- samples
  rownames(resdf2) <- NULL
  # Return TPM frame for each specified group
  resFrameList <- list()
  for (i in 1:length(newList)) {
    samplesNow <- newList[[i]]
    groupNow <- names(newList)[i]
    frameNow <- resdf2[which(resdf2$samples %in% samplesNow),, drop = FALSE]
    resFrameList[[i]] <- frameNow
    names(resFrameList)[i] <- groupNow
  }
  names(resFrameList) <- gsub(x = names(resFrameList),
                              pattern = "repiratory", replacement = "respiratory")

  return(resFrameList)
}


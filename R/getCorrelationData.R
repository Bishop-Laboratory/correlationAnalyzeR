#' Get Gene Correlation Data
#'
#' Queries from MySQL the requested gene correlation data
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#' Either "all", "normal", or "cancer". Can be a single value for all genes,
#' or a vector where each entry corresponds to a gene of interest.
#' @param Tissue Which tissue type should gene correlations be derived from?
#' Default = "all". Can be a single value for all genes, or a vector where each
#' entry corresponds to a gene of interest.
#' Run getTissueTypes() to see available tissue list.
#' @param geneList Vector of genes for which data will be extracted.
#'
#' @return A correlation data frame object
#'
#' @examples
#'
#' correlationDF <- getCorrelationData(Species = "hsapiens",
#'                                     Sample_Type = "normal",
#'                                     Tissue = "blood",
#'                                     geneList = c("ATM", "BRCA1"))
#'
#' @export
getCorrelationData <- function(Species, Sample_Type, Tissue, geneList) {
  # # Bug Testing
  # Species <- "hsapiens"
  # Sample_Type <- "normal"
  # tissue <- "brain"
  # geneList <- c("BRCA1", "NFE2L2", "ATM", "PARP1")
  # Tissue <- c("brain", "thyroid", "brain", "blood")
  # Sample_Type <- c("normal", "cancer", "cancer", "normal")
  # dbListTables(con)
  # # Got this part from the Server
  # load("data/hsapiens_geneNames.RData")
  # hsapiens_corrSmall_geneNames <- geneNames
  # usethis::use_data(hsapiens_corrSmall_geneNames, overwrite = T)
  # load("data/mmusculus_geneNames.RData")
  # mmusculus_corrSmall_geneNames <- geneNames
  # usethis::use_data(mmusculus_corrSmall_geneNames, overwrite = T)

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
      stop("Number of valid genes not equal
           to length of supplied Tissue vector.
           Tissue vector should include 1 entry or the number of
           entries equal to the number of genesOfInterest.")
    }
  }
  if (length(Sample_Type) == 1) {
    Sample_Type <- rep(Sample_Type, length(geneList))
  } else if (length(Sample_Type) > 1) {
    if (length(Sample_Type) != length(geneList)) {
      stop("Number of valid genes not equal
           to length of supplied Sample_Type vector.
           Sample_Type vector should include 1 entry or the number of
           entries equal to the number of genesOfInterest.")
    }
  }

  con <- DBI::dbConnect(RMySQL::MySQL(), user = "public-rds-user", port = 3306,
                        dbname="bishoplabdb",
                        password='public-user-password',
                        host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
  if (length(unique(Tissue)) == 1 & length(unique(Sample_Type)) == 1) {
    sql <- paste0("SELECT * FROM correlations_",
                  Species, "_",
                  tolower(unique(Sample_Type)), "_",
                  tolower(unique(Tissue)),
                  " WHERE row_names IN ('",
                  paste(geneList, collapse = "','"), "')")
    res <- DBI::dbSendQuery(conn = con, statement = sql)
    resdf <- DBI::dbFetch(res, n=-1)
    DBI::dbClearResult(res)
    resdf2 <- sapply(resdf$values, strsplit, ",")
    names(resdf2) <- resdf$row_names
    resdf2 <- lapply(resdf2, as.numeric)
    resdf2 <- as.data.frame(resdf2)
    rownames(resdf2) <- geneNames
  } else {
    resDfList <- list()
    for ( i in 1:length(geneList) ) {
      geneName <- geneList[i]
      TissueNow <- Tissue[i]
      Sample_TypeNow <- Sample_Type[i]
      sql <- paste0("SELECT * FROM correlations_",
                    Species, "_",
                    tolower(Sample_TypeNow), "_", tolower(TissueNow),
                    " WHERE row_names IN ('",
                    geneName, "')")
      res <- DBI::dbSendQuery(conn = con, statement = sql)
      resdf <- DBI::dbFetch(res, n=-1)
      DBI::dbClearResult(res)
      resdf2 <- sapply(resdf$values, strsplit, ",")
      names(resdf2) <- resdf$row_names
      resdf2 <- lapply(resdf2, as.numeric)
      resdf2 <- as.data.frame(resdf2)
      rownames(resdf2) <- geneNames
      resDfList[[i]] <- resdf2
    }
    resdf2 <- do.call(cbind, resDfList)
  }
  DBI::dbDisconnect(conn = con)
  return(resdf2)
}




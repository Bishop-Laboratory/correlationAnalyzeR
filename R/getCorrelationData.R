#' Get Gene Correlation Data
#'
#' Queries from MySQL the requested gene correlation data
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#'     Either "All", "Normal_Tissues", or "Tumor_Tissues".
#' @param geneList Vector of genes for which data will be extracted.
#'
#' @return A correlation data frame object
#'
#' @examples
#' correlationDF <- getCorrelationData(Species = "hsapiens", Sample_Type = "Normal_Tissues",
#'                                     geneList = c("ATM", "BRCA1"))
#'
#' @export
getCorrelationData <- function(Species, Sample_Type, geneList) {
  # # Bug Testing
  # Species <- "hsapiens"
  # Sample_Type <- "Normal_Tissues"
  # geneList <- c("BRCA1", "EZH2", "PARP1", "TOP2A", "NFE2L2")
  # geneList <- intGenes

  # # Got this part from the Server
  # load("data/hsapiens_corrSmall_geneNames.RData")
  # hsapiens_corrSmall_geneNames <- geneNames
  # usethis::use_data(hsapiens_corrSmall_geneNames)

  if (Species == "hsapiens") {
    geneNames <- hsapiens_corrSmall_geneNames
  } else {
    geneNames <- mmusculus_corrSmall_geneNames
  }


  con <- DBI::dbConnect(RMySQL::MySQL(), user = "public-rds-user", port = 3306, dbname="bishoplabdb",
                        password='public-user-password', host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
  sql <- paste0("SELECT * FROM ",
                Species, "_",
                tolower(Sample_Type),
                "_corrsmall WHERE row_names IN ('",
                paste(geneList, collapse = "','"), "')")
  res <- DBI::dbSendQuery(conn = con, statement = sql)
  resdf <- DBI::dbFetch(res, n=-1)
  DBI::dbClearResult(res)
  DBI::dbDisconnect(conn = con)
  resdf2 <- sapply(resdf$values, strsplit, ",")
  names(resdf2) <- resdf$row_names
  resdf2 <- lapply(resdf2, as.numeric)
  resdf2 <- as.data.frame(resdf2)
  rownames(resdf2) <- geneNames
  return(resdf2)
}




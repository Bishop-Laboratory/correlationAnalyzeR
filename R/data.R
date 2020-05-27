# # Get datasets
# load("data/MSIGDB_Geneset_Names.rda")
# usethis::use_data(MSIGDB_Geneset_Names, internal = F, overwrite = T)
#
# load("data/blackListHuman.rda")
# usethis::use_data(blackListHuman, internal = F, overwrite = T)
#
# load("data/blackListMouse.rda")
# usethis::use_data(blackListMouse, internal = F, overwrite = T)
#
#
# load("data/hsapiens_corrSmall_geneNames.rda")
# usethis::use_data(hsapiens_corrSmall_geneNames, internal = F, overwrite = T)
#
# load("data/mmusculus_corrSmall_geneNames.rda")
# usethis::use_data(mmusculus_corrSmall_geneNames, internal = F, overwrite = T)
#
# load("data/humanGenesVST.rda")
# usethis::use_data(humanGenesVST, internal = F, overwrite = T)
#
# load("data/mouseGenesVST.rda")
# usethis::use_data(mouseGenesVST, internal = F, overwrite = T)
#
# load("data/sampleVSTOrderHuman.rda")
# usethis::use_data(sampleVSTOrderHuman, internal = F, overwrite = T)
#
# load("data/sampleVSTOrderMouse.rda")
# usethis::use_data(sampleVSTOrderMouse, internal = F, overwrite = T)
#
# load("data/MSIGDB_Geneset_Small_Names.rda")
# usethis::use_data(MSIGDB_Geneset_Small_Names, internal = F, overwrite = T)

# load("misc/groupList_human_raw.RData")
# human_coldata <- colData_human_raw
# usethis::use_data(human_coldata)
# human_grouplist <- groupList_human_raw
# usethis::use_data(human_grouplist)

# tabDF <- as.data.frame(table( paste0(human_coldata$Tissue, " - ", human_coldata$disease)), stringsAsFactors = F)
# keepSamps <- unique(tabDF$Var1[tabDF$Freq > 30])
# keepSamps <- tolower(keepSamps)
# newList <- unlist(correlationAnalyzeR::human_grouplist, recursive = F)
# groupNow <- gsub(names(newList), pattern = "\\.", replacement = " - ")
# blackListHuman <- groupNow[! groupNow %in% keepSamps]
# blackListHuman <- blackListHuman[grep(blackListHuman, pattern = " all", invert = T)]
# usethis::use_data(blackListHuman, overwrite = T)


#' A data frame of sample info from GEO (human)
#'@docType data
#'@keywords data
"human_coldata"

#' A list of sample categorizations (human)
#'@docType data
#'@keywords data
"human_grouplist"

#' A vector of valid MSIGDB geneset names
#'@source msigdbr()
#'@docType data
#'@keywords data
"MSIGDB_Geneset_Names"

#' A vector of valid MSIGDB geneset names with fewer then 500 genes associated with them
#'@source msigdbr()
#'@docType data
#'@keywords data
"MSIGDB_Geneset_Small_Names"

#' A vector of blacklisted tissue-disease categories for human samples
#' @docType data
#' @keywords data
"blackListHuman"

#' A vector of valid human genes to extract correlations from
#' @docType data
#' @keywords data
"hsapiens_corrSmall_geneNames"


#' A vector of valid human genes with VST data available
#' @docType data
#' @keywords data
"humanGenesVST"


#' A vector containing the order of human samples in the sample-tissue SQL table
#' @docType data
#' @keywords data
"sampleVSTOrderHuman"


#' A vector containing the names of valid TERM2GENE categories for GSEA_Type or pathwayType input.
#' @docType data
#' @keywords data
"pathwayCategories"

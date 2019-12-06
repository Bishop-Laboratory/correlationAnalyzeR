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
# load("data/humanGenesTPM.rda")
# usethis::use_data(humanGenesTPM, internal = F, overwrite = T)
#
# load("data/mouseGenesTPM.rda")
# usethis::use_data(mouseGenesTPM, internal = F, overwrite = T)
#
# load("data/sampleTPMOrderHuman.rda")
# usethis::use_data(sampleTPMOrderHuman, internal = F, overwrite = T)
#
# load("data/sampleTPMOrderMouse.rda")
# usethis::use_data(sampleTPMOrderMouse, internal = F, overwrite = T)

#' A vector of valid MSIGDB geneset names
#'@source msigdbr()
#'@docType data
#'@keywords data
"MSIGDB_Geneset_Names"

#' A vector of blacklisted tissue-disease categories for human samples
#' @docType data
#' @keywords data
"blackListHuman"

#' A vector of blacklisted tissue-disease categories for mouse samples
#' @docType data
#' @keywords data
"blackListMouse"

#' A vector of valid human genes to extract correlations from
#' @docType data
#' @keywords data
"hsapiens_corrSmall_geneNames"

#' A vector of valid mouse genes to extract correlations from
#' @docType data
#' @keywords data
"mmusculus_corrSmall_geneNames"

#' A vector of valid human genes with TPM data available
#' @docType data
#' @keywords data
"humanGenesTPM"

#' A vector of valid mouse genes with TPM data available
#' @docType data
#' @keywords data
"mouseGenesTPM"

#' A vector containing the order of human samples in the sample-tissue SQL table
#' @docType data
#' @keywords data
"sampleTPMOrderHuman"

#' A vector containing the order of mouse samples in the sample-tissue SQL table
#' @docType data
#' @keywords data
"sampleTPMOrderMouse"

#' A vector containing the names of valid TERM2GENE categories for GSEA_Type or pathwayType input.
#' @docType data
#' @keywords data
"pathwayCategories"

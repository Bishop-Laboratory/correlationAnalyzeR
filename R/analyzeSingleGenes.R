#' Analyze Single Genes
#'
#' Obtains correlation info and GSEA results for each gene of interest
#'
#' @param genesOfInterest A vector of genes to analyze.
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @param GSEA_Type Whether GSEA should consider all msigdb annotations,
#'     or just those in the most popular categories. Should be one of either
#'     'simple' or 'complex'.
#'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#'     Either "All", "Normal_Tissues", or "Tumor_Tissues".
#'
#' @param outputPrefix Prefix for saved files. Should include directory info.
#'
#' @param runGSEA If TRUE will run GSEA using gene correlation values.
#'
#' @param returnDataOnly if TRUE will return only a dataframe of correlation
#'    values and will not generate any folders or files.
#'
#' @return A dataframe of correlation values for each gene of interest.
#'
#' @examples
#' genesOfInterest <- c("ATM", "SLC7A11")
#' Result <- analyzeSingleGenes(genesOfInterest = genesOfInterest,
#'                               Species = "hsapiens",
#'                               GSEA_Type = "simple",
#'                               Sample_Type = "Normal_Tissues")
#'
#' @export
analyzeSingleGenes <- function(genesOfInterest,
                               Species = c("hsapiens", "mmusculus"),
                               GSEA_Type = c("simple", "complex"),
                               Sample_Type = c(#"All",
                                          "Normal_Tissues",
                                          "Tumor_Tissues"),
                               outputPrefix = "CorrelationAnalyzeR_Output",
                               runGSEA = T, returnDataOnly = F) {

  # # Debug/Test
  # genesOfInterest <- c("ATM", "SLC3A2", "SON", "BRCA1")
  # GSEA_Type = "simple"
  # runGSEA = F
  # outputPrefix = "tests/CorrelationAnalyzeR_Output"
  # Species = "hsapiens"
  # Sample_Type = "Normal_Tissues"
  # returnDataOnly = F

  # Parse arguments
  genesOfInterest <- unique(genesOfInterest)
  geneString <- paste(genesOfInterest, collapse = ", ")
  # Fix possible parameter conflict
  if (returnDataOnly) {
    runGSEA <- F
  }

  # Create output folder
  if (! dir.exists(outputPrefix) & ! returnDataOnly) {
    dir.create(outputPrefix)
  }
  # Load appropriate TERM2GENE file built from msigdbr()
  if (Species[1] %in% c("hsapiens", "mmusculus")) {
    if (runGSEA) {
      if (! GSEA_Type %in% c("simple", "complex")) {
        stop("\nPlease enter either 'simple' or 'complex' for GSEA_Type\n")
      } else if (GSEA_Type[1] == "simple") {
        if (Species[1] == "hsapiens") {
          data("hsapiens_simple_TERM2GENE")
          TERM2GENE <- hsapiens_simple_TERM2GENE
        } else {
          data("mmusculus_simple_TERM2GENE")
          TERM2GENE <- mmusculus_simple_TERM2GENE
        }
      } else {
        if (Species[1] == "hsapiens") {
          data("hsapiens_complex_TERM2GENE")
          TERM2GENE <- hsapiens_complex_TERM2GENE
        } else {
          data("mmusculus_complex_TERM2GENE")
          TERM2GENE <- mmusculus_complex_TERM2GENE
        }
      }
    }
  } else {
    stop("\ncorrelationAnalyzeR currently supports only Human and Mouse data.
         Please select either 'hsapiens' or 'mmusculus' for Species parameter.
         \n")
  }
  # Check genes to make sure they exist
  avGenes <- getAvailableGenes(Species = Species)
  avGenes <- as.character(avGenes$geneName)
  badGenes <- genesOfInterest[which(! genesOfInterest %in% avGenes)]
  if (length(badGenes) > 0) {
    stop(paste0("\n\t\t\t'", paste(badGenes, collapse = ", "), "' not found
                in correlation data.
                Please check available data with getAvailableGenes().
                Your gene(s) of interest may have an updated name or
                have a species-specific identifier.\n"))
  }
  cat("\nRetrieving any missing correlation data...\n")
  # Call downloadData to get all required files
  corrDF <- getCorrelationData(Species = Species[1],
                               Sample_Type = Sample_Type[1],
                               geneList = genesOfInterest)
  # Main code
  for (i in 1:length(colnames(corrDF))) {
    gene <- colnames(corrDF)[i]
    cat(paste0("\n", gene))
    # Create output folder for gene
    geneOutDir <- file.path(outputPrefix, gene)
    if (! dir.exists(geneOutDir) & ! returnDataOnly) {
      dir.create(geneOutDir)
    }
    if (! returnDataOnly) {
      # Remove gene from correlation values to build normal distribution
      vec <- corrDF[,i]
      names(vec) <- rownames(corrDF)
      vec <- vec[order(vec, decreasing = T)]
      vec <- vec[c(-1)]
      # Make a histogram of gene correlations
      png(file.path(geneOutDir, paste0(gene, ".png")))
      hist(vec,
           main = paste0(gene, " gene correlations"),
           breaks = 100, xlab = "Correlation Value")
      dev.off()
      if (runGSEA) {
        # Perform GSEA
        cat("\nGSEA\n")
        resGSEA <- suppressWarnings(myGSEA(ranks = vec, TERM2GENE = TERM2GENE,
                                           plotFile = paste0(gene, ".corrPathways"),
                                           outDir = geneOutDir,
                                           Condition = paste0(gene, " Correlated Genes")))
      }
    }
  }
  if (! returnDataOnly) {
    # Save results from all queried genes
    if (length(colnames(corrDF)) < 25 & length(colnames(corrDF)) > 1) {
      corrDF <- cbind(rownames(corrDF), corrDF)
      colnames(corrDF)[1] <- "geneName"
      write.csv(corrDF, file = file.path(outputPrefix, "geneCorrelationData.csv"), row.names = F)
      save(corrDF, file = file.path(outputPrefix, "geneCorrelationData.RData"))
    } else if (length(colnames(corrDF)) > 1) {
      save(corrDF, file = file.path(outputPrefix, "geneCorrelationData.RData"))
    }
  }
  # Return results
  return(corrDF)
}

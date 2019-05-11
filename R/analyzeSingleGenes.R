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
                               Sample_Type = c("All",
                                          "Normal_Tissues",
                                          "Tumor_Tissues"),
                               outputPrefix = "CorrelationAnalyzeR_Output",
                               runGSEA = T) {

  # # Debug/Test
  # genesOfInterest <- c("ATM", "SLC3A2", "SON", "BRCA1")
  # GSEA_Type = "simple"
  # runGSEA = T
  # outputPrefix = "tests/CorrelationAnalyzeR_Output"
  # Species = "hsapiens"
  # Sample_Type = "Normal_Tissues"

  # Parse arguments
  genesOfInterest <- unique(genesOfInterest)
  geneString <- paste(genesOfInterest, collapse = ", ")
  # Create output folder
  if (! dir.exists(outputPrefix)) {
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
  downloadData(Species = Species[1],
               Sample_Type = Sample_Type[1],
               geneList = genesOfInterest)
  downloadFolder <- system.file("data", package = "correlationAnalyzeR")
  dataDir <- file.path(downloadFolder, "Correlation_Data",
                              Species[1], Sample_Type[1])
  # Main code
  for (i in 1:length(genesOfInterest)) {
    gene <- genesOfInterest[i]
    cat(paste0("\n", gene))
    # Get the data files for genesOfInterest
    geneFile <- paste0(gene, ".RData")
    file <- file.path(dataDir, geneFile)
    load(file)
    # Create output folder for gene
    geneOutDir <- file.path(outputPrefix, gene)
    if (! dir.exists(geneOutDir)) {
      dir.create(geneOutDir)
    }
    # Initialize results frame
    if (i == 1) {
      resultsFrame <- data.frame(geneName = names(vec))
    }
    # Build dataframe object from correlation data
    corrdf <- as.data.frame(vec)
    colnames(corrdf)[1] <- gene
    corrdf$geneName <- rownames(corrdf)
    # Add correlation data to final dataframe
    resultsFrame <- merge(x = resultsFrame, y = corrdf, by = "geneName")
    corrdf <- corrdf[,c(2, 1)]
    write.csv(corrdf, file = file.path(geneOutDir, paste0(gene, ".csv")),
              row.names = F)
    # Remove gene from correlation values to build normal distribution
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
  # Save results from all queried genes
  if (length(genesOfInterest) < 25 & length(genesOfInterest) > 1) {
    write.csv(resultsFrame, file = file.path(outputPrefix, "geneCorrelationData.csv"), row.names = F)
    save(resultsFrame, file = file.path(outputPrefix, "geneCorrelationData.RData"))
  } else if (length(genesOfInterest) > 1) {
    save(resultsFrame, file = file.path(outputPrefix, "geneCorrelationData.RData"))
  }
  # Return results
  return(resultsFrame)
}

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
#' Either "all", "normal", or "cancer". Can be a single value for all genes,
#' or a vector where each entry corresponds to a gene of interest.
#'
#' @param Tissue Which tissue type should gene correlations be derived from?
#' Default = "all". Can be a single value for all genes, or a vector where each
#' entry corresponds to a gene of interest.
#'  Run getTissueTypes() to see available tissue list.
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
                                          "normal",
                                          "cancer"),
                               Tissue = "all",
                               outputPrefix = "CorrelationAnalyzeR_Output",
                               runGSEA = T, topPlots = T, returnDataOnly = F) {

  # # Debug/Test
  # genesOfInterest <- c("BRCA1", "ATM", "PARP1")
  # GSEA_Type = "simple"
  # runGSEA = T
  # topPlots = F
  # outputPrefix = "tests/CorrelationAnalyzeR_Output"
  # Species = "hsapiens"
  # Sample_Type = c("cancer", "normal", "cancer")
  # Tissue = "brain"
  # returnDataOnly = F

  # Parse arguments
  geneString <- paste(genesOfInterest, collapse = ", ")

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
          TERM2GENE <- correlationAnalyzeR::hsapiens_simple_TERM2GENE
        } else {
          data("mmusculus_simple_TERM2GENE")
          TERM2GENE <- correlationAnalyzeR::mmusculus_simple_TERM2GENE
        }
      } else {
        if (Species[1] == "hsapiens") {
          data("hsapiens_complex_TERM2GENE")
          TERM2GENE <- correlationAnalyzeR::hsapiens_complex_TERM2GENE
        } else {
          data("mmusculus_complex_TERM2GENE")
          TERM2GENE <- correlationAnalyzeR::mmusculus_complex_TERM2GENE
        }
      }
    }
  } else {
    stop("\ncorrelationAnalyzeR currently supports only Human and Mouse data.
         Please select either 'hsapiens' or 'mmusculus' for Species parameter.
         \n")
  }
  # Check genes to make sure they exist
  avGenes <- correlationAnalyzeR::getAvailableGenes(Species = Species)
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

  if (length(Tissue) == 1) {
    Tissue <- rep(Tissue, length(genesOfInterest))
  } else if (length(Tissue) > 1) {
    if (length(Tissue) != length(genesOfInterest)) {
      stop("Number of valid genes not equal
           to length of supplied Tissue vector.
           Tissue vector should include 1 entry or the number of
           entries equal to the number of genesOfInterest.")
    }
  }
  if (length(Sample_Type) == 1) {
    Sample_Type <- rep(Sample_Type, length(genesOfInterest))
  } else if (length(Sample_Type) > 1) {
    if (length(Sample_Type) != length(genesOfInterest)) {
      stop("Number of valid genes not equal
           to length of supplied Sample_Type vector.
           Sample_Type vector should include 1 entry or the number of
           entries equal to the number of genesOfInterest.")
    }
  }
  # Call downloadData to get all required files
  corrDF <- correlationAnalyzeR::getCorrelationData(Species = Species,
                                                    Tissue = Tissue,
                                                    Sample_Type = Sample_Type,
                                                    geneList = genesOfInterest)
  resList <- list()
  # Main code
  for (i in 1:length(colnames(corrDF))) {
    gene <- colnames(corrDF)[i]
    cat(paste0("\n", gene))
    # Create output folder for gene
    geneOutDir <- file.path(outputPrefix, gene)
    if (! dir.exists(geneOutDir) & ! returnDataOnly) {
      dir.create(geneOutDir)
    }
    # Remove gene from correlation values to build normal distribution
    vec <- corrDF[,i]
    names(vec) <- rownames(corrDF)
    vec <- vec[order(vec, decreasing = T)]
    vec <- vec[c(-1)]
    # Make a histogram of gene correlations
    corrDF2 <- corrDF[which(rownames(corrDF) != gene),i, drop = F]
    geneOne <- genesOfInterest[i]
    tissueOne <- Tissue[i]
    sampleOne <- Sample_Type[i]
    geneOneTitle <- paste0(geneOne, ", ",
                           tools::toTitleCase(tissueOne),
                           " - ",
                           tools::toTitleCase(sampleOne))

    p <- ggpubr::gghistogram(data = corrDF2, x = gene, y = "..count..",
                             bins = 100, ylab = "Frequency\n",
                             title = gene,
                             caption = paste0(tools::toTitleCase(tissueOne),
                                              " - ",
                                              tools::toTitleCase(sampleOne)),
                             xlab = paste0(gene, " correlation values"))

    resList[[i]] <- list()
    names(resList)[i] <- geneOneTitle
    geneOneTitleFile <- gsub(pattern = ",| |-", replacement = "", x = geneOneTitle)
    resList[[geneOneTitle]][["corrHist"]] <- p
    if (! returnDataOnly) {
      ggplot2::ggsave(filename = file.path(geneOutDir,
                                           paste0(geneOneTitleFile, ".png")),
                      plot = p)
    }


    if (runGSEA) {
      # Perform GSEA
      cat("\nGSEA\n")

      resGSEA <- suppressWarnings(
        correlationAnalyzeR::myGSEA(ranks = vec, TERM2GENE = TERM2GENE,
                                    plotFile = paste0(geneOneTitleFile,
                                                      ".corrPathways"),
                                    outDir = geneOutDir,
                                    topPlots = topPlots,
                                    returnDataOnly = returnDataOnly,
                                    Condition = paste0(geneOneTitle,
                                                       ": Correlated Genes"))
      )
      resList[[geneOneTitle]][["GSEA"]] <- resGSEA
    }

  }
  resList[["correlations"]] <- corrDF
  # Return results
  return(resList)
}

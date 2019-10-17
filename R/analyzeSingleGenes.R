#' Analyze Single Genes
#'
#' Obtains correlations and corGSEA results for each gene of interest
#'
#' @param genesOfInterest A vector of genes to analyze.
#'
#' @param Species Species to obtain gene names for.
#' Either 'hsapiens' or 'mmusculus'
#'
#' @param GSEA_Type Whether GSEA should consider all msigdb annotations,
#' or just those in the most popular categories. Should be one of either
#' 'simple' or 'complex'.
#'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#' Either "all", "normal", or "cancer". Can be a single value for all genes,
#' or a vector corresponding to genesOfInterest.
#'
#' @param Tissue Which tissue type should gene correlations be derived from?
#' Default = "all". Can be a single value for all genes,
#' or a vector corresponding to genesOfInterest.
#' Run getTissueTypes() to see available tissues.
#'
#' @param crossCompareMode Instead of normal single gene analysis, will generate
#' correlations for single gene across all tissue-disease groups. GSEA will not be run.
#' Will only consider user input for returnDataOnly, whichCompareGroups,
#' outputPrefix, Species, and genesOfInterest.
#'
#' @param whichCompareGroups For crossCompareMode, select "all", "normal", or "cancer"
#' to analyze correlations from the corresponding groups.
#'
#' @param outputPrefix Prefix for saved files. Should include directory info.
#'
#' @param runGSEA If TRUE will run GSEA using gene correlation values.
#'
#' @param returnDataOnly if TRUE will return only a dataframe of correlation
#'    values and will not generate any folders or files.
#'
#' @param topPlots Logical. If TRUE, myGSEA() will build gsea plots for top correlated genesets.
#'
#' @return A named list of correlation values, corGSEA results,
#' and visualizations for each gene of interest.
#'
#' @examples
#' genesOfInterest <- c("ATM", "SLC7A11")
#' correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = genesOfInterest,
#'                               Species = "hsapiens", returnDataOnly = TRUE,
#'                               GSEA_Type = "simple",
#'                               Sample_Type = c("normal", "cancer"),
#'                               Tissue = c("respiratory", "pancreas"))
#'
#' genesOfInterest <- c("Brca1")
#' correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = genesOfInterest,
#'                               Species = "mmusculus",
#'                               GSEA_Type = "simple", returnDataOnly = TRUE,
#'                               crossCompareMode = TRUE,
#'                               whichCompareGroups = "normal")
#'
#' @export
analyzeSingleGenes <- function(genesOfInterest,
                               Species = c("hsapiens", "mmusculus"),
                               GSEA_Type = c("simple", "complex"),
                               Sample_Type = "normal",
                               Tissue = "all",
                               crossCompareMode = FALSE,
                               whichCompareGroups = c("all", "normal", "cancer"),
                               outputPrefix = "CorrelationAnalyzeR_Output",
                               runGSEA = TRUE, topPlots = TRUE, returnDataOnly = FALSE) {

  # # Bug testing
  # genesOfInterest <- c("FLI1", "EWSR1", "STAG2")
  # Species = c("hsapiens", "mmusculus")
  # GSEA_Type = c("simple", "complex")
  # Sample_Type = "normal"
  # Tissue = "all"
  # crossCompareMode = TRUE
  # whichCompareGroups = c("all", "normal", "cancer")
  # outputPrefix = "CorrelationAnalyzeR_Output"
  # runGSEA = FALSE
  # topPlots = FALSE
  # returnDataOnly = TRUE


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
          TERM2GENE <- correlationAnalyzeR::hsapiens_simple_TERM2GENE
        } else {
          TERM2GENE <- correlationAnalyzeR::mmusculus_simple_TERM2GENE
        }
      } else {
        if (Species[1] == "hsapiens") {
          TERM2GENE <- correlationAnalyzeR::hsapiens_complex_TERM2GENE
        } else {
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

  if (crossCompareMode) {
    cat("\nRunning cross comparison mode ... \n")
    if (Species == "mmusculus") {
      whichCompareGroups <- "normal"
      cat("\nOnly normal tissue comparisons available for mouse",
          " due to black-listing of cancer groups.\n")
      cat("\nContinuing with normal tissues ... \n")
    }
    availTissue <- correlationAnalyzeR::getTissueTypes(Species = Species,
                                                       useBlackList = TRUE)
    runGSEA <- F

    if(whichCompareGroups != "all") {
      availTissue <- availTissue[grep(availTissue,
                                      pattern = whichCompareGroups)]
    }
    availTissue <- strsplit(availTissue, split = " - ")
    Tissue <- sapply(availTissue, "[[", 1)
    genesVec <- rep(genesOfInterest, each = length(Tissue))
    Tissue <- rep(Tissue, length(genesOfInterest))
    Sample_Type <- sapply(availTissue, "[[", 2)
    Sample_Type <- rep(Sample_Type, length(genesOfInterest))
    genesOfInterest <- genesVec

  }

  # Call downloadData to get all required files
  corrDF <- correlationAnalyzeR::getCorrelationData(Species = Species,
                                                    Tissue = Tissue,
                                                    Sample_Type = Sample_Type,
                                                    geneList = genesOfInterest)
  if(crossCompareMode) {
    tissue2 <- gsub(Tissue, pattern = "0", replacement = " ")
    tissue2 <- stringr::str_to_title(tissue2)
    if(whichCompareGroups != "all") {
      namesVec <- tissue2
    } else {
      namesVec <- paste0(tissue2, " - ", stringr::str_to_title(Sample_Type))
    }
    topName <- paste0(genesOfInterest, "_", Tissue, "_", Sample_Type)
    resList <- list()
    geneList <- unique(genesOfInterest)
    for (i in 1:length(geneList)) {
      geneNow <- geneList[i]
      resList[[i]] <- list()
      names(resList)[i] <- geneNow
      inds <- which(genesOfInterest == geneNow)
      newDF <- corrDF[,inds]
      topNameNow <- topName[inds]
      colnames(newDF) <- topNameNow
      namesVecNow <- namesVec[inds]
      resList[[i]][["correlations"]] <- newDF
      newDFNorm <- preprocessCore::normalize.quantiles(as.matrix(newDF))
      newDFNorm <- as.data.frame(newDFNorm)
      rownames(newDFNorm) <- rownames(newDF)
      colnames(newDFNorm) <- colnames(newDF)
      newDFNorm$Variance <- matrixStats::rowVars(as.matrix(newDFNorm))
      newDFNorm <- newDFNorm[order(newDFNorm$Variance, decreasing = TRUE),]
      n <- length(colnames(newDFNorm))
      pMat <- newDFNorm[c(1:30),c(-n)]
      pMatBig <- newDFNorm[c(1:500),c(-n)]
      titleName <- ifelse(whichCompareGroups == "all", paste0(geneNow,
                                                              " correlations across conditions"),
                          no = ifelse(whichCompareGroups == "normal",
                                      paste0(geneNow,
                                             " correlations across tissue types"),
                                      no = paste0(geneNow,
                                                  " correlations across tumor tissues")))

      phSmall <- pheatmap::pheatmap(pMat, silent = TRUE,
                                    angle_col = 45, #main = titleName,
                                    labels_col = namesVecNow)
      phBig <- pheatmap::pheatmap(pMatBig, silent = TRUE, angle_col = 45,
                                  show_rownames = FALSE, #main = titleName,
                                  labels_col = namesVecNow)

      # Get TPM for gene
      geneTPMList <- correlationAnalyzeR::getTissueTPM(genesOfInterest = geneNow,
                                                       Species = Species,
                                                       Tissues = "all",
                                                       Sample_Type = whichCompareGroups[1],
                                                       useBlackList = TRUE)
      # Make TPM plot
      geneTPMDF <- data.table::rbindlist(geneTPMList, idcol = "group")
      colnames(geneTPMDF)[3] <- "value"
      geneTPMDF$value <- log2(geneTPMDF$value + 1)
      rawGroup <- geneTPMDF$group
      rawGroup1 <- gsub(rawGroup, pattern = "_.*", replacement = "")
      rawGroup2 <- gsub(rawGroup, pattern = ".*_", replacement = "")
      rawGroup2 <- gsub(rawGroup2, pattern = "0", replacement = " ")
      geneTPMDF$group <- paste0(rawGroup2, " - ", rawGroup1)
      geneTPMDF$group <- stringr::str_to_title(geneTPMDF$group)
      availTissue <- correlationAnalyzeR::getTissueTypes(Species = Species,
                                                         useBlackList = TRUE)
      availTissue <- gsub(availTissue, pattern = "0", replacement = " ")
      availTissue <- stringr::str_to_title(availTissue)
      # all(geneTPMDF$group %in% availTissue) -- should be TRUE
      geneTPMDF <- geneTPMDF[order(match(geneTPMDF$group, availTissue)),]
      TPMBP <- ggpubr::ggboxplot(data = geneTPMDF,
                                 x = "group",
                                 title = paste0(geneNow, " expression across groups"),
                                 ylab = "log2(TPM + 1)",
                                 fill = "group",
                                 y = "value") +
        ggpubr::rotate_x_text() +
        ggpubr::rremove("legend") +
        ggpubr::rremove("xlab")
      colnames(geneTPMDF)[3] <- paste0(geneNow, "_log2TPM")
      resList[[i]][["TPM_DF"]] <- geneTPMDF
      resList[[i]][["TPM_boxPlot"]] <- TPMBP
      resList[[i]][["heatmapSmall"]] <- phSmall
      resList[[i]][["heatmapSmallData"]] <- pMat
      resList[[i]][["heatmapBig"]] <- phBig
      resList[[i]][["heatmapBigData"]] <- pMatBig
    }
    return(resList)
  }

  resList <- list()
  # Main code
  for (i in 1:length(colnames(corrDF))) {
    gene <- colnames(corrDF)[i]
    cat(paste0("\nAnalyzing: ", gene))
    # Create output folder for gene
    geneOutDir <- file.path(outputPrefix, gene)
    if (! dir.exists(geneOutDir) & ! returnDataOnly) {
      dir.create(geneOutDir)
    }
    # Remove gene from correlation values to build normal distribution
    vec <- corrDF[,i]
    names(vec) <- rownames(corrDF)
    vec <- vec[order(vec, decreasing = TRUE)]
    vec <- vec[c(-1)]
    # Make a histogram of gene correlations
    corrDF2 <- corrDF[which(rownames(corrDF) != gene),i, drop = F]
    geneOne <- genesOfInterest[i]
    tissueOne <- Tissue[i]
    tissueOne <- gsub(tissueOne, pattern = "0", replacement = " ")

    sampleOne <- Sample_Type[i]
    geneOneTitle <- paste0(geneOne, ", ",
                           stringr::str_to_title(tissueOne),
                           " - ",
                           stringr::str_to_title(sampleOne))

    p <- ggpubr::gghistogram(data = corrDF2, x = gene, y = "..count..",
                             bins = 100, ylab = "Frequency\n",
                             title = gene,
                             caption = paste0(stringr::str_to_title(tissueOne),
                                              " - ",
                                              stringr::str_to_title(sampleOne)),
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
                                    nperm = 2000,
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

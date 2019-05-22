#' Analyze Gene Pairs
#'
#' Explores how a list of secondary genes relates to a primary gene of interest
#'
#' @param pairedGenesList A list, named with primary genes of interest with
#'     vectors of secondary genes to test against OR a string containing the
#'     official MSIGDB name for a gene set of interest.
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#'     Either "All", "Normal_Tissues", or "Tumor_Tissues".
#'
#' @param outputPrefix Prefix for saved files. Should include directory info.
#'
#' @param sigTest Should the results be compared against random genes?
#'
#' @param nPerm Number of bootstrap sampling events to run during sigTest.
#'
#' @param plotMaxMinCorr If TRUE, the top correlated and anti-correlated genes
#'     will be plotted alongside the selected secondary genes.
#'
#' @param plotLabels If TRUE, correlation histograms will contain labeled lines showing
#'     secondary genes and their correlation values.
#'     If list of secondary genes is large, set this to FALSE or onlyTop to TRUE
#'     to avoid cluttering the plot.
#'
#' @param onlyTop For larger secondary gene lists -- This will filter the
#'     number of secondary genes which are plotted to avoid clutter if plotLabels = TRUE.
#'
#' @param topCutoff The value used for filtering if 'onlyTop' is 'TRUE'
#'
#' @param autoRug If the size of a secondary gene list > 50, plot lines will be replaced
#'     by an auto-generated rug. To disable this behavior, set to FALSE.
#'
#' @param returnDataOnly if TRUE will only return a list containing correlations
#'     and significance testing results if applicable.
#'
#' @return A list containing correlation values and signficance testing results
#'
#' @examples
#' pairedGenesList <- list("TP53" = c("BRCA1", "CDK12", "PARP1"),
#'                         "SON" = c("AURKB", "SFPQ", "DHX9"))
#'
#' pairedGenesList <- list("ATM" = "PUJANA_BRCA1_PCC_NETWORK")
#'
#' Result <- pairedGenesAnalyzeR(pairedGenesList = pairedGenesList,
#'                               Species = "hsapiens",
#'                               Sample_Type = "Normal_Tissues")
#'
#' @export
pairedGenesAnalyzeR <- function(pairedGenesList,
                                Species = c("hsapiens", "mmusculus"),
                                Sample_Type = c("Normal_Tissues",
                                                "Tumor_Tissues"),
                                outputPrefix = "CorrelationAnalyzeR_Output",
                                plotLabels = T,
                                sigTest = T, nPerm = 2000,
                                plotMaxMinCorr = T,
                                onlyTop = F,
                                topCutoff = .5,
                                autoRug = T,
                                returnDataOnly = F) {

  # load("data/hsapiens_complex_TERM2GENE.rda")
  # load("data/mmusculus_complex_TERM2GENE.rda")
  # hs_names <- unique(hsapiens_complex_TERM2GENE$gs_name)
  # mm_names <- unique(mmusculus_complex_TERM2GENE$gs_name)
  # MSIGDB_Geneset_Names <- hs_names[order(hs_names)]
  # # devtools::use_data(MSIGDB_Geneset_Names)
  # load("data/MSIGDB_Geneset_Names.rda")
  #
  # # Bug checking
  # Species = c("hsapiens", "mmusculus")
  # Sample_Type = "Normal_Tissues"
  # returnDataOnly <- F
  # outputPrefix = "tests/pairedOut"
  # # pairedGenesList <- list("ATM" = "MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_UP",
  # #                         "SON" = "TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_UP",
  # #                         "BRCA1" = "BILD_E2F3_ONCOGENIC_SIGNATURE")
  # # pairedGenesList <- list("ATM" = c("TP53", "NFE2L2", "BRCA2"))
  # pairedGenesList <- list("BRCA1" = "PUJANA_BRCA1_PCC_NETWORK")
  # library(correlationAnalyzeR)
  # onlyTop <- F
  # topCutoff <- .5
  # plotLabels <- T
  # sigTest <- T
  # autoRug <- T
  # nPerm <- 2000
  # outputPrefix = "tests/pairedTestFour"
  # plotMaxMinCorr <- T

  # Create output folder
  if (! dir.exists(outputPrefix)) {
    dir.create(outputPrefix)
  }
  # Ensure input data is correct
  if (typeof(pairedGenesList) != "list" | is.null(names(pairedGenesList))) {
    stop("\tFormat pairedGenesList as a named list in which:

        1) Names are the primary genes of interest.
        2) List values contain vectors of secondary genes to compare against -- OR
           a string with the official MSIGDB name of a geneset to compare against.\n
        e.g. pairedGenesList = list('TP53' = c('BRCA1', 'CDK12', 'PARP1'),
                                    'SON' = c('DHX9'),
                                    'MCM2' = c('PCNA', 'STAG2'))  \n")
  }

  # Initialize results object
  resList <- list()

  # Check primary genes to make sure they exist
  avGenes <- getAvailableGenes(Species = Species)
  avGenes <- as.character(avGenes$geneName)
  intGenes <- names(pairedGenesList)
  badGenes <- intGenes[which(! intGenes %in% avGenes)]
  if (length(badGenes) > 0) {
    stop(paste0("\n\t\t\t'", paste(badGenes, collapse = ", "), "' not found
                      in correlation data.
                      Please check available data with getAvailableGenes().
                      Your gene(s) of interest may have an updated name or
                      have a species-specific identifier.\n"))
  }

  # Check secondary genes to make sure they exist -- only a warning
  intGenes_secondary <- unlist(pairedGenesList, use.names = T)
  badGenes_secondary <- intGenes_secondary[which(! intGenes_secondary %in% avGenes &
                                                   ! intGenes_secondary %in% MSIGDB_Geneset_Names)]

  if (length(badGenes_secondary) > 0) {
    warning(paste0("\n\t\t\t'", paste(badGenes, collapse = ", "), "'
                      not found in correlation data and is not an official MSIGDB name.
                      Please check available gene data with getAvailableGenes().
                      Your gene(s) of interest may have an updated name or
                      have a species-specific identifier. Find offical MSIGDB
                      names by examining the MSIGDB_Geneset_Names object.\n
                      Continuing without this/these gene(s)..."))
  }

  # Make list of terms inputted by the user
  termGenes_secondary <- intGenes_secondary[which(! intGenes_secondary %in% avGenes &
                                                    intGenes_secondary %in% MSIGDB_Geneset_Names)]
  if (length(termGenes_secondary > 0)) {
    if (Species[1] == "hsapiens") {
      TERM2GENE <- hsapiens_complex_TERM2GENE
    } else {
      TERM2GENE <- mmusculus_complex_TERM2GENE
    }
    for(i in 1:length(termGenes_secondary)) {
      term <- termGenes_secondary[i]
      print(term)
      nameStr <- names(term)
      termGenes <- TERM2GENE$human_gene_symbol[which(TERM2GENE$gs_name == term)]
      termGenes <- termGenes[which(termGenes %in% avGenes)] # Ensure actionable genes
      pairedGenesList[[nameStr]] <- termGenes
    }
  }


  cat("\nRetrieving any missing correlation data...\n")
  # Call downloadData to get all required files
  downloadData(Species = Species[1],
               Sample_Type = Sample_Type[1],
               geneList = intGenes)
  downloadFolder <- system.file("data", package = "correlationAnalyzeR")
  downloadFolder <- file.path(downloadFolder, "Correlation_Data",
                              Species[1], Sample_Type[1])
  # Main code
  for (i in 1:length(pairedGenesList)) {
    gene <- names(pairedGenesList)[i]
    cat(paste0("\n", gene))
    # Load gene data
    geneFile <- paste0(gene, ".RData")
    file <- file.path(downloadFolder, geneFile)
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
    # Create correlation dataframe object
    corrDF <- as.data.frame(vec)
    corrDF$geneName <- names(vec)
    colnames(corrDF)[1] <- "correlationValue"
    # Filter for secondary genes
    secondaryGenes <- pairedGenesList[[i]]
    secondaryGenes <- secondaryGenes[which(secondaryGenes != gene)]
    corrDF <- corrDF[which(corrDF$geneName != gene),]
    corrDF$secondaryGene <- F
    corrDF$secondaryGene[which(corrDF$geneName %in% secondaryGenes)] <- T
    res <- numeric()

    # Check to make sure secondary genes exist within data
    for (j in 1:length(secondaryGenes)) {
      if (! secondaryGenes[j] %in% corrDF$geneName) {
        warning(paste0("\n\t\t\t'", secondaryGenes[j], "' not found
                      in correlation data.
                      Please check available data with getAvailableGenes().
                      Your gene(s) of interest may have an updated name or
                      have a species-specific identifier.\n"))
        secondaryGenes <- secondaryGenes[c(-j)]

      }
    }

    # Create histogram of gene correlations
    his <- ggplot2::ggplot(data = corrDF,
                           mapping = ggplot2::aes(x = correlationValue)) +
      ggplot2::geom_histogram(bins = 100, color = "black", fill = "white") +
      ggplot2::labs(title = paste0("Histogram of ",
                                   gene,
                                   " correlations")) +
      ggplot2::ylab("Frequency") +
      ggplot2::scale_y_continuous(limits = c(0, 3500), expand = c(0,0)) +
      ggplot2::theme_classic()
    # Add min-max vals
    if (plotMaxMinCorr) {
      maxGene <- corrDF$geneName[which.max(corrDF$correlationValue)]
      maxVal <- max(corrDF$correlationValue)
      minGene <- corrDF$geneName[which.min(corrDF$correlationValue)]
      minVal <- min(corrDF$correlationValue)
      maxLab <- paste0(maxGene, " (", round(maxVal, 2), ")\n")
      minLab <- paste0("\n", minGene, " (", round(minVal, 2), ")")
      # Add in max line
      his <- his + ggplot2::geom_vline(xintercept = maxVal, linetype = 3)
      if (plotLabels) {
        his <- his + ggplot2::annotate(geom = "text", x = maxVal, y = 2000,
                                       label = maxLab,
                                       color = "black", angle = 90, size = 2)
      }
      # Add in min line
      his <- his + ggplot2::geom_vline(xintercept = minVal, linetype = 3)

      if (plotLabels) {
        his <- his + ggplot2::annotate(geom = "text", x = minVal, y = 2000,
                                       label = minLab,
                                       color = "black", angle = 90, size = 2)
      }
    }

    # Harmonize original corrDF frame
    corrDF_small <- corrDF[which(corrDF$geneName %in% secondaryGenes),]
    # Choose top genes
    if (onlyTop) {
      corrDF_small <- corrDF_small[order(abs(corrDF_small$correlationValue), decreasing = T),]
      takeNum <- round(length(corrDF_small$geneName) * (1-topCutoff))
      corrDF_small <- corrDF_small[c(1:takeNum),]
    }
    secondaryGenes <- corrDF_small$geneName
    errorCounter <- 0

    # Warn if plot is going to be cluttered
    if (! onlyTop & plotLabels & length(secondaryGenes) > 8) {
      warning(paste0("\n\t\t\tAttempting to plot '", length(secondaryGenes), "' secondary genes, along with
                        text labels on one histogram. This will probably
                        produce a cluttered plot. Consider setting onlyTop = TRUE,
                        plotLabels = FALSE, or choosing fewer secondary genes to plot.\n"))
    }


    # Code to analyze how histogram should be formatted
    if (length(secondaryGenes) > 50 & autoRug) {
      warning("Due to large secondary gene size -- plot lines will be implemented as a rug.
              You may disable this behavior by setting autoRug = F.")
      # Create rug by subsetting original data
      his <- his + ggplot2::geom_rug(data = corrDF[which(corrDF$secondaryGene),], color = "red")

    }

    for (j in 1:length(secondaryGenes)) {
      secondaryGene <- secondaryGenes[j]
      corrVal <- corrDF$correlationValue[which(
        corrDF$geneName == secondaryGene
        )]
      res[j] <- corrVal
      if (length(secondaryGenes) > 50 & autoRug) {
        next
      }
      if (! plotLabels) {
        his <- his + ggplot2::geom_vline(xintercept = corrVal, linetype = 2, color = "red")
      } else {
        # List of previous values excluding current one.
        res2 <- res[c(-j)]
        # Are any previous values similar to the current one?
        corrPrev <- which(abs(round(res2, 2) - round(corrVal,2)) < .02)
        if (length(corrPrev) == 1) {
          corrPrevVal <- res[corrPrev]
          if (corrPrevVal < corrVal ) {
            # Code to change direction of label if new corr val is greater
            ind <- (corrPrev-errorCounter)*2 + 1
            his$layers[[ind]]$aes_params$label <- paste0(
              substr(his$layers[[ind]]$aes_params$label, 2, 100
              ), "\n")
            label <- paste0("\n", secondaryGene, " (", round(corrVal, 2), ")")
          } else {
            label <- paste0(secondaryGene, " (", round(corrVal, 2), ")\n")
          }
        } else if (length(corrPrev) > 1) {
          warning(paste0("Difficulty plotting with labels for ", gene,
                         " -- more than two values are similar"))
          corrPrevVal <- res[corrPrev[1]]
          if (corrPrevVal < corrVal ) {
            # Code to change direction of label if new corr val is greater
            ind <- corrPrev*2 + 1
            ind <- ind[1]
            his$layers[[ind]]$aes_params$label <- paste0(
              substr(his$layers[[ind]]$aes_params$label, 2, 100
              ), "\n")
            label <- paste0("\n", secondaryGene, " (", round(corrVal, 2), ")")
          } else {
            label <- paste0(secondaryGene, " (", round(corrVal, 2), ")\n")
          }
          label <- paste0("\n", secondaryGene, " (", round(corrVal, 2), ")")
        } else {
          label <- paste0("\n", secondaryGene, " (", round(corrVal, 2), ")")
        }
        if (abs(corrVal) < .15) {
          his <- his +
            ggplot2::geom_vline(xintercept = corrVal, linetype = 3) +
            ggplot2::annotate(geom = "text", x = corrVal,
                              y = sample(c(1800:2600), 1),
                              label = label,
                              color = "red", angle = 90, size = 2)
        } else {
          his <- his +
            ggplot2::geom_vline(xintercept = corrVal, linetype = 3) +
            ggplot2::annotate(geom = "text", x = corrVal,
                              y = sample(c(1500:2200), 1),
                              label = label,
                              color = "red", angle = 90, size = 2)
        }
      }
    }
    # Finalizes histogram
    his <- his + ggplot2::theme(text = ggplot2::element_text(size = 10))
    if (! returnDataOnly) {
      ggplot2::ggsave(plot = his, width = 8, height = 5, filename = file.path(geneOutDir,
                                                       paste0(gene, ".png")))
    }

    names(res) <- secondaryGenes

    # Significance testing
    if (sigTest) {
      # Get correlation values for selected genes
      origSamples <- c(pairedGenesList[[i]])
      selectVec <- vec[origSamples]

      # Construct distribution for random at same input size
      n <- length(selectVec)
      # Function to boot
      meanBoot <- function(data, indices) {
        d <- data[indices] # allows boot to select sample
        d <- sample(d, size = n)
        ttest <- t.test(x = selectVec, y = d)

        c(mean(d), median(d), ttest$p.value)
      }

      # Reproducibility
      set.seed(1)

      # bootstrapping with 1000 replications
      results <- boot::boot(data=vec, statistic=meanBoot,
                            R=nPerm)

      # Get value for input data
      Mean <- mean(selectVec)
      Median <- median(selectVec)

      # Get bootstrapped values
      meanVec <- results$t[,1]
      medianVec <- results$t[,2]
      ttestVec <- results$t[,3]

      # Build plots
      tdf <- as.data.frame(results$t)
      colnames(tdf) <- c("Means", "Medians", "TTest_pVals")
      # Means
      pMeans <- ggpubr::gghistogram(data = tdf, bins = 40,
                                    x = "Means", xlab = "Bootstrapped means",
                                    title = paste0(gene, " correlation means distribution")) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = Mean, color = "Selected Genes Mean")) +
        ggplot2::scale_colour_manual("",
                                     breaks = c("Selected Genes Mean"),
                                     values = c("red"))
      # Means
      pMedians <- ggpubr::gghistogram(data = tdf, bins = 40,
                                    x = "Medians", xlab = "Bootstrapped medians",
                                    title = paste0(gene, " correlation medians distribution")) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = Mean, color = "Selected Genes Median")) +
        ggplot2::scale_colour_manual("",
                                     breaks = c("Selected Genes Median"),
                                     values = c("red"))

      # TTest pVals
      Threshold <- .05
      pTtest <- ggpubr::gghistogram(data = tdf, bins = 40,
                                      x = "TTest_pVals", xlab = "Bootstrapped T-Test pVals (selected vs random)",
                                      title = paste0(gene, " correlations (selected vs random) pval distribution")) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = Threshold, color = "Threshold (p < .05)")) +
        ggplot2::scale_colour_manual("",
                            breaks = c("Threshold (p < .05)"),
                            values = c("grey"))

      # Plot if not returnDataOnly
      if (! returnDataOnly) {
        ggplot2::ggsave(plot = pMeans, filename = file.path(geneOutDir,
                                                      paste0(gene, ".bootstrapped_means.png")),
                        width = 8, height = 5)
        ggplot2::ggsave(plot = pMedians, filename = file.path(geneOutDir,
                                                            paste0(gene, ".bootstrapped_medians.png")),
                        width =8, height = 5)
        ggplot2::ggsave(plot = pTtest, filename = file.path(geneOutDir,
                                                            paste0(gene, ".bootstrapped_TTest_PVals.png")),
                        width = 8, height = 5)
      }
    }

    # Return results based on whether sigTest is true
    if (! sigTest) {
      resList[[i]] <- list("Correlation_Values" = res,
                           "Correlation_histogram" = his)
    } else {
      resList[[i]] <- list("Correlation_Values" = res,
                           "Correlation_histogram" = his,
                           "sigTest" =  list("means" = meanVec,
                                             "meansPlot" = pMeans,
                                             "medians" = medianVec,
                                             "mediansPlot" = pMedians,
                                             "tTest_pvals" = ttestVec,
                                             "tTest_pvalsPlot" = pTtest))
    }
    names(resList)[i] <- gene

  }
  return(resList)
}



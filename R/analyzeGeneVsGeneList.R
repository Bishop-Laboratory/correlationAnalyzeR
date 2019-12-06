#' Analyze gene vs gene list relationship
#'
#' Explores how a list of secondary genes relates to a primary gene of interest
#'
#' @param pairedGenesList A list, named with primary genes of interest with
#' vectors of secondary genes to test against OR a string containing the
#' official MSIGDB name for a gene set of interest.
#'
#' @param Species Species to obtain gene names for.
#' Either 'hsapiens' or 'mmusculus'
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
#' @param outputPrefix Prefix for saved files. Should include directory info.
#'
#' @param sigTest Should the results be compared against random? Default: TRUE.
#'
#' @param nPerm Number of bootstrap sampling events to run during sigTest. Default: 2000.
#'
#' @param plotMaxMinCorr If TRUE, the top correlated and anti-correlated genes
#' will be plotted alongside the selected secondary genes. Default: TRUE.
#'
#' @param plotLabels If TRUE, correlation histograms will contain labeled lines showing
#' secondary genes and their correlation values.
#' If list of secondary genes is large, set this to FALSE or onlyTop to TRUE
#' to avoid cluttering the plot. Default: TRUE.
#'
#' @param plotTitle Logical. If TRUE, plot title will be added to visualizations. Default: TRUE.
#'
#' @param onlyTop For larger secondary gene lists -- This will filter the
#' number of secondary genes which are plotted to avoid clutter if plotLabels = TRUE.
#' Default: FALSE.
#'
#' @param topCutoff The value used for filtering if 'onlyTop' is 'TRUE'. Default: .5
#'
#' @param autoRug If the size of a secondary gene list > 50, plot lines will be replaced
#' by an auto-generated rug. Default: TRUE.
#'
#' @param returnDataOnly if TRUE will only return a list containing correlations
#' and significance testing results if applicable. Default: TRUE.
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied.
#' @return A list containing correlation values and signficance testing results
#'
#' @examples
#' pairedGenesList <- list("TP53" = c("BRCA1", "CDK12", "PARP1"),
#'                         "SON" = c("AURKB", "SFPQ", "DHX9"))
#'
#'correlationAnalyzeR::geneVsGeneListAnalyze(pairedGenesList = pairedGenesList,
#'                               Species = "hsapiens", returnDataOnly = TRUE,
#'                               Sample_Type = "normal",
#'                               Tissue = "brain")
#'
#' @export
geneVsGeneListAnalyze <- function(pairedGenesList,
                                  Species = c("hsapiens", "mmusculus"),
                                  Sample_Type = c("normal", "cancer"),
                                  Tissue = "all",
                                  outputPrefix = "CorrelationAnalyzeR_Output",
                                  plotLabels = TRUE,
                                  sigTest = TRUE, nPerm = 2000,
                                  plotMaxMinCorr = TRUE,
                                  onlyTop = FALSE,
                                  topCutoff = .5,
                                  autoRug = TRUE,
                                  plotTitle = TRUE,
                                  returnDataOnly = TRUE,
                                  pool = NULL) {

  # pairedGenesList = list("ATM" = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
  # Species = c("hsapiens", "mmusculus")
  # Sample_Type = c("normal", "cancer")
  # Tissue = "all"
  # outputPrefix = "CorrelationAnalyzeR_Output"
  # plotLabels = TRUE
  # sigTest = TRUE
  # nPerm = 2000
  # plotMaxMinCorr = TRUE
  # onlyTop = FALSE
  # topCutoff = .5
  # autoRug = TRUE
  # plotTitle = TRUE
  # returnDataOnly = TRUE

  if (is.null(pool)) {
    retryCounter <- 1
    cat("\nEstablishing connection to database ... \n")
    while(is.null(pool)) {
      pool <- try(silent = T, eval({
        pool::dbPool(
          drv = RMySQL::MySQL(),
          user = "public-rds-user", port = 3306,
          dbname="bishoplabdb",
          password='public-user-password',
          host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com"
        )
      }))
      if ("try-error" %in% class(pool)) {
        if (retryCounter == 3) {
          stop("Unable to connect to database. Check internet connection and please contanct",
               " package maintainer if you believe this is an error.")
        }
        warning(paste0("Failed to establish connection to database ... retrying now ... ",
                       (4-retryCounter), " attempts left."),
                immediate. = T)
        pool <- NULL
        retryCounter <- retryCounter + 1
      }
    }

    on.exit(function() {
      pool::poolClose(pool)
    })
  }


  # Create output folder
  if (! dir.exists(outputPrefix) & ! returnDataOnly) {
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
  avGenes <- correlationAnalyzeR::getAvailableGenes(Species = Species, pool = pool)
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
  intGenes_secondary <- unlist(pairedGenesList, use.names = TRUE)
  badGenes_secondary <- intGenes_secondary[
    which(! intGenes_secondary %in% avGenes &
          ! intGenes_secondary %in% correlationAnalyzeR::MSIGDB_Geneset_Names)
    ]

  if (length(badGenes_secondary) > 0) {
    warning(paste0("\n\t\t\t'", paste(badGenes_secondary, collapse = ", "), "'
                      not found in correlation data and is not an official MSIGDB name.
                      Please check available gene data with getAvailableGenes().
                      Your gene(s) of interest may have an updated name or
                      have a species-specific identifier. Find offical MSIGDB
                      names by examining the MSIGDB_Geneset_Names object.\n
                      Continuing without this/these gene(s)..."))
  }

  # Make list of terms inputted by the user
  termGenes_secondary <- intGenes_secondary[
    which(! intGenes_secondary %in% avGenes &
            intGenes_secondary %in% correlationAnalyzeR::MSIGDB_Geneset_Names)
    ]
  if (length(termGenes_secondary > 0)) {
    TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "complex",
                                                   Species = Species)
    for(i in 1:length(termGenes_secondary)) {
      term <- termGenes_secondary[i]
      print(term)
      nameStr <- names(term)
      termGenes <- TERM2GENE$gene_symbol[which(TERM2GENE$gs_name == term)]
      termGenes <- termGenes[which(termGenes %in% avGenes)] # Ensure actionable genes
      pairedGenesList[[nameStr]] <- termGenes
    }
  }


  # Call getCorrelationData to get all required files
  corrDFFull <- correlationAnalyzeR::getCorrelationData(Species = Species[1],
                                                        Tissue = Tissue[1],
                                                        Sample_Type = Sample_Type[1],
                                                        geneList = intGenes,
                                                        pool = pool)
  colLengths <- ifelse(Species[1] == "hsapiens", yes = 28685, no = 24924)
  pVals <- apply(as.matrix(corrDFFull), MARGIN = 1:2, n = colLengths, FUN = function(x, n) {
    dt(abs(x)/sqrt((1-x^2)/(n-2)), df = 2)
  })
  pVals <- as.data.frame(pVals)

  # Main code
  for (i in 1:length(colnames(corrDFFull))) {
    gene <- colnames(corrDFFull)[i]
    cat(paste0("\n", gene))
    # Create output folder for gene
    geneOutDir <- ifelse(length(colnames(corrDFFull)) == 1,
                         yes = outputPrefix,
                         no = file.path(outputPrefix, gene))
    if (! dir.exists(geneOutDir) & ! returnDataOnly) {
      dir.create(geneOutDir)
    }
    # Filter for secondary genes
    secondaryGenes <- pairedGenesList[[gene]]
    secondaryGenes <- secondaryGenes[which(secondaryGenes != gene)] # Incase they added the original gene in
    corrDF <- data.frame(geneName = rownames(corrDFFull),
                         correlationValue = corrDFFull[,gene],
                         stringsAsFactors = FALSE)
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
                           mapping = ggplot2::aes_string(x = "correlationValue")) +
      ggplot2::geom_histogram(bins = 100, color = "black", fill = "white") +
      ggplot2::ylab("Frequency\n") +
      ggplot2::scale_y_continuous(limits = c(0, 3500), expand = c(0,0)) +
      ggplot2::theme_classic()
    if (plotTitle) {
      his <- his + ggplot2::labs(title = paste0("Histogram of ",
                                                gene,
                                                " correlations\n"))
    }
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

    # Cut the corrDF frame to show only the secondary genes
    corrDF_small <- corrDF[which(corrDF$secondaryGene),]
    # Choose top genes
    if (onlyTop) {
      corrDF_small <- corrDF_small[order(abs(corrDF_small$correlationValue), decreasing = TRUE),]
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
      vec <- corrDF$correlationValue
      names(vec) <- corrDF$geneName
      vec <- vec[which(! is.na(vec))]
      selectVec <- vec[secondaryGenes]

      # Construct distribution for random at same input size
      n <- length(selectVec)
      # Function to boot
      meanBoot <- function(data, indices) {
        d <- data[indices] # allows boot to select sample
        d <- sample(d, size = n)
        ttest <- stats::t.test(x = abs(selectVec), y = abs(d)) # Considers absolute value comparison
        c(mean(d), stats::median(d), ttest$p.value)
      }
      # Function to give significance stars
      getSigStars <- function(pVal) {
        pMap <- c(.05, .01, .001, .0001)
        names(pMap) <- c("*", "**", "***", "****")
        sig <- ifelse(test = pVal < pMap[4],
                      yes = names(pMap)[4],
                      no = ifelse(test = pVal < pMap[3],
                                  yes = names(pMap)[3],
                                  no = ifelse(pVal < pMap[2],
                                              yes = names(pMap)[2],
                                              no = ifelse(test = pVal < pMap[1],
                                                          yes = names(pMap)[1],
                                                          no = "n.s."))))
      }

      # bootstrapping with # replications
      results <- boot::boot(data=vec, statistic=meanBoot,
                            R=nPerm)

      # Get value for input data
      Meanabs <- mean(abs(selectVec), na.rm = TRUE) # For pvalue calc
      Medianabs <- stats::median(abs(selectVec), na.rm = TRUE) # For pvalue calc
      Mean <- mean((selectVec), na.rm = TRUE) # For plotting
      Median <- stats::median((selectVec), na.rm = TRUE) # For plotting

      # Get bootstrapped values
      meanVec <- results$t[,1]
      medianVec <- results$t[,2]
      ttestVec <- results$t[,3]

      # Build plots
      tdf <- as.data.frame(results$t)
      colnames(tdf) <- c("Means", "Medians", "TTest_pVals")
      # Means

      # Determine if the abs mean of selected genes
      pValMeans <- sum(abs(tdf$Means) > Meanabs)/sum(! is.na(tdf$Means))
      meanStr <- paste0("Selected Genes Mean, p = ", pValMeans, " [", getSigStars(pValMeans), "]")
      pMeans <- ggpubr::gghistogram(data = tdf, bins = 60, ylab = "Frequency\n",
                                    title = switch(plotTitle + 1, NULL,
                                                   paste0(gene, " correlation means distribution")),
                                    x = "Means", xlab = "Bootstrapped means") +
        ggplot2::geom_vline(ggplot2::aes(xintercept = Mean, color = meanStr)) +
        ggplot2::scale_colour_manual("",
                                     breaks = meanStr,
                                     values = c("red")) +
        ggplot2::scale_y_continuous(expand = c(0,0))

      # Means
      pValMedians <- sum(abs(tdf$Medians) > Medianabs)/sum(! is.na(tdf$Medians))
      medianStr <- paste0("Selected Genes Median, p = ", pValMedians, " [", getSigStars(pValMedians), "]")
      pMedians <- ggpubr::gghistogram(data = tdf, bins = 60, ylab = "Frequency\n",
                                    x = "Medians", xlab = "Bootstrapped medians",
                                    title = switch(plotTitle + 1, NULL,
                                                   paste0(gene, " correlation medians distribution"))) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = Mean, color = medianStr)) +
        ggplot2::scale_colour_manual("",
                                     breaks = medianStr,
                                     values = c("blue")) +
        ggplot2::scale_y_continuous(expand = c(0,0))

      # TTest pVals
      title = paste0(gene, " correlations (selected vs random) pval distribution")
      denY <- which.max(stats::density(tdf$TTest_pVals)$y)
      pValTTest <- stats::density(tdf$TTest_pVals)$x[denY]
      pValStr <- paste0("P(summit) = ", round(pValTTest, digits = 3), " [", getSigStars(pValTTest), "]")
      pTtest <- ggpubr::ggdensity(data = tdf, ylab = "P value density\n",
                                      x = "TTest_pVals",
                                  xlab = "Bootstrapped T-Test pVals (selected vs random)",
                                  title = switch(plotTitle + 1, NULL,
                                                 paste0(gene, " permutation t-test p-val distribution"))) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = pValTTest, color = pValStr)) +
        ggplot2::scale_colour_manual("",
                            breaks = pValStr,
                            values = c("grey")) +
        ggplot2::scale_y_continuous(expand = c(0,0))

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
                           "P values" = pVals,
                           "Correlation_histogram" = his)
    } else {
      resList[[i]] <- list("Correlation_Values" = res,
                           "P values" = pVals,
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
  resList[["P values"]] <- pVals
  resList[["correlations"]] <- corrDFFull
  return(resList)
}



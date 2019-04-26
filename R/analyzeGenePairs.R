require(ggplot2)
require(ggpubr)
# Main function
pairedGenesAnalyzeR <- function(pairedGenesList,
                                Species = c("hsapiens", "mmusculus"),
                                Sample_Type = c("All",
                                                "Normal_Tissues",
                                                "Tumor_Tissues"),
                                outputPrefix = "CorrelationAnalyzeR_Output",
                                sigTest = T,
                                plotMaxMinCorr = T,
                                onlyTop = F,
                                topCutoff = .5) {
  # Bug checking
  Species = c("hsapiens", "mmusculus")
  Sample_Type = "Normal_Tissues"
  pairedGenesList <- list("TP53" = c("BRCA1", "CDK12", "PARP1"),
                          "SON" = c("AURKB", "SFPQ", "DHX9"),
                          "MCM2" = c("POLE", "PCNA", "STAG2"),
                          "REST" = c("ERCC3", "MAPT"),
                          "NFE2L2" = c("GSTP1", "KEAP1", "HMOX1", "GCLM"),
                          "ATM" = c("SLC7A11", "SLC3A2", "CHEK2", "HMOX1", "GCLM"),
                          "AURKB" = c("ATR", "CENPI", "DHX9", "RNASEH2A"),
                          "SETX" = c("SON", "PARP1", "DHX9", "BRCA1", "TOP2A"))
  # load("Data/geneCorrelationDataFiles/A1BG.RData")
  # geneNames <- names(corr)
  # set.seed(4)
  # pairedGenesList <- list(sample(geneNames, 3),
  #                         sample(geneNames, 2),
  #                         sample(geneNames, 4))
  # names(pairedGenesList) <- sample(geneNames, 3)
  outputPrefix = "tests/pairedTestOne"
  plotMaxMinCorr <- T

  # Create output folder
  if (! dir.exists(outputPrefix)) {
    dir.create(outputPrefix)
  }
  # Ensure input data is correct
  if (typeof(pairedGenesList) != "list" | is.null(names(pairedGenesList))) {
    stop("\tFormat pairedGenesList as a named list in which:

        1) Names are the primary genes of interest.
        2) List values contain vectors of secondary genes to compare against.\n
        e.g. pairedGenesList = list('TP53' = c('BRCA1', 'CDK12', 'PARP1'),
                                    'SON' = c('DHX9'),
                                    'MCM2' = c('PCNA', 'STAG2'))  \n")
  }

  # Initialize results object
  resList <- list()

  # Check genes to make sure they exist
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
  cat("\nRetrieving any missing correlation data...\n")
  # Call downloadData to get all required files
  downloadData(Species = Species,
               Sample_Type = Sample_Type,
               geneList = intGenes)
  downloadFolder <- system.file("data", package = "correlationAnalyzeR")
  downloadFolder <- file.path(downloadFolder, "Correlation_Data",
                              Species, Sample_Type)
  # Main code
  for (i in 1:length(pairedGenesList)) {
    gene <- names(pairedGenesList)[i]
    cat(paste0("\n", gene))
    # Load gene data
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
    # Create correlation dataframe object
    corrDF <- as.data.frame(vec)
    corrDF$geneName <- row.names(vec)
    colnames(corrDF)[1] <- "correlationValue"
    # Filter for secondary genes
    secondaryGenes <- pairedGenesList[[i]]
    corrDF <- corrDF[which(corrDF$geneName != primaryGene),]
    res <- numeric()
    # Create histogram of gene correlations
    his <- ggplot(corrDF, aes(x = correlationValue)) +
      geom_histogram(bins = 100, color = "black", fill = "white") +
      labs(title = paste0("Histogram of ", primaryGene, " correlations")) +
      ylab("Frequency") +
      scale_y_continuous(limits = c(0, 7500), expand = c(0,0)) +
      theme_classic()
    # Add min-max vals
    if (plotMaxMinCorr) {
      maxGene <- corrDF$geneName[which.max(corrDF$correlationValue)]
      maxVal <- max(corrDF$correlationValue)
      minGene <- corrDF$geneName[which.min(corrDF$correlationValue)]
      minVal <- min(corrDF$correlationValue)
      maxLab <- paste0(maxGene, " (", round(maxVal, 2), ")\n")
      minLab <- paste0("\n", minGene, " (", round(minVal, 2), ")")
      # Add in max line
      his <- his + geom_vline(xintercept = maxVal, linetype = 3) +
        annotate(geom = "text", x = maxVal, y = 4000,
                 label = maxLab,
                 color = "black", angle = 90, size = 2)
      # Add in min line
      his <- his + geom_vline(xintercept = minVal, linetype = 3) +
        annotate(geom = "text", x = minVal, y = 4000,
                 label = minLab,
                 color = "black", angle = 90, size = 2)
    }
    # Check to make sure secondary genes exist within data
    for (j in 1:length(secondaryGenes)) {
      if (! secondaryGenes[j] %in% corrDF$geneName) {
        secondaryGenes <- secondaryGenes[c(-j)]
        warning(paste0("\n\t\t\t'", secondaryGenes[j], "' not found in correlation data.
                      Please check available data with getAvailableGenes().
                      Your gene(s) of interest may have an updated name or
                      have a species-specific identifier.\n"))
      }
    }
    # Harmonize original corrDF frame
    corrDF <- corrDF[which(corrDF$geneName %in% secondaryGenes),]
    # Choose top genes
    if (onlyTop) {
      corrDF <- corrDF[order(abs(corrDF$correlationValue), decreasing = T),]
      takeNum <- round(length(corrDF$geneName) * (1-topCutoff))
      corrDF <- corrDF[c(1:takeNum),]
    }
    secondaryGenes <- corrDF$geneName
    errorCounter <- 0
    # Code to analyze how histogram should be formatted
    for (j in 1:length(secondaryGenes)) {
      secondaryGene <- secondaryGenes[j]
      corrVal <- corrDF$correlationValue[which(corrDF$geneName == secondaryGene)]
      corrPrev <- which(abs(round(corrVal,2) - round(res, 2)) < .02)
      if (length(corrPrev) == 1) {
        corrPrevVal <- res[corrPrev]
        if (corrPrevVal < corrVal ) {
          # Code to change direction of label if new corr val is greater
          ind <- (corrPrev-errorCounter)*2 + 1
          his$layers[[ind]]$aes_params$label <- paste0(substr(his$layers[[ind]]$aes_params$label, 2, 100), "\n")
          label <- paste0("\n", secondaryGene, " (", round(corrVal, 2), ")")
        } else {
          label <- paste0(secondaryGene, " (", round(corrVal, 2), ")\n")
        }
      } else if (length(corrPrev) > 1) {
        warning(paste0("Difficulty plotting with labels for ", primaryGene, " -- more than two values are similar"))
        corrPrevVal <- res[corrPrev[1]]
        if (corrPrevVal < corrVal ) {
          # Code to change direction of label if new corr val is greater
          ind <- corrPrev*2 + 1
          ind <- ind[1]
          his$layers[[ind]]$aes_params$label <- paste0(substr(his$layers[[ind]]$aes_params$label, 2, 100), "\n")
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
          geom_vline(xintercept = corrVal, linetype = 3) +
          annotate(geom = "text", x = corrVal, y = sample(c(4000:6500), 1),
                   label = label,
                   color = "red", angle = 90, size = 2)
      } else {
        his <- his +
          geom_vline(xintercept = corrVal, linetype = 3) +
          annotate(geom = "text", x = corrVal, y = sample(c(1500:2700), 1),
                   label = label,
                   color = "red", angle = 90, size = 2)
      }
      res[j] <- corrVal
    }
    # Finalizes histogram
    his <- his + theme(text = element_text(size = 10))
    ggsave(plot = his, filename = file.path(geneOutDir, paste0(primaryGene, ".png")))
    names(res) <- secondaryGenes
    resList[[i]] <- res
    names(resList)[i] <- primaryGene
    # Significance testing
    if (sigTest) {
      origSamples <- c(primaryGene, pairedGenesList[[i]])
      # Get random sample of the same size as user input
      randGenesVec <- sample(names(corr)[which(! names(corr) %in% origSamples)],
                             length(res))
      # Convert to dataframe
      corrRand <- data.frame("geneName" = randGenesVec,
                             "correlationValue" = corr[randGenesVec])
      corrRand$Group <- "randomGenes"
      # Build dataframe for result values
      res <- data.frame("geneName" = names(res),
                        "correlationValue" = as.numeric(res))
      res$Group <- "testGenes"
      corrTest <- rbind(res, corrRand)
      corrTest$correlationValue <- abs(corrTest$correlationValue)
      bp <- ggboxplot(data = corrTest, x = "Group", y = "correlationValue",
                color = "Group", palette = "jco",
                ylab = "Absolute Correlation Value\n",
                title = paste0(primaryGene,
                               " Gene Correlation Significance Test\n")) +
        stat_compare_means(method = "t.test",
                           method.args = list(alternative = "greater"),
                           comparisons = list(c("testGenes", "randomGenes"))) +
        rremove("legend") + rremove('xlab')
      ggsave(plot = bp, filename = file.path(geneOutDir,
                                             paste0(primaryGene,
                                                    ".sigTest.png")))
      resSave <- paste0(primaryGene, "_T.Test")
      corrTest$Group <- factor(corrTest$Group)
      tt <- t.test(x = corrTest$correlationValue[which(
        corrTest$Group == "testGenes"
        )],
                   alternative = "greater",
             y = corrTest$correlationValue[which(
               corrTest$Group == "randomGenes"
               )], paired = F)
      saveVals <- list("t.test" = tt,
                       "correlationTestValues" = corrTest)
      resList[[resSave]] <- saveVals
    }
  }
  return(resList)
}



#' Analyze Gene List Topology
#'
#' Analyzes the topology of a gene list using gene correlation data and dimension-reduction techniques.
#'
#' @param genesOfInterest A vector of genes to analyze or the name of an official MSIGDB term.
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#'     Either "All", "Normal_Tissues", or "Tumor_Tissues".
#'
#' @param outputPrefix Prefix for saved files. Should include directory info.
#'
#' @param returnDataOnly if TRUE will return only a list of analysis results.
#'
#' @param setComparisonCutoff Only relevant for cocorrelation analysis -- the number of genes which
#'      must aggree for a gene to be considered co-correlative within the input gene list.
#'
#' @param numTopGenesToPlot When creating a heatmap of the top co-correlative or top variant genes,
#'      how many genes should be plotted on the y axis?
#'
#' @param crossComparisonType The type of topology tests to run. (see details)
#'
#' @return A list of analysis results and plotting data.
#'
#' @examples
#' genesOfInterest <- c("CDK12", "AURKB", "SFPQ", "PARP1", "BRCC3", "BRCA2", "PARP1",
#'                      "DHX9", "SON", "AURKA", "SETX", "BRCA1", "ATMIN")
#' Result <- analyzeGenesetTopology(genesOfInterest = genesOfInterest,
#'                                  Species = "hsapiens",
#'                                  Sample_Type = "Normal_Tissues",
#'                                  crossComparisonType = c("variantGenes", "PCA"))
#'
#' @details
#'
#' Cross Comparison Types:
#' - variantGenes: These are the genes which best explain variation between genes within the input list.
#'                 These genes can divide a list into functional groups.
#' - coCorrelativeGenes: These are the genes which best explain similarities between all genes in the input list.
#'                       These genes can explain what biological processes unify the input genes.
#' - PCA: This is a dimensionality reduction technique for exploring the topology of a gene list.
#'        The PCA analyses here employes hclust to divide the gene list into functional clusters.
#'        If the input list is > 100 genes, RTsne will be used for visualization.
#'
#' @export
analyzeGenesetTopology <-  function(genesOfInterest,
                                    Species = c("hsapiens", "mmusculus"),
                                    Sample_Type = c("Normal_Tissues",
                                                    "Tumor_Tissues"),
                                    crossComparisonType = c("variantGenes",
                                                            "coCorrelativeGenes",
                                                            # "hierarchicalClustering",
                                                            "PCA"),
                                    setComparisonCutoff = "Auto",
                                    numTopGenesToPlot = "Auto",
                                    outputPrefix = "CorrelationAnalyzeR_Output",
                                    returnDataOnly = F) {

  # # Test parameters for debugging
  # genesOfInterest <- c("CDK12", "AURKB", "SFPQ", "PARP1", "BRCC3", "BRCA2", "PARP1",
  #                      "DHX9", "SON", "AURKA", "SETX", "BRCA1", "ATMIN")
  # genesOfInterest <- "PUJANA_BRCA1_PCC_NETWORK"
  # outputPrefix = "tests/topologyOutput1"
  # setComparisonCutoff = "Auto"
  # numTopGenesToPlot = "Auto"
  # Species = "hsapiens"
  # Sample_Type <- "Normal_Tissues"
  # alternativeTSNE <- T
  # numClusters = "Auto"
  # returnDataOnly <- F
  # crossComparisonType = c("variantGenes",
  #                         "coCorrelativeGenes",
  #                         "hierarchicalClustering",
  #                         "PCA")
  # setComparisonCutoff = "Auto"

  # Create output folder
  if (! dir.exists(outputPrefix) & ! returnDataOnly) {
    dir.create(outputPrefix)
  }

  # Initialize results object
  resList <- list()

  # Get available gene names
  avGenes <- getAvailableGenes(Species)

  # Check genes to make sure they exist -- only a warning
  intGenes <- genesOfInterest
  badGenes <- intGenes[which(! intGenes %in% avGenes$geneName &
                                                   ! intGenes %in% MSIGDB_Geneset_Names)]
  if (length(badGenes) > 0) {
    warning(paste0("\n\t\t\t'", paste(badGenes, collapse = ", "), "'
                      not found in correlation data and is not an official MSIGDB name.
                      Please check available gene data with getAvailableGenes().
                      Your gene(s) of interest may have an updated name or
                      have a species-specific identifier. Find offical MSIGDB
                      names by examining the MSIGDB_Geneset_Names object.\n
                      Continuing without this/these gene(s)..."))
  }

  # Make list of terms inputted by the user
  termlist <- intGenes[which(! intGenes %in% avGenes$geneName &
                               intGenes %in% MSIGDB_Geneset_Names)]

  if (length(termlist > 0)) {
    if (Species[1] == "hsapiens") {
      TERM2GENE <- hsapiens_complex_TERM2GENE
    } else {
      TERM2GENE <- mmusculus_complex_TERM2GENE
    }
    for(i in 1:length(termlist)) {
      term <- termlist[i]
      print(term)
      nameStr <- names(term)
      termGenes <- TERM2GENE$human_gene_symbol[which(TERM2GENE$gs_name == term)]
      termGenes <- termGenes[which(termGenes %in% avGenes$geneName)] # Ensure actionable genes
      intGenes <- unique(c(intGenes, termGenes)) # Append to intgenes vector
    }
    intGenes <- intGenes[which(! intGenes %in% termlist)]
  }




  cat("\nRetrieving correlation data...\n")
  # Call downloadData to get all required files
  timestamp()
  corrDF <- getCorrelationData(Species = Species[1],
                               Sample_Type = Sample_Type[1],
                               geneList = intGenes)
  timestamp()

  resList[["Correlation_Data"]] <- corrDF


  # Comparison Cutoff scaled to gene list size
  if (setComparisonCutoff == "Auto") {
    setComparisonCutoff <- round(length(genesOfInterest)/2.5)
  }


  # Cross comparison code block
  cat("\nStarting cross comparison\n")
  # prepare for PCA/variance analysis
  resultsMat <- corrDF
  rownames(resultsMat) <- rownames(corrDF)
  resultsMat <- as.matrix(resultsMat)
  # Remove duplicates
  resultsMat <- resultsMat[which(! duplicated(resultsMat)),]
  # Remove identity genes
  genesToRemove <- colnames(resultsMat)
  resultsMat <- resultsMat[which(! rownames(resultsMat) %in% genesToRemove),]


  height <- floor((length(intGenes) - 10)/8) + 10
  if (height > 20) {
    height <- 20
  } else if (height < 10) {
    height <- 10
  }
  width <- round(length(intGenes)/8)
  if (width < 5) {
    width <- 5
  }
  compList_up <- list()
  compList_dn <- list()
  # Get the top 10% of values within each dataset
  for ( i in 1:length(colnames(resultsMat))) {
    gene <- colnames(resultsMat)[i]
    vals <- resultsMat[,i]
    top <- vals[which(vals > quantile(vals, prob = .95))]
    bottom <- vals[which(vals < quantile(vals, prob = .05))]
    compList_up[[i]] <- names(top)
    names(compList_up)[i] <- gene
    compList_dn[[i]] <- names(bottom)
    names(compList_dn)[i] <- gene
  }
  # Get the variant gene results
  if ("variantGenes" %in% crossComparisonType) {
    # select the genes with top variance -- they explain differences between genes of interest
    rv <- metaMA::rowVars(resultsMat)
    select <- order(rv, decreasing=TRUE)[seq_len(min(1500, length(rv)))]
    topVarMat <- resultsMat[select,]
    topVarHeat <- pheatmap::pheatmap(topVarMat, color = gplots::greenred(100), show_rownames = F,
                                     main = "Variable Genes",silent = T,
                                     width = width, height = height)

    resList[["variantGenesHeatmap"]] <- topVarHeat
    resList[["variantGenesHeatmap_MAT"]] <- topVarMat

    varGenes <- rownames(topVarMat)
    # Perform pathway enrichment with Co-Correlative genes
    VarGenesEGMT <- clusterProfiler::enricher(gene = varGenes, TERM2GENE = TERM2GENE)
    eres <- as.data.frame(VarGenesEGMT)
    resList[["variantGenes_pathways"]] <- VarGenesEGMT

    # Modify gene set names to fit plotting window
    VarGenesEGMT@result$Description[which(nchar(VarGenesEGMT@result$Description) > 40)] <- paste0(substr(VarGenesEGMT@result$Description[which(nchar(VarGenesEGMT@result$Description) > 40)], 1, 40), "...")
    dp <- enrichplot::dotplot(VarGenesEGMT)
    dp <- dp + ggplot2::labs(title = "Variant Genes Pathway Enrichment")
    resList[["variantGenes_pathways_dotplot"]] <- dp

    if (numTopGenesToPlot == "Auto") {
      numTopGenesToPlot <- 50
    }

    select <- order(rv, decreasing=TRUE)[seq_len(min(numTopGenesToPlot, length(rv)))]
    topVarMat <- resultsMat[select,]
    topVarHeattop <- pheatmap::pheatmap(topVarMat, color = gplots::greenred(100), show_rownames = T,
             main = "Top Variable Genes", silent = T,
             width = width, height = height)

    resList[["variantGenesHeatmap_Top"]] <- topVarHeattop
    resList[["variantGenesHeatmap_Top_MAT"]] <- topVarMat

    if (! returnDataOnly) {
      png(filename = file.path(outputPrefix, "VarGeneHeatmap.png"),
          height = height, width = width, units = "in", res = 300)
      print(topVarHeat)
      dev.off()
      png(filename = file.path(outputPrefix, "VarGeneHeatmap_top.png"),
          height = height, width = width, units = "in", res = 300)
      print(topVarHeattop)
      dev.off()

      write.csv(eres,
                file = file.path(outputPrefix, "variantGenes.Pathway.Analysis.csv"),
                row.names = F)
      ggplot2::ggsave(plot = dp,
                      filename = file.path(outputPrefix, "variantGenes.Pathway.Analysis.png"),
                      height = 7, width = 10)

    }


  }


  if ("coCorrelativeGenes" %in% crossComparisonType) {
    # Select genes that best correlate together -- they explain similarities between genes of interest

    # Use exact test to find co-correlative genes
    ie_up <- SuperExactTest::intersectElements(compList_up)
    barCode_up <- strsplit(as.character(ie_up$barcode), "")
    df_up <- data.frame(matrix(unlist(barCode_up), nrow=length(barCode_up), byrow=T))
    rownames(df_up) <- ie_up$Entry
    df_up <- apply(df_up, 1:2, as.numeric)
    ie_dn <- SuperExactTest::intersectElements(compList_dn)
    barCode_dn <- strsplit(as.character(ie_dn$barcode), "")
    df_dn <- data.frame(matrix(unlist(barCode_dn), nrow=length(barCode_dn), byrow=T))
    rownames(df_dn) <- ie_dn$Entry
    df_dn <- apply(df_dn, 1:2, as.numeric)

    olGenes_determine <- function(df, setComparisonCutoff) {
      # Get genes with at least n overlaps, determined by setComparisonCutoff
      olGenes <- rownames(df)[which(rowSums(df) > setComparisonCutoff)]
      if (length(olGenes) < 10) {
        warning("Comparison cutoff value too high -- adjusted down for set analysis")
        olGenes <- rownames(df)[which(rowSums(df) > setComparisonCutoff-1)]
      }
      if (length(olGenes) < 10) {
        olGenes <- rownames(df)[which(rowSums(df) > setComparisonCutoff-2)]
      }
      if (length(olGenes) < 10) {
        olGenes <- rownames(df)[which(rowSums(df) > setComparisonCutoff-3)]
      }
      return(olGenes)
    }

    olGenes_up <- olGenes_determine(df = df_up, setComparisonCutoff = setComparisonCutoff)
    olGenes_dn <- olGenes_determine(df = df_dn, setComparisonCutoff = setComparisonCutoff)

    olGenes <- unique(c(olGenes_up, olGenes_dn))

    select <- which(rownames(resultsMat) %in% olGenes)

    olMat <- resultsMat[select,]

    cocorheatmap <- pheatmap::pheatmap(olMat, color = gplots::greenred(100), show_rownames = F,
             main = "Co-Correlative Genes", silent = T,
             width = width, height = height)

    resList[["cocorrelativeGenesHeatmap"]] <- cocorheatmap
    resList[["cocorrelativeGenesHeatmap_MAT"]] <- olMat

    # # Test with a smaller subset
    # olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-1)]
    # # If too few genes overlapping
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-2)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-3)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-4)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-5)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-6)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-8)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-10)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-14)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-18)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-24)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-32)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-42)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-50)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-60)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-75)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-95)]
    # }
    # if (length(olGenesSmall) < 10) {
    #   olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-120)]
    # }
    # # If too many genes overlapping
    # if (length(olGenesSmall) > numTopGenesToPlot) {
    #   select <- which(rownames(resultsMat) %in% olGenesSmall)
    #   olMat <- topGenes[select,]
    #   # Get the top rowSums
    #   olRows <- rowSums(olMat)
    #   olRows <- olRows[order(olRows, decreasing = T)]
    #   olGenesSmall <- names(olRows)[c(1:numTopGenesToPlot)]
    # }
    # select <- which(rownames(resultsMat) %in% olGenesSmall)
    # olMat <- resultsMat[select,]
    # pheatmap(olMat, color = greenred(100), show_rownames = T,
    #          main = "Heatmap of Genes of Interest by Top Co-Correlative Genes",
    #          filename = file.path(crossDir, "TopCoCorrelativeGenesHeatmap.png"),
    #          width = width, height = height)



    # Perform pathway enrichment with Co-Correlative genes
    CCGenesEGMT <- clusterProfiler::enricher(gene = olGenes, TERM2GENE = TERM2GENE)
    eres <- as.data.frame(CCGenesEGMT)
    resList[["coCorrelativeGenes_pathways"]] <- eres

    # Modify gene set names to fit plotting window
    CCGenesEGMT@result$Description[which(nchar(CCGenesEGMT@result$Description) > 40)] <- paste0(substr(CCGenesEGMT@result$Description[which(nchar(CCGenesEGMT@result$Description) > 40)], 1, 40), "...")
    dp <- enrichplot::dotplot(CCGenesEGMT)
    dp <- dp + ggplot2::labs(title = "Co-Correlated Genes Pathway Enrichment")
    resList[["coCorrelativeGenes_pathways_dotplot"]] <- dp

    if (! returnDataOnly) {
      write.csv(eres,
                file = file.path(outputPrefix, "coCorrelativeGenes.Pathway.Analysis.csv"),
                row.names = F)
      ggplot2::ggsave(plot = dp,
             filename = file.path(outputPrefix, "coCorrelativeGenes.Pathway.Analysis.png"),
             height = 7, width = 10)

      png(filename = file.path(outputPrefix, "coCorrelativeGeneHeatmap.png"),
          height = height, width = width, units = "in", res = 300)
      print(cocorheatmap)
      dev.off()

    }


  }
  # Set up for topological distance analysis
  if (numClusters == "Auto") {
    numClusters <- ceiling(length(intGenes) / 25)
  }
  if (numClusters > 20) {
    numClusters <- 20
  }
  if (numClusters < 2) {
    numClusters <- 2
  }
  # Begin topological distance analysis
  if ("PCA" %in% crossComparisonType & length(intGenes) <= 10) {
    # Principle component analysis
    pca <- prcomp(resultsMat)
    dd <- data.frame(summary(pca)$importance)
    percentVar <- as.numeric(round(100 * dd[2,]))
    percentVar <- percentVar[1:2]
    pcaData <- as.data.frame(pca$rotation)
    pcaData <- pcaData[,c(1:2)]
    pcaData$Gene <- rownames(pcaData)
    plt1 <- ggplot2::ggplot(pcaData, ggplot2::aes(PC1, PC2, color=Gene)) +
      ggplot2::geom_point(size = 5) +
      ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggplot2::labs(title = "PCA of Gene Correlation Values") + ggplot2::theme_classic()

    resList[["PCA_plot"]] <- plt1
    resList[["PCA_data"]] <- pcaData

    if (! returnDataOnly) {
      ggplot2::ggsave(plt1, file = file.path(outputPrefix, "geneCorrelationData.PCA.png"),
             height = 6, width = 7)
    }

  } else if ("PCA" %in% crossComparisonType & (length(intGenes) <= 100 || ! alternativeTSNE)) {
    # Principle component analysis without colors
    pca <- prcomp(resultsMat)
    dd <- data.frame(summary(pca)$importance)
    percentVar <- as.numeric(round(100 * dd[2,]))
    percentVar <- percentVar[1:2]
    pcaData <- as.data.frame(pca$rotation)
    pcaData <- pcaData[,c(1:2)]
    pcaData$Gene <- rownames(pcaData)
    lenPCA <- 5 + (length(genesOfInterest) - 5)/10
    pointSize <- 2
    if (lenPCA > 12) {
      lenPCA <- 12
      cexPCA <- 9/(3 + (.01 * length(genesOfInterest)))
      pointSize <- 1
      if (cexPCA < 1.4) {
        cexPCA <- 1.4
      }
    } else {
      cexPCA <- 3
      pointsize <- 4
    }
    plt1 <- ggplot2::ggplot(pcaData, ggplot2::aes(PC1, PC2)) +
      ggplot2::geom_point(size = pointSize, color = "red") +
      ggrepel::geom_text_repel(ggplot2::aes(label=Gene),
                      size = cexPCA, color = "black",
                      min.segment.length = 0.02, segment.alpha = 1,
                      box.padding = 0.1) +
      ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggplot2::labs(title = "PCA of Gene Correlation Values") + ggplot2::theme_classic()

    resList[["PCA_plot"]] <- plt1
    resList[["PCA_data"]] <- pcaData

    if (! returnDataOnly) {
      ggplot2::ggsave(plt1, file = file.path(outputPrefix, "geneCorrelationData.PCA.png"),
                      height = 6, width = 7)
    }


  } else if (length(intGenes) > 100 & alternativeTSNE) {
    cat("\nUsing TSNE instead of PCA for large sample sizes.
To disable this behavior, set 'alternativeTSNE' to FALSE")
    rv <- metaMA::rowVars(resultsMat)
    select <- order(rv, decreasing=TRUE)[seq_len(min(2500, length(rv)))]
    topVarMat <- resultsMat[select,]
    pca <- prcomp(topVarMat)
    # TSNE after PCA
    pcaData <- as.data.frame(pca$rotation)
    pcaData <- as.matrix(pcaData)
    tsne_out <- Rtsne::Rtsne(X = pcaData, pca = F,
                      verbose = T, max_iter = 5000, perplexity = 30, exaggeration_factor = 16)
    # Cluster by TSNE distance
    require(dplyr)
    hc.norm <- hclust(dist(tsne_out$Y))
    info.norm <- data.frame(geneNames = rownames(pcaData))
    info.norm <- info.norm %>% mutate(tsne1 = tsne_out$Y[, 1], tsne2 = tsne_out$Y[,2])
    info.norm$hclust <- factor(cutree(hc.norm, numClusters))
    hc.norm.cent <- info.norm %>% group_by(hclust) %>% select(tsne1,
                                                             tsne2) %>% summarize_all(mean)
    # Plot TSNE with clusters
    gp <- ggplot2::ggplot(info.norm, ggplot2::aes(x = tsne1, y = tsne2, colour = hclust)) +
      ggplot2::geom_point(alpha = 0.3) + ggplot2::theme_bw() +
      ggrepel::geom_label_repel(ggplot2::aes(label = hclust), data = hc.norm.cent) +
      ggplot2::guides(colour = FALSE) +
      ggplot2::ggtitle("Genes of Interest TSNE")

    resList[["TSNE_plot"]] <- gp
    resList[["TSNE_data"]] <- info.norm

    if (! returnDataOnly) {
      ggplot2::ggsave(plot = gp, filename = file.path(outputPrefix, "geneCorrelationData.TSNE.png"))
      # Save clustering info
      write.csv(info.norm, file = file.path(outputPrefix, "geneClusterData.TSNE.csv"), row.names = F)

    }

  } else {
    warning("Cannot perform PCA -- List too large. Please allow for TSNE to substitute.")
  }

  return(resList)

#   # Hierarchical clustering for distance-based analysis
#   if (! file.exists()) {
#     cat("\nCalculating distances between genes of interest.
# If more than 200 genes defined, this step may take a while.
# Set numDistGenes to estimate distance quicker.")
#     sampleDists <- dist(t(resultsMat))
#     save(sampleDists, file = file.path(outDir, "sampleDists.RData"))
#   } else {
#     load(file.path(outDir, "sampleDists.RData"))
#   }
#
#   sampleDistMatrix <- as.matrix(sampleDists)
#   colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#   # Distance heatmap
#   pheatmap(sampleDistMatrix,
#            clustering_distance_rows=sampleDists,
#            clustering_distance_cols=sampleDists,
#            col=colors,
#            main = "Heatmap of Gene-To-Gene Distances",
#            filename = file.path(crossDir, "geneDistances.png"),
#            width = width, height = width)
#   # Dendrogram
#   if (n)
#   colors <- c("red", "blue", "green", "black",
#              "darkgreen", "purple", "brown", "deeppink",
#              "darkorange2", "dimgray", "firebrick", "goldenrod",
#              "darkred", "darkturquoise", "indianred3", "navy",
#              "orangered", "peru", "tan", "yellowgreen")
#   colors <- sample(colors, numClusters)
#
#   clus <- cutree(hc, numClusters)
#   png(filename = file.path(crossDir, "geneDistanceTree.png"),
#       width = 8, height = 8, units = c("in"), res = 600)
#   plot(as.phylo(hc), type = "fan", tip.color = colors[clus],
#        label.offset = .15, cex = 0.22, main = "Distance plot for genes of interest")
#   dev.off()
#   ctDF <- as.data.frame(ct)
#   ctDF$geneName <- rownames(ctDF)
#   ctDF <- ctDF[,c(2,1)]
#   colnames(ctDF)[2] <- "Cluster"
#   write.csv(ctDF, file = file.path(crossDir, "geneClusterData.hClust.csv"), row.names = F)
#
#   # If there are only 1-2 genes, set analysis is not useful
#   if (! setAnalysis | length(compList) < 3) {
#     cat("\nDONE -- No set analysis performed.\n\nSet analysis is not informative with only 1 or 2 samples. \nSet 'setAnalysis' variable to TRUE and specify 3 or more 'genesOfInterest'.")
#     return()
#   } else {
#     cat("\nPerforming set analysis.\n")
#   }
#   # Set analysis
#   # If too many sets to compare -- exclude those which do not cluster together
#   setNum <- length(genesOfInterest)
#   if (length(setNum) >= 10 & is.null(genesToCompare)) {
#     clust <- kmeans(x = resultsMat, 3)
#     cm <- unique(which.max(rowMeans(clust$centers)))
#     clustSamps <- clust$centers
#     clustSamps <- colnames(clustSamps)[which(clustSamps[cm,] > 0)]
#     if (length(clustSamps) <= 10) {
#       resultsMat <- resultsMat[which(colnames(resultsMat) %in% clustSamps),]
#       cat("\nSome samples removed from comparative analysis via K-Means Clustering")
#       cat("\nSet comparisons should be limited to 10 samples. These can be manually defined\nby setting the genesToCompare variable.\n")
#     } else if (length(clustSamps) > 10) {
#       cat("\nSome samples removed from comparative analysis via K-Means Clustering")
#       cat("\nSet comparisons should be limited to 10 samples. These can be manually defined\nby setting the genesToCompare variable.\n")
#       stop("Still too many samples after clustering. Re-run with manually defined genesToCompare variable with 10 or less genes")
#     }
#   }
#
#   cat("\nStarting Exact Test for set comparison\n")
#   cat("\nIf you are attempting to compare the correlations of more than 5 genes, this step may take a while.")
#   total <- length(resultsFrame$geneName)
#   res <- SuperExactTest::supertest(x = compList, n = total)
#   # Relationship between samples
#   top <- length(compList)
#   bot <- round(length(compList)/2.5)
#   if (bot < 3) {
#     bot <- 2
#   }
#   png(file.path(crossDir, "coCorrelativeGenes.setAnalysis.Spiral.png"),
#       height = 9, width = 11, units = c("in"), res = 300)
#   plot.new()
#   title("Correlated genes (top 95%) between multiple sets")
#   plot(res, sort.by="size", margin=c(2,2,2,2),
#        color.scale.pos=c(0.15,0.02), legend.pos=c(.85,0.02),
#        max.intersection.size = 500, degree = c(bot:top), new.gridPage = F)
#   dev.off()
#   png(file.path(crossDir, "coCorrelativeGenes.setAnalysis.barChart.png"),
#       height = 9, width = 11, units = c("in"), res = 300)
#   plot.new()
#   title("Correlated genes (top 95%) between multiple sets")
#   plot(res, Layout="landscape",
#        sort.by="size", margin=c(0.5,5,1,2), degree = c(bot:top),
#        max.intersection.size = 500, yfrac = .6, show.overlap.size = F,
#        new.gridPage = F)
#   dev.off()
#   cat("\n\nDONE\n")
#   timestamp()
}


# # Packages required
require(clusterProfiler)
# require(ggplot2)
# require(pheatmap)
# require(ggrepel)
# require(gplots)
# require(ape)
# require(Rtsne)
# require(dplyr)
# require(RColorBrewer)
# require(metaMA)
# require(SuperExactTest)
#
# # Function main code
# GeneCorrelationAnalyzeR <- function(genesOfInterest,
#                                     pathToDataDir = file.path(getwd(), "Data/geneCorrelationDataFiles"),
#                                     pathToResultsDir = file.path(getwd(), "Results"),
#                                     outputPrefix = "GeneCorrelationAnalyzeR_Output",
#                                     runGSEA = T, complexGSEA = F, setAnalysis = T,
#                                     genesToCompare = NULL, setComparisonCutoff = "Auto",
#                                     alternativeTSNE = T, numClusters = "Auto",
#                                     numDistGenes = "Auto",
#                                     TERM2GENE = NULL,
#                                     crossComparisonType = c("variantGenes",
#                                                             "coCorrelativeGenes",
#                                                             "hierarchicalClustering",
#                                                             "PCA"),
#                                     numTopGenesToPlot = "Auto",
#                                     singleAnalysis = T) {
#
#   # # Test parameters for debugging
#   # # genesOfInterest <- c("CDK12", "AURKB", "SFPQ", "PARP1", "BRCC3")
#   # # # Check against RNAi
#   # # RNAi <- read.csv("Data/All.genes.clean.csv")
#   # # genesOfInterest <- as.character(RNAi$geneName)[c(1:600)]
#   # pathToDataDir = file.path(getwd(), "Data/geneCorrelationDataFiles")
#   # pathToResultsDir = file.path(getwd(), "Results")
#   # # outputPrefix = "RNAi_Genes_Top600"
#   # genesOfInterest <- load("Data/TC32MSC.Genes2check.RData")
#   # genesOfInterest <- genesToCheck
#   # outputPrefix = "genesTC32MSCDGE2"
#   # runGSEA = F
#   # complexGSEA = F
#   # setAnalysis = T
#   # singleAnalysis = F
#   # genesToCompare = NULL
#   # setComparisonCutoff = "Auto"
#   # numTopGenesToPlot = "Auto"
#   # crossComparisonType = c("variantGenes",
#   #                         "coCorrelativeGenes",
#   #                         "hierarchicalClustering",
#   #                         "PCA")
#
#   # Parse arguments
#   genesOfInterest <- unique(genesOfInterest)
#   geneString <- paste(genesOfInterest, collapse = ", ")
#   outDir <- file.path(pathToResultsDir, outputPrefix)
#   dir.create(outDir)
#   cat("\n\nInitializing Gene Correlation Analyzer \n\n")
#   timestamp()
#   # Load helper functions into memory
#   # Load human TERM2GENE file built from msigdbr
#   if (is.null(TERM2GENE)) {
#     if (complexGSEA) {
#       load(file.path(getwd(), "Data/GSEA_Data/hsapiens/TERM2GENE.RData"))
#     } else if (! complexGSEA) {
#       load(file.path(getwd(), "Data/GSEA_Data/hsapiens/TERM2GENE.SIMPLE.RData"))
#     }
#   }
#
#   # Comparison Cutoff scaled to gene list size
#   if (setComparisonCutoff == "Auto") {
#     setComparisonCutoff <- round(length(genesOfInterest)/2.5)
#   }
#   # Scale number of top genes to plot by size of gene list
#   if (numTopGenesToPlot == "Auto") {
#     if (length(genesOfInterest) <= 30) {
#       numTopGenesToPlot <- 30
#     } else {
#       numTopGenesToPlot <- floor((length(genesOfInterest) - 30)/3) + 30
#       if (numTopGenesToPlot > 80) {
#         numTopGenesToPlot <- 80
#       }
#     }
#   }
#   # Argument read-out for user
#   cat("\n\nArguments: \n")
#   cat("\n########################")
#   cat(paste0("\nGenes of Interest: ", geneString))
#   cat(paste0("\nData Directory: ", pathToDataDir))
#   cat(paste0("\nOutput Directory: ", outDir))
#   cat(paste0("\nPerform Single Gene Analyses: ", singleAnalysis))
#   cat(paste0("\nCross comparison method(s): ", paste(crossComparisonType, collapse = ", ")))
#   cat(paste0("\nComparison Cutoff Setting: ", setComparisonCutoff))
#   cat(paste0("\nNumber of Top Genes to Plot: ", numTopGenesToPlot))
#   cat(paste0("\nPerform Set Analysis: ", setAnalysis))
#   cat(paste0("\nrunGSEA: ", runGSEA))
#   if (runGSEA & ! complexGSEA) {
#     cat("\nGSEA Annotations: Simple")
#   } else if (runGSEA & complexGSEA) {
#     cat("\nGSEA Annotations: Complex")
#   }
#   cat("\n########################\n")
#   cat("\nBeginning Run\n")
#   # Main code
#   if (! file.exists(file.path(outDir, "geneCorrelationData.RData"))) {
#     load(file.path(pathToDataDir, "A1BG.RData"))
#     resultsFrame <- data.frame(geneName = names(corr))
#     for (i in 1:length(genesOfInterest)) {
#       gene <- genesOfInterest[i]
#       cat(paste0("\n", gene))
#       geneFile <- paste0(gene, ".RData")
#       file <- file.path(pathToDataDir, geneFile)
#       if(! file.exists(file)) {
#         cat ("\n", gene, " is not available -- please check to make sure gene name is correct.
#              To view list of availble genes, run 'getAvailableGenes()'")
#         next
#       }
#       er <- FALSE
#       er <- tryCatch(load(file), error = function(e){cat("ERROR: cannot load ", file); return(TRUE)})
#       if (er == TRUE) {
#         cat("Skipping to next gene")
#         next
#       }
#       # Build dataframe object from correlation data
#       corrdf <- as.data.frame(corr)
#       colnames(corrdf)[1] <- gene
#       corrdf$geneName <- rownames(corrdf)
#       # Add correlation data to final dataframe
#       resultsFrame <- merge(x = resultsFrame, y = corrdf, by = "geneName")
#       if (singleAnalysis) {
#         singleAnalysisDir <- file.path(outDir, "singleGeneAnalysis")
#         dir.create(singleAnalysisDir)
#         geneOutDir <- file.path(outDir, "singleGeneAnalysis", gene)
#         dir.create(geneOutDir)
#
#         write.csv(corrdf, file = file.path(geneOutDir, paste0(gene, ".csv")), row.names = F)
#         # Remove gene from correlation values to build normal distribution
#         corr <- corr[order(corr, decreasing = T)]
#         corr <- corr[c(-1)]
#         # Make a histogram of REST gene correlations
#         png(file.path(geneOutDir, paste0(gene, ".png")))
#         hist(corr, main = paste0(gene, " gene correlations"), breaks = 100)
#         dev.off()
#         if (runGSEA) {
#           # Perform GSEA
#           cat("\nStarting GSEA\n")
#           resGSEA <- myGSEA(ranks = corr, TERM2GENE = TERM2GENE,
#                             plotName = paste0(gene, ".corrPathways"), outDir = geneOutDir,
#                             Condition = paste0(gene, " Correlated Genes"))
#           # # Perform pathway analysis with these data
#           # cat("\n\nPlotting Results\n\n")
#           # ep <- emapplot(resGSEA, color = "NES")
#           # ep <- ep + labs(title = paste0(gene, " Correlated Genes GSEA Plot"))
#           # ggsave(plot = ep, filename = paste0(geneOutDir, "/", gene, ".EMAP.png"), height = 11, width = 11)
#         }
#       }
#     }
#     # Many genes may have been lost -- we may need to redo our genesOfInterest
#     if (length(genesOfInterest) != (length(resultsFrame)-1)) {
#       cat("\nDue to missing gene data -- reconfiguring run info")
#       genesOfInterest <- colnames(resultsFrame)[c(-1)]
#       cat("\nNumber of genes remaining for analysis: ", length(genesOfInterest))
#       # Comparison Cutoff scaled to gene list size
#       setComparisonCutoff <- round(length(genesOfInterest)/2.5)
#       # Scale number of top genes to plot by size of gene list
#       if (length(genesOfInterest) <= 30) {
#         numTopGenesToPlot <- 30
#       } else {
#         numTopGenesToPlot <- floor((length(genesOfInterest) - 30)/3) + 30
#         if (numTopGenesToPlot > 80) {
#           numTopGenesToPlot <- 80
#         }
#       }
#       # Argument read-out for user
#       cat("\n\nArguments: \n")
#       cat("\n########################")
#       cat(paste0("\nGenes of Interest: ", geneString))
#       cat(paste0("\nData Directory: ", pathToDataDir))
#       cat(paste0("\nOutput Directory: ", outDir))
#       cat(paste0("\nPerform Single Gene Analyses: ", singleAnalysis))
#       cat(paste0("\nCross comparison method(s): ", paste(crossComparisonType, collapse = ", ")))
#       cat(paste0("\nComparison Cutoff Setting: ", setComparisonCutoff))
#       cat(paste0("\nNumber of Top Genes to Plot: ", numTopGenesToPlot))
#       cat(paste0("\nPerform Set Analysis: ", setAnalysis))
#       cat(paste0("\nrunGSEA: ", runGSEA))
#       if (runGSEA & ! complexGSEA) {
#         cat("\nGSEA Annotations: Simple")
#       } else if (runGSEA & complexGSEA) {
#         cat("\nGSEA Annotations: Complex")
#       }
#       cat("\n########################\n")
#       cat("\nResumming Run\n")
#     }
#     if (length(genesOfInterest) < 25) {
#       write.csv(resultsFrame, file = file.path(outDir, "geneCorrelationData.csv"), row.names = F)
#       save(resultsFrame, file = file.path(outDir, "geneCorrelationData.RData"))
#     } else {
#       save(resultsFrame, file = file.path(outDir, "geneCorrelationData.RData"))
#     }
#     } else {
#       load(file.path(outDir, "geneCorrelationData.RData"))
#     }
#
#   # Cross comparison code block
#   cat("\nStarting cross comparison\n")
#   crossDir <- file.path(outDir, "crossComparison")
#   dir.create(crossDir)
#   # prepare for PCA/variance analysis
#   resultsMat <- resultsFrame[,c(-1)]
#   rownames(resultsMat) <- resultsFrame[,c(1)]
#   resultsMat <- as.matrix(resultsMat)
#   # Remove duplicates
#   resultsMat <- resultsMat[which(! duplicated(resultsMat)),]
#   # Remove identity genes
#   genesToRemove <- colnames(resultsMat)
#   resultsMat <- resultsMat[which(! rownames(resultsMat) %in% genesToRemove),]
#   id <- which(apply(resultsMat, 1, function (x) all(abs(x) <= .9)))
#   resultsMat <- resultsMat[id, ]
#   height <- floor((length(genesOfInterest) - 10)/8) + 10
#   if (height > 20) {
#     height <- 20
#   } else if (height < 10) {
#     height <- 10
#   }
#   width <- round(length(genesOfInterest)/8)
#   if (width < 5) {
#     width <- 5
#   }
#   topGenes <- abs(resultsMat)
#   compList <- list()
#   for ( i in 1:length(colnames(topGenes))) {
#     gene <- colnames(topGenes)[i]
#     vals <- topGenes[,i]
#     top <- vals[vals > quantile(vals, prob = .90)]
#     compList[[i]] <- names(top)
#     names(compList)[i] <- gene
#   }
#   if ("variantGenes" %in% crossComparisonType) {
#     # select the genes with top variance -- they explain differences between genes of interest
#     rv <- rowVars(resultsMat)
#     select <- order(rv, decreasing=TRUE)[seq_len(min(1500, length(rv)))]
#     topVarMat <- resultsMat[select,]
#     pheatmap(topVarMat, color = redgreen(100), show_rownames = F,
#              main = "Heatmap of Genes of Interest by Variable Genes",
#              filename = file.path(crossDir, "varGeneHeatmap.png"),
#              width = width, height = height)
#     select <- order(rv, decreasing=TRUE)[seq_len(min(numTopGenesToPlot, length(rv)))]
#     topVarMat <- resultsMat[select,]
#     pheatmap(topVarMat, color = redgreen(100), show_rownames = T,
#              main = "Heatmap of Genes of Interest by Top Variable Genes",
#              filename = file.path(crossDir, "topVarGeneHeatmap.png"),
#              width = width, height = height)
#   }
#   if ("coCorrelativeGenes" %in% crossComparisonType) {
#     # Select genes that best correlate together -- they explain similarities between genes of interest
#
#     # Use exact test to find co-correlative genes
#     ie <- SuperExactTest::intersectElements(compList)
#     barCode <- strsplit(as.character(ie$barcode), "")
#     df <- data.frame(matrix(unlist(barCode), nrow=length(barCode), byrow=T))
#     rownames(df) <- ie$Entry
#     df <- apply(df, 1:2, as.numeric)
#     # Get genes with at least n overlaps, determined by setComparisonCutoff
#     olGenes <- rownames(df)[which(rowSums(df) > setComparisonCutoff)]
#     if (length(olGenes) < 10) {
#       warning("Comparison cutoff value too high -- adjusted down for set analysis")
#       olGenes <- rownames(df)[which(rowSums(df) > setComparisonCutoff-1)]
#     }
#     if (length(olGenes) < 10) {
#       olGenes <- rownames(df)[which(rowSums(df) > setComparisonCutoff-2)]
#     }
#     if (length(olGenes) < 10) {
#       olGenes <- rownames(df)[which(rowSums(df) > setComparisonCutoff-3)]
#     }
#     select <- which(rownames(resultsMat) %in% olGenes)
#     olMat <- resultsMat[select,]
#     pheatmap(olMat, color = redgreen(100), show_rownames = F,
#              main = "Heatmap of Genes of Interest by Co-Correlative Genes",
#              filename = file.path(crossDir, "CoCorrelativeGenesHeatmap.png"),
#              width = width, height = height)
#     # Test with a smaller subset
#     olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-1)]
#     # If too few genes overlapping
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-2)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-3)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-4)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-5)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-6)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-8)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-10)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-14)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-18)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-24)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-32)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-42)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-50)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-60)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-75)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-95)]
#     }
#     if (length(olGenesSmall) < 10) {
#       olGenesSmall <- rownames(df)[which(rowSums(df) > length(compList)-120)]
#     }
#     # If too many genes overlapping
#     if (length(olGenesSmall) > numTopGenesToPlot) {
#       select <- which(rownames(resultsMat) %in% olGenesSmall)
#       olMat <- topGenes[select,]
#       # Get the top rowSums
#       olRows <- rowSums(olMat)
#       olRows <- olRows[order(olRows, decreasing = T)]
#       olGenesSmall <- names(olRows)[c(1:numTopGenesToPlot)]
#     }
#     select <- which(rownames(resultsMat) %in% olGenesSmall)
#     olMat <- resultsMat[select,]
#     pheatmap(olMat, color = redgreen(100), show_rownames = T,
#              main = "Heatmap of Genes of Interest by Top Co-Correlative Genes",
#              filename = file.path(crossDir, "TopCoCorrelativeGenesHeatmap.png"),
#              width = width, height = height)
#     # Perform pathway enrichment with Co-Correlative genes
#     CCGenesEGMT <- enricher(gene = olGenes, TERM2GENE = TERM2GENE)
#     eres <- as.data.frame(CCGenesEGMT)
#     write.csv(eres,
#               file = file.path(crossDir, "coCorrelativeGenes.Pathway.Analysis.csv"),
#               row.names = F)
#     # Modify gene set names to fit plotting window
#     CCGenesEGMT@result$Description[which(nchar(CCGenesEGMT@result$Description) > 40)] <- paste0(substr(CCGenesEGMT@result$Description[which(nchar(CCGenesEGMT@result$Description) > 40)], 1, 40), "...")
#     dp <- dotplot(CCGenesEGMT)
#     dp <- dp + labs(title = "Top Co-Correlated Genes Pathway Enrichment")
#     ggsave(plot = dp,
#            filename = file.path(crossDir, "coCorrelativeGenes.Pathway.Analysis.png"),
#            height = 7, width = 10)
#   }
#   # Set up for topological distance analysis
#   if (numClusters == "Auto") {
#     numClusters <- ceiling(length(genesOfInterest) / 25)
#   }
#   if (numClusters > 20) {
#     numClusters <- 20
#   }
#   if (numClusters < 2) {
#     numClusters <- 2
#   }
#   # Begin topological distance analysis
#   lenPCA <- 5 + (length(genesOfInterest) - 5)/10
#   pointSize <- 2
#   cexPCA <- 3
#   if (lenPCA > 12) {
#     lenPCA <- 12
#     cexPCA <- 9/(3 + (.01 * length(genesOfInterest)))
#     pointSize <- 1
#     if (cexPCA < 1.4) {
#       cexPCA <- 1.4
#     }
#   }
#   if ("PCA" %in% crossComparisonType & length(genesOfInterest) <= 10) {
#     # Principle component analysis
#     pca <- prcomp(resultsMat)
#     dd <- data.frame(summary(pca)$importance)
#     percentVar <- as.numeric(round(100 * dd[2,]))
#     percentVar <- percentVar[1:2]
#     pcaData <- as.data.frame(pca$rotation)
#     pcaData <- pcaData[,c(1:2)]
#     pcaData$Gene <- rownames(pcaData)
#     plt1 <- ggplot(pcaData, aes(PC1, PC2, color=Gene)) +
#       geom_point(size = 3) +
#       xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#       ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#       labs(title = "PCA of Gene Correlation Values") + theme_classic()
#     ggsave(plt1, file = file.path(crossDir, "geneCorrelationData.PCA.png"),
#            height = lenPCA, width = lenPCA)
#   } else if ("PCA" %in% crossComparisonType) {
#     # Principle component analysis without colors
#     pca <- prcomp(resultsMat)
#     dd <- data.frame(summary(pca)$importance)
#     percentVar <- as.numeric(round(100 * dd[2,]))
#     percentVar <- percentVar[1:2]
#     pcaData <- as.data.frame(pca$rotation)
#     pcaData <- pcaData[,c(1:2)]
#     pcaData$Gene <- rownames(pcaData)
#
#     plt1 <- ggplot(pcaData, aes(PC1, PC2)) +
#       geom_point(size = pointSize, color = "red") +
#       geom_text_repel(aes(label=Gene),
#                       size = cexPCA, color = "black",
#                       min.segment.length = 0.01, segment.alpha = 1,
#                       box.padding = 0) +
#       xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#       ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#       labs(title = "PCA of Gene Correlation Values") + theme_classic()
#     ggsave(plt1, file = file.path(crossDir, "geneCorrelationData.PCA.png"),
#            height = lenPCA, width = lenPCA)
#     ggsave(plt1, file = file.path(crossDir, "geneCorrelationData.PCA.pdf"),
#            height = lenPCA, width = lenPCA, device = "pdf")
#     if (length(genesOfInterest) > 100) {
#       cat("\nUsing TSNE in addition to PCA for large sample sizes.
#           To disable this behavior, set 'alternativeTSNE' to FALSE")
#       # TSNE after PCA
#       pcaData <- as.data.frame(pca$rotation)
#       pcaData <- as.matrix(pcaData)
#       tsne_out <- Rtsne(X = pcaData, pca = F,
#                         verbose = T, max_iter = 5000, perplexity = 30, exaggeration_factor = 16)
#       # Cluster by TSNE distance
#       hc.norm = hclust(dist(tsne_out$Y))
#       info.norm <- data.frame(geneNames = rownames(pcaData))
#       info.norm %<>% mutate(tsne1 = tsne_out$Y[, 1], tsne2 = tsne_out$Y[,2])
#       info.norm$hclust = factor(cutree(hc.norm, numClusters))
#       hc.norm.cent = info.norm %>% group_by(hclust) %>% select(tsne1,
#                                                                tsne2) %>% summarize_all(mean)
#       # Plot TSNE with clusters
#       gp <- ggplot(info.norm, aes(x = tsne1, y = tsne2, colour = hclust)) +
#         geom_point(alpha = 0.3) + theme_bw() +
#         geom_label_repel(aes(label = hclust), data = hc.norm.cent) + guides(colour = FALSE) +
#         ggtitle("Genes of Interest TSNE")
#       ggsave(plot = gp, filename = file.path(crossDir, "geneCorrelationData.TSNE.png"))
#       # Save clustering info
#       write.csv(info.norm, file = file.path(crossDir, "geneClusterData.TSNE.csv"), row.names = F)
#
#     }
#   }
#
#   # Hierarchical clustering for distance-based analysis
#   if (! file.exists(file.path(outDir, "sampleDists.RData"))) {
#     cat("\nCalculating distances between genes of interest.
#         If more than 200 genes defined, this step may take a while.
#         Set numDistGenes to estimate distance quicker.")
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
#   colors <- c("red", "blue", "green", "black",
#               "darkgreen", "purple", "brown", "deeppink",
#               "darkorange2", "dimgray", "firebrick", "goldenrod",
#               "darkred", "darkturquoise", "indianred3", "navy",
#               "orangered", "peru", "tan", "yellowgreen")
#   colors <- sample(colors, numClusters)
#   hc <- hclust(sampleDists)
#   clus <- cutree(hc, numClusters)
#   png(filename = file.path(crossDir, "geneDistanceTree.png"),
#       width = 8, height = 8, units = c("in"), res = 600)
#   plot(as.phylo(hc), type = "fan", tip.color = colors[clus],
#        label.offset = .15, cex = 0.22, main = "Distance plot for genes of interest")
#   dev.off()
#   ct <- clus
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
#   if (length(setNum) >= 6 & is.null(genesToCompare)) {
#     clust <- kmeans(x = resultsMat, 3)
#     cm <- unique(which.max(rowMeans(clust$centers)))
#     clustSamps <- clust$centers
#     clustSamps <- colnames(clustSamps)[which(clustSamps[cm,] > 0)]
#     if (length(clustSamps) <= 6) {
#       resultsMat <- resultsMat[which(colnames(resultsMat) %in% clustSamps),]
#       cat("\nSome samples removed from comparative analysis via K-Means Clustering")
#       cat("\nSet comparisons should be limited to 6 samples. These can be manually defined\nby setting the genesToCompare variable.\n")
#     } else if (length(clustSamps) > 6) {
#       cat("\nSome samples removed from comparative analysis via K-Means Clustering")
#       cat("\nSet comparisons should be limited to 6 samples. These can be manually defined\nby setting the genesToCompare variable.\n")
#       stop("Still too many samples after clustering. Re-run with manually defined genesToCompare variable with 6 or less genes")
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
# }
#

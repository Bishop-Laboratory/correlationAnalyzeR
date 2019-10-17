#' Analyze Gene List Topology
#'
#' Analyzes the topology of a gene list using gene correlation data and dimension-reduction techniques.
#'
#' @param genesOfInterest A vector of genes to analyze or the name of an official MSIGDB term.
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'
#'
#' @param pathwayType Whether pathway enrichment should consider all msigdb annotations
#'     or just those in the most popular categories. Should be one of either
#'     'simple' or 'complex'.
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
#' @param returnDataOnly if TRUE will return only a list of analysis results. Default: FALSE.
#'
#' @param setComparisonCutoff Only relevant for co-correlation analysis -- the number of genes which
#' must aggree for a gene to be considered co-correlative within the input gene list.
#'
#' @param numTopGenesToPlot When creating a heatmap of the top co-correlative or top variant genes,
#' how many genes should be plotted on the y axis? Default: "Auto"
#'
#' @param numClusters The number of clusters to create with hclust or TSNE analysis. Default: "Auto"
#'
#' @param alternativeTSNE Logical. If TRUE, then a TSNE will be run as an alternative to PCA for visualizing
#' large input gene lists. This is highly recommended as 100+ member gene lists cannot be visualized otherwise.
#'
#' @param pathwayEnrichment Logic. If TRUE, pathway enrichment will be performed on variant genes --
#' if 'variantGenes' selected -- and/or on co-correlative genes -- if "coCorrelativeGenes" selected.
#'
#' @param pValueCutoff Numeric. The p value cutoff applied when running all pathway enrichment tests.
#'
#' @param crossComparisonType The type of topology tests to run. (see details)
#'
#' @return A list of correlations for input genes, and the results of chosen analysis + visualizations.
#'
#' @examples
#' genesOfInterest <- c("CDK12", "AURKB", "SFPQ", "NFKB1", "BRCC3", "BRCA2", "PARP1",
#'                      "DHX9", "SON", "AURKA", "SETX", "BRCA1", "ATMIN")
#' correlationAnalyzeR::analyzeGenesetTopology(genesOfInterest = genesOfInterest,
#'                                  Species = "hsapiens",
#'                                  Sample_Type = "cancer", returnDataOnly = TRUE,
#'                                  Tissue = "brain",
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
#' - pathwayEnrich: Cluster profiler's enricher function will be run on the input gene list.
#' @importFrom rlang .data
#' @import dplyr
#' @export
analyzeGenesetTopology <-  function(genesOfInterest,
                                    Species = c("hsapiens", "mmusculus"),
                                    Sample_Type = c("normal", "cancer"),
                                    Tissue = "all",
                                    crossComparisonType = c("PCA",
                                                            "variantGenes",
                                                            "coCorrelativeGenes",
                                                            "pathwayEnrich"),
                                    pathwayType = c("simple", "complex"),
                                    setComparisonCutoff = "Auto",
                                    pathwayEnrichment = FALSE,
                                    pValueCutoff = .05,
                                    numTopGenesToPlot = "Auto",
                                    alternativeTSNE = TRUE,
                                    numClusters = "Auto",
                                    outputPrefix = "CorrelationAnalyzeR_Output",
                                    returnDataOnly = FALSE) {


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
  badGenes <- intGenes[
    which(! intGenes %in% avGenes &
            ! intGenes %in% correlationAnalyzeR::MSIGDB_Geneset_Names)
    ]
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
  termlist <- intGenes[which(! intGenes %in% avGenes &
                               intGenes %in% correlationAnalyzeR::MSIGDB_Geneset_Names)]

  # Load appropriate TERM2GENE file built from msigdbr()
  if (Species[1] %in% c("hsapiens", "mmusculus")) {
    TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = pathwayType,
                                                   Species = Species)
  } else {
    stop("\ncorrelationAnalyzeR currently supports only Human and Mouse data.
         Please select either 'hsapiens' or 'mmusculus' for Species parameter.
         \n")
  }

  if (length(termlist > 0)) {

    for(i in 1:length(termlist)) {
      term <- termlist[i]
      print(term)
      nameStr <- names(term)
      termGenes <- TERM2GENE$gene_symbol[which(TERM2GENE$gs_name == term)]
      termGenes <- termGenes[which(termGenes %in% avGenes)] # Ensure actionable genes
      intGenes <- unique(c(intGenes, termGenes)) # Append to intgenes vector
    }
    intGenes <- intGenes[which(! intGenes %in% termlist)]
  }

  cat("\nRetrieving correlation data...\n")
  # Call downloadData to get all required files
  # timestamp()
  corrDF <- correlationAnalyzeR::getCorrelationData(Species = Species,
                                                    Tissue = Tissue,
                                                    Sample_Type = Sample_Type,
                                                    geneList = genesOfInterest)
  # timestamp()

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
    top <- vals[which(vals > stats::quantile(vals, prob = .95))]
    bottom <- vals[which(vals < stats::quantile(vals, prob = .05))]
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
    topVarHeat <- pheatmap::pheatmap(topVarMat, color = gplots::greenred(100), show_rownames = FALSE,
                                     main = "Variable Genes",silent = TRUE,
                                     width = width, height = height)

    resList[["variantGenesHeatmap"]] <- topVarHeat
    resList[["variantGenesHeatmap_MAT"]] <- topVarMat

    varGenes <- rownames(topVarMat)
    # Perform pathway enrichment with Co-Correlative genes
    if (pathwayEnrichment) {
      VarGenesEGMT <- clusterProfiler::enricher(gene = varGenes, TERM2GENE = TERM2GENE,
                                                pvalueCutoff = pValueCutoff)
      eres <- as.data.frame(VarGenesEGMT)
      resList[["variantGenes_pathways"]] <- VarGenesEGMT
      # Modify gene set names to fit plotting window
      VarGenesEGMT@result$Description[which(nchar(VarGenesEGMT@result$Description) > 40)] <- paste0(substr(VarGenesEGMT@result$Description[which(nchar(VarGenesEGMT@result$Description) > 40)], 1, 40), "...")
      dp <- clusterProfiler::dotplot(VarGenesEGMT)
      dp <- dp + ggplot2::labs(title = "Variant Genes Pathway Enrichment")
      resList[["variantGenes_pathways_dotplot"]] <- dp
    }

    if (numTopGenesToPlot == "Auto") {
      numTopGenesToPlot <- 50
    }

    select <- order(rv, decreasing=TRUE)[seq_len(min(numTopGenesToPlot, length(rv)))]
    topVarMat <- resultsMat[select,]
    topVarHeattop <- pheatmap::pheatmap(topVarMat, color = gplots::greenred(100), show_rownames = TRUE,
             main = "Top Variable Genes", silent = TRUE,
             width = width, height = height)

    resList[["variantGenesHeatmap_Top"]] <- topVarHeattop
    resList[["variantGenesHeatmap_Top_MAT"]] <- topVarMat

    if (! returnDataOnly) {
      grDevices::png(filename = file.path(outputPrefix, "VarGeneHeatmap.png"),
          height = height, width = width, units = "in", res = 300)
      print(topVarHeat)
      grDevices::dev.off()
      grDevices::png(filename = file.path(outputPrefix, "VarGeneHeatmap_top.png"),
          height = height, width = width, units = "in", res = 300)
      print(topVarHeattop)
      grDevices::dev.off()
      if (pathwayEnrichment) {
        utils::write.csv(eres,
                  file = file.path(outputPrefix, "variantGenes.Pathway.Analysis.csv"),
                  row.names = FALSE)
        ggplot2::ggsave(plot = dp,
                        filename = file.path(outputPrefix, "variantGenes.Pathway.Analysis.png"),
                        height = 7, width = 10)
      }


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

    cocorheatmap <- pheatmap::pheatmap(olMat, color = gplots::greenred(100),
                                       show_rownames = FALSE,
             main = "Co-Correlative Genes", silent = TRUE,
             width = width, height = height)

    resList[["cocorrelativeGenesHeatmap"]] <- cocorheatmap
    resList[["cocorrelativeGenesHeatmap_MAT"]] <- olMat
    if (pathwayEnrichment) {
      # Perform pathway enrichment with Co-Correlative genes
      CCGenesEGMT <- clusterProfiler::enricher(gene = olGenes, TERM2GENE = TERM2GENE,
                                               pvalueCutoff = pValueCutoff)
      eres <- as.data.frame(CCGenesEGMT)
      resList[["coCorrelativeGenes_pathways"]] <- eres
      # Modify gene set names to fit plotting window
      CCGenesEGMT@result$Description[
        which(nchar(CCGenesEGMT@result$Description) > 40)
        ] <- paste0(substr(CCGenesEGMT@result$Description[
          which(nchar(CCGenesEGMT@result$Description) > 40)
          ], 1, 40), "...")
      dp <- clusterProfiler::dotplot(CCGenesEGMT)
      dp <- dp + ggplot2::labs(title = "Co-Correlated Genes Pathway Enrichment")
      resList[["coCorrelativeGenes_pathways_dotplot"]] <- dp
    }
    if (! returnDataOnly) {
      if (pathwayEnrichment) {
        utils::write.csv(eres,
                  file = file.path(outputPrefix, "coCorrelativeGenes.Pathway.Analysis.csv"),
                  row.names = FALSE)
        ggplot2::ggsave(plot = dp,
                        filename = file.path(outputPrefix, "coCorrelativeGenes.Pathway.Analysis.png"),
                        height = 7, width = 10)
      }
      grDevices::png(filename = file.path(outputPrefix, "coCorrelativeGeneHeatmap.png"),
          height = height, width = width, units = "in", res = 300)
      print(cocorheatmap)
      grDevices::dev.off()
    }
  }
  # Set up for topological distance analysis
  if (numClusters == "Auto") {
    numClusters <- ceiling(length(intGenes) / 10)
  }
  if (numClusters > 15) {
    numClusters <- 15
  }
  if (numClusters < 3) {
    numClusters <- 3
  }
  # Begin topological distance analysis
  if ("PCA" %in% crossComparisonType & length(intGenes) <= 10) {
    # Principle component analysis
    pca <- stats::prcomp(resultsMat)
    dd <- data.frame(summary(pca)$importance)
    percentVar <- as.numeric(round(100 * dd[2,]))
    percentVar <- percentVar[1:2]
    pcaData <- as.data.frame(pca$rotation)
    pcaData <- pcaData[,c(1:2)]
    pcaData$Gene <- rownames(pcaData)
    plt1 <- ggplot2::ggplot(pcaData, ggplot2::aes_string("PC1", "PC2", color="Gene")) +
      ggplot2::geom_point(size = 5) +
      ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggplot2::labs(title = "PCA of Gene Correlation Values") + ggplot2::theme_classic()

    resList[["PCA_plot"]] <- plt1
    resList[["PCA_data"]] <- pcaData
    resList[["clustered"]] <- F

    if (! returnDataOnly) {
      ggplot2::ggsave(plt1, file = file.path(outputPrefix, "geneCorrelationData.PCA.png"),
             height = 6, width = 7)
    }

  } else if ("PCA" %in% crossComparisonType & (length(intGenes) <= 100 || ! alternativeTSNE)) {
    # Principle component analysis without colors
    pca <- stats::prcomp(resultsMat)
    dd <- data.frame(summary(pca)$importance)
    percentVar <- as.numeric(round(100 * dd[2,]))
    percentVar <- percentVar[1:2]
    pcaData <- as.data.frame(pca$rotation)
    pcaDataSmall <- pcaData[,c(1:2)]
    pcaData$Gene <- rownames(pcaDataSmall)

    # HClust

    sdat <- t(scale(t(resultsMat)))
    pr.dis <- stats::dist(t(sdat), method = "euclidean")
    hc.norm <- stats::hclust(pr.dis)

    info.norm <- data.frame(geneNames = rownames(pca$rotation))
    info.norm <- info.norm %>% mutate(PC1 = pca$rotation[, 1], PC2 = pca$rotation[,2])
    info.norm$clusters <- factor(stats::cutree(hc.norm, numClusters))


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
    plt1 <- ggplot2::ggplot(info.norm, ggplot2::aes_string("PC1", "PC2")) +
      ggplot2::geom_point(size = pointSize, color = info.norm$clusters) +
      ggrepel::geom_text_repel(ggplot2::aes_string(label="geneNames"),
                      size = cexPCA, color = "black",
                      min.segment.length = 0.02, segment.alpha = 1,
                      box.padding = 0.1) +
      ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggplot2::labs(title = "PCA of Gene Correlation Values") + ggplot2::theme_classic()

    resList[["PCA_plot"]] <- plt1
    resList[["PCA_data"]] <- info.norm
    resList[["clustered"]] <- T

    if (! returnDataOnly) {
      ggplot2::ggsave(plt1, file = file.path(outputPrefix, "geneCorrelationData.PCA.png"),
                      height = 6, width = 7)
    }


  } else if (length(intGenes) > 100 & alternativeTSNE & "PCA" %in% crossComparisonType) {
    cat("\nUsing TSNE instead of PCA for large sample sizes.
To disable this behavior, set 'alternativeTSNE' to FALSE")
    rv <- metaMA::rowVars(resultsMat)
    select <- order(rv, decreasing=TRUE)[seq_len(min(2500, length(rv)))]
    topVarMat <- resultsMat[select,]
    pca <- stats::prcomp(topVarMat)
    # TSNE after PCA
    pcaData <- as.data.frame(pca$rotation)
    pcaData <- as.matrix(pcaData)
    # library(umap)
    # umap::umap(pcaData)
    tsne_out <- Rtsne::Rtsne(X = pcaData, pca = FALSE,
                      verbose = TRUE, max_iter = 5000, perplexity = 30, exaggeration_factor = 16)
    # Cluster by TSNE distance
    hc.norm <- stats::hclust(stats::dist(tsne_out$Y))
    info.norm <- data.frame(geneNames = rownames(pcaData))
    info.norm <- info.norm %>% mutate(tsne1 = tsne_out$Y[, 1], tsne2 = tsne_out$Y[,2])
    info.norm$hclust <- factor(stats::cutree(hc.norm, numClusters))
    hc.norm.cent <- info.norm %>% group_by(.data$hclust) %>% select(.data$tsne1,
                                                             .data$tsne2) %>% summarize_all(mean)
    # Plot TSNE with clusters
    gp <- ggplot2::ggplot(info.norm, ggplot2::aes_string(x = "tsne1", y = "tsne2", colour = "hclust")) +
      ggplot2::geom_point(alpha = 0.3) + ggplot2::theme_bw() +
      ggrepel::geom_label_repel(ggplot2::aes_string(label = "hclust"), data = hc.norm.cent) +
      ggplot2::guides(colour = FALSE) +
      ggplot2::ggtitle("Genes of Interest TSNE")

    resList[["TSNE_plot"]] <- gp
    resList[["TSNE_data"]] <- info.norm

    if (! returnDataOnly) {
      ggplot2::ggsave(plot = gp, filename = file.path(outputPrefix, "geneCorrelationData.TSNE.png"))
      # Save clustering info
      utils::write.csv(info.norm, file = file.path(outputPrefix, "geneClusterData.TSNE.csv"), row.names = FALSE)

    }

  } else if ("PCA" %in% crossComparisonType) {
    warning("Cannot perform PCA -- List too large. Please allow for TSNE to substitute.")
  }

  if ("pathwayEnrich" %in% crossComparisonType) {
    if (length(intGenes) < 10) {
      warning("\nPathway enrichment is recommended with at least 10 genes, otherwise results may not be informative.\n")
    }
    EGMT <- clusterProfiler::enricher(gene = intGenes, TERM2GENE = TERM2GENE,
                                      pvalueCutoff = pValueCutoff)
    eres <- as.data.frame(EGMT)
    resList[["inputGenes_pathwayEnrich"]] <- EGMT
    resList[["inputGenes_pathwayEnrich_data"]] <- eres
    # Modify gene set names to fit plotting window
    EGMT@result$Description[which(nchar(EGMT@result$Description) > 40)] <- paste0(substr(EGMT@result$Description[which(nchar(EGMT@result$Description) > 40)], 1, 40), "...")
    dp <- enrichplot::dotplot(EGMT)
    dp <- dp + ggplot2::labs(title = "Input Genes Pathway Enrichment")
    resList[["inputGenes_pathwayEnrich_dotplot"]] <- dp
    if (! returnDataOnly) {
      utils::write.csv(eres,
                file = file.path(outputPrefix, "inputGenes.Pathway.Analysis.csv"),
                row.names = FALSE)
      ggplot2::ggsave(plot = dp,
                      filename = file.path(outputPrefix, "inputGenes.Pathway.Analysis.png"),
                      height = 7, width = 10)
    }
  }
  return(resList)
}


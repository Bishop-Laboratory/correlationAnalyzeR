#' Analyze Gene Pairs
#'
#' Comaprison of two genes of interest using correlation values.
#' This can be 2 different genes in the same tissue or sample type. Or
#' the same gene accross two sample or tissue types.
#'
#' @param genesOfInterest A vector with two genes to compare.
#'
#' @param Sample_Type A vector of length-2 corresponding to the genesOfInterest.
#' Indicates the type of sample for each; "all", "normal", or "cancer".
#' Default: c("normal", "normal")
#'
#' @param Tissue A vector of length-2 corresponding to the genesOfInterest.
#' Indicates the type of tissue for each; Run getTissueTypes() to see available list.
#' Default: c("all", "all")
#'
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'. Default: "hsapiens".
#'
#' @param GSEA_Type Whether GSEA should consider all msigdb annotations,
#'     or just those in the most popular categories. Should be one of either
#'     'simple' or 'complex'.
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
analyzeGenePairs <- function(genesOfInterest,
                             Sample_Type = c("normal", "normal"),
                             Tissue = c("all", "all"),
                             Species = c("hsapiens", "mmusculus"),
                             GSEA_Type = c("simple", "complex"),
                             outputPrefix = "CorrelationAnalyzeR_Output_Paired",
                             runGSEA = T, topPlots = F, returnDataOnly = F) {
  # # Debug/Test
  # genesOfInterest <- c("BRCA1", "BRCA1")
  # GSEA_Type = "simple"
  # # set.seed(1)
  # runGSEA = T
  # outputPrefix = "tests/CorrelationAnalyzeR_Output_paired"
  # Species = "hsapiens"
  # Sample_Type = c("normal", "normal")
  # Tissue = c("all", "female0reproductive")
  # returnDataOnly = T
  # topPlots=F


  require(dplyr)
  require(tibble)

  if (length(genesOfInterest) == 2 ) {
    pairRes <- correlationAnalyzeR::analyzeSingleGenes(
      genesOfInterest = genesOfInterest,
      returnDataOnly = returnDataOnly, topPlots = topPlots,
      outputPrefix = outputPrefix, runGSEA = runGSEA,
      Sample_Type = Sample_Type, Tissue = Tissue,
      Species = Species, GSEA_Type = GSEA_Type
    )

  } else {
    stop("Please enter only 2 genes to compare")
  }

  # Compare correlations -- scatter plot
  correlations <- pairRes$correlations
  geneOne <- genesOfInterest[1]
  tissueOneRaw <- Tissue[1]
  tissueOne <- gsub(tissueOneRaw, pattern = "0", replacement = " ")

  sampleOne <- Sample_Type[1]
  geneOneTitle <- paste0(geneOne, ", ",
                         stringr::str_to_title(tissueOne),
                         " - ",
                         stringr::str_to_title(sampleOne))
  geneTwo <- genesOfInterest[2]
  tissueTwoRaw <- Tissue[2]
  tissueTwo <- gsub(tissueTwoRaw, pattern = "0", replacement = " ")

  sampleTwo <- Sample_Type[2]
  geneTwoTitle <- paste0(geneTwo, ", ",
                         stringr::str_to_title(tissueTwo),
                         " - ",
                         stringr::str_to_title(sampleTwo))
  longName <- ifelse((tissueOne != tissueTwo | sampleOne != sampleTwo),
                     yes = T, no = F)


  pairRes[["compared"]] <- list()
  # Variance heat map
  correlations$average <- rowMeans(correlations)
  correlations$variance <- matrixStats::rowVars(as.matrix(correlations[,c(1,2)]))
  correlations <- correlations[order(correlations$variance,
                                     decreasing = T),]
  pairRes[["compared"]][["correlations"]] <- correlations
  cn <- colnames(correlations)
  correlations2 <- correlations[which(! rownames(correlations) %in%
                                        colnames(correlations)),]

  corHeatOne <- correlations2 %>%
    rownames_to_column('gene') %>%
    dplyr::filter(eval(parse(text = cn[1])) > 0) %>%
    top_n(15, variance) %>%
    column_to_rownames('gene')
  corHeatTwo <- correlations2 %>%
    rownames_to_column('gene') %>%
    dplyr::filter(eval(parse(text = cn[1])) < 0) %>%
    top_n(15, variance) %>%
    column_to_rownames('gene')
  corHeat <- rbind(corHeatOne, corHeatTwo)
  corHeat <- corHeat[, c(1, 2)]
  colnames(corHeat) <- c(geneOne, geneTwo)

  # GSEA compare -- heatmap
  compPaths <- merge( x= pairRes[[geneOneTitle]][["GSEA"]][["eres"]],
                      y = pairRes[[geneTwoTitle]][["GSEA"]][["eres"]],
                      by = c("ID", "Description"))
  compPaths <- compPaths[,c(1, 5, 6, 7, 11, 14, 15, 16, 20)]
  if(longName) {
    colnames(compPaths) <- gsub(x = colnames(compPaths),pattern =  ".x",
                                replacement =  paste0("_", geneOne,
                                                      "_", gsub(tissueOne, pattern = " ", replacement = "_"),
                                                      "_", sampleOne))
    colnames(compPaths) <- gsub(x = colnames(compPaths),pattern =  ".y",
                                replacement =  paste0("_", geneTwo,
                                                      "_", gsub(tissueTwo, pattern = " ", replacement = "_"),
                                                      "_", sampleTwo))
  } else {
    colnames(compPaths) <- gsub(x = colnames(compPaths),pattern =  ".x",
                                replacement =  paste0("_", geneOne))
    colnames(compPaths) <- gsub(x = colnames(compPaths),pattern =  ".y",
                                replacement =  paste0("_", geneTwo))
  }

  compPaths$NES_average <- rowMeans(compPaths[,c(2, 6)])
  compPaths$NES_variance <- matrixStats::rowVars(as.matrix(compPaths[,c(2, 6)]))
  compPaths <- compPaths[order(compPaths$NES_variance, decreasing = T),]

  cn <- colnames(compPaths)
  cnes <- cn[grep(x = cn, pattern = "NES")]
  compHeatOne <- compPaths %>%
    dplyr::filter(eval(parse(text = cnes[1])) > 0) %>%
    top_n(15, NES_variance)  %>% slice(1:15)
  compHeatTwo <- compPaths %>%
    dplyr::filter(eval(parse(text = cnes[2])) > 0 & ! ID %in% compHeatOne$ID) %>%
    top_n(15, NES_variance) %>% slice(1:15)
  compHeat <- unique(rbind(compHeatOne, compHeatTwo))
  titleID <- compHeat$ID
  titleID <- correlationAnalyzeR::fixStrings(titleID)
  titleID[which(nchar(titleID) > 40)] <- paste0(substr(titleID[which(nchar(titleID) > 40)],
                                                       1, 40), "...")
  dups <- which(duplicated(titleID))
  if (length(dups)) {
    ends <- substr(compPaths$ID[dups],
                   nchar(compPaths$ID[dups])-3,
                   nchar(compPaths$ID[dups]))
    titleID[dups] <- paste0(substr(titleID[dups], 1, 34), "...", tolower(ends) )
  }

  dups <- which(duplicated(titleID))
  if (length(dups)) {
    titleID[dups] <- paste0(substr(titleID[dups], 1, 34),
                            "...", replicate(expr = paste0(sample(letters, 3),
                                                           collapse = ""),
                                             n = length(dups)))
  }
  compHeat <- compHeat[,c(2, 6)]

  rownames(compHeat) <- titleID[1:30]


  if (longName) {
    correlationsScatter <- correlations
    colnames(correlationsScatter)[c(1, 2)] <- c("x", "y")
    gs <- ggpubr::ggscatter(correlationsScatter,
                            title = paste0(geneOne, " vs. ",
                                           geneTwo, " Correlations"),
                            x = "x", y = "y",
                            ylab = geneTwoTitle,
                            xlab = geneOneTitle,
                            add = "reg.line", size = .1,
                            # cor.method = "spearman",
                            cor.coef = T, conf.int = T)
    pairRes[["compared"]][["correlationPlot"]] <- gs

    phCor <- pheatmap::pheatmap(corHeat, silent = T, angle_col = 0,
                                main = "Differentially correlated genes",
                                labels_col = c(geneOneTitle, geneTwoTitle),
                                cluster_rows = T, cluster_cols = F)
    pairRes[["compared"]][["correlationVarianceHeatmap"]] <- phCor
    phGSEA <- pheatmap::pheatmap(compHeat, silent = T, angle_col = 0,
                                 main = "Differentially correlated pathways",
                                 labels_col = c(geneOneTitle, geneTwoTitle),
                                 cluster_rows = T,
                                 cluster_cols = F)
    pairRes[["compared"]][["correlatedPathwaysHeatmap"]] <- phGSEA
    if (! returnDataOnly) {
      ggplot2::ggsave(phGSEA, height = 7.5, width = 6,
                      filename = file.path(outputPrefix,
                                           "GSEA_compared_heatmap.png"))
      ggplot2::ggsave(phCor, height = 7.5, width = 4.5,
                      filename = file.path(outputPrefix,
                                           "correlations_compared_heatmap.png"))
      ggplot2::ggsave(gs, filename = file.path(outputPrefix,
                                               "correlationScatterCompare.png"))
    }

  } else {
    gs <- ggpubr::ggscatter(correlations,
                            title = paste0(geneOne, " vs. ",
                                           geneTwo, " Correlations"),
                            caption = paste0(tissueOne, " - ", sampleOne),
                            x = geneOne, y = geneTwo,
                            ylab = geneTwo,
                            xlab = geneOne,
                            add = "reg.line", size = .1,
                            cor.coef = T, conf.int = T)
    pairRes[["compared"]][["correlationPlot"]] <- gs

    phCor <- pheatmap::pheatmap(corHeat,angle_col = 0,
                                silent = T,
                                cluster_rows = T, cluster_cols = F)
    pairRes[["compared"]][["correlationVarianceHeatmap"]] <- phCor
    phGSEA <- pheatmap::pheatmap(compHeat, silent = T,
                                 cluster_rows = T, angle_col = 0,
                                 labels_col = c(geneOne, geneTwo),
                                 cluster_cols = F)
    pairRes[["compared"]][["correlatedPathwaysHeatmap"]] <- phGSEA
    if (! returnDataOnly) {
      ggplot2::ggsave(phGSEA, 7.5, width = 6,
                      filename = file.path(outputPrefix,
                                           "GSEA_compared_heatmap.png"))
      ggplot2::ggsave(phCor, 7.5, width = 4.5,
                      filename = file.path(outputPrefix,
                                           "correlations_compared_heatmap.png"))
      ggplot2::ggsave(gs, filename = file.path(outputPrefix,
                                               "correlationScatterCompare.png"))
    }
  }
  pairRes[["compared"]][["correlatedPathwaysDataFrame"]] <- compPaths
  # Return results
  return(pairRes)

}






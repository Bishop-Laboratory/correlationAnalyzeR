#' Wrapper for clusterProlifer's GSEA()
#'
#' Runs GSEA() from clusterProfiler and creates useful
#' visualizations.
#'
#' @param ranks Numeric of gene 'scores' ordered by decreasing value and
#'     named with gene symbols.
#' @param TERM2GENE Data frame with two columns: gene set identifiers and
#' gene symbols. Can be generated using correlationAnalyzeR::getTERM2GENE()
#' @param plotFile prefix to use for naming output files.
#' @param outDir output directory.
#' @param Condition Name to use for titles of plots. Default = "GSEA Results".
#' @param nperm Number of permutations to run. Default is 2000
#' @param padjustedCutoff Value to use as a cutoff for returned gene sets.
#' @param returnDataOnly Should GSEA data/plots be saved to file? Default: TRUE
#' @param topPlots Should top GSEA pathways be plotted? Default: FALSE
#'
#' @return Named list containing GSEA() output, GSEA data frame, and visualizations.
#'
#' @examples
#' corrDF <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = c("BRCA1"),
#'  returnDataOnly = TRUE, runGSEA = FALSE, Sample_Type = "normal")
#' ranks <- corrDF$correlations[,1]
#' names(ranks) <- rownames(corrDF$correlations)
#' TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "simple",
#' Species = "hsapiens")
#' correlationAnalyzeR::myGSEA(ranks = ranks,
#' TERM2GENE = TERM2GENE,
#'        plotFile = "GSEA_out", outDir = getwd(),
#'        topPlots = FALSE, returnDataOnly=TRUE, Condition = "GSEA Results")
#'
#' @import dplyr
#' @import clusterProfiler
#'
#' @export
myGSEA <- function(ranks,
                   TERM2GENE,
                   padjustedCutoff = .05,
                   returnDataOnly = TRUE,
                   nperm = 2000,
                   topPlots = FALSE,
                   outDir,
                   Condition = "GSEA Results",
                   plotFile = "GSEA_results") {


  # # Bug testing
  # padjustedCutoff = .05
  # topPlots = FALSE
  # nperm = 2000
  # corrDF <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = c("ATM"),
  #                                                   returnDataOnly = TRUE,
  #                                                   runGSEA = FALSE,
  #                                                   Sample_Type = "normal")
  # ranks <- corrDF$correlations[,1]
  # names(ranks) <- rownames(corrDF$correlations)
  # TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "simple",
  #                                                Species = "hsapiens")

  resList <- list()
  ranks <- ranks[which(! duplicated(names(ranks)))]
  ranks <- ranks[which(! is.na(ranks))]
  ranks <- ranks[order(ranks, decreasing = TRUE)]

  EGMT <- GSEA2(TERM2GENE = TERM2GENE, ranks = ranks, nproc = 1,
                nperm = nperm, pvalueCutoff = padjustedCutoff)

  resGSEA <- as.data.frame(EGMT)

  resList[["EGMT"]] <- EGMT

  if (length(resGSEA$ID) < 10){
    warning(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                   padjustedCutoff, ". Rerunning with higher pValue."))
    EGMT <- GSEA2(TERM2GENE = TERM2GENE, ranks = ranks, nproc = 1,
                  nperm = nperm, pvalueCutoff = padjustedCutoff + .15)
    resGSEA <- as.data.frame(EGMT)
    resList[["EGMT"]] <- EGMT
  }
  if (length(resGSEA$ID) < 10){
    warning(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                   padjustedCutoff, ". Rerunning with higher pValue."))
    EGMT <- GSEA2(TERM2GENE = TERM2GENE, ranks = ranks, nproc = 1,
                  nperm = nperm, pvalueCutoff = padjustedCutoff + .45)
    resGSEA <- as.data.frame(EGMT)
    resList[["EGMT"]] <- EGMT
  }
  if (length(resGSEA$ID) < 10){
    stop(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                padjustedCutoff, ". Please check your data. If you believe this ",
                "behavior is a bug, please contact the package maintainer."))
  }
  if (topPlots) {
    resGSEA <- resGSEA[order(resGSEA$NES, decreasing = TRUE),]
    topUP <- resGSEA$ID[1:10]
    resGSEA <- resGSEA[order(resGSEA$NES, decreasing = FALSE),]
    topDOWN <- resGSEA$ID[1:10]
    resGSEA <- resGSEA[order(resGSEA$pvalue),]
    plUP <- list()
    plDOWN <- list()
    for ( i in 1:6 ) {
      pathway <- topUP[i]
      if (nchar(pathway) > 35) {
        pathTitle <- paste0(substr(pathway, 1, 30), "...")
      } else {
        pathTitle <- pathway
      }
      gp <- clusterProfiler::gseaplot(EGMT, pathway, title = NULL)
      gp <- gp + ggplot2::labs(title = pathTitle,
                               subtitle = paste0("Enrichment score: ",
                                                 round(resGSEA$NES[which(
                                                   resGSEA$ID == pathway
                                                 )], 3)))
      gg <- ggplot2::theme_classic()
      gp <-  gp + ggplot2::theme(plot.title = gg[["plot.title"]],
                                 plot.subtitle = gg[["plot.subtitle"]],
                                 plot.margin = gg[["plot.margin"]])
      if ( i == 1) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 45))
      } else if (i == 2) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 20))
      } else if (i == 3) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 45, 20, 20))
      } else if (i == 4) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 45))
      } else if (i == 5) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 20))
      } else if (i == 6) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 45, 45, 20))
      }
      plUP[[i]] <- gp

      pathway <- topDOWN[i]
      if (nchar(pathway) > 35) {
        pathTitle <- paste0(substr(pathway, 1, 30), "...")
      } else {
        pathTitle <- pathway
      }
      gp <- clusterProfiler::gseaplot(EGMT, pathway, title = NULL)
      gp <- gp + ggplot2::labs(title = pathTitle,
                               subtitle = paste0("Enrichment score: ",
                                                 round(resGSEA$NES[which(
                                                   resGSEA$ID == pathway
                                                 )], 3)))
      gg <- ggplot2::theme_classic()
      gp <-  gp + ggplot2::theme(plot.title = gg[["plot.title"]],
                                 plot.subtitle = gg[["plot.subtitle"]],
                                 plot.margin = gg[["plot.margin"]])
      if (i == 1) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 45))
      } else if (i == 2) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 20))
      } else if (i == 3) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 45, 20, 20))
      } else if (i == 4) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 45))
      } else if (i == 5) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 20))
      } else if (i == 6) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 45, 45, 20))
      }
      plDOWN[[i]] <- gp

    }

    gaUP <- ggpubr::ggarrange(plotlist = plUP, nrow = 2, ncol = 3)
    gaUP <- ggpubr::annotate_figure(gaUP,
                                    top = ggpubr::text_grob(paste0(
                                      "Top Over-Expressed Pathways in ", Condition
                                    ),
                                    size = 35)
    )
    resList[["GSEA_up"]] <- gaUP
    if (! returnDataOnly) {
      ggplot2::ggsave(plot = gaUP,
                      filename = file.path(outDir, paste0(plotFile, "_topPathwaysUP.png")),
                      height = 14, width = 20)
    }


    gaDOWN <- ggpubr::ggarrange(plotlist = plDOWN, nrow = 2, ncol = 3)
    gaDOWN <- ggpubr::annotate_figure(gaDOWN,
                                      top = ggpubr::text_grob(paste0(
                                        "Top Under-Expressed Pathways in ", Condition
                                      ),
                                      size = 35)
    )

    if (! returnDataOnly) {
      ggplot2::ggsave(plot = gaDOWN,
                      filename = file.path(outDir, paste0(plotFile, "_topPathwaysDOWN.png")),
                      height = 14, width = 20)
    }
    resList[["GSEA_down"]] <- gaDOWN

  }

  if (! returnDataOnly) {
    data.table::fwrite(x = resGSEA, file = file.path(outDir,
                                                     paste0(plotFile,
                                                            "_GSEA.csv")))
  }

  resList[["eres"]] <- resGSEA

  cat("\nReturning ... ", names(resList), "\n")

  return(resList)
}



# Modified GSEA from clusterProfiler. Reduces computational time.
GSEA2 <- function(TERM2GENE, ranks,
                  nperm = 2000, nproc = "auto",
                  pvalueCutoff = .05) {
  if (nproc == "auto") {
    nproc = parallel::detectCores()
  }
  TERMList <- TERM2GENE %>% split(x = TERM2GENE$gene_symbol, f = TERM2GENE$gs_name)
  EGMT <- suppressWarnings(fgsea::fgsea(pathways = TERMList, nproc = nproc,
                       maxSize = 500,
                       minSize = 15,
                       stats = ranks, nperm = nperm))
  res <- data.frame(
    ID = as.character(EGMT$pathway),
    Description = as.character(EGMT$pathway),
    setSize = EGMT$size,
    enrichmentScore = EGMT$ES,
    NES = EGMT$NES,
    pvalue = EGMT$pval,
    p.adjust = EGMT$padj,
    core_enrichment = vapply(EGMT$leadingEdge, FUN.VALUE = "char",
                             paste0, collapse='/'),
    stringsAsFactors = FALSE
  )
  res <- res[!is.na(res$pvalue),]
  res <- res[ res$pvalue <= pvalueCutoff, ]
  res <- res[ res$p.adjust <= pvalueCutoff, ]
  idx <- order(res$pvalue, decreasing = FALSE)
  res <- res[idx, ]
  params <- list(pvalueCutoff = pvalueCutoff,
                 nPerm = nperm,
                 pAdjustMethod = "BH",
                 exponent = 1,
                 minGSSize = 15,
                 maxGSSize = 500
  )
  row.names(res) <- res$ID
  EGMT <- new("gseaResult",
              result     = res,
              geneSets   = TERMList,
              geneList   = ranks,
              params     = params,
              readable   = FALSE
  )
  EGMT@organism <- "UNKNOWN"
  EGMT@setType <- "UNKNOWN"
  EGMT@keytype <- "UNKNOWN"
  return(EGMT)
}

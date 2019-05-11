#' Wrapper for clusterProlifer's GSEA function
#'
#' Run's the GSEA function from clusterProfiler and creates useful
#'     visualizations. It also outputs GSEA results to csv.
#'
#' @param ranks Numeric of gene 'scores' ordered by decreasing value and
#'     named with gene symbols.
#' @param TERM2GENE Data frame with two columns: gene set identifiers and
#'     gene symbols.
#' @param plotFile prefix to use for naming output files.
#' @param outDir output directory.
#' @param Condition Name to use for titles of plots.
#' @param padjustedCutoff Value to use as a cutoff for returned gene sets.
#'
#' @return Output of GSEA function in clusterProfiler.
#'
#' @examples
#' myGSEA(ranks = DGE_results$GSEA, TERM2GENE = msigdbr()[,c(1,5)],
#'        plotFile = "GSEA_out", outDir = getwd(), Condition = "GSEA Results")
#'
#' @export
myGSEA <- function(ranks, TERM2GENE, plotFile,
                   outDir, Condition = Condition,
                   padjustedCutoff = .05) {
  ranks <- ranks[which(! is.na(ranks))]
  ranks <- ranks[order(ranks, decreasing = T)]
  y <- clusterProfiler::GSEA(ranks, TERM2GENE=TERM2GENE,
                             nPerm = 1000, pvalueCutoff = padjustedCutoff)
  resGSEA <- as.data.frame(y)
  if (length(resGSEA$ID) < 10){
    warning(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                   padjustedCutoff, ". Rerunning with higher pValue."))
    y <- clusterProfiler::GSEA(ranks, TERM2GENE=TERM2GENE,
                               nPerm = 1000, pvalueCutoff = .2)
    resGSEA <- as.data.frame(y)
  }
  if (length(resGSEA$ID) < 10){
    warning(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                   padjustedCutoff, ". Rerunning with higher pValue."))
    y <- clusterProfiler::GSEA(ranks, TERM2GENE=TERM2GENE,
                               nPerm = 1000, pvalueCutoff = .6)
    resGSEA <- as.data.frame(y)
  }
  resGSEA <- resGSEA[order(resGSEA$NES, decreasing = T),]
  topUP <- resGSEA$ID[1:10]
  resGSEA <- resGSEA[order(resGSEA$NES, decreasing = F),]
  topDOWN <- resGSEA$ID[1:10]
  plUP <- list()
  plDOWN <- list()
  for ( i in 1:6 ) {
    pathway <- topUP[i]
    if (nchar(pathway) > 35) {
      pathTitle <- paste0(substr(pathway, 1, 30), "...")
    } else {
      pathTitle <- pathway
    }
    gp <- clusterProfiler::gseaplot(y, pathway, title = NULL)
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
    gp <- clusterProfiler::gseaplot(y, pathway, title = NULL)
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
  ggplot2::ggsave(plot = gaUP,
         filename = file.path(outDir, paste0(plotFile, "_topPathwaysUP.png")),
         height = 14, width = 20)


  gaDOWN <- ggpubr::ggarrange(plotlist = plDOWN, nrow = 2, ncol = 3)
  gaDOWN <- ggpubr::annotate_figure(gaDOWN,
                            top = ggpubr::text_grob(paste0(
                              "Top Under-Expressed Pathways in ", Condition
                              ),
                                            size = 35)
  )
  ggplot2::ggsave(plot = gaDOWN,
         filename = file.path(outDir, paste0(plotFile, "_topPathwaysDOWN.png")),
         height = 14, width = 20)

  resGSEA <- resGSEA[order(resGSEA$NES, decreasing = T),]
  data.table::fwrite(x = resGSEA, file = file.path(outDir,
                                                   paste0(plotFile,
                                                          "_GSEA.csv")))
  return(y)
}

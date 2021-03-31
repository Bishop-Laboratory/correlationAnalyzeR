#' Analyze Gene Pairs
#'
#' Comaprison of two genes of interest using correlation values.
#' This can be 2 different genes in the same tissue or sample type or
#' the same gene accross two sample or tissue types. Alternatively, specify 'crossCompareMode'
#' to view compared correrlations across all available tissue types.
#' @param genesOfInterest A length-two vector with genes to compare.
#' @param Sample_Type A length-two vector of sample types corresponding to genesOfInterest.
#' Choose "all", "normal", or "cancer". Default: c("normal", "normal")
#' @param Tissue A length-two vector of tissue types corresponding to genesOfInterest.
#' Run getTissueTypes() to see available list. Default: c("all", "all")
#' @param GSEA_Type Which GSEA annotations should be considered? Options listed in
#' correlationAnalyzeR::pathwayCategories -- See details of getTERM2GENE for more info.
#' @param nperm Number of permutations to run in GSEA. Default is 2000
#' @param sampler If TRUE, will only return 100,000 random genesets. Useful for reducing GSEA computational burden.
#' @param crossCompareMode Use this mode to generate comparisons
#' across all tissue and disease types. If both genes for genesOfInterest are the
#' same -- will compare normal vs cancer for that gene in each available tissue. Else, will
#' perform comparison of two different genes in all tissue-disease groups.
#' Will only consider user input for returnDataOnly, outputPrefix, and genesOfInterest.
#' @param outputPrefix Prefix for saved files -- the directory name to store output files in. If
#' folder does not exist, it will be created.
#' @param runGSEA If TRUE will run GSEA using gene correlation values.
#' @param TERM2GENE Mapping of geneset IDs to gene names. If NULL, it will be
#' generated automatically. Only applicable if GSEA is to be run.
#' @param returnDataOnly if TRUE will return result list object
#' and will not generate any folders or files.
#' @param topPlots Logical. If TRUE, myGSEA() will build gsea plots for top correlated genesets.
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied.
#' @param makePool Logical. Should a pool be created if one is not supplied? Default: FALSE.
#' @return A named list containing visualizations and correlation data from paired analysis.
#' @examples
#' genesOfInterest <- c("ATM", "SLC7A11")
#' correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#'                               GSEA_Type = "simple", returnDataOnly = TRUE,
#'                               Sample_Type = c("normal", "normal"),
#'                               Tissue = c("brain", "brain"))
#' genesOfInterest <- c("BRCA1", "BRCA1")
#' correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#'                               GSEA_Type = "simple", returnDataOnly = TRUE,
#'                               Sample_Type = c("normal", "cancer"),
#'                               Tissue = c("respiratory", "respiratory"))
#' genesOfInterest <- c("NFKB1", "SOX10")
#' correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#'                               returnDataOnly = TRUE,
#'                               crossCompareMode = TRUE)
#' @importFrom rlang .data
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @export
#'
analyzeGenePairs <- function(genesOfInterest,
                             Sample_Type = c("normal", "normal"),
                             Tissue = c("all", "all"),
                             #Species = c("hsapiens", "mmusculus"),
                             GSEA_Type = c("simple"),
                             outputPrefix = "CorrelationAnalyzeR_Output_Paired",
                             crossCompareMode = FALSE,
                             runGSEA = TRUE,
                             TERM2GENE = NULL,
                             nperm = 2000, sampler = FALSE,
                             topPlots = FALSE, returnDataOnly = TRUE,
                             pool = NULL,
                             makePool = FALSE) {

  # genesOfInterest <- c("BRCA1", "NQO1")
  # Species <- "hsapiens"
  # crossCompareMode = FALSE
  # returnDataOnly = TRUE
  # outputPrefix = "CorrelationAnalyzeR_Output_Paired"
  # runGSEA = TRUE
  # topPlots = FALSE
  # Sample_Type = c("normal", "normal")
  # Tissue = c("all", "all")
  # TERM2GENE = NULL
  # GSEA_Type = c("simple")
  # pool = NULL
  # sampler = FALSE
  # nperm = 2000
  # makePool = FALSE
  # Sample_Type = c("cancer", "cancer")
  # Tissue = c("bone", "bone")

  # Validate inputs
  unGene <- unique(genesOfInterest)
  unTissue <- unique(Tissue)
  unSample <- unique(Sample_Type)
  if (length(unGene) == 1 && length(unTissue) == 1 && length(unSample) == 1) {
    stop("Genes, Tissues, or Samples must be different in gene vs gene mode!")
  }


  getPhBreaks <- function(mat, palette = NULL) {
    # From https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
    if (is.null(palette)) {
      palette <- grDevices::colorRampPalette(rev(
        RColorBrewer::brewer.pal(n = 7, name =
                     "RdYlBu")))(100)
    }
    n <- length(palette)
    breaks <- c(seq(min(mat), 0, length.out=ceiling(n/2) + 1),
                  seq(max(mat)/n, max(mat), length.out=floor(n/2)))
    return(list(palette, breaks))
  }

  if (! is.null(pool)) {
    if (! pool$valid) {
      pool <- NULL
    }
  }
  if (is.null(pool)) {
    if (makePool) {
      retryCounter <- 1
      # cat("\nEstablishing connection to database ... \n")
      while(is.null(pool)) {
        pool <- try(silent = T, eval({
          pool::dbPool(
            drv = RMySQL::MySQL(),
            user = "public-rds-user@m2600az-db01p.mysql.database.azure.com", port = 3306,
            dbname="correlation_analyzer",
            password='public-user-password',
            host="m2600az-db01p.mysql.database.azure.com"
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

      on.exit(pool::poolClose(pool))
    }
  }

  lm_eqn <- function(df){
    m <- stats::lm(eval(parse(text = colnames(df)[2])) ~ eval(parse(text = colnames(df)[1])), df)
    r <- sqrt(summary(m)$r.squared) * sign(unname(stats::coef(m)[2]))
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)~"="~r,
                     list(a = format(unname(stats::coef(m)[1]), digits = 2),
                          b = format(unname(stats::coef(m)[2]), digits = 2),
                          r = format(r, digits = 2)))
    as.character(as.expression(eq));
  }

  # If running in cross-comparison mode...
  if (crossCompareMode) {

    cat("\nRunning cross comparison mode ... \n")
    availTissue <- correlationAnalyzeR::getTissueTypes(#Species = Species,
                                                       pool = pool,
                                                       useBlackList = TRUE)
    runGSEA <- F
    if (genesOfInterest[1] == genesOfInterest[2]) {
      cat("\nGene one is the same as gene two ... \n")
      # if (Species == "mmusculus") {
      #   stop("\nOnly normal tissues available for mouse",
      #        " due to black-listing of cancer groups for low quality.",
      #        "\nEnter two different genes or choose human ... \n")
      # }
      cat("\nWill perform normal vs cancer comparison on",
          genesOfInterest[1], "... \n")
      mode <- "cross_normalVsCancer"
      df <- as.data.frame(table(gsub(availTissue,
                                     pattern = " - .*",
                                     replacement = "")), stringsAsFactors = FALSE)
      goodTissues <- df$Var1[which(df$Freq == 3)]
      Tissue <- rep(goodTissues, each = 2)
      Sample_Type <- rep(c("normal", "cancer"), length(goodTissues))
      genesVec <- rep(genesOfInterest[1], length(Sample_Type))
    } else {
      geneOne <- genesOfInterest[1]
      geneTwo <- genesOfInterest[2]
      cat("\nGene one is not the same as gene two ... \n")
      cat("\nWill perform comparison of",
          geneOne, "and",
          geneTwo, "across all available tissue-disease conditions... \n")
      mode <- "cross_geneVsGene"
      availTissue <- availTissue[grep(availTissue, pattern = "all", invert = T)]
      availTissue <- strsplit(availTissue, split = " - ")
      TissueO <- vapply(availTissue, FUN = "[[", FUN.VALUE = "character", 1)
      Tissue <- rep(TissueO, 2)
      Sample_TypeO <- vapply(availTissue, FUN = "[[", FUN.VALUE = "character", 2)
      Sample_Type <- rep(Sample_TypeO, 2)
      genesVec <- rep(genesOfInterest, each =  length(TissueO))

    }
    # Get VST for each gene
    geneUnique <- unique(genesVec)
    geneVSTList <- correlationAnalyzeR::getTissueVST(genesOfInterest = geneUnique,
                                                     #Species = Species,
                                                     Tissues = "all",
                                                     pool = pool,
                                                     Sample_Type = "all",
                                                     useBlackList = TRUE)
    geneVSTDF <- data.table::rbindlist(geneVSTList, idcol = "group")
    geneVSTDF$Tissue <- stringr::str_to_title(gsub(geneVSTDF$group, pattern = "(.+) - (.+)", replacement = "\\1"))
    geneVSTDF$sampleType <-  stringr::str_to_title(gsub(geneVSTDF$group, pattern = "(.+) - (.+)", replacement = "\\2"))
    if (mode == "cross_geneVsGene") {
      geneVSTDF <- geneVSTDF[, c(1, 5, 6, 2, 3, 4)]
    } else {
      geneVSTDF <- geneVSTDF[, c(1, 4, 5, 2, 3)]
    }
    colnames(geneVSTDF)[c(1:4)] <- c("Group", "Tissue", "sampleType", "Samples")

    # Make VST plot
    crossCompareResVST <- list()
    if (mode == "cross_geneVsGene") {
      geneVSTDF <- geneVSTDF[geneVSTDF$sampleType != "All",]
      geneVSTDFToPlot <- geneVSTDF[,c(1,2, 3, 5, 6)]
      geneVSTDFToPlot <- geneVSTDFToPlot %>%
        gather("Gene", "VST", -.data$Group, -.data$Tissue, -.data$sampleType)
      geneOne <- geneUnique[1]
      geneTwo <- geneUnique[2]
      fillStr <- "Group"
      geneVSTDFToPlot1 <- geneVSTDFToPlot[which(geneVSTDFToPlot$Gene == geneOne),]
      groups <- unique(geneVSTDFToPlot1$Group)
      groups <- groups[order(groups)]
      plotOne <- ggpubr::ggboxplot(data = geneVSTDFToPlot1, order = groups,
                                   x = fillStr, #facet.by = "Gene",
                                   title = paste0(geneOne,
                                                  " expression across tissue groups"),
                                   ylab = "Expression (VST counts)",
                                   fill = fillStr,
                                   y = "VST") +
        ggpubr::rotate_x_text(angle = 45) +
        ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 20)) +
        ggpubr::rremove("legend") +
        ggpubr::rremove("xlab")
      geneVSTDFToPlot2 <- geneVSTDFToPlot[which(geneVSTDFToPlot$Gene == geneTwo),]
      plotTwo <- ggpubr::ggboxplot(data = geneVSTDFToPlot2, order = groups,
                                   x = fillStr, #facet.by = "Gene",
                                   title = paste0(geneTwo,
                                                  " expression across conditions"),
                                   ylab = "Expression (VST counts)",
                                   fill = fillStr,
                                   y = "VST") +
        ggpubr::rotate_x_text(angle = 45) +
        ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 20)) +
        ggpubr::rremove("legend") +
        ggpubr::rremove("xlab")
      colnames(geneVSTDF)[c(5:6)] <- paste0(colnames(geneVSTDF)[c(5:6)], "_VST")
      crossCompareResVST[["VST_boxPlotOne"]] <- plotOne
      crossCompareResVST[["VST_boxPlotTwo"]] <- plotTwo

    } else {
      geneVSTDFToPlot <- geneVSTDF
      colnames(geneVSTDFToPlot)[length(colnames(geneVSTDFToPlot))] <- "VST"

      goodTiss <- unique(geneVSTDFToPlot$Tissue[which(geneVSTDFToPlot$sampleType == "Cancer")])
      goodTiss2 <- unique(geneVSTDFToPlot$Tissue[which(geneVSTDFToPlot$sampleType == "Normal")])
      goodTissFinal <- goodTiss[which(goodTiss %in% goodTiss2)]
      geneVSTDF2 <- geneVSTDFToPlot[which(geneVSTDFToPlot$Tissue %in% goodTissFinal &
                                            geneVSTDFToPlot$sampleType != "All"),]
      maxHeight <- max(geneVSTDF2$VST)
      meds <- sapply(unique(geneVSTDF2$Tissue), FUN = function(groupNow) {
        stats::median(geneVSTDF2$VST[geneVSTDF2$Tissue == groupNow &
                                  geneVSTDF2$sampleType == "Cancer"]) -
          stats::median(geneVSTDF2$VST[geneVSTDF2$Tissue == groupNow &
                                    geneVSTDF2$sampleType == "Normal"])
      })
      meds <- meds[order(abs(meds))]
      VSTBPproto <- ggpubr::ggboxplot(data = geneVSTDF2,
                                      x = "Tissue", order = names(meds),
                                      title = paste0(geneUnique[1],
                                                     " expression across tissues"),
                                      ylab = "Expression (VST counts)",
                                      fill = "sampleType", legend = "right",
                                      y = "VST") +
        ggpubr::stat_compare_means(ggplot2::aes_string(group = "sampleType"),
                                   label.y = (maxHeight * 1.15),
                                   hide.ns = TRUE, label = "p.signif")
      plot <- VSTBPproto +
        ggpubr::rotate_x_text(angle = 45) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_text(size = 16),
          title = ggplot2::element_text(size = 22),
          plot.margin = ggplot2::margin(10, 25, 10, 25)
        ) +
        ggpubr::rremove("xlab") + ggpubr::rremove("legend.title")
      colnames(geneVSTDF)[c(5)] <- paste0(colnames(geneVSTDF)[c(5)], "_VST")
      crossCompareResVST[["VST_boxPlot"]] <- plot
    }

    crossCompareResVST[["VST_DF"]] <- geneVSTDF

    # Do paired to get correlation data
    pairRes <- correlationAnalyzeR::analyzeSingleGenes(
      genesOfInterest = genesVec, pool = pool, TERM2GENE = TERM2GENE,
      returnDataOnly = returnDataOnly, topPlots = topPlots,
      outputPrefix = outputPrefix, runGSEA = FALSE,
      Sample_Type = Sample_Type, Tissue = Tissue,
      #Species = Species,
      GSEA_Type = GSEA_Type
    )
    n <- length(names(pairRes))
    correlations <- pairRes$correlations
    pVals <- pairRes[["P values"]]
    if (mode != 'cross_geneVsGene') {
      newNames <- gsub(names(pairRes), pattern = ".*, |- ", replacement = "")
      oldNames <- gsub(names(pairRes), pattern = ".*, ", replacement = "")
      newNames <- gsub(newNames, pattern = " ", replacement = "_")
      newNames <- newNames[! newNames %in% c("correlations", "P_values", "P values")]
      oldNames <- oldNames[! oldNames %in% c("correlations", "P_values", "P values")]
    } else {
      newNames <- gsub(names(pairRes), pattern = ", | - ", replacement = "_")
      tissueFull <- gsub(names(pairRes), pattern = ".*, |- ", replacement = "")
      tissueFull <- gsub(tissueFull, pattern = "_", replacement = " ")
      oldNames <- gsub(names(pairRes), pattern = ",.*", replacement = "")
      tissueFull <- tissueFull[! tissueFull %in% c("correlations", "P_values", "P values")]
      newNames <- newNames[! newNames %in% c("correlations", "P_values", "P values")]
      oldNames <- oldNames[! oldNames %in% c("correlations", "P_values", "P values")]
      newOrder <- order(tissueFull)
      tissueFull <- tissueFull[newOrder]
      newNames <- newNames[newOrder]
      oldNames <- oldNames[newOrder]
      correlations <- correlations[,newOrder]
      pVals <- pVals[,newOrder]
    }

    colnames(correlations) <- newNames
    colnames(pVals) <- newNames
    resList <- list()
    resList[["pairResList"]] <- list()
    for (i in 1:length(colnames(correlations))) {
      samp <- colnames(correlations)[i]
      if (i %% 2 == 0) {
        dfRaw <- correlations[,c((i-1), i)]
        df <- dfRaw
        tempList <- list()
        sampleType <- stringr::str_match(colnames(df)[1], pattern = "^.+_([a-zA-Z]+$)")[,2]
        tissueSmall <- stringr::str_match(colnames(df)[1], pattern = "^(.+)_[a-zA-Z]+$")[,2]
        tissueSmall <- gsub(tissueSmall, pattern = "_", replacement = " ")
        tissue <- paste0(tissueSmall, " - ", sampleType)
        if (genesOfInterest[1] != genesOfInterest[2]) {
          tissue <- tissueFull[i]
        }
        df$Gene <- F
        df$Gene[which(rownames(df) %in% unique(genesVec))] <- T
        titleStr <- ifelse(genesOfInterest[1] == genesOfInterest[2], yes = tissueSmall,
                           no = tissue)
        labb <- lm_eqn(df)
        # xtex <- eval(parse(text = colnames(df)[1]))
        # ytex <- eval(parse(text = colnames(df)[2]))
        gp <- ggplot2::ggplot(data = df,
                              mapping = ggplot2::aes_string(x = colnames(df)[1],
                                                            y = colnames(df)[2])) +
          ggplot2::stat_bin2d(bins = 150) +
          ggplot2::geom_smooth(colour="black", size = 1.25,
                               method='lm') +
          ggplot2::labs(title = titleStr) +
          ggplot2::ylab(oldNames[(i-1)]) +
          ggplot2::xlab(oldNames[i]) +
          ggplot2::annotate("text", x = 0, y = 1.2,
                            label = labb,
                            parse = TRUE) +
          ggpubr::theme_pubr() +
          ggplot2::theme(axis.title = ggplot2::element_text(size = 20),
                         title = ggplot2::element_text(size = 22),
                         plot.margin = ggplot2::margin(10, 10, 10, 10),
                         axis.text = ggplot2::element_text(size = 16)) +
          ggplot2::theme(legend.position = "none")
        tempList[["scatterPlot"]] <- gp
        dfRaw$Variance <- matrixStats::rowVars(as.matrix(dfRaw))
        tempList[["correlations"]] <- dfRaw
        dfRaw <- dfRaw[which(! rownames(dfRaw) %in% genesOfInterest),]
        dfRawUp <- dfRaw[dfRaw[,1] > 0,]
        dfRawUp <- dfRawUp[order(dfRawUp$Variance, decreasing = TRUE),]
        dfRawUpSmall <- dfRawUp[c(1:15),]
        dfRawDn <- dfRaw[dfRaw[,1] < 0,]
        dfRawDn <- dfRawDn[order(dfRawDn$Variance, decreasing = TRUE),]
        dfRawDnSmall <- dfRawDn[c(1:15),]
        dfPh <- rbind(dfRawUpSmall, dfRawDnSmall)
        dfPh <- dfPh[,c(-3)]
        if (mode == "cross_geneVsGene") {
          breaks <- getPhBreaks(dfPh)
          ph <- pheatmap::pheatmap(dfPh, cluster_cols = FALSE,
                                   breaks = breaks[[2]],# fontsize = 18, fontsize_row = 16, fontsize_col = 18,
                                   silent = TRUE, angle_col = 0, main = titleStr,
                                   labels_col = c(genesOfInterest[1],
                                                  genesOfInterest[2]))
        } else {
          breaks <- getPhBreaks(dfPh)
          ph <- pheatmap::pheatmap(dfPh, cluster_cols = FALSE,
                                   breaks = breaks[[2]], # fontsize = 18, fontsize_row = 16, fontsize_col = 18,
                                   silent = TRUE, angle_col = 0, main = titleStr,
                                   labels_col = c("Normal", "Cancer"))
        }

        tempList[["heatMap"]] <- ggplotify::as.ggplot(ph)

        resList[["pairResList"]][[i/2]] <- tempList
        names(resList[["pairResList"]])[i/2] <- tissue
      }
    }

    correlations$average <- rowMeans(correlations)
    correlations$variance <- matrixStats::rowVars(as.matrix(correlations))
    correlations <- correlations[order(correlations$variance,
                                       decreasing = TRUE),]
    resList[["Correlations"]] <- correlations
    resList[["P values"]] <- pVals
    resList[["crossCompareVST"]] <- crossCompareResVST
    resList[["mode"]] <- mode
    return(resList)

  }

  # If running in normal mode ...
  if (length(genesOfInterest) == 2 ) {
    pairRes <- correlationAnalyzeR::analyzeSingleGenes(
      genesOfInterest = genesOfInterest, pool = pool,
      returnDataOnly = returnDataOnly, topPlots = topPlots,
      outputPrefix = outputPrefix, TERM2GENE = TERM2GENE,
      runGSEA = runGSEA, nperm = nperm, sampler = sampler,
      Sample_Type = Sample_Type, Tissue = Tissue,
      #Species = Species,
      GSEA_Type = GSEA_Type
    )

  } else {
    stop("Please enter only 2 genes to compare")
  }

  # Compare correlations -- scatter plot
  correlations <- pairRes$correlations
  pVals <- pairRes[["P values"]]
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
                     yes = TRUE, no = FALSE)
  if (geneOne != geneTwo) {
    geneOneTitleGP <- geneOneTitle
    geneTwoTitleGP <- geneTwoTitle
  } else {
    geneOneTitleGP <- paste0("Gene correlation (", stringr::str_to_title(tissueOne),
                            " - ",
                            stringr::str_to_title(sampleOne),
                            ")")
    geneTwoTitleGP <- paste0("Gene correlation (", stringr::str_to_title(tissueTwo),
                             " - ",
                             stringr::str_to_title(sampleTwo),
                             ")")
  }



  pairRes[["compared"]] <- list()


  # Variance heat map
  correlations$average <- rowMeans(correlations)
  correlations$variance <- matrixStats::rowVars(as.matrix(correlations[,c(1,2)]))
  correlations <- correlations[order(correlations$variance,
                                     decreasing = TRUE),]
  pairRes[["compared"]][["correlations"]] <- correlations

  pairRes[["compared"]][["P values"]] <- pVals

  cn <- colnames(correlations)
  correlations2 <- correlations[which(! rownames(correlations) %in%
                                        colnames(correlations)),]

  # Divergence of gene correlations
  corHeatOne <- correlations2 %>%
    rownames_to_column('gene') %>%
    dplyr::filter(eval(parse(text = cn[1])) > 0) %>%
    top_n(15, .data$variance) %>%
    column_to_rownames('gene')
  corHeatTwo <- correlations2 %>%
    rownames_to_column('gene') %>%
    dplyr::filter(eval(parse(text = cn[1])) < 0) %>%
    top_n(15, .data$variance) %>%
    column_to_rownames('gene')
  corHeat <- rbind(corHeatOne, corHeatTwo)
  corHeatVar <- corHeat[, c(1, 2)]
  colnames(corHeatVar) <- c(geneOne, geneTwo)

  # Similarity of gene correlations
  # First get the top 1000 genes by both columns
  # ( We want similar and meaninful genes)
  correlations21 <- correlations2[order(abs(correlations2[,1]),
                                       decreasing = TRUE),]
  correlations21 <- correlations21[c(1:1000),]
  correlations22 <- correlations2[order(abs(correlations2[,2]),
                                        decreasing = TRUE),]
  correlations22 <- correlations22[c(1:1000),]
  correlations2 <- unique(rbind(correlations21, correlations22))

  corHeatOne <- correlations2 %>%
    rownames_to_column('gene') %>%
    dplyr::filter(eval(parse(text = cn[1])) > 0) %>%
    top_n(15, -.data$variance) %>%
    column_to_rownames('gene')
  corHeatTwo <- correlations2 %>%
    rownames_to_column('gene') %>%
    dplyr::filter(eval(parse(text = cn[1])) < 0) %>%
    top_n(15, -.data$variance) %>%
    column_to_rownames('gene')
  corHeat <- rbind(corHeatOne, corHeatTwo)
  corHeatSim <- corHeat[, c(1, 2)]
  colnames(corHeatSim) <- c(geneOne, geneTwo)

  if (runGSEA) {
    # GSEA compare -- heatmap
    compPaths <- merge( x= pairRes[[geneOneTitle]][["GSEA"]][["eres"]],
                        y = pairRes[[geneTwoTitle]][["GSEA"]][["eres"]],
                        by = c("ID", "Description"))

    compPaths <- compPaths[,c(1, 5, 6, 7, 8, 11, 12, 13, 14)]
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
    # Divergence of pathways
    compPaths$NES_average <- rowMeans(compPaths[,c(2, 6)])
    compPaths$NES_variance <- matrixStats::rowVars(as.matrix(compPaths[,c(2, 6)]))
    compPaths <- compPaths[order(compPaths$NES_variance, decreasing = TRUE),]
    cn <- colnames(compPaths)
    cnes <- cn[grep(x = cn, pattern = "NES")]

    compHeatOne <- compPaths %>%
      dplyr::filter(eval(parse(text = cnes[1])) > 0) %>%
      top_n(15, .data$NES_variance)  %>% slice(1:15)
    compHeatTwo <- compPaths %>%
      dplyr::filter(eval(parse(text = cnes[2])) > 0 & ! .data$ID %in% compHeatOne$ID) %>%
      top_n(15, .data$NES_variance) %>% slice(1:15)
    compHeat <- unique(rbind(compHeatOne, compHeatTwo))
    if (! length(compHeat$ID)) {
      compHeatOne <- compPaths %>%
        dplyr::filter(eval(parse(text = cnes[1])) < 0) %>%
        top_n(15, .data$NES_variance)  %>% slice(1:15)
      compHeatTwo <- compPaths %>%
        dplyr::filter(eval(parse(text = cnes[2])) < 0 & ! .data$ID %in% compHeatOne$ID) %>%
        top_n(15, .data$NES_variance) %>% slice(1:15)
      compHeat <- unique(rbind(compHeatOne, compHeatTwo))
    }
    if (! length(compHeat$ID)) {
      compHeat <- compPaths
    }
    n <- length(compHeat$ID)
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
    compHeatVar <- compHeat[,c(2, 6)]
    rownames(compHeatVar) <- titleID[1:n]


    # Similarity of pathways
    compHeatOne <- compPaths %>%
      dplyr::filter(eval(parse(text = cnes[1])) > 0) %>%
      top_n(15, -.data$NES_variance)  %>% slice(1:15)
    compHeatTwo <- compPaths %>%
      dplyr::filter(eval(parse(text = cnes[2])) > 0 & ! .data$ID %in% compHeatOne$ID) %>%
      top_n(15, -.data$NES_variance) %>% slice(1:15)
    compHeat <- unique(rbind(compHeatOne, compHeatTwo))
    if (! length(compHeat$ID)) {
      compHeatOne <- compPaths %>%
        dplyr::filter(eval(parse(text = cnes[1])) < 0) %>%
        top_n(15, .data$NES_variance)  %>% slice(1:15)
      compHeatTwo <- compPaths %>%
        dplyr::filter(eval(parse(text = cnes[2])) < 0 & ! .data$ID %in% compHeatOne$ID) %>%
        top_n(15, .data$NES_variance) %>% slice(1:15)
      compHeat <- unique(rbind(compHeatOne, compHeatTwo))
    }
    if (! length(compHeat$ID)) {
      compHeat <- compPaths
    }

    n <- length(compHeat$ID)
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
    compHeatSim <- compHeat[,c(2, 6)]
    rownames(compHeatSim) <- titleID[1:n]
  }



  # Get VST
  tissueOneNow <- gsub(tissueOne, pattern = " ", replacement = "0")
  tissueTwoNow <- gsub(tissueTwo, pattern = " ", replacement = "0")
  VSTList <- correlationAnalyzeR::getTissueVST(genesOfInterest = c(geneOne, geneTwo),
                                               #Species = Species,
                                               pool = pool,
                                               Tissues = c(tissueOneNow, tissueTwoNow),
                                               Sample_Type = c(sampleOne, sampleTwo),
                                               useBlackList = FALSE)
  VSTDF <- data.table::rbindlist(VSTList, idcol = "Group")
  if (tissueOne == "all" & tissueTwo == "all" &
      geneOne == geneTwo) {
    colnames(VSTDF)[3] <- "value"
    VSTDF$tissue <- gsub(VSTDF$Group,
                         pattern = "(.*) - (.*)",
                         replacement = "\\1")
    VSTDF$sample <- gsub(VSTDF$Group,
                         pattern = "(.*) - (.*)",
                         replacement = "\\2")
    goodTiss <- unique(VSTDF$tissue[which(VSTDF$sample == "cancer")])
    goodTiss2 <- unique(VSTDF$tissue[which(VSTDF$sample == "normal")])
    goodTissFinal <- goodTiss[which(goodTiss %in% goodTiss2)]
    VSTDF2 <- VSTDF[which(VSTDF$tissue %in% goodTissFinal),]
    VSTDF2 <- VSTDF2[VSTDF2$sample != "all",]
    maxHeight <- max(VSTDF2$value)
    meds <- sapply(unique(VSTDF2$tissue), FUN = function(groupNow) {
      stats::median(VSTDF2$value[VSTDF2$tissue == groupNow &
                                   VSTDF2$sample == "cancer"]) -
        stats::median(VSTDF2$value[VSTDF2$tissue == groupNow &
                                     VSTDF2$sample == "normal"])
    })
    meds <- meds[order(abs(meds))]
    VSTBPproto <- ggpubr::ggboxplot(data = VSTDF2,
                                    x = "tissue", order = names(meds),
                                    title = paste0(geneOne, " expression"),
                                    ylab = "Expression (VST counts)",
                                    fill = "sample", legend = "right",
                                    y = "value") +
      ggpubr::stat_compare_means(ggplot2::aes_string(group = "sample"),
                                 label.y = (maxHeight * 1.15),
                                 hide.ns = TRUE, label = "p.signif")
    VSTplot <- VSTBPproto +
      ggpubr::rotate_x_text(angle = 45) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_text(size = 16),
        title = ggplot2::element_text(size = 22),
        plot.margin = ggplot2::margin(10, 25, 10, 25)
      ) +
      ggpubr::rremove("xlab") + ggpubr::rremove("legend.title")
    VSTDFFinal <- VSTDF2
  } else {
    if (tissueOne != "all") {
      VSTDFOne <- VSTDF[grep(VSTDF$Group,
                             pattern = gsub(tissueOne,
                                            pattern = "_",
                                            replacement = " ")),]
    } else {
      VSTDFOne <- VSTDF
    }

    VSTDFOne <- VSTDFOne[grep(VSTDFOne$Group,
                              pattern = sampleOne),]

    if (tissueTwo != "all") {
      VSTDFTwo <- VSTDF[grep(VSTDF$Group,
                             pattern = gsub(tissueTwo,
                                            pattern = "_",
                                            replacement = " ")),]
    } else {
      VSTDFTwo <- VSTDF
    }
    VSTDFTwo <- VSTDFTwo[grep(VSTDFTwo$Group,
                              pattern = sampleTwo),]

    if (geneOne != geneTwo) {
      uiNameOne <- paste0(
        stringr::str_to_title(gsub(tissueOne,
                                   pattern = "_",
                                   replacement = " ")), "-",
        stringr::str_to_title(sampleOne))
      uiNameTwo <- paste0(
        stringr::str_to_title(gsub(tissueTwo,
                                   pattern = "_",
                                   replacement = " ")), "-",
        stringr::str_to_title(sampleTwo))
      VSTDFOne <- VSTDFOne[,c(-4)]
      VSTDFOne$Gene <- colnames(VSTDFOne)[3]
      VSTDFTwo <- VSTDFTwo[,c(-3)]
      VSTDFTwo$Gene <- colnames(VSTDFTwo)[3]
      colnames(VSTDFTwo)[3] <- "VST"
      colnames(VSTDFOne)[3] <- "VST"
      titleStr <- paste0(geneOne, " vs. ",
                         geneTwo, " expression")
      if (uiNameOne != uiNameTwo) {
        capStr <- paste0("In ", tolower(uiNameOne), " vs ", tolower(uiNameTwo), " samples.\n(*VST may be inaccurate for comparing different genes. TPM will be used in future versions.)")
      } else {
        capStr <- paste0("In ", tolower(uiNameOne), " samples.\n(*VST may be inaccurate for comparing different genes. TPM will be used in future versions.)")
      }
      fillStr <- "Gene"
      VSTDFFinal <- rbind(VSTDFOne, VSTDFTwo)
      VSTDFFinal$Group <- stringr::str_to_title(gsub(VSTDFFinal$Group, pattern = "_",
                                                     replacement = " - "))
      VSTplot <- ggpubr::ggboxplot(data = VSTDFFinal,
                                   x = fillStr, #facet.by = "Gene",
                                   title = titleStr,
                                   caption = capStr,
                                   ylab = "Expression (VST counts)",
                                   fill = fillStr,
                                   y = "VST") +
        ggpubr::rremove("legend") +
        ggpubr::rremove("xlab") +
        ggplot2::scale_y_continuous(limits = c(.95*min(VSTDFFinal$VST), 1.1*max(VSTDFFinal$VST))) +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = 20),
                       title = ggplot2::element_text(size = 22),
                       plot.margin = ggplot2::margin(10, 10, 10, 10),
                       axis.text.x = ggplot2::element_text(size = 16))
    } else {
      capStr <- NULL
      VSTDFOne$Gene <- paste0(colnames(VSTDFOne)[3], "_", VSTDFOne$Group)
      VSTDFTwo$Gene <- paste0(colnames(VSTDFTwo)[3], "_", VSTDFTwo$Group)
      colnames(VSTDFTwo)[3] <- "VST"
      colnames(VSTDFOne)[3] <- "VST"
      titleStr <- paste0(geneOne, " expression")
      fillStr <- "Group"
      VSTDFFinal <- unique(rbind(VSTDFOne, VSTDFTwo))
      VSTDFFinal$Group <- stringr::str_to_title(gsub(VSTDFFinal$Group, pattern = "_",
                                                     replacement = " - "))
      groupOrder <- unique(VSTDFFinal$Group[order(VSTDFFinal$Group)])
      VSTplot <- ggpubr::ggboxplot(data = VSTDFFinal,
                                   x = fillStr, #facet.by = "Gene",
                                   title = titleStr,
                                   caption = capStr, order = groupOrder,
                                   ylab = "Expression (VST counts)",
                                   fill = fillStr,
                                   y = "VST") +
        ggpubr::rremove("legend") +
        ggpubr::rremove("xlab") +
        ggplot2::scale_y_continuous(limits = c(.95*min(VSTDFFinal$VST), 1.1*max(VSTDFFinal$VST))) +
        ggpubr::rotate_x_text(45) +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = 20),
                       title = ggplot2::element_text(size = 22),
                       plot.margin = ggplot2::margin(10, 10, 10, 25),
                       axis.text.x = ggplot2::element_text(size = 16)) +
        ggpubr::stat_compare_means(comparisons = list(unique(as.data.frame(VSTDFFinal)[, which(colnames(VSTDFFinal) == fillStr)])),
                                   size = 5, method = "t.test", label = "p.signif")

    }

  }


  pairRes[["compared"]][["VST_boxPlot"]] <- VSTplot
  pairRes[["compared"]][["VST_Data"]] <- VSTDFFinal



  correlationsScatter <- correlations
  correlationsScatter$Gene <- F
  correlationsScatter$Gene[which(rownames(correlationsScatter) %in%
                                   c(geneOne, geneTwo))] <- T
  colnames(correlationsScatter)[c(1, 2)] <- c("x", "y")
  if (geneOne != geneTwo) {
    titleStr <- paste0(geneOne, " vs. ",
                       geneTwo, " correlations")
    geneOneTitleHeat <- geneOne
    geneTwoTitleHeat <- geneTwo
  } else {
    titleStr <- paste0(geneOne, " correlations")
    geneOneTitleHeat <- geneOneTitle
    geneTwoTitleHeat <- geneTwoTitle
  }
  gs <- ggpubr::ggscatter(correlationsScatter,
                          title = titleStr,
                          x = "x", y = "y",
                          ylab = geneTwoTitleGP,
                          xlab = geneOneTitleGP,
                          label.rectangle = FALSE,
                          size = .5,
                          repel = TRUE, label = rownames(correlationsScatter),
                          font.label = c(12, "bold", "black"),
                          label.select = unique(c(geneOne, geneTwo)),
                          add = "reg.line") +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 20),
                   title = ggplot2::element_text(size = 22),
                   plot.margin = ggplot2::margin(10, 10, 10, 10),
                   axis.text = ggplot2::element_text(size = 16)) +
    ggpubr::stat_cor(label.y = 1.2, size = 5)
  pairRes[["compared"]][["correlationPlot"]] <- gs

  df <- correlationsScatter
  df <- unique(df[,c(1,2)])
  # df <- pairResNow$correlations
  lm_eqn <- function(df){
    m <- stats::lm(eval(parse(text = colnames(df)[2])) ~ eval(parse(text = colnames(df)[1])), df)
    r <- sqrt(summary(m)$r.squared) * sign(unname(stats::coef(m)[2]))
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)~"="~r,
                     list(a = format(unname(stats::coef(m)[1]), digits = 2),
                          b = format(unname(stats::coef(m)[2]), digits = 2),
                          r = format(r, digits = 2)))
    as.character(as.expression(eq));
  }
  labb <- lm_eqn(df)
  gp <- ggplot2::ggplot(data = df,
                        mapping = ggplot2::aes_string(x = colnames(df)[1],
                                                      y = colnames(df)[2])) +
    ggplot2::stat_bin2d(bins = 150) +
    ggplot2::geom_smooth(colour="black", size = 1.25,
                         method='lm') +
    ggplot2::labs(title = titleStr) +
    ggplot2::ylab(geneTwoTitleGP) +
    ggplot2::xlab(geneOneTitleGP) +
    ggplot2::annotate("text", x = 0, y = 1.2,
                      label = labb, size = 5,
                      parse = TRUE) +
    ggpubr::theme_pubr() +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 20),
                   title = ggplot2::element_text(size = 22),
                   plot.margin = ggplot2::margin(10, 10, 10, 10),
                   axis.text = ggplot2::element_text(size = 16)) +
    ggplot2::theme(legend.position = "none")
  pairRes[["compared"]][["correlationPlotBin"]] <- gp



  breaks <- getPhBreaks(corHeatVar)
  phCorVar <- pheatmap::pheatmap(corHeatVar, silent = TRUE, angle_col = 45,
                                 breaks = breaks[[2]], # fontsize = 18, fontsize_row = 16, fontsize_col = 18,
                                 main = "Top differentially correlated genes",
                                 labels_col = c(geneOneTitleHeat, geneTwoTitleHeat),
                                 cluster_rows = TRUE, cluster_cols = FALSE)
  pairRes[["compared"]][["correlationVarianceHeatmap"]] <- ggplotify::as.ggplot(phCorVar)
  if (runGSEA) {
    breaks <- getPhBreaks(compHeatVar)
    phGSEAVar <- pheatmap::pheatmap(compHeatVar, silent = TRUE, angle_col = 45,
                                    breaks = breaks[[2]], # fontsize = 18, fontsize_row = 16, fontsize_col = 18,
                                    main = "Top differentially correlated pathways",
                                    labels_col = c(geneOneTitleHeat, geneTwoTitleHeat),
                                    cluster_rows = TRUE,
                                    cluster_cols = FALSE)
    pairRes[["compared"]][["pathwayVarianceHeatmap"]] <- ggplotify::as.ggplot(phGSEAVar)
  }

  breaks <- getPhBreaks(corHeatSim)
  phCorSim <- pheatmap::pheatmap(corHeatSim, silent = TRUE, angle_col = 45,
                                 breaks = breaks[[2]], # fontsize = 18, fontsize_row = 16, fontsize_col = 18,
                                 main = "Top similarly correlated genes",
                                 labels_col = c(geneOneTitleHeat, geneTwoTitleHeat),
                                 cluster_rows = TRUE, cluster_cols = FALSE)
  pairRes[["compared"]][["correlationSimilarityHeatmap"]] <- ggplotify::as.ggplot(phCorSim)
  if (runGSEA) {
    breaks <- getPhBreaks(compHeatSim)
    phGSEASim <- pheatmap::pheatmap(compHeatSim, silent = TRUE, angle_col = 45,
                                    breaks = breaks[[2]], # fontsize = 18, fontsize_row = 16, fontsize_col = 18,
                                    main = "Top similarly correlated pathways",
                                    labels_col = c(geneOneTitleHeat, geneTwoTitleHeat),
                                    cluster_rows = TRUE,
                                    cluster_cols = FALSE)
    pairRes[["compared"]][["pathwaySimilarityHeatmap"]] <- ggplotify::as.ggplot(phGSEASim)
    pairRes[["compared"]][["correlatedPathwaysDataFrame"]] <- compPaths

  }

  # Get the corrPlot
  to_sample <- ifelse(length(VSTDFFinal$samples) > 10000, 10000, length(VSTDFFinal$samples))
  VSTWide <- tidyr::pivot_wider(dplyr::filter(VSTDFFinal, samples %in% sample(samples, to_sample)),
                                values_from = VST, names_from = Gene)
  VSTWide <- dplyr::inner_join(VSTWide, y = correlationAnalyzeR::human_coldata, by = "samples")
  geneOne_Corr <- correlations[, geneOne, drop = FALSE]
  pDF <- apply(geneOne_Corr, MARGIN = 1:2, n = length(geneOne_Corr[,1]), FUN = function(x, n) {
    stats::dt(abs(x)/sqrt((1-x^2)/(n-2)), df = 2)
  })
  padj <- p.adjust(pDF[,1], method = "BH")
  Rval <- geneOne_Corr[geneTwo,]
  Padj <- padj[geneTwo]
  plt1 <- ggplot2::ggplot(VSTWide, aes_string(x = geneOne, y = geneTwo,
                                              group="Group",
                                              text = "samples")) +
    ggplot2::geom_point(alpha = .5) +
    ggplot2::labs(title = titleStr,
                  subtitle = paste0("Pearson's R = ", round(Rval, 3),
                                    " (padj = ", signif(Padj, 3), ")")) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::xlab(paste0(geneOne, " Expression (VST)")) +
    ggplot2::ylab(paste0(geneTwo, " Expression (VST)"))
  plt2 <- ggplot2::ggplot(VSTWide, aes_string(x = geneOne, y = geneTwo,
                                              group="Group", color = "Tissue",
                                              text = "samples")) +
    ggplot2::geom_point(alpha = .5) +
    ggplot2::labs(title = titleStr,
                  subtitle = paste0("Pearson's R = ", round(Rval, 3),
                                    " (padj = ", signif(Padj, 3), ")")) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::xlab(paste0(geneOne, " Expression (VST)")) +
    ggplot2::ylab(paste0(geneTwo, " Expression (VST)"))
  plt3 <- ggplot2::ggplot(VSTWide, aes_string(x = geneOne, y = geneTwo,
                                              group="Group", color = "disease",
                                              text = "samples")) +
    ggplot2::geom_point(alpha = .5) +
    ggplot2::labs(title = titleStr,
                  subtitle = paste0("Pearson's R = ", round(Rval, 3),
                                    " (padj = ", signif(Padj, 3), ")")) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::scale_color_manual(name = "Disease", values = c("Cancer" = "firebrick",
                                                             "Normal" = "forestgreen")) +
    ggplot2::xlab(paste0(geneOne, " Expression (VST)")) +
    ggplot2::ylab(paste0(geneTwo, " Expression (VST)"))


  pairRes[["compared"]][["VST_corrPlot"]] <- list(
    "corrPlot_all" = plt1,
    "corrPlot_tissue" = plt2,
    "corrPlot_disease" = plt3,
    "Rval" = Rval,
    "Padj" = Padj,
    "corrPlot_VST_data" = VSTWide
  )

  # Return results
  return(pairRes)
}






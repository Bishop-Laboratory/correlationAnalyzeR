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
#' @param Species Species to obtain gene names for.
#'     Either 'hsapiens' or 'mmusculus'. Default: "hsapiens".
#' @param GSEA_Type Whether GSEA should consider all msigdb annotations,
#'     or just those in the most popular categories. Should be one of either
#'     'simple' or 'complex'.
#' @param nperm Number of permutations to run in GSEA. Default is 2000
#' @param sampler If TRUE, will only return 100,000 random genesets from either
#' simple or complex TERM2GENEs. Useful for reducing GSEA computational burden.
#' @param crossCompareMode Use this mode to generate comparisons
#' across all tissue and disease types. If both genes for genesOfInterest are the
#' same -- will compare normal vs cancer for that gene in each available tissue. Else, will
#' perform comparison of two different genes in all tissue-disease groups.
#' Will only consider user input for returnDataOnly, outputPrefix, Species, and genesOfInterest.
#' @param outputPrefix Prefix for saved files -- the directory name to store output files in. If
#' folder does not exist, it will be created.
#' @param runGSEA If TRUE will run GSEA using gene correlation values.
#' @param returnDataOnly if TRUE will return result list object
#' and will not generate any folders or files.
#' @param topPlots Logical. If TRUE, myGSEA() will build gsea plots for top correlated genesets.
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied.
#' @return A named list containing visualizations and correlation data from paired analysis.
#' @examples
#' genesOfInterest <- c("ATM", "SLC7A11")
#' correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#'                               Species = "hsapiens",
#'                               GSEA_Type = "simple", returnDataOnly = TRUE,
#'                               Sample_Type = c("normal", "normal"),
#'                               Tissue = c("brain", "brain"))
#' genesOfInterest <- c("BRCA1", "BRCA1")
#' correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#'                               Species = "hsapiens",
#'                               GSEA_Type = "simple", returnDataOnly = TRUE,
#'                               Sample_Type = c("normal", "cancer"),
#'                               Tissue = c("respiratory", "respiratory"))
#' genesOfInterest <- c("NFKB1", "SOX10")
#' correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#'                               Species = "hsapiens", returnDataOnly = TRUE,
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
                             Species = c("hsapiens", "mmusculus"),
                             GSEA_Type = c("simple", "complex"),
                             outputPrefix = "CorrelationAnalyzeR_Output_Paired",
                             crossCompareMode = FALSE,
                             runGSEA = TRUE, nperm = 2000, sampler = FALSE,
                             topPlots = FALSE, returnDataOnly = FALSE, pool = NULL) {

  # genesOfInterest <- c("BRCA1", "BRCA1")
  # Species <- "hsapiens"
  # crossCompareMode = FALSE
  # returnDataOnly = TRUE
  # returnDataOnly = TRUE
  # outputPrefix = "CorrelationAnalyzeR_Output_Paired"
  # runGSEA = TRUE
  # topPlots = FALSE
  # Sample_Type = c("normal", "normal")
  # Tissue = c("respiratory", "female0reproductive")
  # GSEA_Type = c("simple", "complex")
  # pool = NULL

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
    availTissue <- correlationAnalyzeR::getTissueTypes(Species = Species,
                                                       pool = pool,
                                                       useBlackList = TRUE)
    runGSEA <- F
    if (genesOfInterest[1] == genesOfInterest[2]) {
      cat("\nGene one is the same as gene two ... \n")
      if (Species == "mmusculus") {
        stop("\nOnly normal tissues available for mouse",
             " due to black-listing of cancer groups for low quality.",
             "\nEnter two different genes or choose human ... \n")
      }
      cat("\nWill perform normal vs cancer comparison on",
          genesOfInterest[1], "... \n")
      mode <- "cross_normalVsCancer"
      df <- as.data.frame(table(gsub(availTissue,
                                     pattern = " - .*",
                                     replacement = "")), stringsAsFactors = FALSE)
      goodTissues <- df$Var1[which(df$Freq == 2)]
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
      availTissue <- strsplit(availTissue, split = " - ")
      Tissue <- vapply(availTissue, FUN = "[[", FUN.VALUE = "character", 1)
      genesVec <- rep(genesOfInterest, length(Tissue))
      Tissue <- rep(Tissue, each = 2)
      Sample_Type <- vapply(availTissue, FUN = "[[", FUN.VALUE = "character", 2)
      Sample_Type <- rep(Sample_Type, each = 2)
    }
    # Get TPM for each gene
    geneUnique <- unique(genesVec)
    geneTPMList <- correlationAnalyzeR::getTissueTPM(genesOfInterest = geneUnique,
                                                     Species = Species,
                                                     Tissues = "all",
                                                     pool = pool,
                                                     Sample_Type = "all",
                                                     useBlackList = TRUE)
    geneTPMDF <- data.table::rbindlist(geneTPMList, idcol = "group")
    rawGroup <- geneTPMDF$group
    rawGroup1 <- stringr::str_to_title(gsub(rawGroup,
                                            pattern = "_.*",
                                            replacement = ""))
    rawGroup2 <- gsub(rawGroup, pattern = ".*_", replacement = "")
    rawGroup2 <- stringr::str_to_title(gsub(rawGroup2,
                                            pattern = "0",
                                            replacement = " "))
    geneTPMDF <- cbind(paste0(rawGroup2,
                              " - ",
                              rawGroup1),
                       rawGroup2,
                       rawGroup1,
                       geneTPMDF[,c(-1)])
    colnames(geneTPMDF)[c(1:4)] <- c("Group", "Tissue", "sampleType", "Samples")
    availTissue <- correlationAnalyzeR::getTissueTypes(Species = Species,
                                                       pool = pool,
                                                       useBlackList = TRUE)
    availTissue <- gsub(availTissue, pattern = "0", replacement = " ")
    availTissue <- stringr::str_to_title(availTissue)
    # all(geneTPMDF$Group %in% availTissue) #-- should be TRUE
    geneTPMDF <- geneTPMDF[order(match(geneTPMDF$Group, availTissue)),]
    # Make TPM plot
    crossCompareResTPM <- list()
    if (mode == "cross_geneVsGene") {
      geneTPMDF[,5] <- log2(geneTPMDF[,5] + 1)
      geneTPMDF[,6] <- log2(geneTPMDF[,6] + 1)
      geneTPMDFToPlot <- geneTPMDF[,c(1,2, 3, 5, 6)]
      geneTPMDFToPlot <- geneTPMDFToPlot %>%
        gather("Gene", "TPM", -.data$Group, -.data$Tissue, -.data$sampleType)
      geneOne <- geneUnique[1]
      geneTwo <- geneUnique[2]
      if (Species == "mmusculus") {
        fillStr <- "Tissue"
      } else {
        fillStr <- "Group"
      }
      geneTPMDFToPlot1 <- geneTPMDFToPlot[which(geneTPMDFToPlot$Gene == geneOne),]
      plotOne <- ggpubr::ggboxplot(data = geneTPMDFToPlot1,
                                   x = fillStr, #facet.by = "Gene",
                                   title = paste0(geneOne,
                                                  " expression across tissues"),
                                   ylab = "log2(TPM + 1)",
                                   fill = fillStr,
                                   y = "TPM") +
        ggpubr::rotate_x_text() +
        ggpubr::rremove("legend") +
        ggpubr::rremove("xlab")
      geneTPMDFToPlot2 <- geneTPMDFToPlot[which(geneTPMDFToPlot$Gene == geneTwo),]
      plotTwo <- ggpubr::ggboxplot(data = geneTPMDFToPlot2,
                                   x = fillStr, #facet.by = "Gene",
                                   title = paste0(geneTwo,
                                                  " expression across tissues"),
                                   ylab = "log2(TPM + 1)",
                                   fill = fillStr,
                                   y = "TPM") +
        ggpubr::rotate_x_text() +
        ggpubr::rremove("legend") +
        ggpubr::rremove("xlab")
      colnames(geneTPMDF)[c(5:6)] <- paste0(colnames(geneTPMDF)[c(5:6)], "_log2TPM")
      crossCompareResTPM[["TPM_boxPlotOne"]] <- plotOne
      crossCompareResTPM[["TPM_boxPlotTwo"]] <- plotTwo

    } else {
      geneTPMDF[,5] <- log2(geneTPMDF[,5] + 1)
      geneTPMDFToPlot <- geneTPMDF
      colnames(geneTPMDFToPlot)[length(colnames(geneTPMDFToPlot))] <- "TPM"
      plot <- ggpubr::ggboxplot(data = geneTPMDFToPlot,
                                   x = "Group", #facet.by = "Gene",
                                   title = paste0(geneUnique[1],
                                                  " expression across tissues"),
                                   ylab = "log2(TPM + 1)",
                                   fill = "Group",
                                   y = "TPM") +
        ggpubr::rotate_x_text() +
        ggpubr::rremove("legend") +
        ggpubr::rremove("xlab")
      colnames(geneTPMDF)[c(5)] <- paste0(colnames(geneTPMDF)[c(5)], "_log2TPM")
      crossCompareResTPM[["TPM_boxPlot"]] <- plot
    }
    crossCompareResTPM[["TPM_DF"]] <- geneTPMDF

    # Do paired to get correlation data
    pairRes <- correlationAnalyzeR::analyzeSingleGenes(
      genesOfInterest = genesVec, pool = pool,
      returnDataOnly = returnDataOnly, topPlots = topPlots,
      outputPrefix = outputPrefix, runGSEA = FALSE,
      Sample_Type = Sample_Type, Tissue = Tissue,
      Species = Species, GSEA_Type = GSEA_Type
    )
    n <- length(names(pairRes))
    oldNames <- names(pairRes)[1:(n-1)]
    newNames <- gsub(oldNames, pattern = ", ", replacement = "_")
    newNames <- gsub(newNames, pattern = " - ", replacement = "_")
    newNames <- gsub(newNames, pattern = " ", replacement = "0")
    correlations <- pairRes$correlations
    colnames(correlations) <- newNames
    resList <- list()
    resList[["pairResList"]] <- list()
    for (i in 1:length(colnames(correlations))) {
      samp <- colnames(correlations)[i]
      if (i %% 2 == 0) {
        dfRaw <- correlations[,c((i-1), i)]
        df <- dfRaw
        tempList <- list()
        tissue <- stringr::str_match(colnames(df)[1], pattern = "_(.*)")[,2]
        tissue <- gsub(tissue, pattern = "_", replacement = " - ")
        tissue <- gsub(tissue, pattern = "0", replacement = " ")
        tissueSmall <- gsub(tissue, pattern = " - .*", replacement = "")
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
          ggplot2::annotate("text", x = 0, y = 1.1,
                            label = labb,
                            parse = TRUE) +
          ggpubr::theme_pubr() +
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
          ph <- pheatmap::pheatmap(dfPh, cluster_cols = FALSE,
                                   silent = TRUE, angle_col = 0, main = titleStr,
                                   labels_col = c(genesOfInterest[1],
                                                  genesOfInterest[2]))
        } else {
          ph <- pheatmap::pheatmap(dfPh, cluster_cols = FALSE,
                                   silent = TRUE, angle_col = 0, main = titleStr,
                                   labels_col = c("Normal", "Cancer"))
        }

        tempList[["heatMap"]] <- ph

        resList[["pairResList"]][[i/2]] <- tempList
        names(resList[["pairResList"]])[i/2] <- tissue
      }
    }

    correlations$average <- rowMeans(correlations)
    correlations$variance <- matrixStats::rowVars(as.matrix(correlations))
    correlations <- correlations[order(correlations$variance,
                                       decreasing = TRUE),]
    resList[["Correlations"]] <- correlations
    resList[["crossCompareTPM"]] <- crossCompareResTPM
    resList[["mode"]] <- mode
    return(resList)

  }


  # If running in normal mode ...
  if (length(genesOfInterest) == 2 ) {
    pairRes <- correlationAnalyzeR::analyzeSingleGenes(
      genesOfInterest = genesOfInterest, pool = pool,
      returnDataOnly = returnDataOnly, topPlots = topPlots,
      outputPrefix = outputPrefix,
      runGSEA = runGSEA, nperm = nperm, sampler = sampler,
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
                     yes = TRUE, no = FALSE)


  pairRes[["compared"]] <- list()
  # Variance heat map
  correlations$average <- rowMeans(correlations)
  correlations$variance <- matrixStats::rowVars(as.matrix(correlations[,c(1,2)]))
  correlations <- correlations[order(correlations$variance,
                                     decreasing = TRUE),]
  pairRes[["compared"]][["correlations"]] <- correlations
  cn <- colnames(correlations)
  correlations2 <- correlations[which(! rownames(correlations) %in%
                                        colnames(correlations)),]

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
  # Get TPM
  tissueOneNow <- gsub(tissueOne, pattern = " ", replacement = "0")
  tissueTwoNow <- gsub(tissueTwo, pattern = " ", replacement = "0")
  TPMList <- correlationAnalyzeR::getTissueTPM(genesOfInterest = c(geneOne, geneTwo),
                                               Species = Species, pool = pool,
                                               Tissues = c(tissueOneNow, tissueTwoNow),
                                               Sample_Type = c(sampleOne, sampleTwo),
                                               useBlackList = FALSE)
  TPMDF <- data.table::rbindlist(TPMList, idcol = "Group")
  TPMDF$Group <- gsub(TPMDF$Group, pattern = "0", replacement = " ")
  TPMDFOne <- TPMDF[grep(TPMDF$Group,
                         pattern = paste0("_", tissueOne)),]
  TPMDFOne <- TPMDFOne[grep(TPMDFOne$Group,
                         pattern = paste0(sampleOne, "_")),]

  TPMDFTwo <- TPMDF[grep(TPMDF$Group,
                         pattern = paste0("_", tissueTwo)),]
  TPMDFTwo <- TPMDFTwo[grep(TPMDFTwo$Group,
                            pattern = paste0(sampleTwo, "_")),]
  if (geneOne != geneTwo) {
    uiNameOne <- paste0(
      stringr::str_to_title(tissueOne), "-",
      stringr::str_to_title(sampleOne))
    uiNameTwo <- paste0(
      stringr::str_to_title(tissueTwo), "-",
      stringr::str_to_title(sampleTwo))
    TPMDFOne <- TPMDFOne[,c(-4)]
    TPMDFOne$Gene <- colnames(TPMDFOne)[3]
    TPMDFTwo <- TPMDFTwo[,c(-3)]
    TPMDFTwo$Gene <- colnames(TPMDFTwo)[3]
    TPMDFOne[,3] <- log2(TPMDFOne[,3] + 1)
    TPMDFTwo[,3] <- log2(TPMDFTwo[,3] + 1)
    colnames(TPMDFTwo)[3] <- "Log2TPM"
    colnames(TPMDFOne)[3] <- "Log2TPM"
    titleStr <- paste0(geneOne, " vs. ",
                       geneTwo, " expression")
    if (uiNameOne != uiNameTwo) {
      capStr <- paste0("In ", tolower(uiNameOne), " vs ", tolower(uiNameTwo), " samples")
    } else {
      capStr <- paste0("In ", tolower(uiNameOne), " samples")
    }
    fillStr <- "Gene"
  } else {
    capStr <- NULL
    TPMDFOne$Gene <- paste0(colnames(TPMDFOne)[3], "_", TPMDFOne$Group)
    TPMDFTwo$Gene <- paste0(colnames(TPMDFTwo)[3], "_", TPMDFTwo$Group)
    TPMDFOne[,3] <- log2(TPMDFOne[,3] + 1)
    TPMDFTwo[,3] <- log2(TPMDFTwo[,3] + 1)
    colnames(TPMDFTwo)[3] <- "Log2TPM"
    colnames(TPMDFOne)[3] <- "Log2TPM"
    titleStr <- paste0(geneOne, " expression")
    fillStr <- "Group"
  }
  TPMDFFinal <- rbind(TPMDFOne, TPMDFTwo)
  TPMDFFinal$Group <- stringr::str_to_title(gsub(TPMDFFinal$Group, pattern = "_",
                                                     replacement = " - "))
  TPMplot <- ggpubr::ggboxplot(data = TPMDFFinal,
                               x = fillStr, #facet.by = "Gene",
                               title = titleStr,
                               caption = capStr,
                               ylab = "log2(TPM + 1)",
                               fill = fillStr,
                               y = "Log2TPM") +
    ggpubr::rremove("legend") +
    ggpubr::rremove("xlab")
  pairRes[["compared"]][["TPM_boxPlot"]] <- TPMplot
  pairRes[["compared"]][["TPM_Data"]] <- TPMDFFinal
  if (longName) {
    uiNameOne <- paste0(
                        stringr::str_to_title(tissueOne), " - ",
                        stringr::str_to_title(sampleOne))
    uiNameTwo <- paste0(
                        stringr::str_to_title(tissueTwo), " - ",
                        stringr::str_to_title(sampleTwo))
    correlationsScatter <- correlations
    correlationsScatter$Gene <- F
    correlationsScatter$Gene[which(rownames(correlationsScatter) %in%
                                              c(geneOne, geneTwo))] <- T
    colnames(correlationsScatter)[c(1, 2)] <- c("x", "y")
    if (geneOne != geneTwo) {
      titleStr <- paste0(geneOne, " vs. ",
                         geneTwo, " correlations")
    } else {
      titleStr <- paste0(geneOne, " correlations")
    }
    gs <- ggpubr::ggscatter(correlationsScatter,
                            title = titleStr,
                            x = "x", y = "y",
                            ylab = geneTwoTitle,
                            xlab = geneOneTitle,
                            label.rectangle = FALSE,
                            size = .5,
                            repel = TRUE, label = rownames(correlationsScatter),
                            font.label = c(12, "bold", "black"),
                            label.select = unique(c(geneOne, geneTwo)),
                            add = "reg.line") +
      ggpubr::stat_cor(label.y = 1.1)
    pairRes[["compared"]][["correlationPlot"]] <- gs

    phCor <- pheatmap::pheatmap(corHeat, silent = TRUE, angle_col = 0,
                                main = "Differentially correlated genes",
                                labels_col = c(geneOneTitle, geneTwoTitle),
                                cluster_rows = TRUE, cluster_cols = FALSE)
    pairRes[["compared"]][["correlationVarianceHeatmap"]] <- phCor
    phGSEA <- pheatmap::pheatmap(compHeat, silent = TRUE, angle_col = 0,
                                 main = "Differentially correlated pathways",
                                 labels_col = c(geneOneTitle, geneTwoTitle),
                                 cluster_rows = TRUE,
                                 cluster_cols = FALSE)
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
    correlationsScatter <- correlations
    correlationsScatter$Gene <- F
    correlationsScatter$Gene[which(rownames(correlationsScatter) %in%
                                              c(geneOne, geneTwo))] <- T
    gs <- ggpubr::ggscatter(correlationsScatter,
                            title = paste0(geneOne, " vs. ",
                                           geneTwo, " Correlations"),
                            x = geneOne, y = geneTwo,
                            ylab = geneTwoTitle,
                            xlab = geneOneTitle,
                            label.rectangle = FALSE,
                            size = .5,
                            repel = TRUE, label = rownames(correlationsScatter),
                            font.label = c(12, "bold", "black"),
                            label.select = unique(c(geneOne, geneTwo)),
                            add = "reg.line") +
      ggpubr::stat_cor(label.y = 1.1)
    pairRes[["compared"]][["correlationPlot"]] <- gs

    phCor <- pheatmap::pheatmap(corHeat,angle_col = 0,
                                silent = TRUE,
                                main = "Differentially correlated genes",
                                cluster_rows = TRUE, cluster_cols = FALSE)
    pairRes[["compared"]][["correlationVarianceHeatmap"]] <- phCor
    phGSEA <- pheatmap::pheatmap(compHeat, silent = TRUE,
                                 main = "Differentially correlated pathways",
                                 cluster_rows = TRUE, angle_col = 0,
                                 labels_col = c(geneOne, geneTwo),
                                 cluster_cols = FALSE)
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






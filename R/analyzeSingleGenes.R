#' Analyze Single Genes
#'
#' Obtains correlations and corGSEA results for each gene of interest.
#'
#' @param genesOfInterest A vector of genes to analyze.
#'
#' @param Sample_Type Type of RNA Seq samples used to create correlation data.
#' Either "all", "normal", or "cancer". Can be a single value for all genes,
#' or a vector corresponding to genesOfInterest. Default: "normal"
#'
#' @param Tissue Which tissue type should gene correlations be derived from?
#' Can be a single value for all genes, or a vector corresponding to genesOfInterest.
#' Run getTissueTypes() to see available tissues. Default: "all"
#'
#' @param runGSEA If TRUE will run GSEA using gene correlation values. Default: TRUE.
#'
#' @param GSEA_Type Character vector listing the gene set databases to use.
#' Options are listed in correlationAnalyzeR::pathwayCategories --
#' See details of ?getTERM2GENE for more info.
#'
#' @param TERM2GENE Mapping of geneset IDs to gene names. If not supplied, it will be
#' generated automatically. Only applicable if GSEA is to be run. TERM2GENE objects
#' can be generated manually using the getTERM2GENE() function.
#'
#' @param nperm Number of permutations to run in GSEA. Default: 2000
#'
#' @param crossCompareMode Instead of normal single gene analysis, analyzeSingleGenes() will generate
#' correlations for a single gene across all tissue-disease groups. GSEA will not be run.
#' analyzeSingleGenes() will only consider user input for returnDataOnly, whichCompareGroups,
#' outputPrefix, and genesOfInterest.
#'
#' @param whichCompareGroups For crossCompareMode, select "all", "normal", or "cancer"
#' to analyze correlations from the corresponding groups. Default: "all".
#'
#' @param returnDataOnly if TRUE will return result list object
#' and will not generate any folders or files. Default: TRUE.
#'
#' @param topPlots Logical. If TRUE, myGSEA() will build gsea plots for top correlated genesets.
#' Default: TRUE.
#'
#' @param corrMat A custom correlation matrix generated by generateCorrelations()
#' to use instead of pre-supplied databases. If supplied, "Tissue" and "Sample_Type" are ignored.
#'
#' @param corrMat_label If corrMat is provided, this label will be used for plotting. Default: "User-Supplied".
#'
#' @param pool an object created by pool::dbPool to accessing SQL database.
#' It will be created if not supplied as long as makePool is TRUE.
#'
#' @param makePool Logical. Should a database pool be created if one is not supplied? Default: FALSE.
#'
#' @param outputPrefix Prefix for saved files -- the directory name to store output files in.
#' This is ignored unless returnDataOnly is FALSE. Default: "CorrelationAnalyzeR_Output"
#'
#' @param sampler Logical. If TRUE, will only return 100,000 random genesets from either
#' simple or complex TERM2GENEs. Useful for reducing GSEA computational burden.
#'
#' @return A named list of correlation values, corGSEA results,
#' and visualizations for each gene of interest.
#'
#' @details analyzeSingleGenes() performs most of the core tasks for analyzing gene function
#' via co-expression correlations. Please view the vignette for more detail about this function,
#' including the structure of the ouput data list.
#'
#' @examples
#' genesOfInterest <- c("ATM", "SLC7A11")
#' res <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = genesOfInterest,
#'                               returnDataOnly = TRUE,
#'                               GSEA_Type = "simple",
#'                               Sample_Type = c("normal", "cancer"),
#'                               Tissue = c("respiratory", "pancreas"))
#'
#' genesOfInterest <- c("BRCA1")
#' res <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = genesOfInterest,
#'                               GSEA_Type = "simple", returnDataOnly = TRUE,
#'                               crossCompareMode = TRUE,
#'                               whichCompareGroups = "normal")
#'
#' @export
analyzeSingleGenes <- function(genesOfInterest,
                               # Species = c("hsapiens", "mmusculus"),
                               GSEA_Type = c("simple"),
                               Sample_Type = "normal",
                               Tissue = "all",
                               crossCompareMode = FALSE,
                               nperm = 2000,
                               TERM2GENE = NULL,
                               whichCompareGroups = c("all", "normal", "cancer"),
                               outputPrefix = "CorrelationAnalyzeR_Output",
                               sampler = FALSE, runGSEA = TRUE,
                               topPlots = TRUE, returnDataOnly = TRUE,
                               pool = NULL,
                               corrMat = NULL,
                               corrMat_label = "User-Supplied",
                               makePool = FALSE) {

  # # Bug testing
  # genesOfInterest <- c("SRY")
  # GSEA_Type = c("simple")
  # Sample_Type = "normal"
  # Tissue = "all"
  # crossCompareMode = FALSE
  # nperm = 2000
  # TERM2GENE = NULL
  # whichCompareGroups = c("all", "normal", "cancer")
  # outputPrefix = "CorrelationAnalyzeR_Output"
  # sampler = FALSE
  # runGSEA = TRUE
  # topPlots = TRUE
  # returnDataOnly = TRUE
  # pool = NULL
  # makePool = FALSE
  # corrMat = corrMat
  # corrMat_label = "User-Supplied"


  if (! is.null(corrMat) & crossCompareMode) {
    stop(paste0("Currently analyzeSingleGenes() does not support custom ",
                "datasets via the corrMat argument in crossCompareMode.",
                " Please set crossCompareMode = FALSE or do not specify the",
                " corrMat argument."))
  }

  Species = "hsapiens"

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


  # Parse arguments
  geneString <- paste(genesOfInterest, collapse = ", ")

  # Create output folder
  if (! dir.exists(outputPrefix) & ! returnDataOnly) {
    dir.create(outputPrefix)
  }
  # Load appropriate TERM2GENE file built from msigdbr()
  if (is.null(TERM2GENE)) {
    if (Species[1] %in% c("hsapiens", "mmusculus")) {
      if (runGSEA) {
        cat("\nRetrieving GSEA annotations...")
        TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = GSEA_Type,
                                                       sampler = sampler,
                                                       Species = Species)
      }
    } else {
      stop("\ncorrelationAnalyzeR currently supports only Human and Mouse data.
         Please select either 'hsapiens' or 'mmusculus' for Species parameter.
         \n")
    }
  }

  # Check genes to make sure they exist
  avGenes <- correlationAnalyzeR::getAvailableGenes(#Species = Species,
                                                    pool = pool)

  badGenes <- genesOfInterest[which(! genesOfInterest %in% avGenes)]
  if (length(badGenes) > 0) {
    stop(paste0("\n\t\t\t'", paste(badGenes, collapse = ", "), "' not found
                in correlation data.
                Please check available data with getAvailableGenes().
                Your gene(s) of interest may have an updated name or
                have a species-specific identifier.\n"))
  }
  cat("\nRetrieving any missing correlation data...\n")

  if (length(Tissue) == 1) {
    Tissue <- rep(Tissue, length(genesOfInterest))
  } else if (length(Tissue) > 1) {
    if (length(Tissue) != length(genesOfInterest)) {
      stop("Number of valid genes not equal
           to length of supplied Tissue vector.
           Tissue vector should include 1 entry or the number of
           entries equal to the number of genesOfInterest.")
    }
  }
  if (length(Sample_Type) == 1) {
    Sample_Type <- rep(Sample_Type, length(genesOfInterest))
  } else if (length(Sample_Type) > 1) {
    if (length(Sample_Type) != length(genesOfInterest)) {
      stop("Number of valid genes not equal
           to length of supplied Sample_Type vector.
           Sample_Type vector should include 1 entry or the number of
           entries equal to the number of genesOfInterest.")
    }
  }

  if (crossCompareMode) {
    cat("\nRunning cross comparison mode ... \n")
    if (Species[1] == "mmusculus") {
      whichCompareGroups <- "normal"
      cat("\nOnly normal tissue comparisons available for mouse",
          " due to black-listing of cancer groups.\n")
      cat("\nContinuing with normal tissues ... \n")
    }
    availTissue <- correlationAnalyzeR::getTissueTypes(#Species = Species,
                                                       pool = pool,
                                                       useBlackList = TRUE)
    runGSEA <- F

    whichCompareGroups <- whichCompareGroups[1]
    if(whichCompareGroups != "all") {
      availTissue <- availTissue[grep(availTissue,
                                      pattern = whichCompareGroups)]
    } else {
      availTissue <- availTissue[grep(availTissue,
                                      pattern = "all", invert = T)]
    }
    availTissue <- strsplit(availTissue, split = " - ")
    Tissue <- vapply(availTissue, FUN = "[[", FUN.VALUE = "character", 1)
    genesVec <- rep(genesOfInterest, each = length(Tissue))
    Tissue <- rep(Tissue, length(genesOfInterest))
    Sample_Type <- vapply(availTissue, FUN = "[[", FUN.VALUE = "character", 2)
    Sample_Type <- rep(Sample_Type, length(genesOfInterest))
    genesOfInterest <- genesVec

  }

  # Call downloadData to get all required files
  corrDF <- suppressWarnings(correlationAnalyzeR::getCorrelationData(#Species = Species,
                                                    Tissue = Tissue, pool = pool,
                                                    Sample_Type = Sample_Type,
                                                    corrMat = corrMat,
                                                    geneList = genesOfInterest
                                                    ))
  if (is.null(corrMat)) {
    geneNames <- correlationAnalyzeR::hsapiens_corrSmall_geneNames
    rownames(corrDF) <- geneNames
  }

  if(crossCompareMode) {
    tissue2 <- gsub(Tissue, pattern = "0", replacement = " ")
    tissue2 <- stringr::str_to_title(tissue2)
    if(whichCompareGroups != "all") {
      namesVec <- tissue2
    } else {
      namesVec <- paste0(tissue2, " - ", stringr::str_to_title(Sample_Type))
    }
    topName <- paste0(genesOfInterest, "_", Tissue, "_", Sample_Type)
    resList <- list()
    geneList <- unique(genesOfInterest)
    for (i in 1:length(geneList)) {
      geneNow <- geneList[i]
      resList[[i]] <- list()
      names(resList)[i] <- geneNow
      inds <- which(genesOfInterest == geneNow)
      newDF <- corrDF[,inds]
      newDF[,which(is.na(newDF[1,]))] <- 0
      topNameNow <- topName[inds]
      colnames(newDF) <- topNameNow
      namesVecNow <- namesVec[inds]
      resList[[i]][["correlations"]] <- newDF
      newDF <- newDF[which(rownames(newDF) != geneNow),]
      newDFNorm <- preprocessCore::normalize.quantiles(as.matrix(newDF))
      newDFNorm <- as.data.frame(newDFNorm)
      rownames(newDFNorm) <- rownames(newDF)
      colnames(newDFNorm) <- colnames(newDF)
      newDFNorm$Variance <- matrixStats::rowVars(as.matrix(newDFNorm))
      n <- length(colnames(newDFNorm))
      newDFNorm$Average <- rowMeans(newDFNorm[,c(1:(n-1))])
      newDFNorm21 <- newDFNorm[order(newDFNorm$Average, decreasing = TRUE),]
      newDFNorm21 <- newDFNorm21[c(1:1500),]
      newDFNorm21 <- newDFNorm21[order(newDFNorm21$Variance, decreasing = FALSE),]
      n <- length(colnames(newDFNorm21))
      pMatCo1 <- newDFNorm21[c(1:15),c(1:(n-2))]
      pMatCoBig1 <- newDFNorm21[c(1:250),c(1:(n-2))]
      newDFNorm22 <- newDFNorm[order(newDFNorm$Average, decreasing = FALSE),]
      newDFNorm22 <- newDFNorm22[c(1:1500),]
      newDFNorm22 <- newDFNorm22[order(newDFNorm22$Variance, decreasing = FALSE),]
      n <- length(colnames(newDFNorm22))
      pMatCo2 <- newDFNorm22[c(1:15),c(1:(n-2))]
      pMatCoBig2 <- newDFNorm22[c(1:250),c(1:(n-2))]
      pMatCo <- unique(rbind(pMatCo1, pMatCo2))
      pMatCoBig <- unique(rbind(pMatCoBig1, pMatCoBig2))
      newDFNorm <- newDFNorm[order(newDFNorm$Variance, decreasing = TRUE),]
      pMatVar <- newDFNorm[c(1:30),c(1:(n-2))]
      pMatVarBig <- newDFNorm[c(1:500),c(1:(n-2))]
      titleName <- ifelse(whichCompareGroups == "all", paste0(geneNow,
                                                              " correlations across conditions"),
                          no = ifelse(whichCompareGroups == "normal",
                                      paste0(geneNow,
                                             " correlations across normal tissues"),
                                      no = paste0(geneNow,
                                                  " correlations across tumor tissues")))
      titleNameExp <- ifelse(whichCompareGroups == "all", paste0(geneNow,
                                                              " expression across conditions"),
                          no = ifelse(whichCompareGroups == "normal",
                                      paste0(geneNow,
                                             " expression in normal tissues"),
                                      no = paste0(geneNow,
                                                  " expression in tumor tissues")))
      breaks <- getPhBreaks(pMatVar)
      phSmallVar <- pheatmap::pheatmap(pMatVar, silent = TRUE,
                                    angle_col = 45, breaks = breaks[[2]],
                                    labels_col = namesVecNow)
      breaks <- getPhBreaks(pMatCo)
      phSmallCo <- pheatmap::pheatmap(pMatCo, silent = TRUE,
                                       angle_col = 45, breaks = breaks[[2]],
                                       labels_col = namesVecNow)
      breaks <- getPhBreaks(pMatVarBig)
      phBigVar <- pheatmap::pheatmap(pMatVarBig, silent = TRUE, angle_col = 45,
                                  show_rownames = FALSE, breaks = breaks[[2]],
                                  labels_col = namesVecNow)
      breaks <- getPhBreaks(pMatCoBig)
      phBigCo <- pheatmap::pheatmap(pMatCoBig, silent = TRUE, angle_col = 45,
                                     show_rownames = FALSE, breaks = breaks[[2]],
                                     labels_col = namesVecNow)

      # Get VST for gene
      geneVSTList <- correlationAnalyzeR::getTissueVST(genesOfInterest = geneNow,
                                                       # Species = Species,
                                                       Tissues = "all", pool = pool,
                                                       Sample_Type = whichCompareGroups[1],
                                                       useBlackList = TRUE)
      # Make VST plot
      geneVSTDF <- data.table::rbindlist(geneVSTList, idcol = "group")
      colnames(geneVSTDF)[3] <- "value"
      if(whichCompareGroups != "all") {
        geneVSTDF <- geneVSTDF[grep(geneVSTDF$group,
                                        pattern = whichCompareGroups),]
      }
      geneVSTDF$group <- stringr::str_to_title(geneVSTDF$group)

      availTissue <- correlationAnalyzeR::getTissueTypes(#Species = Species,
                                                         pool = pool,
                                                         useBlackList = TRUE)
      availTissue <- stringr::str_to_title(availTissue)
      if (whichCompareGroups != "all") {
        if (whichCompareGroups == "normal") {
          geneVSTDF$group <- gsub(geneVSTDF$group,
                                  pattern = "(.*) - (.*)",
                                  replacement = "\\1")
        } else {
          geneVSTDF$group <- gsub(geneVSTDF$group,
                                  pattern = "(.*) - (.*)",
                                  replacement = "\\1")
        }
        meds <- sapply(unique(geneVSTDF$group), FUN = function(groupNow) {
          stats::median(geneVSTDF$value[geneVSTDF$group == groupNow])
        })
        meds <- meds[order(meds)]
        VSTBPproto <- ggpubr::ggboxplot(data = geneVSTDF,
                                        x = "group", order = names(meds),
                                        title = titleNameExp,
                                        ylab = "Expression (VST counts)",
                                        fill = "group",
                                        y = "value")
        VSTBP <- VSTBPproto +
          ggplot2::scale_fill_manual(values = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(length(meds))) +
          ggpubr::rotate_x_text(angle = 45) +
          ggplot2::theme(
            axis.title.y = ggplot2::element_text(size = 16),
            title = ggplot2::element_text(size = 22),
            plot.margin = ggplot2::margin(10, 25, 10, 25)
          ) +
          ggpubr::rremove("legend") +
          ggpubr::rremove("xlab")
      } else {
        geneVSTDF$tissue <- gsub(geneVSTDF$group,
                                pattern = "(.*) - (.*)",
                                replacement = "\\1")
        geneVSTDF$sample <- gsub(geneVSTDF$group,
                                pattern = "(.*) - (.*)",
                                replacement = "\\2")
        goodTiss <- unique(geneVSTDF$tissue[which(geneVSTDF$sample == "Cancer")])
        goodTiss2 <- unique(geneVSTDF$tissue[which(geneVSTDF$sample == "Normal")])
        goodTissFinal <- goodTiss[which(goodTiss %in% goodTiss2)]
        geneVSTDF2 <- geneVSTDF[which(geneVSTDF$tissue %in% goodTissFinal),]
        geneVSTDF2 <- geneVSTDF2[geneVSTDF2$sample != "All",]
        maxHeight <- max(geneVSTDF2$value)
        meds <- sapply(unique(geneVSTDF2$tissue), FUN = function(groupNow) {
          stats::median(geneVSTDF2$value[geneVSTDF2$tissue == groupNow &
                                    geneVSTDF2$sample == "Cancer"]) -
            stats::median(geneVSTDF2$value[geneVSTDF2$tissue == groupNow &
                                      geneVSTDF2$sample == "Normal"])
        })
        meds <- meds[order(abs(meds))]
        VSTBPproto <- ggpubr::ggboxplot(data = geneVSTDF2,
                                        x = "tissue", order = names(meds),
                                        title = titleNameExp[1],
                                        ylab = "Expression (VST counts)",
                                        fill = "sample", legend = "right",
                                        y = "value") +
          ggpubr::stat_compare_means(ggplot2::aes_string(group = "sample"),
                                     label.y = (maxHeight * 1.15),
                                     hide.ns = TRUE, label = "p.signif")
        VSTBP <- VSTBPproto +
          ggpubr::rotate_x_text(angle = 45) +
          ggplot2::theme(
            axis.title.y = ggplot2::element_text(size = 16),
            title = ggplot2::element_text(size = 22),
            plot.margin = ggplot2::margin(10, 25, 10, 25)
          ) +
          ggpubr::rremove("xlab") + ggpubr::rremove("legend.title")

      }


      colnames(geneVSTDF)[3] <- paste0(geneNow, "_VST")
      resList[[i]][["VST_DF"]] <- geneVSTDF
      resList[[i]][["VST_boxPlot"]] <- VSTBP
      resList[[i]][["heatmapSmallCo"]] <- phSmallCo
      resList[[i]][["heatmapSmallDataCo"]] <- pMatCo
      resList[[i]][["heatmapBigCo"]] <- phBigCo
      resList[[i]][["heatmapBigDataCo"]] <- pMatCoBig
      resList[[i]][["heatmapSmallVar"]] <- phSmallVar
      resList[[i]][["heatmapSmallDataVar"]] <- pMatVar
      resList[[i]][["heatmapBigVar"]] <- phBigVar
      resList[[i]][["heatmapBigDataVar"]] <- pMatVarBig
    }
    return(resList)
  }

  resList <- list()
  # Main code
  for (i in 1:length(colnames(corrDF))) {
    gene <- colnames(corrDF)[i]
    cat(paste0("\nAnalyzing: ", gene))
    # Create output folder for gene
    geneOutDir <- file.path(outputPrefix, gene)
    if (! dir.exists(geneOutDir) & ! returnDataOnly) {
      dir.create(geneOutDir)
    }
    # Remove gene from correlation values to build normal distribution
    vec <- corrDF[,i]
    names(vec) <- rownames(corrDF)
    vec <- vec[order(vec, decreasing = TRUE)]
    vec <- vec[c(-1)]
    # Make a histogram of gene correlations
    corrDF2 <- corrDF[which(rownames(corrDF) != gene),i, drop = F]
    geneOne <- genesOfInterest[i]

    if (is.null(corrMat)) {
      tissueOne <- Tissue[i]
      tissueOne <- gsub(tissueOne, pattern = "0", replacement = " ")

      sampleOne <- Sample_Type[i]
      geneOneTitle <- paste0(geneOne, ", ",
                             stringr::str_to_title(tissueOne),
                             " - ",
                             stringr::str_to_title(sampleOne))
      tissueSampleCaption <- paste0(stringr::str_to_title(tissueOne),
                                    " - ",
                                    stringr::str_to_title(sampleOne))
    } else {
      geneOneTitle <- paste0(geneOne, ", ", corrMat_label)
      tissueSampleCaption <- corrMat_label
    }


    p <- ggpubr::gghistogram(data = corrDF2, x = gene, y = "..count..",
                             bins = 100, ylab = "Frequency\n",
                             title = gene,
                             caption = tissueSampleCaption,
                             xlab = paste0(gene, " correlation values"))

    resList[[i]] <- list()
    names(resList)[i] <- geneOneTitle
    geneOneTitleFile <- gsub(pattern = ",| |-", replacement = "", x = geneOneTitle)
    resList[[geneOneTitle]][["corrHist"]] <- p
    if (! returnDataOnly) {
      ggplot2::ggsave(filename = file.path(geneOutDir,
                                           paste0(geneOneTitleFile, ".png")),
                      plot = p)
    }


    if (runGSEA) {
      # Perform GSEA
      cat("\nGSEA\n")

      resGSEA <- suppressWarnings(
        correlationAnalyzeR::myGSEA(ranks = vec, TERM2GENE = TERM2GENE,
                                    plotFile = paste0(geneOneTitleFile,
                                                      ".corrPathways"),
                                    outDir = geneOutDir,
                                    nperm = nperm,
                                    topPlots = topPlots,
                                    returnDataOnly = returnDataOnly,
                                    Condition = paste0(geneOneTitle,
                                                       ": Correlated Genes"))
      )
      resList[[geneOneTitle]][["GSEA"]] <- resGSEA
    }

  }
  resList[["correlations"]] <- corrDF
  pDF <- apply(corrDF, MARGIN = 1:2, n = length(corrDF[,1]), FUN = function(x, n) {

    ## Using R to Z conversion method
    # z <- 0.5 * log((1+x)/(1-x))
    # zse <- 1/sqrt(colLengths-3)
    # p <- min(pnorm(z, sd=zse), pnorm(z, lower.tail=F, sd=zse))*2

    # Using the t statistic method
    stats::dt(abs(x)/sqrt((1-x^2)/(n-2)), df = 2)

  })
  resList[["P values"]] <- as.data.frame(pDF)

  # Return results
  return(resList)
}

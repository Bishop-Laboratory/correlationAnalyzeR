#' Generate a correlation matrix from user-supplied data
#'
#' @param cts a gene count matrix where rownames are genes and colnames are sample IDs.
#'
#' @param transformed Boolean. Indicates whether data is already transformed using
#' VST or a similar approach. If TRUE, VST transformation with DESeq2 will not be
#' performed. Default: FALSE
#'
#' @param cores Numeric. Number of cores to use for calculating size factors.
#' NOTE: cores > 1 does not work on Windows. Detault: 1.
#'
#' @return Matrix with gene co-expression correlations.
#'
#' @details This function performs the same normalization and transformation steps on
#' a user-supplied dataset that were originally used to generate the data provided
#' in the pre-calculated databases used by this package. The resulting correlation
#' matrix can be supplied to analyzeSingleGenes() as an input.
#' NOTE: the resulting matrix is very large and will take up ~8 GB of memory.
#'
#' @examples
#'
#' if (! 'airway' in rownames(install.packages())) {
#'     if (!requireNamespace("BiocManager", quietly = TRUE))
#'         install.packages("BiocManager")
#'     BiocManager::install("airway")
#' }
#'
#' if (! 'EnsDb.Hsapiens.v86' in rownames(install.packages())) {
#'     if (!requireNamespace("BiocManager", quietly = TRUE))
#'         install.packages("BiocManager")
#'     BiocManager::install("EnsDb.Hsapiens.v86")
#' }
#'
#' if (! 'dplyr' in rownames(install.packages())) {
#'     install.packages("dplyr")
#' }
#'
#' data(airway)
#' cts <- assay(airway)
#' ens2gene <- ensembldb::select(EnsDb.Hsapiens.v86, keys = rownames(cts),
#'                               columns = c("SYMBOL"), keytype = "GENEID") %>%
#'   dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
#'   dplyr::inner_join(y = data.frame("GENEID" = rownames(cts)))
#'
#' cts <- cts[ens2gene$GENEID,]
#' rownames(cts) <- ens2gene$SYMBOL
#'
#' corrMat <- generateCorrelations(cts)
#'
#' @export
generateCorrelations <- function(cts,
                                 transformed=FALSE,
                                 cores=1) {

  # library(airway)
  # library(EnsDb.Hsapiens.v86)
  # library(dplyr)
  #
  # data(airway)
  # cts <- assay(airway)
  # ens2gene <- ensembldb::select(EnsDb.Hsapiens.v86, keys = rownames(cts),
  #                               columns = c("SYMBOL"), keytype = "GENEID") %>%
  #   dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
  #   dplyr::inner_join(y = data.frame("GENEID" = rownames(cts)))
  # cts <- cts[ens2gene$GENEID,]
  # rownames(cts) <- ens2gene$SYMBOL
  #
  # label <- 'my_data'
  # transformed <- FALSE
  # cores <- 1

  if (! transformed) {
    # Keep genes expressed in 10%+ of samples
    nonZeroCount <- apply(cts, 1, function(row) {
      return(sum(row > 0))
    })
    keepInd <- which(nonZeroCount > (length(colnames(cts)) * .1))
    cts <- cts[keepInd,]

    # Transform and normalize
    dds <- DESeq2::DESeqDataSetFromMatrix(cts, colData = data.frame(sampleID=colnames(cts)),
                                          design = ~1)
    # from https://support.bioconductor.org/p/62246/#62250
    inds <- rownames(cts)
    geoMeansList <- parallel::mclapply(inds, FUN = function(ind) {
      row <- cts[ind,]
      if (all(row == 0)) {
        0
      } else {
        exp(sum(log(row[row != 0]))/length(row))
      }
    }, mc.cores = cores)
    geoMeans <- unlist(geoMeansList)
    dds <- DESeq2::estimateSizeFactors(dds, geoMeans=geoMeans)
    vsd <- DESeq2::vst(dds)
    cts <- SummarizedExperiment::assay(vsd)
  }

  corrMat <- WGCNA::cor(x = t(cts), verbose = 1, nThreads = cores)

  return(corrMat)
}


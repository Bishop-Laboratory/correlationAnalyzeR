library(rhdf5)
library(DESeq2)
library(uwot)
library(Seurat)
library(jsonlite)

source("Scripts/helpers.R")

# Make directory tree
dir.create("Data", showWarnings = F)
dir.create("Data/ARCHS4_Download", showWarnings = F)

# Download files from ARCHS4
humanExpURL <- "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
mouseExpURL <- "https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5"
downList <- c(humanExpURL, mouseExpURL)
for (i in 1:length(downList)) {
  file <- downList[i]
  outFile <- file.path("Data/ARCHS4_Download", gsub(file, pattern = ".+/([a-zA-Z0-9_]+\\.h5)$",
                                                    replacement = "\\1"))
  if (! file.exists(outFile)) {
    download.file(file, destfile = outFile)
  }
}

fileList <- paste0("Data/ARCHS4_Download/",
                   gsub(downList[c(1)],
                        pattern = "https://s3.amazonaws.com/mssm-seq-matrix/",
                        replacement = ""))

cat("\n", timestamp2(), " Loading & filtering expression data...\n", sep = "")
if (! file.exists("Data/fullRawCountsFiltered.rda")) {
  # Load and filter expression data
  expression <- h5read(dataFile, "data/expression")
  samples <- h5read(dataFile, "meta/Sample_geo_accession")
  genes <- h5read(dataFile, name = "meta/genes")
  rownames(expression) <- genes
  colnames(expression) <- samples
  colDataFinal <- colDataFinal[colDataFinal$samples %in% colnames(expression),]
  expression <- expression[, which(colnames(expression) %in% colDataFinal$samples)]
  colDataFinal <- colDataFinal[order(match(colDataFinal$samples, colnames(expression))),]
  if (! all(colDataFinal$samples == colnames(expression))) {
    stop("ColData samples are not identical to colnames of expression data...",
         " Please email Code/generateFigures_logFile.txt to author if you find this error and/or submit issue on github.")
  }
  # # Keep genes expressed in 10%+ of samples
  # nonZeroCount <- apply(expression, 1, nonZeroSamps)
  # expression <- expression[which(nonZeroCount > (length(colnames(expression)) * .01)),]
  timestamp()
  cat("\nDone. Saving expression data...\n")
  save(expression, file = "Data/fullRawCountsFiltered.rda")
}


if (! file.exists("Data/vsd_for_corr.rda")) {
  load("Data/fullRawCountsFiltered.rda")
  colDataFinal <- colDataFinal[colDataFinal$samples %in% colnames(expression),]
  colDataFinal <- colDataFinal[order(match(colDataFinal$samples, colnames(expression))),]
  all(colnames(expression) == colDataFinal$samples)
  rownames(colDataFinal) <- colDataFinal$samples
  # Keep genes expressed in 10%+ of samples
  nonZeroCount <- apply(expression, 1, nonZeroSamps)
  keepInd <- which(nonZeroCount > (length(colnames(expression)) * .1))
  expression <- expression[keepInd,]

  # Transform and normalize
  dds <- DESeqDataSetFromMatrix(expression, colData = colDataFinal,
                                design = ~1)
  # from https://support.bioconductor.org/p/62246/#62250
  inds <- rownames(expression)
  geoMeansList <- mclapply(inds, FUN = function(ind) {
    row <- expression[ind,]
    if (all(row == 0)) {
      0
    } else {
      exp(sum(log(row[row != 0]))/length(row))
    }
  }, mc.cores = 25)
  geoMeans <- unlist(geoMeansList)
  dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
  vsd <- vst(dds)
  timestamp()
  cat("\nDone. Saving vst data...\n")
  save(vsd, file = "Data/vsd_for_corr.rda")
} else {
  cat("\nVSD found -- loading...\n")
  load("Data/vsd_for_corr.rda")
}
























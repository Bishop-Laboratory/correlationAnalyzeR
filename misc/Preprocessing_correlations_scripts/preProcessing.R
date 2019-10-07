library(ontologyIndex)
library(rhdf5)
library(preprocessCore)
library(sva)

# Download files from ARCHS4
humanExpURL <- "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
humanTxURL <- "https://s3.amazonaws.com/mssm-seq-matrix/human_transcript_v7.h5"
mouseExpURL <- "https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5"
mouseTxURL <- "https://s3.amazonaws.com/mssm-seq-matrix/mouse_transcript_v7.h5"
downList <- c(humanExpURL, humanTxURL, mouseExpURL, mouseTxURL)
topdir <- getwd()
archdatDir <- file.path(topdir, "Data/ARCHS4_Download/")
for (i in 1:length(downList)) {
  setwd(archdatDir)
  file <- downList[i]
  print(file)
  system(paste0("wget ", file), ignore.stdout = T, ignore.stderr = T, wait = F)
  setwd(topdir)
}

fileList <- paste0("Data/ARCHS4_Download/", 
                   gsub(downList[c(1, 3)], 
                        pattern = "https://s3.amazonaws.com/mssm-seq-matrix/", 
                        replacement = ""))

# Generate corrected counts for human tissues
dataFile <- fileList[1]
samples <- h5read(dataFile, "meta/Sample_geo_accession")
tissue <- h5read(dataFile, name = "meta/Sample_source_name_ch1")
genes <- h5read(dataFile, name = "meta/genes")
description <- h5read(dataFile, name = "meta/Sample_description")
series <- h5read(dataFile, name = "meta/Sample_series_id")
organism <- h5read(dataFile, name = "meta/Sample_organism_ch1")
institute <- h5read(dataFile, name = "meta/Sample_contact_institute")
extractProtocol <- h5read(dataFile, name = "meta/Sample_extract_protocol_ch1")
reads_aligned <- h5read(dataFile, name = "meta/reads_aligned")


colData <- data.frame(samples, tissue, series, organism, institute,
                      extractProtocol, reads_aligned)
patternSingle <- c("single cell", "single-cell", "smart seq", "smart-seq", "drop-seq", "drop seq",
                   "fluidigm", "tenX", "scRNASeq", "scRNA-Seq", "chromium")
colData <- colData[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                               x = colData$extractProtocol, ignore.case = T)),]
tDFSeries <- as.data.frame(table(colData$series), stringsAsFactors = F)
goodSeries <- tDFSeries$Var1[which(tDFSeries$Freq > 1)]
colData <- colData[which(colData$series %in% goodSeries),]

print("reading in expression")
dir.create("Data/Final_Expression_Matrices")
dir.create("Data/Final_Expression_Matrices/Human")

expression <- h5read(dataFile, "data/expression")
rownames(expression) <- genes
colnames(expression) <- samples
expression <- expression[, which(colnames(expression) %in% colData$samples)]

# Filter expression matrix
nonZeroSamps <- function(row) {
  return(sum(row > 0))
}
print("filtering")
cols <- colSums(expression)
expressionFilt <- expression[, which(cols > 5000000)]
nonZeroCount <- apply(expressionFilt, 1, nonZeroSamps)
expressionFilt <- expressionFilt[which(nonZeroCount > (length(colnames(expressionFilt)) * .10)),]
colData <- colData[which(colData$samples %in% colnames(expressionFilt)),]

print("normalizing")
# Normalize with log2 and quantile
normCounts <- log2(expressionFilt + 1)
normCounts <- normalize.quantiles(normCounts)
colnames(normCounts) <- colnames(expressionFilt)
rownames(normCounts) <- rownames(expressionFilt)
rm(expression)
rm(expressionFilt)
print("Removing bad series and saving")
# Remove bad series again
tDFSeries <- as.data.frame(table(colData$series), stringsAsFactors = F)
goodSeries <- tDFSeries$Var1[which(tDFSeries$Freq > 1)]
colData <- colData[which(colData$series %in% goodSeries),]
normCounts <- normCounts[, which(colnames(normCounts) %in% colData$samples)]
# save(normCounts, file = "Data/Final_Expression_Matrices/Human/normCountsHuman.RData")
save(colData, file = "Data/Final_Expression_Matrices/Human/colData.RData")
save(normCounts, file = "Data/Final_Expression_Matrices/Human/normCountsHuman.RData")


# Generate corrected counts for mouse tissues
dataFile <- fileList[2]
samples <- h5read(dataFile, "meta/Sample_geo_accession")
tissue <- h5read(dataFile, name = "meta/Sample_source_name_ch1")
genes <- h5read(dataFile, name = "meta/genes")
description <- h5read(dataFile, name = "meta/Sample_description")
series <- h5read(dataFile, name = "meta/Sample_series_id")
organism <- h5read(dataFile, name = "meta/Sample_organism_ch1")
institute <- h5read(dataFile, name = "meta/Sample_contact_institute")
extractProtocol <- h5read(dataFile, name = "meta/Sample_extract_protocol_ch1")
colData <- data.frame(samples, tissue, series, organism, institute, extractProtocol)
patternSingle <- c("single cell", "single-cell", "smart seq", "smart-seq", "drop-seq", "drop seq",
                   "fluidigm", "tenX", "scRNASeq", "scRNA-Seq", "chromium")
colData <- colData[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                               x = colData$extractProtocol, ignore.case = T)),]
patternXeno <- c("xeno", "pdx", "cdx", "human", "patient")
colData <- colData[unique(grep(pattern = paste(patternXeno, collapse="|"), invert = T,
                               x = colData$extractProtocol, ignore.case = T)),]
tDFSeries <- as.data.frame(table(colData$series), stringsAsFactors = F)
goodSeries <- tDFSeries$Var1[which(tDFSeries$Freq > 1)]
colData <- colData[which(colData$series %in% goodSeries),]

print("reading in expression")
dir.create("Data/Final_Expression_Matrices/Mouse")

expression <- h5read(dataFile, "data/expression")
rownames(expression) <- genes
colnames(expression) <- samples
expression <- expression[, which(colnames(expression) %in% colData$samples)]

# Filter expression matrix
nonZeroSamps <- function(row) {
  return(sum(row > 0))
}
print("filtering")
cols <- colSums(expression)
expressionFilt <- expression[, which(cols > 5000000)]
nonZeroCount <- apply(expressionFilt, 1, nonZeroSamps)
expressionFilt <- expressionFilt[which(nonZeroCount > (length(colnames(expressionFilt)) * .10)),]
colData <- colData[which(colData$samples %in% colnames(expressionFilt)),]

print("normalizing")
# Normalize with log2 and quantile
normCounts <- log2(expressionFilt + 1)
normCounts <- normalize.quantiles(normCounts)
colnames(normCounts) <- colnames(expressionFilt)
rownames(normCounts) <- rownames(expressionFilt)
rm(expression)
rm(expressionFilt)
print("Removing bad series and saving")
# Remove bad series again
tDFSeries <- as.data.frame(table(colData$series), stringsAsFactors = F)
goodSeries <- tDFSeries$Var1[which(tDFSeries$Freq > 1)]
colData <- colData[which(colData$series %in% goodSeries),]
normCounts <- normCounts[, which(colnames(normCounts) %in% colData$samples)]
save(normCounts, file = "Data/Final_Expression_Matrices/Mouse/normCounts.RData")
save(colData, file = "Data/Final_Expression_Matrices/Mouse/colData.RData")


library(ontologyIndex)
library(GEOmetadb)
library(rhdf5)
library(jsonlite)

# Get REGEX dictionary
tissueDict <- jsonlite::read_json("Scripts/tissueDictionary.json")

dataFile <- "Data/ARCHS4_Download/human_matrix.h5"

# Categorize GEO metadata by tissue type and cancer vs normal
if (! file.exists("Data/groupList_human_raw.RData")) {
  # Get metadata
  samples <- h5read(dataFile, "meta/Sample_geo_accession")
  tissue <- h5read(dataFile, name = "meta/Sample_source_name_ch1")
  genes <- h5read(dataFile, name = "meta/genes")
  description <- h5read(dataFile, name = "meta/Sample_description")
  series <- h5read(dataFile, name = "meta/Sample_series_id")
  organism <- h5read(dataFile, name = "meta/Sample_organism_ch1")
  institute <- h5read(dataFile, name = "meta/Sample_contact_institute")
  char <- h5read(dataFile, "meta/Sample_characteristics_ch1")
  extractProtocol <- h5read(dataFile, name = "meta/Sample_extract_protocol_ch1")
  reads_aligned <- h5read(dataFile, name = "meta/reads_aligned")
  colData <- data.frame(samples, tissue, series, institute, char,
                        extractProtocol, description, reads_aligned)

  # Reads aligned > 5 million
  colData <- colData[colData$reads_aligned > 5e6,]

  # Remove single-cell
  patternSingle <- unlist(tissueDict$SingleCell)
  colData <- colData[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                                 x = colData$extractProtocol, ignore.case = T, perl = TRUE)),]
  colData <- colData[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                                 x = colData$tissue, ignore.case = T, perl = TRUE)),]
  colData <- colData[unique(grep(pattern = paste(patternSingle, collapse="|"), invert = T,
                                 x = colData$description, ignore.case = T, perl = TRUE)),]
  # -- because ARCHS4 cut off the text in descriptions, we need to download the data from GEOmetadb
  if (! file.exists("Data/hiddenSingleCell.rda")) {
    con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
    metaTable <- dbGetQuery(con,'select gsm,extract_protocol_ch1,extract_protocol_ch2,description,data_processing from gsm')
    dbDisconnect(con)
    metaTableNow <- metaTable[metaTable$gsm %in% colData$samples,]
    hiddenSingleCell <- metaTableNow$gsm[grep(x = paste(metaTableNow$extract_protocol_ch1,
                                                        metaTableNow$extract_protocol_ch2,
                                                        metaTableNow$description,
                                                        metaTableNow$data_processing),
                                              pattern = paste(patternSingle, collapse="|"), invert = F,
                                              ignore.case = T, perl = TRUE)]
    save(hiddenSingleCell, file = "Data/hiddenSingleCell.rda")
  } else {
    load("Data/hiddenSingleCell.rda")
  }
  # -- Final filter for single cell samples
  colData <- colData[! colData$samples %in% hiddenSingleCell,]

  # Keep samples that can be batch corrected
  tDFSeries <- as.data.frame(table(colData$series), stringsAsFactors = F)
  goodSeries <- tDFSeries$Var1[which(tDFSeries$Freq > 1)]
  colData <- colData[which(colData$series %in% goodSeries),]

  # Prep for grep
  colData$char <- gsub(colData$char, pattern = "Xx-xX", replacement = " ", perl = T, ignore.case = F)
  newTissue <- paste0(colData$tissue, " ", colData$char)

  # Re-format tissue names so we can grep cell types
  newTissue = gsub(pattern = "\\.", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "'", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "\\(", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "\\)", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "\\[", replacement = " ", x = newTissue)
  newTissue = gsub(pattern = "\\]", replacement = " ", x = newTissue)
  newTissue = tolower(newTissue)
  colData$tissue <- newTissue

  ## Manual assignments ##
  # First, manually assign disease and tissue types
  colData$disease <- "Unknown"

  # Remove known problematic studies
  toRemove <- c("GSE46665", "GSE45326", "GSE50957", "GSE68086", "GSE72420Xx-xXGSE72509",
                "GSE88850", "GSE100127", "GSE83402", "GSE68799", "GSE47774Xx-xXGSE47792")
  colData <- colData[! colData$series %in% toRemove,]

  # Remove 4su seq
  colData <- colData[unique(grep(pattern = "\\b4su", invert = T,
                                 x = colData$tissue, ignore.case = T, perl = TRUE)),]

  # Remove human refrence
  colData <- colData[unique(grep(pattern = "\\buhrr\\b|\\bbhrr\\b|\\bhbrr\\b|\\burr\\b|Universal", invert = T,
                                 x = colData$tissue, ignore.case = T, perl = TRUE)),]

  # Remove other unwanted types
  colData <- colData[unique(grep(pattern = "thrombocyt|extracellular|platelet|exosome", invert = T,
                                 x = colData$tissue, ignore.case = T, perl = TRUE)),]

  # Spot fix
  colData$tissue[colData$series == "GSE48865"] <- "glioma"
  colData$tissue[colData$series == "GSE72790"] <- "leukemia"
  colData$tissue[colData$series %in% c("GSE49155", "GSE39121", "GSE62118")] <- "lung cancer"
  colData$tissue[colData$series %in% c("GSE65185Xx-xXGSE65186", "GSE78220")] <- "melanoma"
  colData$tissue[colData$series %in% c("GSE79492", "GSE56066")] <- "breast cancer"
  colData$tissue[colData$series == "GSE60052" & colData$tissue == "tumor"] <- "lung cancer"
  colData$tissue[colData$series %in% c("GSE47944Xx-xXGSE47965")] <- "skin"
  colData$tissue[colData$series == "GSE66301"] <- "breast cancer"
  colData$tissue[colData$series == "GSE66306"] <- "pbmc"

  # Categorize cancer samples
  cancerString <- unlist(tissueDict$tumors[[1]])
  normalString <- unlist(tissueDict$tumors[[2]])
  cancerGrep <- grep(colData$tissue, pattern = paste0(cancerString,
                                                      collapse = "|"),
                     perl = TRUE, ignore.case = TRUE)
  colData$disease[cancerGrep] <- "Cancer"
  normalGrep <- grep(colData$tissue, pattern = paste0(normalString,
                                                      collapse = "|"),
                     perl = TRUE, ignore.case = TRUE)
  colData$disease[normalGrep] <- "Unknown"

  # hidden tumor samps -- poorly annotated...
  hiddenSamps <- paste0("GSM", c(
    1185643:1185648
  ))
  colData$disease[colData$samples %in% hiddenSamps] <- "Cancer"

  # Categorize samples by tissue type
  tissueDictNow <- tissueDict[c(1:27)]
  colDataCat <- categorizeMetaData(metadata = colData,
                                   cols = "tissue",
                                   dictionary = tissueDictNow)
  rSums <- rowSums(colDataCat[,c(10:36)])
  colDataMain <- unique(colDataCat[rSums == 1,]) # Keep only unique assigments
  colDataMain$category <- "category"
  possibles <- colnames(colDataMain[,c(10:36)])
  resVec <- c()
  for (i in 1:length(colDataMain$samples)) {
    resVec <- c(resVec, possibles[which(colDataMain[i,c(10:36)] == 1)])
  }
  colDataMain$category <- resVec

  # Samples which received no classification -- use Cellosaurus for secondary classification
  colDataNoClass <- unique(colDataCat[rSums == 0,])
  # -- Get list of cancer cells and cancer types from cellosaurus
  if (! file.exists("Data/cellDisease.rda")) {
    lines <- readLines("Data/cellosaurus.txt")
    # Get list of all cancer cell lines
    cancerCellVec <- c()
    diseaseVec <- c()
    for (i in 1:length(lines)) {
      line <- lines[i]
      grepRes1 <- grep(line, pattern = "^ID")
      grepRes2 <- grep(line, pattern = "^SY")
      grepRes3 <- grep(line, pattern = "^CA")
      grepRes4 <- grep(line, pattern = "^DI")
      if (length(grepRes1)) {
        tumors <- FALSE
        disease <- NULL
        cellsNow <- gsub(x = line, pattern = "^ID\\s+",
                         replacement = "", perl = TRUE)
      }
      if (length(grepRes2)) {
        cellsNow2 <- gsub(x = line, pattern = "^SY\\s+",
                          replacement = "", perl = TRUE)
        cellsNow <- c(cellsNow, unlist(strsplit(cellsNow2, split = "; ")))
      }
      if (length(grepRes4)) {
        disease <- gsub(x = line, pattern = "^DI.+;.+; ",
                        replacement = "", perl = TRUE)
      }
      if (length(grepRes3)) {
        if (line == "CA   Cancer cell line") {
          if (is.null(disease)) {
            next
          }
          cancerCellVec <- c(cancerCellVec, cellsNow)
          diseaseVec <- c(diseaseVec, rep(disease, length(c(cellsNow))))
          if (! length(cancerCellVec) == length(diseaseVec)) {
            stop()
          }
        }
      }
    }

    cellDisease <- data.frame(
      "cell" = cancerCellVec,
      "disease" = diseaseVec
    )
    save(cellDisease, file = "Data/cellDisease.rda")

  } else {
    load("Data/cellDisease.rda")
  }
  # -- grep cell types and assign disease
  if (! file.exists("Data/tissueDisease.rda")) {
    tissueUnique <- unique(colDataNoClass$tissue)
    # Clean for grep
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\+", replacement = "\\\\+")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\[", replacement = " ")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\]", replacement = " ")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\(", replacement = " ")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\)", replacement = " ")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\.", replacement = "\\\\.")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\^", replacement = "\\\\^")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\$", replacement = "\\\\$")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\?", replacement = "\\\\?")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\|", replacement = "\\\\|")
    cellDisease$cell <- gsub(cellDisease$cell, pattern = "\\*", replacement = "\\\\*")
    cellDisease$cell2 <- gsub(cellDisease$cell, pattern = "\\#", replacement = "")

    # Prepare for vec
    tissueUnique <- unique(colDataNoClass$tissue)
    # tissueUnique <- gsub(tissueUnique, pattern = "h sts", replacement = "hsts")

    # Pattern to match cell types with or without add info separated by - or _
    # This basically says, there can be something such that "-cellName-",
    # "_cellName_" or " cellName " and all the possible combinations.
    cellDisease$cell1 <- sapply(cellDisease$cell, FUN = function(x) {
      paste0("\\b", x, "[_-]+|\\b", x, "\\b|[_-]+", x, "\\b|[_-]+",x,"[_-]+")
    })
    cellDisease$nchar <- sapply(cellDisease$cell2, FUN = nchar)
    cellDisease <- cellDisease[cellDisease$nchar > 2,]

    tableDF <- as.data.frame(table(cellDisease$disease))
    cellDiseaseList <- split(cellDisease, f = cellDisease$disease)
    calculateGrep <- function(wordVec, patternVec) {
      # wordVec <-  tissueUnique[1:100]
      # patternVec <- cancerCellVec1[30344:3344]
      return(sapply(wordVec, FUN = grep, ignore.case = TRUE,
                    pattern = paste0(as.character(patternVec), collapse = "|"), perl = TRUE))
    }

    # Split the tissueUnique -- from https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
    chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
    tissueUniqueList <- chunk2(tissueUnique, 100)

    # apply to each chunk -- map GEO tissue annotations to cellosaurus diseases
    tissueDiseaseList <- parallel::mclapply(tissueUniqueList,
                                            mc.cores = length(tissueUniqueList),
                                            FUN = function(tissueNow) {
                                              namesList <- c()
                                              diseaseList <- c()
                                              for (i in 1:length(cellDiseaseList)) {
                                                print(i)
                                                dataNow <- cellDiseaseList[[i]]
                                                diseaseNow <- names(cellDiseaseList)[i]
                                                cancerVecNow <- unique(dataNow$cell1)
                                                cancerVecNowSplit <- split(cancerVecNow, ceiling(seq_along(cancerVecNow)/100)) # Chunks of 500
                                                names(cancerVecNowSplit) <- paste0(names(cancerVecNowSplit), "_xxx_")
                                                grepRes <- lapply(cancerVecNowSplit, FUN = function(splitNow) {
                                                  calculateGrep(tissueNow, patternVec = splitNow)
                                                })
                                                grepRes <- unlist(grepRes, recursive = T, use.names = T)
                                                namesNow <- gsub(names(grepRes), pattern = ".+_xxx_\\.", replacement = "")
                                                namesList <- c(namesList, namesNow)
                                                diseaseList <- c(diseaseList, rep(diseaseNow, length(namesNow)))
                                                if (! length(namesList) == length(diseaseList)) {
                                                  stop()
                                                }
                                              }

                                              tissueDisease <- data.frame(
                                                "tissue" = namesList,
                                                "disease" = diseaseList
                                              )
                                              tissueDisease
                                            })
    tissueDisease <- data.table::rbindlist(tissueDiseaseList)
    save(tissueDisease, file = "Data/tissueDisease.rda")
  } else {
    load("Data/tissueDisease.rda")
  }

  # Remove non-human classes
  tissueDisease <- tissueDisease[grep(tissueDisease$disease, pattern = "\\bhamster\\b|\\brat\\b|\\bmouse\\b|\\bcanine\\b|\\bfeline\\b|Li-Fraumeni|Mycosis|\\bgoldfish\\b",
                                      ignore.case = T, invert = T),]

  # Clean up duplicates
  tissueDisease <- tissueDisease[which(! duplicated(tissueDisease$tissue)),]

  # Wrangle
  colDataCell <- colDataNoClass[colDataNoClass$tissue %in% tissueDisease$tissue,]
  colDataCell <- merge(x = colDataCell[,c(1:8)], y = tissueDisease, by = "tissue")
  colDataCell2 <- colDataCell[,c(-9)]
  colDataCell2$tissue_orig <- colDataCell2$tissue
  colDataCell2$tissue <- colDataCell$disease
  colDataCell2$disease <- "Cancer"

  # Make final tissue assignments
  colDataCellCat <- categorizeMetaData(metadata = colDataCell2,
                                       cols = "tissue",
                                       dictionary = tissueDictNow)
  rSums <- rowSums(colDataCellCat[,c(11:37)])
  colDataCellCatMain <- unique(colDataCellCat[rSums == 1,]) # Keep only unique assigments
  colDataCellCatMain$category <- "category"
  possibles <- colnames(colDataCellCatMain[,c(11:37)])
  resVec <- c()
  for (i in 1:length(colDataCellCatMain$samples)) {
    resVec <- c(resVec, possibles[which(colDataCellCatMain[i,c(11:37)] == 1)])
  }
  colDataCellCatMain$category <- resVec

  # Combine data
  colDataFinal <-  data.frame(
    "samples" = c(colDataMain$samples, colDataCellCatMain$samples),
    "series" = c(colDataMain$series, colDataCellCatMain$series),
    "orig_tissue" = c(colDataMain$tissue, colDataCellCatMain$tissue_orig),
    "readsAligned" = c(colDataMain$reads_aligned, colDataCellCatMain$reads_aligned),
    "disease" = c(colDataMain$disease, colDataCellCatMain$disease),
    "Tissue" = c(colDataMain$category, colDataCellCatMain$category), stringsAsFactors = F
  )
  colDataFinal <- colDataFinal[! duplicated(colDataFinal$sample),]
  table(colDataFinal$Tissue, colDataFinal$disease)
  # Final clean-up
  colDataFinal <- colDataFinal[! colDataFinal$Tissue %in% c(
    "thymus", "spleen", "unknown"
  ),]
  colDataFinal$disease[colDataFinal$disease == "Unknown"] <- "Normal"
  colDataFinal$Tissue[colDataFinal$series %in% c("GSE140442")] <- "stem-like"
  colDataFinal$disease[colDataFinal$series %in% c("GSE140442")] <- "Normal"
  colDataFinal <- colDataFinal[! is.na(colDataFinal$Tissue),]

  # Build group list
  groupList <- lapply(unique(colDataFinal$Tissue), FUN = function(tissue) {
    list(
      "all" = colDataFinal$samples[colDataFinal$Tissue == tissue],
      "cancer" = colDataFinal$samples[colDataFinal$Tissue == tissue &
                                        colDataFinal$disease == "Cancer"],
      "normal" = colDataFinal$samples[colDataFinal$Tissue == tissue &
                                        colDataFinal$disease == "Normal"]
    )
  })
  names(groupList) <- unique(colDataFinal$Tissue)

  # This constitutes the colData prior to dim reduction assisted assignment
  groupList_human_raw <- groupList
  colData_human_raw <- colDataFinal
  save(groupList_human_raw, colData_human_raw, file = "Data/groupList_human_raw.RData")
} else {
  load("Data/groupList_human_raw.RData")
  colDataFinal <- colData_human_raw
}

# ## Dimensionality reduction assisted assignments ##
# cat("\n", timestamp2(), " Loading & filtering expression data...\n", sep = "")
# if (! file.exists("Data/fullRawCountsFiltered.rda")) {
#   # Load and filter expression data
#   expression <- h5read(dataFile, "data/expression")
#   samples <- h5read(dataFile, "meta/Sample_geo_accession")
#   genes <- h5read(dataFile, name = "meta/genes")
#   rownames(expression) <- genes
#   colnames(expression) <- samples
#   colDataFinal <- colDataFinal[colDataFinal$samples %in% colnames(expression),]
#   expression <- expression[, which(colnames(expression) %in% colDataFinal$samples)]
#   colDataFinal <- colDataFinal[order(match(colDataFinal$samples, colnames(expression))),]
#   if (! all(colDataFinal$samples == colnames(expression))) {
#     stop("ColData samples are not identical to colnames of expression data...",
#          " Please email Code/generateFigures_logFile.txt to author if you find this error and/or submit issue on github.")
#   }
#   # # Keep genes expressed in 10%+ of samples
#   # nonZeroCount <- apply(expression, 1, nonZeroSamps)
#   # expression <- expression[which(nonZeroCount > (length(colnames(expression)) * .01)),]
#   timestamp()
#   cat("\nDone. Saving expression data...\n")
#   save(expression, file = "Data/fullRawCountsFiltered.rda")
# } else if (file.exists("Data/fullVSTCounts.rda")) {
#   cat("\nVST-transformed expression data found... Skipping this step...\n")
# } else {
#   cat("\nFiltered expression data found! Loading and continuing to next step...\n")
#   load("Data/fullRawCountsFiltered.rda")
#   colDataFinal <- colDataFinal[colDataFinal$samples %in% colnames(expression),]
#   colDataFinal <- colDataFinal[order(match(colDataFinal$samples, colnames(expression))),]
#   all(colnames(expression) == colDataFinal$samples)
#   rownames(colDataFinal) <- colDataFinal$samples
# }

# cat("\n", timestamp2(), " Applying Variance Stabilizing Transform to dataset...\n", sep = "")
# UMAP <- FALSE
# if (! file.exists("Data/fullVSTCounts.rda")) {
#   # Make data Homoscedastic with VST
#   dds <- DESeqDataSetFromMatrix(expression, colData = colDataFinal,
#                                 design = ~1)
#   # from https://support.bioconductor.org/p/62246/#62250
#   inds <- rownames(expression)
#   geoMeansList <- mclapply(inds, FUN = function(ind) {
#     row <- expression[ind,]
#     if (all(row == 0)) {
#       0
#     } else {
#       exp(sum(log(row[row != 0]))/length(row))
#     }
#   }, mc.cores = 60)
#   geoMeans <- unlist(geoMeansList)
#   dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
#   vsd <- vst(dds)
#   timestamp()
#   cat("\nDone. Saving vst data...\n")
#   save(vsd, file = "Data/fullVSTCounts.rda")
# } else if (file.exists("Data/fullUMAPData.rda")) {
#   UMAP <- TRUE
#   cat("\nUMAP results found... Skipping this step...\n")
# } else {
#   cat("\nVST-transformed expression data found! Loading and continuing to next step...\n")
#   load("Data/fullVSTCounts.rda")
# }
#
# if (! file.exists("Data/fullUMAPData.rda") & ! UMAP) {
#   cat("\n", timestamp2(), " Calculating UMAP...\n", sep = "")
#   # Calculate highly-variable genes (HVGs)
#   rv <- matrixStats::rowVars(assay(vsd))
#   names(rv) <- rownames(assay(vsd))
#   rv <- rv[order(rv, decreasing = T)]
#   hvgs <- rv[c(1:10000)]
#   vsdHVG <- assay(vsd)[names(hvgs),]
#   # rm(vsd)
#   # gc()
#   colDataFinal <- colDataFinal[colDataFinal$samples %in% colnames(vsdHVG),]
#   colDataFinal <- colDataFinal[order(match(colDataFinal$samples, colnames(vsdHVG))),]
#   if (! all(colDataFinal$samples == colnames(vsdHVG))) {
#     stop("ColData samples are not identical to colnames of expression data...",
#          " Please email Code/generateFigures_logFile.txt to author if you find this error and/or submit issue on github.")
#   }
#   # Calculate PCA
#   if (! file.exists( "Data/fullPCData.rda")) {
#     pcDataFull <- prcomp(t(vsdHVG))
#     pcData <- pcDataFull$x
#     pcData <- pcData[,c(1:500)]
#     save(pcData, file = "Data/fullPCData.rda")
#   } else {
#     load("Data/fullPCData.rda")
#   }
#   colDataFinal$PC1 <- pcData[,c(1)]
#   colDataFinal$PC2 <- pcData[,c(2)]
#
#   # Calculate Neighbors
#   neighbors <- FindNeighbors(pcData[,c(1:200)], k.param = 200)
#   # Calculate clustering
#   clusters <- FindClusters(neighbors$snn)
#   colDataFinal$cluster <- clusters[,c(1)]
#   # Calculate UMAP embedding
#   nnNow <- floor(length(colnames(vsdHVG)) * .15)
#   umapData <- umap(t(vsdHVG), verbose = T, min_dist = 1,
#                    n_neighbors = nnNow, n_threads = 120,
#                    pca = 100)
#   colDataFinal$samples == colnames(vsdHVG)
#
#   colDataFinal$UMAP_1 <- umapData[,c(1)]
#   colDataFinal$UMAP_2 <- umapData[,c(2)]
#   colDataUMAP <- colDataFinal
#   timestamp()
#   save(colDataUMAP, file = "Data/fullUMAPData.rda")
# } else {
#   if (UMAP) {
#     cat("\nUMAP data found! Loading and continuing to next step...\n")
#   }
#   load("Data/bulkRNASeq/fullUMAPData.rda")
# }
#
# library(ggpubr)
# g1 <- ggplot(data = colDataUMAP, mapping = aes(x = UMAP_1,
#                                                y = UMAP_2,
#                                                color = Tissue)) +
#   geom_point(size = .4) +
#   theme_pubr(border = T, base_size = 22) +
#   guides(colour = guide_legend(override.aes = list(size=3))) +
#   theme(legend.position="right") +
#   theme(#border = element_line(colour = "black", size = 1),
#     panel.border = element_rect(colour = "black", fill=NA, size=1))
# ggsave(g1, filename = "umap_tissue.png", height = 20, width = 40)
#
#
#
#
















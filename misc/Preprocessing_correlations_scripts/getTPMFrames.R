library(foreach)
doParallel::registerDoParallel()

makeUploadTable <- function(tx2gene, geneNames, transcripts, 
                            dataFile, samplesInd, samples) {
  resList <- list()
  # geneNames <- geneNames[c(1:2)]
  for (j in 1:length(geneNames)) {
    gene <- geneNames[j]
    print(j)
    # gene <- "AC012454.4"
    print(gene)
    tx2geneNow <- tx2gene[which(tx2gene$gene == gene),]
    transcriptsInd <- which(transcripts %in% tx2geneNow$tx)
    if (! length(transcriptsInd)) {
      resList[[j]] <- NULL
      next
    }
    H5close()
    expression <- h5read(dataFile, "data/expression", 
                         index = list(transcriptsInd, samplesInd))
    H5close()
    colnames(expression) <- samples[samplesInd]
    rownames(expression) <- transcripts[transcriptsInd]
    geneId <- rep(gene, length(rownames(expression)))
    expression <- rowsum(expression, geneId)
    vec <- signif(as.numeric(expression),digits = 3)
    vecStr <- paste(vec, collapse = ",")
    resList[[j]] <- data.frame(geneName = gene, values = vecStr)
  }
  resDF <- data.table::rbindlist(resList)
  return(resDF)
}

# Download TPM files
humanTxURL <- "https://s3.amazonaws.com/mssm-seq-matrix/human_tpm_v7.h5"
mouseTxURL <- "https://s3.amazonaws.com/mssm-seq-matrix/mouse_tpm_v7.h5"
topdir <- getwd()
downList <- c(humanTxURL, mouseTxURL)
archdatDir <- file.path(topdir, "Data/ARCHS4_Download/")
# for (i in 1:length(downList)) {
#   setwd(archdatDir)
#   file <- downList[i]
#   print(file)
#   system(paste0("wget ", file), ignore.stdout = T, ignore.stderr = T, wait = F)
#   setwd(topdir)
# }

# Get tx2gene objects

library(biomaRt)
print("Getting biomaRt")
# Get v87 for human
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl",
                   host = "http://dec2016.archive.ensembl.org")
tx2geneHuman <- getBM(attributes = c("ensembl_transcript_id", 
                                     "external_gene_name"),
                      mart = ensembl)
# Get v88 for mouse
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl",
                   host = "http://mar2017.archive.ensembl.org")
tx2geneMouse <- getBM(attributes = c("ensembl_transcript_id", 
                                     "external_gene_name"),
                      mart = ensembl)

## Human
# Aggregate transcript TPM to gene level
library(tximport)
library(rhdf5)
fileList <- paste0("Data/ARCHS4_Download/", 
                   gsub(downList[c(1, 2)], 
                        pattern = "https://s3.amazonaws.com/mssm-seq-matrix/", 
                        replacement = ""))

for (species in c("Human", "Mouse")) {
  print(species)
  speciesScien <- ifelse(test = species == "Human", yes = "hsapiens", no = "mmusculus")
  cat("Loading data ... \n")
  dataDir <- "Data/Final_Expression_Matrices/"
  load(file.path(dataDir, species, "colData.RData"))
  load(file.path(dataDir, species, "groupList.RData"))
  load(file.path(dataDir, species, "geneNames.RData"))
  load(file.path(dataDir, species, "finalSamples.RData"))
  sampsKeep <- sampList
  if (species == "Human") {
    dataFile <- fileList[1]
    tx2gene <- tx2geneHuman
  } else {
    dataFile <- fileList[2]
    tx2gene <- tx2geneMouse
  }
  samples <- h5read(dataFile, "meta/Sample_geo_accession")
  
  colnames(tx2gene) <- c("tx", "gene")
  samples <- h5read(dataFile, "meta/Sample_geo_accession")
  transcripts <- h5read(dataFile, "meta/transcripts")
  transcripts <- gsub(transcripts, pattern = "\\..*", replacement = "")
  samplesInd <- which(samples %in% sampsKeep)
  if (species == "Human") {
    sampleTPMOrderHuman <- samples[samplesInd]
    
    save(sampleTPMOrderHuman, file = file.path(dataDir, species,
                                     "sampleTPMOrderHuman.RData"))
  } else {
    sampleTPMOrderMouse <- samples[samplesInd]
    
    save(sampleTPMOrderMouse, file = file.path(dataDir, species,
                                               "sampleTPMOrderMouse.RData"))
  }
  
  # Set up foreach
  cores <- 60
  nr <- length(geneNames)
  n <- ceiling(nr/cores)
  forUploadList <- split(geneNames, rep(1:ceiling(nr/n), each=n, length.out=nr))
  k <- length(forUploadList)
  print("running foreach")
  start <- proc.time()
  resList <- foreach(j = 1:k, 
                     .verbose = T,
                     .export = c("k", "forUploadList")) %dopar%
    makeUploadTable(tx2gene = tx2gene, geneNames = forUploadList[[j]], 
                    transcripts = transcripts, dataFile = dataFile, 
                    samplesInd = samplesInd, samples = samples)
  finalDF <- data.table::rbindlist(resList)
  ans <- proc.time() - start
  print(ans)
  print("Saving")
  save(finalDF, file = file.path(dataDir, species, "geneTPM_forUpload.RData"))
}




# load(file.path(dataDir, species, "geneTPM_forUpload.RData"))
# humanGenesTPM <- unique(as.character(finalDF$geneName))
# save(humanGenesTPM,file = file.path(dataDir, species, "geneTPM_forUpload_geneList.RData"))
# 
# species <- "Mouse"
# load(file.path(dataDir, species, "geneTPM_forUpload.RData"))
# mouseGenesTPM <- unique(as.character(finalDF$geneName))
# save(mouseGenesTPM, file = file.path(dataDir, species, "geneTPM_forUpload_geneList.RData"))
# 

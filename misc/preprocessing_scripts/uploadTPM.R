library(foreach)
doParallel::registerDoParallel()
credFile <- "Data/credFile.txt"
source("Scripts/helpers.R")

# Get TPM data
if (! file.exists("Data/ARCHS4_Download/human_tpm_v8.h5")){
  # Download TPM files
  humanTxURL <- "https://s3.amazonaws.com/mssm-seq-matrix/human_tpm_v8.h5"
  download.file(humanTxURL, destfile = "Data/ARCHS4_Download/human_tpm_v8.h5")
}

# Process for upload
if (! file.exists("Data/geneTPM_forUpload_geneList.RData")) {
  # Get tx2gene objects
  library(biomaRt)
  print("Getting biomaRt")
  # Get v87 for human
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl",
                     host = "http://dec2016.archive.ensembl.org")
  tx2geneHuman <- getBM(attributes = c("ensembl_transcript_id", 
                                       "external_gene_name"),
                        mart = ensembl)
  # # Get v88 for mouse
  # ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl",
  #                    host = "http://mar2017.archive.ensembl.org")
  # tx2geneMouse <- getBM(attributes = c("ensembl_transcript_id", 
  #                                      "external_gene_name"),
  #                       mart = ensembl)
  
  ## Human
  # Aggregate transcript TPM to gene level
  library(tximport)
  library(rhdf5)
  speciesScien <- "hsapiens"
  cat("Loading data ... \n")
  dataFile <- "Data/ARCHS4_Download/human_tpm_v8.h5"
  load("Data/groupList_human_raw.RData")
  samples <- h5read(dataFile, "meta/Sample_geo_accession")
  tx2gene <- tx2geneHuman
  colnames(tx2gene) <- c("tx", "gene")
  load("Data/hsapiens_corrSmall_geneNames.rda")
  transcripts <- h5read(dataFile, "meta/transcripts")
  transcripts <- gsub(transcripts, pattern = "\\..*", replacement = "")
  samplesInd <- which(samples %in% colData_human_raw$samples)
  names(samplesInd) <- samples[samplesInd]
  sampleTPMOrderHuman <- names(samplesInd)
  save(sampleTPMOrderHuman, file = file.path("Data/sampleTPMOrderHuman.RData"))
  
  geneNames <- hsapiens_corrSmall_geneNames
  
  # Set up foreach
  cores <- 120
  nr <- length(samplesInd)
  n <- ceiling(nr/cores)
  forUploadList <- split(samplesInd, rep(1:ceiling(nr/n), each=n, length.out=nr))
  k <- length(forUploadList)
  print("running foreach")
  start <- proc.time()
  resList <- foreach(j = 1:k, 
                     .verbose = T,
                     .export = c("k", "forUploadList")) %dopar%
    makeUploadTable(tx2gene = tx2gene, samplesInd = forUploadList[[j]], 
                    transcripts = transcripts, dataFile = dataFile)
  
  finalDF <- data.table::rbindlist(resList)
  ans <- proc.time() - start
  print(ans)
  print("Saving")
  save(finalDF, file = file.path("Data/geneTPM_forUpload.RData"))
  humanSamplesTPM <- unique(as.character(finalDF$sampleName))
  save(humanSamplesTPM,file = "Data/geneTPM_forUpload_sampleList.RData")
}

# Upload
finalDF2 <- finalDF[,c(-1)]
rownames(finalDF2) <- finalDF$sampleName
checkBool <- uploadToAzure(finalDF2, credentials = credFile,
                            tableName = "tpm_hsapiens", check = T)
if (checkBool & ! doKey) {
  # Upload to AWS
  print("uploading to Azure ... ")
  uploadToAzure(finalDF2, credentials = credFile,
                tableName = "tpm_hsapiens", check = F)
  # sql <- paste0("ALTER TABLE ", tableName, " ADD PRIMARY KEY (row_names);")
  # dbExecute(conn = conn, statement = sql)
  # tryCatch(expr = {uploadToAzure(finalDF2, tableName)},
  #          error = function(e) {cat("AWS error -- probably not enough space")})
  
} else {
  print("Table already uploaded ... ")
}
if (doKey) {
  credentials <- "Data/credFile.txt"
  credentials <- suppressWarnings(read.delim(credentials, sep = ";",
                                             header = FALSE, stringsAsFactors = F))
  uName <- credentials$V1
  pWord <- credentials$V2
  tableName = "tpm_hsapiens"
  conn <- dbConnect(drv = RMySQL::MySQL(), user = uName, 
                    port = 3306, dbname="correlation_analyzer",
                    password=pWord,
                    host="m2600az-db01p.mysql.database.azure.com")
  sql <- paste0("ALTER TABLE ", tableName, " ADD PRIMARY KEY (row_names);")
  done <- 0
  while (done  == 0) {
    done <- tryCatch(expr = {
      dbExecute(conn = conn, statement = sql)
      dbDisconnect(conn)
      rm(conn)
      gc()
      1
    },
    error = function(e) {
      print(e$message)
      retCall <- "could not run statement: Multiple primary key defined"
      if (e$message == retCall) {
        cat("Already defined key")
        dbDisconnect(conn)
        rm(conn)
        gc()
        1
      } else {
        cat("Fail -- retry with new connection")
        rm(conn)
        gc()
        conn <- dbConnect(drv = RMySQL::MySQL(), user = uName,
                          port = 3306, dbname="correlation_analyzer",
                          password=pWord,
                          host="m2600az-db01p.mysql.database.azure.com")
        0
      }
    })
    print(done)
  }
}
lapply(dbListConnections(drv = RMySQL::MySQL()), dbDisconnect)
Sys.sleep(3)
cat("\nNext sample ... \n")
rm(finalDF2)
gc()





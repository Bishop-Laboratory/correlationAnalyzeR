load("Data/vsd_for_corr.rda")

library(foreach)
doParallel::registerDoParallel()
credFile <- "Data/credFile.txt"
source("Scripts/helpers.R")

# Process for upload
if (! file.exists("Data/geneVSD_forUpload_geneList.RData")) {
  # Get tx2gene objects
  # Aggregate transcript TPM to gene level
  cores <- 120
  genes <- rownames(vsd)
  nr <- length(genes)
  n <- ceiling(nr/cores)
  dfVsd <- as.data.frame(assay(vsd))
  forUploadList <- split(dfVsd, rep(1:ceiling(nr/n), each=n, length.out=nr))
  samplesVSD <- colnames(dfVsd)    
  save(samplesVSD, file = "Data/samplesVSD.rda")
  k <- length(forUploadList)
  print("running foreach")
  start <- proc.time()
  resList <- mclapply(forUploadList, FUN = makeUploadVSD, mc.cores = 10)
  
  finalDF <- data.table::rbindlist(resList)
  ans <- proc.time() - start
  print(ans)
  print("Saving")
  save(finalDF, file = file.path("Data/geneVSD_forUpload.RData"))
  humanGenesTPM <- unique(as.character(finalDF$geneName))
  save(humanGenesTPM,file = "Data/geneVSD_forUpload_geneList.RData")
}

# Upload
finalDF2 <- finalDF[,c(-1)]
rownames(finalDF2) <- finalDF$geneName
checkBool <- uploadToAzureGeneKey(finalDF2, credentials = credFile,
                                  tableName = "vsd_hsapiens", check = T)
if (checkBool & ! doKey) {
  # Upload to AWS
  print("uploading to Azure ... ")
  uploadToAzureGeneKey(finalDF2, credentials = credFile,
                tableName = "vsd_hsapiens", check = F)
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
  tableName = "vsd_hsapiens"
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









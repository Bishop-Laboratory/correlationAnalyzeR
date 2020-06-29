library(sva)
library(BiocParallel)
library(foreach)
library(data.table)
library(DBI)
#library(future)
library(DESeq2)
#options(future.globals.maxSize= 5e10)
#plan(multicore(workers = 10))
library(WGCNA)
WGCNA::enableWGCNAThreads()
doParallel::registerDoParallel()

source("Scripts/helpers.R")

set.seed(42)

# Should we just use the primary key logic?
doKey <- TRUE

species <- "Human"
cat("Species is ", species)
speciesScien <- ifelse(test = species == "Human", yes = "hsapiens", no = "mmusculus")
cat("Loading data ... \n")
dataDir <- "Data/"
load("Data/groupList_human_raw.RData")
credFile <- "Data/credFile.txt"
colData <- colData_human_raw
groupList <- groupList_human_raw

# Add the all category
cancerList <- c()
normalList <- c()
for (i in 1:length(groupList)) {
  if (length(groupList[[i]]$cancer)) {
    cancerList <- c(cancerList, sample(x = groupList[[i]]$cancer, 100, replace = T))
  }
  if (length(groupList[[i]]$normal)) {
    normalList <- c(normalList, sample(x = groupList[[i]]$normal, 100, replace = T))
  }
}

groupList[["all"]] <- list("all" = unique(c(cancerList, normalList)),
                           "cancer" = unique(cancerList),
                           "normal" = unique(normalList))

dir.create("Data/corr_mats")

cat("Starting analysis ... \n")
for (i in 1:length(groupList)) {
  groups <- groupList[[i]]
  name <- names(groupList)[i]
  name <- gsub(name, pattern = " ", replacement = "_")
  dir.create(file.path("Data/corr_mats", name))
  cat("Current group: ", name, "\n")
  for (type in c("all", "cancer", "normal")) {
    #  samples
    Samps <- groups[[type]]
    n <- length(Samps)
    cat(paste0("\n", type, " - ", n), "\n")
    if (n < 30) {
      warning("Not enough samples")
    } else {
      Path <- file.path("Data/corr_mats", name, type)
      dir.create(path = Path, showWarnings = F)
      if (! file.exists(file.path(Path, "corMatForUpload.RData"))) {
        load("Data/vsd_for_corr.rda")
        expNow <- assay(vsd)[,colnames(assay(vsd)) %in% Samps]
        rm(vsd)
        gc()
        colDataNow <- colData[which(colData$samples %in% colnames(expNow)),]
        colDataNow <- colDataNow[order(match(colDataNow$samples, colnames(expNow))),]
        all(colDataNow$samples %in% colnames(expNow))
        Samps <- colDataNow$samples
        corMat <- WGCNA::cor(x = t(expNow), verbose = 1, nThreads = 100)
        forUpload <- cbind(as.data.frame(rownames(corMat)), corMat)
        cores <- 20
        colnames(forUpload)[1] <- "geneName"
        nr <- nrow(forUpload)
        n <- ceiling(nr/cores)
        forUploadList <- split(forUpload, rep(1:ceiling(nr/n), each=n, length.out=nr))
        k <- length(forUploadList)
        print("running foreach")
        start <- proc.time()
        resList <- mclapply(forUploadList, convertForUpload, mc.cores = cores)
        finalDF <- data.table::rbindlist(resList)
        ans <- proc.time() - start
        print(ans)
        print("Saving")
        finalDF2 <- finalDF[,c(-1), drop = F]
        rownames(finalDF2) <- finalDF$geneName
        save(finalDF2, file = file.path(Path, "corMatForUpload.RData"))
      } else {
        # Check whether upload necessary before loading
        tableName <- paste0("correlations_", speciesScien, "_", type, "_", name)
        tableName <- gsub(tableName, pattern = "-", replacement = "_")
        checkBool <- uploadToAzure(tableName = tableName, check = T,
                                   credentials = credFile)
        if (checkBool & ! doKey) {
          print("Loading matrix for upload ... ")
          load(file.path(Path, "corMatForUpload.RData"))
          
        }
      }
      tableName <- paste0("correlations_", speciesScien, "_", type, "_", name)
      tableName <- gsub(tableName, pattern = "-", replacement = "_")
      checkBool <- uploadToAzure(tableName = tableName, credentials = credFile,
                                 check = T)
      if (checkBool & ! doKey) {
        # Upload to AWS
        print("uploading to Azure ... ")
        uploadToAzure(finalDF2, credentials = credFile,
                      tableName, check = F)
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
    } 
  }
}

# 
# hsapiens_corrSmall_geneNames <- rownames(finalDF2)
# save(hsapiens_corrSmall_geneNames, file = "Data/hsapiens_corrSmall_geneNames.rda")
# 






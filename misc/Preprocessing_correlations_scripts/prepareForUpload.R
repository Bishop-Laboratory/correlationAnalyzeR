library(sva)
library(BiocParallel)
library(foreach)
library(data.table)
library(DBI)
library(future)

doParallel::registerDoParallel()


for (species in c("Human", "Mouse")) {
  print(species)
  speciesScien <- ifelse(test = species == "Human", yes = "hsapiens", no = "mmusculus")
  print("Loading data")
  dataDir <- "Data/Final_Expression_Matrices/"
  load(file.path(dataDir, species, "colData.RData"))
  load(file.path(dataDir, species, "groupList.RData"))
  load(file.path(dataDir, species, "geneNames.RData"))
  
  cat("\nStarting analysis ... \n")
  for (i in 1:length(groupList)) {
    groups <- groupList[[i]]
    name <- names(groupList)[i]
    cat("\n", name, "\n")
    namePath <- file.path("Data/Final_Expression_Matrices/", species, name)
    name <- gsub(x = name, pattern = " ", replacement = "0")
    for (type in c("cancer", "normal")) {
      convertForUpload <- function(forUploadListPart) {
        listVec <- list()
        for (i in 1:length(forUploadListPart[,1])) {
          corr_frame_row <- forUploadListPart[i,]
          geneName <- as.character(corr_frame_row[1,1])
          corr_frame_row <- corr_frame_row[1,-1]
          vec <- signif(as.numeric(corr_frame_row),digits = 3)
          vecStr <- paste(vec, collapse = ",")
          listVec[[i]] <- data.frame(geneName = geneName, values = vecStr)
        }
        resDF <- rbindlist(listVec)
        return(resDF)
      }
      
      correctSeries <- function(colData) {
        colData <- colDataNow
        sDF <- as.data.frame(table(colData$series), stringsAsFactors = F)
        badSeries <- sDF$Var1[which(sDF$Freq < 2)]
        if (length(badSeries)) {
          print("series correction")
          colData <- colData[which(! colData$series %in% badSeries),]
          return(colData)
        } else {
          print("No series correction")
          return(colData)
        }
      }
      
      uploadToAWS <- function(finalDF2= NULL, tableName, check = F, geneNames) {
        con <- dbConnect(RMySQL::MySQL(), user = "xxxxxx", 
                         port = 3306, dbname="bishoplabdb",
                         password='xxxxxxx', 
                         host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
        sql <- paste0("SELECT COUNT(*) FROM ", tableName)
        tabs <- DBI::dbListTables(con)
        if(! tableName %in% tabs & check) {
          DBI::dbDisconnect(con)
          return(TRUE)
        } else if (check) {
          print("Before query")
          res <- DBI::dbSendQuery(con, sql)
          res <- DBI::dbFetch(res)
          print("after query")
          countVal <- res[1,1]
          print(countVal)
          n <- length(geneNames)
          print(n)
          if (countVal != n) {
            DBI::dbDisconnect(con)
            return(TRUE)
          } else {
            DBI::dbDisconnect(con)
            return(FALSE)
          }
        } else {
          print("Making table ... ")
          dbRemoveTable(con, tableName)
          n <- length(rownames(finalDF2))
          seqN <- seq(from = 1, to = n, by = 100)
          for ( j in 1:length(seqN)) {
            msg <- paste0("\n", j , " out of ", length(seqN))
            cat(msg)
            start <- seqN[j]
            if (j == length(seqN)) {
              end <- n
            } else {
              end <- seqN[j+1]-1
            }
            uploadDF <- finalDF2[c(start:end), , drop=F]
            rownames(uploadDF) <- rownames(finalDF2)[c(start:end)]
            dbWriteTable(con, tableName, uploadDF, row.names = T,
                         field.types = c("row_names" = "VARCHAR(255)", "values" = "mediumtext"), append = T)
          }
          Sys.sleep(4)
          sql <- paste0("ALTER TABLE ", tableName, " ADD PRIMARY KEY (row_names);")
          DBI::dbExecute(conn = con, statement = sql)
          DBI::dbDisconnect(con)
          
        }
        # DBI::dbDisconnect(con)
        print("DONE")
      }
      print(type)
      #  samples
      Samps <- groups[[type]]
      colDataNow <- colData[which(colData$samples %in% Samps),]
      colDataNow <- correctSeries(colDataNow)
      Samps <- colDataNow$samples
      n <- length(Samps)
      print(paste0(type, " - ", n))
      if (n < 30 | length(unique(colDataNow$series)) < 4) {
        print("Not enough samples")
      } else {
        Path <- file.path(namePath, type)
        dir.create(path = Path)
        
        
        if (! file.exists(file.path(Path, "corMat.RData"))) {
          print("Correlations not present ... ")
          break
        } else {
          print("Correlations detected .. continuing")
          if (! file.exists(file.path(Path, "corMatForUpload.RData"))) {
            print("Loading correlations matrix")
            load(file.path(Path, "corMat.RData"))
          }
          
        }
        
        if (! file.exists(file.path(Path, "corMatForUpload.RData"))) {
          # Convert to compressed character format
          print("Converting for upload")
          
          forUpload <- cbind(as.data.frame(rownames(corMat)), corMat)
          cores <- 20
          colnames(forUpload)[1] <- "geneName"
          nr <- nrow(forUpload)
          n <- ceiling(nr/cores)
          forUploadList <- split(forUpload, rep(1:ceiling(nr/n), each=n, length.out=nr))
          k <- length(forUploadList)
          print("running foreach")
          start <- proc.time()
          resList <- foreach(j = 1:k, .verbose = F, .export = c("k", "forUploadList")) %dopar% convertForUpload(forUploadList[[j]])
          finalDF <- data.table::rbindlist(resList)
          ans <- proc.time() - start
          print(ans)
          print("Saving")
          finalDF2 <- finalDF[,c(-1), drop = F]
          rownames(finalDF2) <- finalDF$geneName
          save(finalDF2, file = file.path(Path, "corMatForUpload.RData"))
        } else {
          print("Correlations detected ... proceeding to upload")
          # Check whether upload necessary before loading
          tableName <- paste0("correlations_", speciesScien, "_", type, "_", name)
          checkBool <- uploadToAWS(tableName = tableName,
                                   check = T, geneNames = geneNames)
          if (checkBool) {
            print("Loading matrix for upload ... ")
            if (file.exists(file.path(Path, "corMatForUpload.RData"))) {
              load(file.path(Path, "corMatForUpload.RData"))
            } else {
              print("Not found -- please run correlation calculations ...  continuing")
              next
            }
            
            
          } 
          
        }
        tableName <- paste0("correlations_", speciesScien, "_", type, "_", name)
        
        checkBool <- uploadToAWS(tableName = tableName, check = T, geneNames = geneNames)
        if (checkBool) {
          # Upload to AWS
          print("uploading to AWS ... ")
          uploadToAWS(finalDF2, tableName, geneNames = geneNames)
          
          
        } else {
          print("Table already uploaded ... ")
        }
        
        
        Sys.sleep(3)
        cat("\nNext sample ... \n")
        
        
      }
      
    }
    
  }
  
}


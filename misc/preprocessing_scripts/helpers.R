# Categorizes Bulk RNA-Seq metadata using custom regex dictionary
categorizeMetaData <- function(metadata, cols, dictionary) {
  for (i in 1:length(names(dictionary))) {
    termNow <- names(dictionary)[i]
    dictNow <- dictionary[[i]]
    posInd <- c()
    negInd <- c()
    metadata$newColNow <- 0
    if ("yes" %in% names(dictNow)) {
      for (j in 1:length(cols)) {
        colNow <- cols[j]
        colIndNow <- which(colnames(metadata) == colNow)
        posInd <- c(posInd, grep(x = metadata[,colIndNow], 
                                 pattern = paste0(dictNow$yes, collapse = "|"), 
                                 perl = T, ignore.case = T))
        if (is.null(dictNow$no)) {
          negInd <- c(negInd)
        } else {
          negInd <- c(negInd, grep(x = metadata[,colIndNow], pattern = paste0(dictNow$no, collapse = "|"), perl = T, ignore.case = T))
        }
      }
      posInd <- unique(posInd[! posInd %in% negInd])
    } else {
      for (j in 1:length(cols)) {
        colNow <- cols[j]
        colIndNow <- which(colnames(metadata) == colNow)
        posInd <- c(posInd, grep(x = metadata[,colIndNow], pattern = paste0(dictNow, collapse = "|"), perl = T, ignore.case = T))
      }
    }
    metadata$newColNow[posInd] <- 1
    colnames(metadata)[which(colnames(metadata) == "newColNow")] <-  termNow
  }
  return(metadata)
}

# Finds number of samples gene is expressed in
nonZeroSamps <- function(row) {
  return(sum(row > 0))
}

# Pretty timestamp
timestamp2 <- function() {
  # Partially from https://stackoverflow.com/questions/1962278/dealing-with-timestamps-in-r
  now <- Sys.time()
  timeList <- unclass(as.POSIXlt(now))
  secNow <- ifelse(round(timeList$sec) < 10, paste0(0, round(timeList$sec)), round(timeList$sec))
  minNow <- ifelse(round(timeList$min) < 10, paste0(0, round(timeList$min)), round(timeList$min))
  paste0("[", timeList$hour,":",minNow,":",secNow,
         " ", (timeList$mon+1), "/", timeList$mday, "/", (timeList$year + 1900), "]", sep = "")
}




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

uploadToAzureSampleKey <- function(finalDF2= NULL, tableName,
                          credentials,
                          check = F) {
  
  # # bug testing
  # check = F
  # finalDF2 <- finalDF
  # tableName = "tpm_hsapiens"
  # credentials <- "Data/credFile.txt"
  
  
  require(DBI)
  
  # credentials <- "Data/credFile.txt"
  
  credentials <- suppressWarnings(read.delim(credentials, sep = ";",
                                             header = FALSE, stringsAsFactors = F))
  uName <- credentials$V1
  pWord <- credentials$V2
  conn <- dbConnect(drv = RMySQL::MySQL(), user = uName, 
                    port = 3306, dbname="correlation_analyzer",
                    password=pWord,
                    host="m2600az-db01p.mysql.database.azure.com")
  
  sql <- paste0("SELECT COUNT(*) FROM ", tableName)
  tabs <- dbListTables(conn)
  if(! tableName %in% tabs & check) {
    return(TRUE)
  } else if (! tableName %in% tabs) {
    countVal <- 0
  } else {
    res <- dbSendQuery(conn, sql)
    res <- dbFetch(res)
    countVal <- res[1,1]
  }
  
  if (! is.null(finalDF2)) {
    n <- nrow(finalDF2)
  } else {
    countVal <- 1.34
  }
  
  if (countVal != n) {
    if (check) {
      return(TRUE)
    }
    print("Making table ... ")
    if (tableName %in% tabs) {
      dbRemoveTable(conn, tableName)
    }
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
      
      # tryCatch({dbExecute(conn, sql)}, function(e) "hey")
      # testQuery <- dbExecute(conn, sql)
      done <- 0
      while (done  == 0) {
        done <- tryCatch(expr = {
          dbWriteTable(conn = conn, tableName, uploadDF, row.names = T,
                       field.types = c("row_names" = "VARCHAR(255)", "values" = "mediumtext"), 
                       append = T)
          1
        },
        error = function(e) {
          cat("Fail -- retry with new connection")
          rm(conn)
          gc()
          conn <- dbConnect(drv = RMySQL::MySQL(), user = uName, 
                            port = 3306, dbname="correlation_analyzer",
                            password=pWord,
                            host="m2600az-db01p.mysql.database.azure.com")
          0
        })
        print(done)
      }
    }
    Sys.sleep(4)
    # sql <- paste0("ALTER TABLE ", tableName, " ADD PRIMARY KEY (row_names);")
    # dbExecute(conn = conn, statement = sql)
  } else {
    if (check) {
      return(FALSE)
    }
    print("Table already filled ... continuing")
  }
  print("DONE")
  if (length(dbListConnections(drv = RMySQL::MySQL())) != 0) {
    lapply(dbListConnections(drv = RMySQL::MySQL()), dbDisconnect)
  }
    
}

# Same but with genes as the key
uploadToAzureGeneKey <- function(finalDF2= NULL, tableName,
                                   credentials,
                                   check = F) {
  
  # # bug testing
  # check = F
  # finalDF2 <- finalDF
  # tableName = "tpm_hsapiens"
  # credentials <- "Data/credFile.txt"
  
  
  require(DBI)
  
  # credentials <- "Data/credFile.txt"
  
  credentials <- suppressWarnings(read.delim(credentials, sep = ";",
                                             header = FALSE, stringsAsFactors = F))
  uName <- credentials$V1
  pWord <- credentials$V2
  conn <- dbConnect(drv = RMySQL::MySQL(), user = uName, 
                    port = 3306, dbname="correlation_analyzer",
                    password=pWord,
                    host="m2600az-db01p.mysql.database.azure.com")
  
  sql <- paste0("SELECT COUNT(*) FROM ", tableName)
  tabs <- dbListTables(conn)
  if(! tableName %in% tabs & check) {
    return(TRUE)
  } else if (! tableName %in% tabs) {
    countVal <- 0
  } else {
    res <- dbSendQuery(conn, sql)
    res <- dbFetch(res)
    countVal <- res[1,1]
  }
  
  if (! is.null(finalDF2)) {
    n <- nrow(finalDF2)
  } else {
    countVal <- 1.34
  }
  
  if (countVal != n) {
    if (check) {
      return(TRUE)
    }
    print("Making table ... ")
    if (tableName %in% tabs) {
      dbRemoveTable(conn, tableName)
    }
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
      
      # tryCatch({dbExecute(conn, sql)}, function(e) "hey")
      # testQuery <- dbExecute(conn, sql)
      done <- 0
      while (done  == 0) {
        done <- tryCatch(expr = {
          dbWriteTable(conn = conn, tableName, uploadDF, row.names = T,
                       field.types = c("row_names" = "VARCHAR(255)", "values" = "mediumtext"), 
                       append = T)
          1
        },
        error = function(e) {
          cat("Fail -- retry with new connection")
          rm(conn)
          gc()
          conn <- dbConnect(drv = RMySQL::MySQL(), user = uName, 
                            port = 3306, dbname="correlation_analyzer",
                            password=pWord,
                            host="m2600az-db01p.mysql.database.azure.com")
          0
        })
        print(done)
      }
    }
    Sys.sleep(4)
    # sql <- paste0("ALTER TABLE ", tableName, " ADD PRIMARY KEY (row_names);")
    # dbExecute(conn = conn, statement = sql)
  } else {
    if (check) {
      return(FALSE)
    }
    print("Table already filled ... continuing")
  }
  print("DONE")
  if (length(dbListConnections(drv = RMySQL::MySQL())) != 0) {
    lapply(dbListConnections(drv = RMySQL::MySQL()), dbDisconnect)
  }
  
}

# TPM upload table
makeUploadTable <- function(tx2gene, transcripts, 
                            dataFile, samplesInd) {
  resList <- list()
  # geneNames <- geneNames[c(1:2)]
  for (j in 1:length(samplesInd)) {
    # gene <- geneNames[j]
    sample <- names(samplesInd)[j]
    sampleInd <- samplesInd[j]
    print(j)
    # gene <- "AC012454.4"
    print(sample)
    tx2geneNow <- tx2gene[which(tx2gene$tx %in% transcripts),]
    transcriptsInd <- which(transcripts %in% tx2geneNow$tx)
    if (! length(transcriptsInd)) {
      resList[[j]] <- NULL
      next
    }
    H5close()
    expression <- as.data.frame(h5read(dataFile, "data/expression", 
                         index = list(transcriptsInd, sampleInd)))
    H5close()
    expression$tx <- transcripts[transcriptsInd]
    expressionMerge <- merge(x = expression, y = tx2gene, by = "tx")
    geneExp <- aggregate(
      expressionMerge$V1, by = list(expressionMerge$gene), FUN = sum
    )
    rownames(geneExp) <- geneExp$Group.1
    vec <- signif(as.numeric(geneExp$x), digits = 4)
    vecStr <- paste(vec, collapse = ",")
    resList[[j]] <- data.frame(sampleName = sample, values = vecStr)
  }
  resDF <- data.table::rbindlist(resList)
  return(resDF)
}


# VSD and TPM upload table
makeUploadTableExp <- function(upNow) {
  geneVec <- rownames(upNow)
  resNow <- apply(upNow, MARGIN = 1, FUN = function(row) {
    row <- signif(as.numeric(row), digits = 4)
    vecStr <- paste(row, collapse = ",")
    vecStr
  })
  resDF <- data.frame(geneName = geneVec, values = resNow)
  return(resDF)
}


library(sva)
library(BiocParallel)
library(foreach)
library(data.table)
library(DBI)
library(future)

uploadToAWS <- function(finalDF2= NULL, tableName, check = F, geneNames) {
  con <- dbConnect(RMySQL::MySQL(), user = "...", 
                   port = 3306, dbname="bishoplabdb",
                   password="...", 
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
                   field.types = c("row_names" = "VARCHAR(255)",
                                   "values" = "mediumtext"), append = T)
    }
    Sys.sleep(4)
    sql <- paste0("ALTER TABLE ", tableName, " ADD PRIMARY KEY (row_names);")
    DBI::dbExecute(conn = con, statement = sql)
    DBI::dbDisconnect(con)
    
  }
  # DBI::dbDisconnect(con)
  print("DONE")
}



for (species in c("Human", "Mouse")) {
  print(species)
  speciesScien <- ifelse(test = species == "Human", yes = "hsapiens", no = "mmusculus")
  print("Loading data")
  dataDir <- "Data/Final_Expression_Matrices/"
  load(file.path(dataDir, species, "geneTPM_forUpload_geneList.RData"))
  load(file.path(dataDir, species, "geneTPM_forUpload.RData"))  
  if (species == "Human") {
    genesTPM <- humanGenesTPM
  } else {
    genesTPM <- mouseGenesTPM
  }
  tableName <- paste0("TPM_", speciesScien)
  checkBool <- uploadToAWS(tableName = tableName, check = T, geneNames = genesTPM)
  
  if (checkBool) {
    # Upload to AWS
    print("uploading to AWS ... ")
    finalDF2 <- finalDF[,c(-1)]
    rownames(finalDF2) <- finalDF$geneName
    uploadToAWS(finalDF2 = finalDF2, tableName, geneNames = geneNames)
  } else {
    print("Table already uploaded ... ")
  }
  Sys.sleep(3)
  cat("\nNext sample ... \n")
  
}


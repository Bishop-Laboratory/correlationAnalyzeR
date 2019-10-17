library(ontologyIndex)
library(rhdf5)
library(preprocessCore)
library(sva)
library(jsonlite)


tissueDict <- list("brain" = list("yes" = c("cortex", "brain", "lobe", "hippoc",
                                            "alzheim", "frontal", "nerv", "gangli", "bipol", "medull",
                                            "hippocamp", "neur", "glia", "amygdala", "oligodendro",
                                            "spine", "spinal", "astrocyt", "cereb"),
                                  "no" = c("liver", "kidney")),
                   "thyroid" = c("thyroid"),
                   "blood" = list("yes" = c("blood", "eryth"),
                                  "no" = c("immune", "macroph", "spleen", "leuk", "marrow", "killer",
                                           "NKT", "thymus",
                                           "b cell", "t cell", "cd[0-9]+", "monocyt", "dendrit",
                                           "lympho", "mononucle", "pbmc", "neutro", "treg")),
                   "repiratory" = c("lung", "airway", "nasal",
                                    "trach", "pleura", "alveol", "bronch"),
                   "skin" = c("skin", "keratin", "derm", "melano", "psoriasis"),
                   "pancreas" = c("pancreas", "pdac", "pancrea",
                                  "islet", "alpha", "beta", "delta", "epsilon"),
                   "kidney" = c("kidney", "nephr", "glomerul", "renal", "clear cell"),
                   "fetal" = c("placent", "fetal", "huvec", "decidua", "germ",
                               "fetus", "embry", "umbil", "cord blood"),
                   "cartilage" = c("cartilag", "chondr", "joint"),
                   "mammary" = c("mammary", "breast", "imec"),
                   "stomach"= c("stomach", "gastric"),
                   "esophagus" = c("esophag"),
                   "intestines" = c("instestine", "colo", "duoden",
                                    "ileum", "gut", "bowel", "jejun"),
                   "muscle" = c("muscle", "myo", "lateralis", "gastrocnemius",
                                "skeletal", "brach", "satellite", "ceps"),
                   "liver" = c("liver", "hepat", "kupffer"),
                   "fibroblasts" = c("fibroblas", "IMR90"),
                   "adipose" = c("adipose", "fat", "adipo"),
                   "stem cells" = list("yes" = c("stem", "progen", "blast", "hesc",
                                                 "npsc", "ipsc", "poetic", "ips cell", "satellite"),
                                       "no" = c("fibroblast")),
                   "cardiac" = c("cardiac", "heart", "atria",
                                 "coron", "aort", "ventric"),
                   "endothelial" = c("endoth", "huvec", "vascul"),
                   "spleen" = c("spleen", "splen"),
                   "bladder" = c("bladder", "urin", "urothe"),
                   "retina" = c("retina", "macular", "retin", "photo"),
                   "thymus" = c("thymus", "thymic"),
                   "male reproductive" = c("testi", "epidid", "gonad", "spermat", "prost"),
                   "female reproductive" = c("ovar", "oocy", "uter", "cervic", "vagi"),
                   "immune" = c("immune", "macroph", "leuk", "marrow", "killer",
                                "NKT",
                                "b cell", "t cell", "cd4", "cd8", "monocyt", "dendrit",
                                "lympho", "mononucle", "pbmc", "neutro", "treg"),
                   "bone" = c("femur", "osteo", "mandible", "bone", "joint", "ewing"))

jsonlite::write_json(x = tissueDict, path = "Code/tissueDictionary.json")
tissueDict <- jsonlite::read_json("Code/tissueDictionary.json")



# Human
load("Data/Final_Expression_Matrices/Human/colData.RData")
tissue <- colData$tissue
# Re-format tissue
newTissue = gsub(pattern = "_", replacement = "", x = tissue)
newTissue = gsub(pattern = "\\.", replacement = "", x = newTissue)
newTissue = gsub(pattern = "'", replacement = "", x = newTissue)
newTissue = gsub(pattern = "-", replacement = "", x = newTissue)
newTissue = tolower(newTissue)
colData$tissue <- newTissue
# Remove scRNA-Seq samples 
patternSingle <- c("single cell", "single-cell", "smart seq", "smart-seq", "drop-seq", "drop seq",
                   "fluidigm", "tenX", "scRNASeq", "scRNA-Seq", "chromium")
singleSamps <- colData$samples[unique(grep(pattern = paste(patternSingle, collapse="|"),
                           x = colData$extractProtocol, ignore.case = T))]

colData <- colData[which(! colData$samples %in% c(singleSamps)),]
tDF <- as.data.frame(table(colData$tissue))
tDF <- tDF[order(tDF$Freq, decreasing = T),]

rmList <- c(patternSingle)
tDF <- tDF[unique(grep(pattern = paste(rmList, collapse="|"),
                       x = tDF$Var1, ignore.case = T, invert = T)),]


resFrame <- data.frame(samples = colData$samples, tissue = colData$tissue)
for (j in 1:length(tissueDict)) {
  terms <- tissueDict[[j]]
  name <- names(tissueDict)[j]
  resFrame$newCol <- 0
  if (! is.null(names(terms))) {
    yes <- terms[["yes"]]
    no <- terms[["no"]]
    print(yes)
    print(no)
    resFrame$newCol[grep(pattern = paste(yes, collapse = "|"),
                         x = resFrame$tissue, ignore.case = T)] <- 1
    resFrame$newCol[grep(pattern = paste(no, collapse = "|"),
                         x = resFrame$tissue, ignore.case = T)] <- 0
  } else {
    resFrame$newCol[grep(pattern = paste(terms, collapse = "|"),
                         x = resFrame$tissue, ignore.case = T)] <- 1
  }

  colnames(resFrame)[colnames(resFrame) == "newCol"] <- name
}
colSums(resFrame[,c(3:30)])
colData <- merge(x = colData, y = resFrame, by = c("samples", "tissue"))

# cancer
cancer = c("hela", "hek", "k562", "reh", "jurkat", "leukemi", "293", "bewo",
             "kras", "mcf", "lncap", "bjab", "gbm", "rko", "ramos", "mel888",
             "vcap", "saos2", "vapc", "nalm6", "set2", "tov21",
             "cancer", "carcin", "sarcom", "metasta", "tumor", "oma ", "oma$")

cl <- ontologyIndex::get_OBO("Code/bto.obo")
cancerLines <- cl$name
cancer <- c(cancer, 
             as.character(cancerLines[grep(x = cancerLines, ignore.case = T,
                                           pattern = paste0(cancer, collapse = "|"))]))
colData$cancer <- 0
colData$cancer[grep(x = colData$tissue, pattern = paste0(cancer, collapse = "|"),
                    ignore.case = T)] <- 1
colData$normal <- 0
colData$normal[which(colData$cancer == 0)] <- 1
groupList <- list()
for (i in 8:35) {
  name <- colnames(colData)[i]
  allsamps <- colData$samples[which(colData[,i] == 1)]
  cancersamps <- colData$samples[which(colData[,i] == 1 & colData[,36] == 1)]
  normSamps <- colData$samples[which(colData[,i] == 1 & colData[,37] == 1)]
  listInd <- i-6
  groupList[[listInd]] <- list("all" = allsamps,
                               "cancer" = cancersamps,
                               "normal" = normSamps)
  names(groupList)[listInd] <- name
}
groupList[["all"]] <- list("all" = colData$samples,
                           "cancer" = colData$samples[which(colData$cancer == 1)],
                           "normal" = colData$samples[which(colData$normal == 1)])
groupList <- groupList[which(! is.na(names(groupList)))]
save(groupList, colData, file = "Data/Final_Expression_Matrices/Human/groupList.RData")

# Make table and upload to AWS
dfList <- list()
for (i in 1:length(groupList)) {
  tissue <- names(groupList)[i]
  listTemp <- groupList[[i]]
  listNew <- lapply(listTemp, FUN = function(x){
    tmp = paste0(x, collapse = ",")
    return(tmp)
  })
  dfList[[i]] <- data.frame(sampleType = names(listNew),
                            samples = unlist(listNew, use.names = F))
  names(dfList)[[i]] <- tissue
  
}
tissueSampleDF <- data.table::rbindlist(dfList)
tissueSampleDF$tissue <- rep(names(dfList), each = 3)
tissueSampleDF$tissue <- gsub(tissueSampleDF$tissue, pattern = " ", replacement = "0")
tissueSampleDF2 <- tissueSampleDF[,c(-1, -3)]
rownames(tissueSampleDF2) <- paste0(tissueSampleDF$sampleType, "_", tissueSampleDF$tissue)
con <- DBI::dbConnect(RMySQL::MySQL(), user = "millerh1", 
                 port = 3306, dbname="bishoplabdb",
                 password='...', 
                 host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
DBI::dbRemoveTable(conn = con, name = "hsapiens_sample_group_key")
DBI::dbWriteTable(con, "hsapiens_sample_group_key", tissueSampleDF2, 
                  row.names = T, overwrite = T,
             field.types = c("row_names" = "VARCHAR(255)", 
                             "samples" = "mediumtext"), append = F)
sql <- paste0("ALTER TABLE ", "hsapiens_sample_group_key", 
              " ADD PRIMARY KEY (row_names);")
DBI::dbExecute(conn = con, statement = sql)
DBI::dbDisconnect(con)



# Mouse
load("Data/Final_Expression_Matrices/Mouse/colData.RData")
tissue <- colData$tissue
# Re-format tissue
newTissue = gsub(pattern = "_", replacement = "", x = tissue)
newTissue = gsub(pattern = "\\.", replacement = "", x = newTissue)
newTissue = gsub(pattern = "'", replacement = "", x = newTissue)
newTissue = gsub(pattern = "-", replacement = "", x = newTissue)
newTissue = tolower(newTissue)
colData$tissue <- newTissue
# Remove scRNA-Seq samples 
patternSingle <- c("single cell", "single-cell", "smart seq", "smart-seq", "drop-seq", "drop seq",
                   "fluidigm", "tenX", "scRNASeq", "scRNA-Seq", "chromium")
singleSamps <- colData$samples[unique(grep(pattern = paste(patternSingle, collapse="|"),
                                           x = colData$extractProtocol, ignore.case = T))]
patternXeno <- c("xeno", "patient", "pdx", "cdx", "graft")
xenoSamps <- colData$samples[unique(grep(pattern = paste(patternXeno, collapse="|"),
                                           x = colData$extractProtocol, ignore.case = T))]
colData <- colData[which(! colData$samples %in% c(singleSamps, xenoSamps)),]
tDF <- as.data.frame(table(colData$tissue))
tDF <- tDF[order(tDF$Freq, decreasing = T),]

rmList <- c(patternSingle, xenoSamps)
tDF <- tDF[unique(grep(pattern = paste(rmList, collapse="|"),
                       x = tDF$Var1, ignore.case = T, invert = T)),]


resFrame <- data.frame(samples = colData$samples, tissue = colData$tissue)
for (j in 1:length(tissueDict)) {
  terms <- tissueDict[[j]]
  name <- names(tissueDict)[j]
  resFrame$newCol <- 0
  if (! is.null(names(terms))) {
    yes <- terms[["yes"]]
    no <- terms[["no"]]
    print(yes)
    print(no)
    resFrame$newCol[grep(pattern = paste(yes, collapse = "|"),
                         x = resFrame$tissue, ignore.case = T)] <- 1
    resFrame$newCol[grep(pattern = paste(no, collapse = "|"),
                         x = resFrame$tissue, ignore.case = T)] <- 0
  } else {
    resFrame$newCol[grep(pattern = paste(terms, collapse = "|"),
                         x = resFrame$tissue, ignore.case = T)] <- 1
  }
  
  colnames(resFrame)[colnames(resFrame) == "newCol"] <- name
}
colSums(resFrame[,c(3:30)])
colData <- merge(x = colData, y = resFrame, by = c("samples", "tissue"))

# cancer
cancer = c("hela", "hek", "k562", "reh", "jurkat", "leukemi", "293", "bewo",
           "kras", "mcf", "lncap", "bjab", "gbm", "rko", "ramos", "mel888",
           "vcap", "saos2", "vapc", "nalm6", "set2", "tov21",
           "cancer", "carcin", "sarcom", "metasta", "tumor", "oma ", "oma$")

cl <- ontologyIndex::get_OBO("Code/bto.obo")
cancerLines <- cl$name
cancer <- c(cancer, 
            as.character(cancerLines[grep(x = cancerLines, ignore.case = T,
                                          pattern = paste0(cancer, collapse = "|"))]))
colData$cancer <- 0
colData$cancer[grep(x = colData$tissue, pattern = paste0(cancer, collapse = "|"),
                    ignore.case = T)] <- 1
colData$normal <- 0
colData$normal[which(colData$cancer == 0)] <- 1
groupList <- list()
for (i in 7:34) {
  name <- colnames(colData)[i]
  allsamps <- colData$samples[which(colData[,i] == 1)]
  cancersamps <- colData$samples[which(colData[,i] == 1 & colData[,35] == 1)]
  normSamps <- colData$samples[which(colData[,i] == 1 & colData[,36] == 1)]
  listInd <- i-6
  groupList[[listInd]] <- list("all" = allsamps,
                               "cancer" = cancersamps,
                               "normal" = normSamps)
  names(groupList)[listInd] <- name
}
groupList[["all"]] <- list("all" = colData$samples,
                           "cancer" = colData$samples[which(colData$cancer == 1)],
                           "normal" = colData$samples[which(colData$normal == 1)])
groupList <- groupList[which(! is.na(names(groupList)))]
save(groupList, colData, file = "Data/Final_Expression_Matrices/Mouse/groupList.RData")

# Make table and upload to AWS
dfList <- list()
for (i in 1:length(groupList)) {
  tissue <- names(groupList)[i]
  listTemp <- groupList[[i]]
  listNew <- lapply(listTemp, FUN = function(x){
    tmp = paste0(x, collapse = ",")
    return(tmp)
  })
  dfList[[i]] <- data.frame(sampleType = names(listNew),
                            samples = unlist(listNew, use.names = F))
  names(dfList)[[i]] <- tissue
  
}
tissueSampleDF <- data.table::rbindlist(dfList)
tissueSampleDF$tissue <- rep(names(dfList), each = 3)
tissueSampleDF$tissue <- gsub(tissueSampleDF$tissue, pattern = " ", replacement = "0")
tissueSampleDF2 <- tissueSampleDF[,c(-1, -3)]
rownames(tissueSampleDF2) <- paste0(tissueSampleDF$sampleType, "_", tissueSampleDF$tissue)
con <- DBI::dbConnect(RMySQL::MySQL(), user = "millerh1", 
                      port = 3306, dbname="bishoplabdb",
                      password='', 
                      host="bishoplabdb.cyss3bq5juml.us-west-2.rds.amazonaws.com")
DBI::dbRemoveTable(conn = con, name = "mmusculus_sample_group_key")
DBI::dbWriteTable(con, "mmusculus_sample_group_key", tissueSampleDF2, 
                  row.names = T, overwrite = T,
                  field.types = c("row_names" = "VARCHAR(255)", 
                                  "samples" = "mediumtext"), append = F)
sql <- paste0("ALTER TABLE ", "mmusculus_sample_group_key", 
              " ADD PRIMARY KEY (row_names);")
DBI::dbExecute(conn = con, statement = sql)
DBI::dbDisconnect(con)




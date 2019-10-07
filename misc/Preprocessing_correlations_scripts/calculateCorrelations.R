library(sva)
library(BiocParallel)
library(foreach)
library(data.table)
library(DBI)
#library(future)
#options(future.globals.maxSize= 5e10)
#plan(multicore(workers = 10))


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

for (species in c("Human", "Mouse")) {
  cat("Species is ", species)
  speciesScien <- ifelse(test = species == "Human", yes = "hsapiens", no = "mmusculus")
  cat("Loading data ... \n")
  dataDir <- "Data/Final_Expression_Matrices/"
  load(file.path(dataDir, species, "colData.RData"))
  load(file.path(dataDir, species, "groupList.RData"))

  cat("Starting analysis ... \n")
  for (i in 1:length(groupList)) {
    groups <- groupList[[i]]
    name <- names(groupList)[i]
    cat("Current group: ", name, "\n")
    namePath <- file.path("Data/Final_Expression_Matrices/", species, name)
    for (type in c("all", "cancer", "normal")) {
      #  samples
      Samps <- groups[[type]]
      colDataNow <- colData[which(colData$samples %in% Samps),]
      colDataNow <- correctSeries(colDataNow)
      Samps <- colDataNow$samples
      n <- length(Samps)
      cat(paste0("\n", type, " - ", n), "\n")
      if (n < 30 | length(unique(colDataNow$series)) < 4) {
        warning("Not enough samples")
      } else {
        Path <- file.path(namePath, type)
        dir.create(path = Path)

        if (! file.exists(file.path(Path, "corMat.RData"))) {
          if (! file.exists(file.path(Path, "correctedCounts.RData"))) {
            stop("No batch corrected counts available for ", species, " ",
                 name, "-", type, " despite eligibility for analysis.",
                 " \nRe-run batchCorrect.R to supply missing samples")
          } else {
            load(file.path(Path, "correctedCounts.RData"))
          }
          cat("Calculating correlations ... \n")
          start <- proc.time()
          corMat <- WGCNA::cor(x = t(batchCorrected), verbose = 1)
          ans <- proc.time() - start
          cat("Done!\n")
          cat("Time taken: \n", ans, "\n")
          cat("Saving ... \n")
          save(corMat, colDataNow, file = file.path(Path, "corMat.RData"))
#          future(save(corMat, colDataNow, file = file.path(Path, "corMat.RData")))
        } else {
          cat("Correlations detected .. continuing \n")
          
          
        }
        
        
        
      }
      
    }
    
  }
  
}


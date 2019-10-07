library(sva)
library(data.table)


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
  print(species)
  speciesScien <- ifelse(test = species == "Human", yes = "hsapiens", no = "mmusculus")
  cat("Loading data ... \n")
  dataDir <- "Data/Final_Expression_Matrices/"
  load(file.path(dataDir, species, "colData.RData"))
  load(file.path(dataDir, species, "groupList.RData"))
  load(file.path(dataDir, species, "normCounts.RData"))
  geneNames <- rownames(normCounts)
  save(geneNames, 
       file = file.path(file.path("Data/Final_Expression_Matrices/", 
                                  species, "geneNames.RData")))
  cat("Starting analysis ... \n")
  for (i in 1:length(groupList)) {
    groups <- groupList[[i]]
    name <- names(groupList)[i]
    cat("Current group: ", name, "\n")
    namePath <- file.path("Data/Final_Expression_Matrices/", species, name)
    if (! dir.exists(namePath)) {
      dir.create(path = namePath)
    }
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
        if (! dir.exists(Path)) {
          dir.create(path = Path)
        }
        if (! file.exists(file.path(Path, "correctedCounts.RData"))) {
          cat("Batch correcting ... \n")
          countsForCorrect <- normCounts[,which(colnames(normCounts) %in% colDataNow$samples)]
          countsForCorrect[is.na(countsForCorrect)] <- 0
          start <- proc.time()
          batchCorrected <- ComBat(dat = countsForCorrect, batch = colDataNow$series,
                                   par.prior = T, prior.plots = F)
          ans <- proc.time() - start
          cat("Saving ... \n")
          save(colDataNow, batchCorrected, file = file.path(Path, "correctedCounts.RData"))
        } else {
          cat("Detected corrected counts .. continuing \n")
        }
      }
      
    }
  }
}


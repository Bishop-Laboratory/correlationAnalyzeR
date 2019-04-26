getAvailableGenes <- function(Species = c("hsapiens", "mmusculus")) {
  downloadData(Species = Species,
               sampleType = "Normal_Tissue",
               geneList = "A1BG")
  downloadFolder <- system.file("data", package = "correlationAnalyzeR")
  if (Species[1] == "hsapiens") {
    geneFile <- "A1BG.RData"
  } else if (Species[1] == "mmusculus") {
    geneFile <- "A1bg.RData"
  } else {
    stop("\ncorrelationAnalyzeR currently supports only Human and Mouse data.
         Please select either 'hsapiens' or 'mmusculus' for Species parameter.
         \n")
  }
  file <- file.path(downloadFolder, "Correlation_Data",
                    Species, sampleType, geneFile)
  load(file)
  geneNames <- names(vec)
  geneNamesDF <- as.data.frame(geneNames)
  colnames(geneNamesDF)[1] <- "geneName"
  if (Species[1] == "hsapiens") {
    load(file.path(downloadFolder, "biomaRtData", "humanGenes.RData"))
    colnames(humanGenes)[1] <- "geneName"
    geneNamesDF <- merge(x = humanGenes, y = geneNamesDF, by = "geneName",
                         all.y = T)
  } else if (Species[1] == "mmusculus") {
    load(file.path(downloadFolder, "biomaRtData", "mouseGenes.RData"))
    colnames(MouseGenes)[1] <- "geneName"
    geneNamesDF <- merge(x = MouseGenes, y = geneNamesDF, by = "geneName",
                         all.y = T)
  }
  return(geneNamesDF)
}

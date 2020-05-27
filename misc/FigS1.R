
# Get the ARCHS4 correlation data
load("Data/human_correlation.rda")
BRCA1_ARCHS4 <- cc["BRCA1",]
BRCA1_ARCHS4 <- data.frame(geneName = names(BRCA1_ARCHS4),
                           ARCHS4 = BRCA1_ARCHS4)

library(org.Hs.eg.db) 

aa <- listAttributes(ensembl)
entrez2gene <- getBM()
file <- list.files("Data/Hsa-u.v18-12.G26050-S164823.combat_pca_subagging.mrgeo.d/",
                   pattern = "^672$", full.names = T)
BRCA1_COXPR <- read.table(file)
info4gene = select(org.Hs.eg.db, keys = as.character(BRCA1_COXPR$V1), 
                   columns =  c("SYMBOL"))
BRCA1_COXPR <- data.frame(geneName = info4gene$SYMBOL, COXPRESdb = rev(BRCA1_COXPR$V2))


corCompare <- merge(x = BRCA1_ARCHS4, y = BRCA1_COXPR, by = "geneName")
BRCA1_me <- correlationAnalyzeR::getCorrelationData(Species = "hsapiens", Sample_Type = "normal",
                                                    Tissue = "all", geneList = "BRCA1")
BRCA1_me$geneName <- rownames(BRCA1_me)
colnames(BRCA1_me)[1] <- "correlationAnalyzeR"
corCompare <- merge(x = BRCA1_me, y = corCompare, by = "geneName")



library("PerformanceAnalytics")

chart.Correlation(corCompare[,c(-1)],   method = "spearman")

BRCA1_me <- BRCA1_me[order(BRCA1_me$correlationAnalyzeR, decreasing = TRUE),]
BRCA1_ARCHS4 <- BRCA1_ARCHS4[order(BRCA1_ARCHS4$ARCHS4, decreasing = TRUE),]
library(VennDiagram)
topList <- list(
  "correlationAnalyzeR" = as.character(BRCA1_me$geneName[c(1:100)]),
  "ARCHS4" = as.character(BRCA1_ARCHS4$geneName[c(1:100)]),
  "COXPRESdb" = as.character(BRCA1_COXPR$geneName[c(1:100)])
)
vd <- venn.diagram(topList, filename = NULL, fill = c("cornflowerblue", "salmon", "goldenrod"),
                   margin = .05)
grid.draw(vd)





# AURKB

AURKB_ARCHS4 <- cc["AURKB",]
AURKB_ARCHS4 <- data.frame(geneName = names(AURKB_ARCHS4),
                           ARCHS4 = AURKB_ARCHS4)

file <- list.files("Data/Hsa-u.v18-12.G26050-S164823.combat_pca_subagging.mrgeo.d/",
                   pattern = "^9212$", full.names = T)
AURKB_COXPR <- read.table(file)
info4gene = select(org.Hs.eg.db, keys = as.character(AURKB_COXPR$V1), 
                   columns =  c("SYMBOL"))
AURKB_COXPR <- data.frame(geneName = info4gene$SYMBOL, COXPRESdb = rev(AURKB_COXPR$V2))


corCompare <- merge(x = AURKB_ARCHS4, y = AURKB_COXPR, by = "geneName")
AURKB_me <- correlationAnalyzeR::getCorrelationData(Species = "hsapiens", Sample_Type = "normal",
                                                    Tissue = "all", geneList = "AURKB")
AURKB_me$geneName <- rownames(AURKB_me)
colnames(AURKB_me)[1] <- "correlationAnalyzeR"
corCompare <- merge(x = AURKB_me, y = corCompare, by = "geneName")


chart.Correlation(corCompare[,c(-1)],   method = "spearman")

AURKB_me <- AURKB_me[order(AURKB_me$correlationAnalyzeR, decreasing = TRUE),]
AURKB_ARCHS4 <- AURKB_ARCHS4[order(AURKB_ARCHS4$ARCHS4, decreasing = TRUE),]
library(VennDiagram)
topList <- list(
  "correlationAnalyzeR" = as.character(AURKB_me$geneName[c(1:100)]),
  "ARCHS4" = as.character(AURKB_ARCHS4$geneName[c(1:100)]),
  "COXPRESdb" = as.character(AURKB_COXPR$geneName[c(1:100)])
)
vd <- venn.diagram(topList, filename = NULL, fill = c("cornflowerblue", "salmon", "goldenrod"),
                   margin = .05)
grid.draw(vd)





# HSP90AA1

HSP90AA1_ARCHS4 <- cc["HSP90AA1",]
HSP90AA1_ARCHS4 <- data.frame(geneName = names(HSP90AA1_ARCHS4),
                           ARCHS4 = HSP90AA1_ARCHS4)

file <- list.files("Data/Hsa-u.v18-12.G26050-S164823.combat_pca_subagging.mrgeo.d/",
                   pattern = "^3320$", full.names = T)
HSP90AA1_COXPR <- read.table(file)
info4gene = select(org.Hs.eg.db, keys = as.character(HSP90AA1_COXPR$V1), 
                   columns =  c("SYMBOL"))
HSP90AA1_COXPR <- data.frame(geneName = info4gene$SYMBOL, COXPRESdb = rev(HSP90AA1_COXPR$V2))


corCompare <- merge(x = HSP90AA1_ARCHS4, y = HSP90AA1_COXPR, by = "geneName")
HSP90AA1_me <- correlationAnalyzeR::getCorrelationData(Species = "hsapiens", Sample_Type = "normal",
                                                    Tissue = "all", geneList = "HSP90AA1")
HSP90AA1_me$geneName <- rownames(HSP90AA1_me)
colnames(HSP90AA1_me)[1] <- "correlationAnalyzeR"
corCompare <- merge(x = HSP90AA1_me, y = corCompare, by = "geneName")


chart.Correlation(corCompare[,c(-1)],   method = "spearman")

HSP90AA1_me <- HSP90AA1_me[order(HSP90AA1_me$correlationAnalyzeR, decreasing = TRUE),]
HSP90AA1_ARCHS4 <- HSP90AA1_ARCHS4[order(HSP90AA1_ARCHS4$ARCHS4, decreasing = TRUE),]
library(VennDiagram)
topList <- list(
  "correlationAnalyzeR" = as.character(HSP90AA1_me$geneName[c(1:100)]),
  "ARCHS4" = as.character(HSP90AA1_ARCHS4$geneName[c(1:100)]),
  "COXPRESdb" = as.character(HSP90AA1_COXPR$geneName[c(1:100)])
)
vd <- venn.diagram(topList, filename = NULL, fill = c("cornflowerblue", "salmon", "goldenrod"),
                   margin = .05)
grid.draw(vd)
















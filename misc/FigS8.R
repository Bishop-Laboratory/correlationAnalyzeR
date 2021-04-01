# Script used to generate Figure S7
# It assumes you've downloaded the correlation matrix from ARCHS4: https://s3.amazonaws.com/mssm-seq-matrix/human_correlation.rda
# It also assumes you have downloaded and unzipped the correlations from COXPRESSdb: https://coxpresdb.jp/download/Hsa-u.c2-0/coex/
# It also assumes you have downloaded and decompressed the correlations from GeneFriends: "Pearson correlation genes" in https://genefriends.org/RNAseq/about/
# It also assumes you have used BioGRiD on BRCA1, AURKB, and HSP90AA1 and downloaded the interactors as a CSV.

# Get libraries / annotations
library(PerformanceAnalytics)
library(org.Hs.eg.db)
library(tidyverse)
library(VennDiagram)
library(biomaRt)
aa <- listAttributes(ensembl)
entrez2gene <- getBM()
source("https://raw.githubusercontent.com/Bishop-Laboratory/Ewing-sarcoma-paper-Miller-2020/master/helpers_v2.R")

## BRCA1 cor compare ##
# Get the ARCHS4 correlation data
load("Data/human_correlation.rda")
BRCA1_ARCHS4 <- cc["BRCA1",]
BRCA1_ARCHS4 <- data.frame(geneName = names(BRCA1_ARCHS4),
                           ARCHS4 = BRCA1_ARCHS4)

# Get the COXPRESdb data
file <- list.files("Data/Hsa-u.v18-12.G26050-S164823.combat_pca_subagging.mrgeo.d/",
                   pattern = "^672$", full.names = T)
BRCA1_COXPR <- read.table(file)
info4gene = select(org.Hs.eg.db, keys = as.character(BRCA1_COXPR$V1),
                   columns =  c("SYMBOL"))
BRCA1_COXPR <- data.frame(geneName = info4gene$SYMBOL, COXPRESdb = rev(BRCA1_COXPR$V2))

# Get the GeneFriends data
file <- list.files("Data/home/nis/priyanka/work/GF_2019_Uploads/Human/HumanGenes_Correlation/",
                   pattern = "ENSG00000012048", full.names = T)
BRCA1_GF <- read.table(file)
BRCA1_GF$ENSEMBL <- rownames(BRCA1_GF)
info4gene = select(org.Hs.eg.db, keys = as.character(rownames(BRCA1_GF)),
                   columns =  c("SYMBOL"), keytype = "ENSEMBL")
BRCA1_GF <- merge(x = info4gene, y = BRCA1_GF, by = "ENSEMBL")
BRCA1_GF <- unique(BRCA1_GF[which(! is.na(BRCA1_GF$SYMBOL)),c(-1)])
BRCA1_GF <- BRCA1_GF[order(BRCA1_GF[,2], decreasing = TRUE),]
colnames(BRCA1_GF) <- c("geneName", "GeneFriends")

corCompare <- merge(x = BRCA1_ARCHS4, y = BRCA1_COXPR, by = "geneName")
corCompare <- merge(x = corCompare, y = BRCA1_GF, by = "geneName")

BRCA1_me <- correlationAnalyzeR::getCorrelationData(Species = "hsapiens", Sample_Type = "normal",
                                                    Tissue = "all", geneList = "BRCA1")
BRCA1_me$geneName <- rownames(BRCA1_me)
colnames(BRCA1_me)[1] <- "correlationAnalyzeR"
corCompare <- merge(x = BRCA1_me, y = corCompare, by = "geneName")

## Compare BRCA1 top 100 to biogrid
biogrid <- read_tsv("misc/BIOGRID-GENE-107140-4.1.190.tab3.txt")

top50 <- BRCA1_me %>%
  top_n(wt = BRCA1, n = 500)

interacts <- biogrid %>%
  select(`Official Symbol Interactor B`) %>%
  distinct()

vl <- list(
  "Co-expression" = top50$geneName,
  "Interaction" = interacts$`Official Symbol Interactor B`
)
ol <- calculate.overlap(vl)
calculate.overlap.and.pvalue(list1 = vl$`Co-expression`, list2 = vl$Interaction,
                             total.size = unique(length(BRCA1_me$geneName)),
                             lower.tail = FALSE)
vd <- venn.diagram(vl, filename = NULL, fill = c("forestgreen", "firebrick"),
                   margin = .05)
dev.off()
grid.draw(vd)


## AURKB cor compare ##

# ARCHS4
AURKB_ARCHS4 <- cc["AURKB",]
AURKB_ARCHS4 <- data.frame(geneName = names(AURKB_ARCHS4),
                          ARCHS4 = AURKB_ARCHS4)
# COXPRESSdb
file <- list.files("Data/Hsa-u.v18-12.G26050-S164823.combat_pca_subagging.mrgeo.d/",
                   pattern = "^9212$", full.names = T)
AURKB_COXPR <- read.table(file)
info4gene = select(org.Hs.eg.db, keys = as.character(AURKB_COXPR$V1),
                   columns =  c("SYMBOL"))
AURKB_COXPR <- data.frame(geneName = info4gene$SYMBOL, COXPRESdb = rev(AURKB_COXPR$V2))

# GeneFriends
file <- list.files("Data/home/nis/priyanka/work/GF_2019_Uploads/Human/HumanGenes_Correlation/",
                   pattern = "ENSG00000178999", full.names = T)
AURKB_GF <- read.table(file)
AURKB_GF$ENSEMBL <- rownames(AURKB_GF)
info4gene = select(org.Hs.eg.db, keys = as.character(rownames(AURKB_GF)),
                   columns =  c("SYMBOL"), keytype = "ENSEMBL")
AURKB_GF <- merge(x = info4gene, y = AURKB_GF, by = "ENSEMBL")
AURKB_GF <- unique(AURKB_GF[which(! is.na(AURKB_GF$SYMBOL)),c(-1)])
AURKB_GF <- AURKB_GF[order(AURKB_GF[,2], decreasing = TRUE),]
colnames(AURKB_GF) <- c("geneName", "GeneFriends")

corCompare <- merge(x = AURKB_ARCHS4, y = AURKB_COXPR, by = "geneName")
corCompare <- merge(x = corCompare, y = AURKB_GF, by = "geneName")

# correlationAnalyzeR
AURKB_me <- correlationAnalyzeR::getCorrelationData(Species = "hsapiens", Sample_Type = "normal",
                                                   Tissue = "all", geneList = "AURKB")
AURKB_me$geneName <- rownames(AURKB_me)
colnames(AURKB_me)[1] <- "correlationAnalyzeR"
corCompare <- merge(x = AURKB_me, y = corCompare, by = "geneName")

# Plot
chart.Correlation(corCompare[,c(-1)],   method = "spearman")

## Compare AURKB top 100 to biogrid
biogrid <- read_tsv("misc/BIOGRID-GENE-114646-4.1.190.tab3.txt")

top50 <- AURKB_me %>%
  top_n(wt = AURKB, n = 500)

interacts <- biogrid %>%
  select(`Official Symbol Interactor B`) %>%
  distinct()
vl <- list(
  "Co-expression" = top50$geneName,
  "Interaction" = interacts$`Official Symbol Interactor B`
)
ol <- calculate.overlap(vl)
calculate.overlap.and.pvalue(list1 = vl$`Co-expression`, list2 = vl$Interaction,
                             total.size = unique(length(BRCA1_me$geneName)),
                             lower.tail = FALSE)
vd <- venn.diagram(vl, filename = NULL, fill = c("forestgreen", "firebrick"),
                   margin = .05)

dev.off()
grid.draw(vd)


## HSP90AA1 ##
# ARCHS4
HSP90AA1_ARCHS4 <- cc["HSP90AA1",]
HSP90AA1_ARCHS4 <- data.frame(geneName = names(HSP90AA1_ARCHS4),
                              ARCHS4 = HSP90AA1_ARCHS4)
# COXPRESDB
file <- list.files("Data/Hsa-u.v18-12.G26050-S164823.combat_pca_subagging.mrgeo.d/",
                   pattern = "^3320$", full.names = T)
HSP90AA1_COXPR <- read.table(file)
info4gene = select(org.Hs.eg.db, keys = as.character(HSP90AA1_COXPR$V1),
                   columns =  c("SYMBOL"))
HSP90AA1_COXPR <- data.frame(geneName = info4gene$SYMBOL, COXPRESdb = rev(HSP90AA1_COXPR$V2))

# Get the GeneFriends data
file <- list.files("Data/home/nis/priyanka/work/GF_2019_Uploads/Human/HumanGenes_Correlation/",
                   pattern = "ENSG00000080824", full.names = T)
HSP90AA1_GF <- read.table(file)
HSP90AA1_GF$ENSEMBL <- rownames(HSP90AA1_GF)
info4gene = select(org.Hs.eg.db, keys = as.character(rownames(HSP90AA1_GF)),
                   columns =  c("SYMBOL"), keytype = "ENSEMBL")
HSP90AA1_GF <- merge(x = info4gene, y = HSP90AA1_GF, by = "ENSEMBL")
HSP90AA1_GF <- unique(HSP90AA1_GF[which(! is.na(HSP90AA1_GF$SYMBOL)),c(-1)])
HSP90AA1_GF <- HSP90AA1_GF[order(HSP90AA1_GF[,2], decreasing = TRUE),]
colnames(HSP90AA1_GF) <- c("geneName", "GeneFriends")

corCompare <- merge(x = HSP90AA1_ARCHS4, y = HSP90AA1_COXPR, by = "geneName")
corCompare <- merge(x = corCompare, y = HSP90AA1_GF, by = "geneName")

# correlationAnalyzeR
HSP90AA1_me <- correlationAnalyzeR::getCorrelationData(Species = "hsapiens", Sample_Type = "normal",
                                                       Tissue = "all", geneList = "HSP90AA1")
HSP90AA1_me$geneName <- rownames(HSP90AA1_me)
colnames(HSP90AA1_me)[1] <- "correlationAnalyzeR"
corCompare <- merge(x = HSP90AA1_me, y = corCompare, by = "geneName")

# Plot
chart.Correlation(corCompare[,c(-1)],   method = "spearman")


## Compare HSP90AA1 top 100 to biogrid
biogrid <- read_tsv("misc/BIOGRID-GENE-109552-4.1.190.tab3.txt")

top50 <- HSP90AA1_me %>%
  top_n(wt = HSP90AA1, n = 500)

interacts <- biogrid %>%
  select(`Official Symbol Interactor B`) %>%
  distinct()
vl <- list(
  "Co-expression" = top50$geneName,
  "Interaction" = interacts$`Official Symbol Interactor B`
)
ol <- calculate.overlap(vl)
calculate.overlap.and.pvalue(list1 = vl$`Co-expression`, list2 = vl$Interaction,
                             total.size = unique(length(BRCA1_me$geneName)),
                             lower.tail = FALSE)
vd <- venn.diagram(vl, filename = NULL, fill = c("forestgreen", "firebrick"),
                   margin = .05)

dev.off()
grid.draw(vd)











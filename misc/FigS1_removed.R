# Refseq summary
library(biomaRt)
library(openxlsx)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
martRes <- getBM(attributes = c("external_gene_name",
                                "refseq_mrna"),
                 mart = ensembl)
geneInfo <- read.xlsx("../summaryRefseq.xlsx")
colnames(geneInfo) <- c("refseq_mrna", "something", "Function")
geneInfoFinal <- merge(x = martRes, y = geneInfo, by = "refseq_mrna")


geneInfoFinal <- geneInfoFinal[,c(2, 4)]
geneInfoFinal <- unique(geneInfoFinal)
geneWithRef <- geneInfoFinal$external_gene_name[which(geneInfoFinal$Function != "")]


geneInfo$refSeq <- FALSE
geneInfo$refSeq[grep(geneInfo$Function, pattern = "RefSeq")] <- TRUE
toPlot <- as.data.frame(table(geneInfo$refSeq))

# Fig 1A
pdf("fig1Pie.pdf", height = 6, width = 8)
pie(toPlot$Freq, col = c("azure", "coral1"),
    labels = c("No annotation (59.6%)", "RefSeq (40.4%)"),
    init.angle = 90)
dev.off()



# Get the survival tables from GEPIA and determine # with refseq acc

names <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
           "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
           "LUAD", "LUSC", "MESO", "OV")
files <- list.files("data/survivalTables/")
orderNum <- gsub(files, pattern = ".*\\((.+)\\).*", replacement = "\\1")
orderNum[20] <- 0
orderNum <- as.numeric(orderNum)
files <- files[order(orderNum)]
names(files) <- names

nIn <- 100
inVec <- c()
outVec <- c()
for (i in 1:length(files)) {
  file <- files[i]
  group <- names(file)
  print(group)
  df <- read.table(file.path("data/survivalTables", file), stringsAsFactors = FALSE,
                   sep = "\t", header = TRUE)
  df <- unique(df[c(1:nIn),c(1, 3)])
  outVec <- c(outVec, length(df$Gene.Symbol[which(! df$Gene.Symbol %in% geneWithRef)]))
  inVec <- c(inVec, length(df$Gene.Symbol[which(df$Gene.Symbol %in% geneWithRef)]))
}

library(ggpubr)
plotDF <- data.frame(group = rep(names(files), 2), Annotation = c(rep("RefSeq", length(files)),
                                                                  rep("No Annotation", length(files))),
                     number = c(inVec, outVec), stringsAsFactors = FALSE)
plotDF <- plotDF[order(plotDF$Annotation, plotDF$number),]

# Fig 1B
ggbarplot(data = plotDF, x = "group", y = "number", legend = "right",
          font.tickslab = c(15), font.legend = c(15),
          palette = c("azure", "coral1"),
          fill = "Annotation") +
  rotate() + rremove("xlab") +
  rremove("ylab") +
  rremove("legend.title")






# # Make public user
# library(DBI)
# credFile <- "misc/credFile.txt"
# credentials <- suppressWarnings(read.delim(credFile, sep = ";",
#                                            header = FALSE, stringsAsFactors = F))
# uName <- credentials$V1
# pWord <- credentials$V2
# conn <- dbConnect(drv = RMySQL::MySQL(), user = uName,
#                   port = 3306, dbname="correlation_analyzer",
#                   password=pWord,
#                   host="m2600az-db01p.mysql.database.azure.com")
# sql <- "CREATE USER 'public-rds-user'@'%' IDENTIFIED BY 'public-user-password';"
# dbExecute(conn, sql)
# dbDisconnect(conn)
# sql <- "GRANT SELECT ON * TO 'public-rds-user'@'%'"
# dbExecute(conn, sql)















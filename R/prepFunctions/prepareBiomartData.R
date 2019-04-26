# library(biomaRt)
# dir.create("data/biomaRtData")
# # Build human mart dataset
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# aa <- listAttributes(ensembl)
# attributes <- c("external_gene_name",
#                 "description")
# humanGenes <- getBM(attributes = attributes, mart = ensembl)
# save(humanGenes, file = "data/biomaRtData/humanGenes.RData")
# # Build mouse mart dataset
# ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
# aa <- listAttributes(ensembl)
# attributes <- c("external_gene_name",
#                 "description")
# mouseGenes <- getBM(attributes = attributes, mart = ensembl)
# save(mouseGenes, file = "data/biomaRtData/mouseGenes.RData")

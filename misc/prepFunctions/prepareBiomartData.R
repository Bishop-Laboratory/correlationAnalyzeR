library(biomaRt)
library(usethis)
dir.create("data/biomaRtData")
# Build human mart dataset
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
aa <- listAttributes(ensembl)
attributes <- c("external_gene_name",
                "description")
humanGenes <- getBM(attributes = attributes, mart = ensembl)
use_data(humanGenes)
# Build mouse mart dataset
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
aa <- listAttributes(ensembl)
attributes <- c("external_gene_name",
                "description")
mouseGenes <- getBM(attributes = attributes, mart = ensembl)
use_data(mouseGenes)

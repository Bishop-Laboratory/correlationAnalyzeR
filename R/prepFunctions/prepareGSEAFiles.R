# library(msigdbr)
# m_df <- msigdbr("Mus musculus")
# TERM2GENE <- m_df[,c(1, 8)]
# save(TERM2GENE, file = "data/GSEA_Data/mmusculus/TERM2GENE_MM_COMPLEX.RData")
# m_df <- msigdbr("Homo sapiens")
# TERM2GENE <- m_df[,c(1, 5)]
# save(TERM2GENE, file = "data/GSEA_Data/hsapiens/TERM2GENE_HS_COMPLEX.RData")
# m_df <- msigdbr("Mus musculus", category = c("H", "C1", "C2", "C5", "C6"))
# TERM2GENE <- m_df[,c(1, 8)]
# save(TERM2GENE, file = "data/GSEA_Data/mmusculus/TERM2GENE_MM_SIMPLE.RData")
# m_df <- msigdbr("Homo sapiens", category = c("H", "C1", "C2", "C5", "C6"))
# TERM2GENE <- m_df[,c(1, 5)]
# save(TERM2GENE, file = "data/GSEA_Data/hsapiens/TERM2GENE_HS_SIMPLE.RData")

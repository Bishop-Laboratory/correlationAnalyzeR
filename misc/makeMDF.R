# Make new TERM2GENE
mdf1 <- msigdbr::msigdbr()
mdf2 <- msigdbr::msigdbr(species = "Mus musculus")
library(dplyr)
mdf1 <- mdf1 %>% select(gs_name, gs_id, gs_cat, gs_subcat, human_gene_symbol)
mdf2 <- mdf2 %>% select(gs_name, gs_id, gs_cat, gs_subcat, human_gene_symbol, gene_symbol)
mdf <- full_join(mdf1, mdf2, by = colnames(mdf1))
which(is.na(mdf$gene_symbol))
mdf <- mdf %>% select(-gs_id)
colnames(mdf)[c(4:5)] <- c("human_gene_symbol", "mouse_gene_symbol")
# TERM2GENEHuman <- mdf %>%
#   filter(gs_cat %in% c("H")) %>%
#   select(gs_name, human_gene_symbol) %>%
#   distinct()

table(mdf$gs_subcat)
mdf$gs_subcat <- gsub(mdf$gs_subcat, pattern = "", replacement = "None")
mdf$gs_cat <- paste0(mdf$gs_cat, "_", mdf$gs_subcat)
mdf <- mdf %>% select(-gs_subcat)
MDF <- mdf %>% distinct()
# save(MDF, file = "misc/MDF.rda", compression_level = 9)



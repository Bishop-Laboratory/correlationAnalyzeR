library(msigdbr)
library(usethis)
m_df <- msigdbr("Mus musculus")
mmusculus_complex_TERM2GENE <- m_df[,c(1, 8)]
mmusculus_complex_TERM2GENE <- unique(mmusculus_complex_TERM2GENE)
use_data(mmusculus_complex_TERM2GENE, overwrite = T)

m_df <- msigdbr("Homo sapiens")
hsapiens_complex_TERM2GENE <- m_df[,c(1, 5)]
hsapiens_complex_TERM2GENE <- unique(hsapiens_complex_TERM2GENE)
use_data(hsapiens_complex_TERM2GENE, overwrite = T)

categories <- list("H" = c("all"),
                   "C2" = c("all"),
                   "C5" = c("BP"),
                   "C6" = c("all"))
msigListHuman <- list()
msigListMouse <- list()
for (i in 1:length(categories)) {
  category <- names(categories)[i]
  subcategory <- categories[[i]]
  print(category)
  print(subcategory)
  if (category != "C5") {
    msigListHuman[[i]] <- msigdbr("Homo sapiens", category = category)
    msigListMouse[[i]] <- msigdbr("Mus musculus", category = category)
  } else {
    msigListHuman[[i]] <- msigdbr("Homo sapiens",
                                  category = category,
                                  subcategory = subcategory)
    msigListMouse[[i]] <- msigdbr("Mus musculus",
                                  category = category,
                                  subcategory = subcategory)
  }
}
m_dfFinal <- data.table::rbindlist(msigListHuman)
m_dfFinal <- unique(m_dfFinal)
hsapiens_simple_TERM2GENE <- m_dfFinal[,c(1, 5)]
use_data(hsapiens_simple_TERM2GENE, overwrite = T)

m_dfFinal <- data.table::rbindlist(msigListMouse)
m_dfFinal <- unique(m_dfFinal)
mmusculus_simple_TERM2GENE <- m_dfFinal[,c(1, 8)]
mmusculus_simple_TERM2GENE <- unique(mmusculus_simple_TERM2GENE)
use_data(mmusculus_simple_TERM2GENE, overwrite = T)



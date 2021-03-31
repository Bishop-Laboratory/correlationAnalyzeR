### Additional Figure S2 based on requests from reviewers ###

library(tidyverse)

## Part #1: Add in correlation of two genes with different normalization ##
set.seed(42)

gene1 <- "BRCA1"
gene2 <- "NQO1"
Sample_Type <- "cancer"
Tissue <- "bone"

vstB <- getTissueVST(genesOfInterest = c(gene1, gene2),
                     Sample_Type = Sample_Type,
                     Tissues = Tissue) %>%
  bind_rows() %>%
  inner_join(correlationAnalyzeR::human_coldata, by = "samples")
B1_corr <- getCorrelationData(geneList = gene1,
                              Sample_Type = Sample_Type, Tissue = Tissue)
pDF <- apply(B1_corr, MARGIN = 1:2, n = length(B1_corr[,1]), FUN = function(x, n) {
  stats::dt(abs(x)/sqrt((1-x^2)/(n-2)), df = 2)
})
padj <- p.adjust(pDF[,1], method = "BH")
B1B2r <- B1_corr[gene2,]
B1B2p <- pDF[gene2,]
B1B2padj <- padj[gene2]
to_sample <- ifelse(length(vstB$samples) > 5000, 5000, length(vstB$samples))
plt <- vstB %>%
  filter(samples %in% sample(samples, to_sample)) %>%
  ggplot(aes(x = !!sym(gene1), y = !!sym(gene2), color = disease,
             text = paste0("samples", "\n", "disease", "\n", "Tissue"))) +
  geom_point(alpha = .5) +
  labs(title = paste0(gene1, " vs ", gene2),
       subtitle = paste0("Pearson's R = ", round(B1B2r, 3),
                         " (padj = ", signif(B1B2padj, 3), ")")) +
  ggplot2::theme_bw(base_size = 16) +
  scale_color_manual(name = "Disease", values = c("firebrick", "forestgreen")) +
  xlab(paste0(gene1, " Expression (VST)")) +
  ylab(paste0(gene2, " Expression (VST)"))
plt




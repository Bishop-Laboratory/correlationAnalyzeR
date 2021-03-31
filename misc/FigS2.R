### Additional Figure S2 based on requests from reviewers ###

library(tidyverse)

# ## Part #1: Add in correlation of two genes with different normalization ##
# set.seed(42)
#
# gene1 <- "BRCA1"
# gene2 <- "NQO1"
# Sample_Type <- "cancer"
# Tissue <- "bone"
#
# vstB <- getTissueVST(genesOfInterest = c(gene1, gene2),
#                      Sample_Type = Sample_Type,
#                      Tissues = Tissue) %>%
#   bind_rows() %>%
#   inner_join(correlationAnalyzeR::human_coldata, by = "samples")
# B1_corr <- getCorrelationData(geneList = gene1,
#                               Sample_Type = Sample_Type, Tissue = Tissue)
# pDF <- apply(B1_corr, MARGIN = 1:2, n = length(B1_corr[,1]), FUN = function(x, n) {
#   stats::dt(abs(x)/sqrt((1-x^2)/(n-2)), df = 2)
# })
# padj <- p.adjust(pDF[,1], method = "BH")
# B1B2r <- B1_corr[gene2,]
# B1B2p <- pDF[gene2,]
# B1B2padj <- padj[gene2]
# to_sample <- ifelse(length(vstB$samples) > 5000, 5000, length(vstB$samples))
# plt <- vstB %>%
#   filter(samples %in% sample(samples, to_sample)) %>%
#   ggplot(aes(x = !!sym(gene1), y = !!sym(gene2), color = disease,
#              text = paste0("samples", "\n", "disease", "\n", "Tissue"))) +
#   geom_point(alpha = .5) +
#   labs(title = paste0(gene1, " vs ", gene2),
#        subtitle = paste0("Pearson's R = ", round(B1B2r, 3),
#                          " (padj = ", signif(B1B2padj, 3), ")")) +
#   ggplot2::theme_bw(base_size = 16) +
#   scale_color_manual(name = "Disease", values = c("firebrick", "forestgreen")) +
#   xlab(paste0(gene1, " Expression (VST)")) +
#   ylab(paste0(gene2, " Expression (VST)"))
# plt

## Part 2: Compare different correlation types ##
# genes_to_test <- sample(correlationAnalyzeR::humanGenesVST, 1000)
# my_corr <- correlationAnalyzeR::getTissueVST(genesOfInterest = genes_to_test,
#                                              Sample_Type = "all")
# my_corr2 <- bind_rows(my_corr) %>%
#   inner_join(correlationAnalyzeR::human_coldata, by = "samples") %>%
#   distinct(samples, .keep_all = TRUE)
# save(genes_to_test, my_corr2, file = "misc/my_corr2.rda")
load("misc/my_corr2.rda")

res <- lapply(seq(genes_to_test), FUN = function(i) {
  print(i)
  gene_now <- genes_to_test[i]
  resS <- lapply(seq(1000), function(j) {
    my_corr_tmp <- my_corr2  %>%
      filter(samples %in% sample(my_corr2$samples, 25))
    df <- tryCatch({
      ddb1 <- shapiro.test(my_corr_tmp[,gene_now])
      pval <- ddb1$p.value
      w <- ddb1$statistic
      data.frame(pval, w)
      },
                     error=function(cond){
                       return("STOP")
                     })
    return(df)
  })

  if (any(unlist(resS) == "STOP")) {
    return(NULL)
  }
  dd <- bind_rows(resS)
  dd$padj <- p.adjust(dd$pval)
  dd$gene <- gene_now
  return(dd)
})

boot_shapiro <- bind_rows(res)
readr::write_csv(boot_shapiro, "misc/boot_shapiro.csv")
# boot_shapiro <- read_csv("misc/boot_shapiro.csv")
#
# hist(-log10(boot_shapiro$padj), breaks = 100)
#
# quantile(-log10(boot_shapiro$padj))
#
#
# # boot_shapiro %>%
# #   group_by(gene) %>%
# #   summarise(avg = median(padj)) -> dd
# # hist(dd$avg)
#
# #
# #
#
# x <- my_corr2$AC012354.8
# hist(x)
# qqnorm(x); qqline(x)
#
#

res <- analyzeGenePairs(
  genesOfInterest = c("BRCA1", "BRCA2"), runGSEA = F
)

res$compared$VST_corrPlot$corrPlot_all +
  labs(subtitle = NULL)

uuu <- res$compared$VST_corrPlot$corrPlot_VST_data
cor(x = uuu$BRCA1, y = uuu$BRCA2, method = "kendall")
cor(x = uuu$BRCA1, y = uuu$BRCA2, method = "pearson")
cor(x = uuu$BRCA1, y = uuu$BRCA2, method = "spearman")



correlationAnalyzeR::getTissueVST()
randGenes <- sample(correlationAnalyzeR::humanGenesVST )


vstNow <- correlationAnalyzeR::getTissueVST(genesOfInterest = sample(unique(t2g$gene_symbol),
                                                                     500))
vstNow2 <- bind_rows(vstNow) %>%
  inner_join(correlationAnalyzeR::human_coldata, by = "samples") %>%
  distinct(samples, .keep_all = TRUE)
save(vstNow2, file = "misc/vstNow2.rda")
load("misc/vstNow2.rda")
t2g <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "Hallmark")
genes_to_testNow <- colnames(vstNow2)
genes_to_testNow <- genes_to_testNow[which(genes_to_testNow %in% t2g$gene_symbol)]

res <- lapply(seq(5000), function(x) {
  gs <- sample(unique(t2g$gs_name), 1)
  genes <- sample(t2g$gene_symbol[t2g$gs_name == gs &
                                    t2g$gene_symbol %in% genes_to_testNow], 2)
  # cor(x = vstNow2[,genes[1]], y = vstNow2[,genes[2]], method = "kendall")
  pr <- cor(x = vstNow2[,genes[1]], y = vstNow2[,genes[2]], method = "pearson")
  sr <- cor(x = vstNow2[,genes[1]], y = vstNow2[,genes[2]], method = "spearman")
  return(data.frame(pr, sr, gs, gene1=genes[1], gene2=genes[2]))
})

dd <- bind_rows(res) %>%
  mutate(diff = pr - sr) %>%
  distinct(diff, .keep_all = TRUE)
save(dd, file = "misc/simulated_spearman_vs_pearson.rda")

library(correlationAnalyzeR)
ex1 <- analyzeGenePairs(
  genesOfInterest = c("MYL3", "TCAP"), runGSEA = F
)

g11 <- ex1$compared$VST_corrPlot$corrPlot_tissue +
  labs(subtitle = NULL, title = "Top Pearson-Specific Correlation")

ex2 <- analyzeGenePairs(
  genesOfInterest = c("TFF2", "IL12B"), runGSEA = F
)
g22 <- ex2$compared$VST_corrPlot$corrPlot_tissue +
  labs(subtitle = NULL, title = "Top Spearman-Specific Correlation")

ggg2 <- ggarrange(g11, g22, align = "hv", common.legend = TRUE, legend = "bottom") +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare_example_VST.pdf") +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare_example_VST.png",
         height = 7, width = 12)

g1 <- dd %>%
  select(Pearson = pr, Spearman = sr) %>%
  ggplot(aes(x = Pearson, y = Spearman)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1,
              size = 1.5, color = "firebrick") +
  scale_x_continuous(limits = c(-.75, 1)) +
  scale_y_continuous(limits = c(-.75, 1)) +
  theme_bw(base_size = 16) +
  labs(title = "Correlation Methods in Hallmark Collection",
       subtitle = "2,832 Simulations via bootstrapping") +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare_scatter.pdf") +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare__scatter.png")

g2 <- dd %>%
  pivot_longer(cols = c(pr, sr)) %>%
  mutate(name = case_when(name == "pr" ~ "Pearson",
                          TRUE ~ "Spearman")) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(method = "t.test",
                             comparisons = list(c("Pearson", "Spearman")),
                             method.args = list("alternative" = "greater"),
                             label = "p.signif") +
  theme_bw(base_size = 16) +
  ggpubr::rremove("legend") +
  labs(title = "Correlation Methods in Hallmark Collection",
       subtitle = "2,832 Simulations via bootstrapping") +
  ylab("Correlation Coefficient") +
  xlab(NULL) +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare.pdf") +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare.png")

g2 <- dd %>%
  pivot_longer(cols = c(pr, sr)) %>%
  mutate(name = case_when(name == "pr" ~ "Pearson",
                          TRUE ~ "Spearman")) %>%
  ggplot(aes(x = name, y = value, fill = name)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(method = "t.test",
                             comparisons = list(c("Pearson", "Spearman")),
                             method.args = list("alternative" = "greater"),
                             label = "p.signif") +
  theme_bw(base_size = 16) +
  ggpubr::rremove("legend") +
  labs(title = "Correlation Methods in Hallmark Collection",
       subtitle = "2,832 Simulations via bootstrapping") +
  ylab("Correlation Coefficient") +
  xlab(NULL) +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare.pdf") +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare.png")


library(cowplot)
library(ggpubr)

ggg <- ggarrange(g1, g2, align = "hv")
annotate_figure(ggg, top = text_grob("Correlation comparison using 'Hallmark' collection\n(2,832 Simulations)",
                                face = "bold", size = 18)) +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare_comb.pdf",
         height = 5, width = 10) +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare_comb.png",
         height = 5, width = 10)


gggg <- ggarrange(ggg,ggg2,  align = "hv", nrow = 2)
gggg




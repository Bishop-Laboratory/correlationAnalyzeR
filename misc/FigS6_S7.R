### Additional Figure S5 and S6 based on requests from reviewers ###
library(tidyverse)
library(ggpubr)
library(correlationAnalyzeR)

# Hallmark geneset
t2g <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "Hallmark")

# Get VST for 500 genes in the Hallmark gene set and wrangle
vstNow <- correlationAnalyzeR::getTissueVST(genesOfInterest = sample(unique(t2g$gene_symbol),
                                                                     500))
vstNow2 <- bind_rows(vstNow) %>%
  inner_join(correlationAnalyzeR::human_coldata, by = "samples") %>%
  distinct(samples, .keep_all = TRUE)
save(vstNow2, file = "misc/vstNow2.rda")

## Boostrap Pearson and Spearman ##
load("misc/vstNow2.rda")
genes_to_testNow <- colnames(vstNow2)
genes_to_testNow <- genes_to_testNow[which(genes_to_testNow %in% t2g$gene_symbol)]

# Bootstrap 5000 simulations
res <- lapply(seq(5000), function(x) {
  gs <- sample(unique(t2g$gs_name), 1)
  genes <- sample(t2g$gene_symbol[t2g$gs_name == gs &
                                    t2g$gene_symbol %in% genes_to_testNow], 2)
  # cor(x = vstNow2[,genes[1]], y = vstNow2[,genes[2]], method = "kendall")
  pr <- cor(x = vstNow2[,genes[1]], y = vstNow2[,genes[2]], method = "pearson")
  sr <- cor(x = vstNow2[,genes[1]], y = vstNow2[,genes[2]], method = "spearman")
  return(data.frame(pr, sr, gs, gene1=genes[1], gene2=genes[2]))
})

# Find difference between Pearson and Spearman AND remove duplicate gene pairs
dd <- bind_rows(res) %>%
  mutate(diff = pr - sr) %>%
  distinct(diff, .keep_all = TRUE)
save(dd, file = "misc/simulated_spearman_vs_pearson.rda")
load("misc/simulated_spearman_vs_pearson.rda")

## Generate Figures ##

# Fig S5
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
       subtitle = "2,832 Simulations via bootstrapping")
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
  xlab(NULL)

ggg <- ggarrange(g1, g2, align = "hv")
annotate_figure(ggg, top = text_grob("Correlation comparison using 'Hallmark' collection\n(2,832 Simulations)",
                                     face = "bold", size = 18)) +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare_comb.pdf",
         height = 5, width = 10) +
  ggsave(filename = "../Manuscript/Assets/Correlation_Methods_Compare_comb.png",
         height = 5, width = 10)


# Fig S6
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






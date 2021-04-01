library(tidyverse)
library(ggpubr)
library(correlationAnalyzeR)
res <- correlationAnalyzeR::analyzeGenePairs(
  genesOfInterest = c("BRCA1", "BRCA2"),
  Sample_Type = "all",
  runGSEA = FALSE
)
g1 <- res$compared$VST_corrPlot$corrPlot_disease +
  labs(title = "All samples")

resNorm <- correlationAnalyzeR::analyzeGenePairs(
  genesOfInterest = c("BRCA1", "BRCA2"),
  Sample_Type = "normal",
  runGSEA = FALSE
)
g2 <- resNorm$compared$VST_corrPlot$corrPlot_disease +
  labs(title = "Normal samples")

resCancer <- correlationAnalyzeR::analyzeGenePairs(
  genesOfInterest = c("BRCA1", "BRCA2"),
  Sample_Type = "cancer",
  runGSEA = FALSE
)
g3 <- resCancer$compared$VST_corrPlot$corrPlot_disease +
  labs(title = "Cancer samples")


ggarrange(g1, g2, g3, nrow = 1, align = "hv") +
  ggsave(filename = "../Manuscript/FinalAssets/FigureS2_raw.png",
         height = 5, width = 18)

cd <- correlationAnalyzeR::human_coldata
write_csv(cd, file = "misc/colData.csv")

### IL1B; IL1RN

resRev <- analyzeGenePairs(genesOfInterest = c("IL1B", "IL1RN"),
                           runGSEA = F)
resRev$compared$VST_corrPlot$corrPlot_tissue



g1 <- res$compared$VST_corrPlot$corrPlot_disease +
  labs(title = "All samples")

resNorm <- correlationAnalyzeR::analyzeGenePairs(
  genesOfInterest = c("BRCA1", "NQO1"),
  Sample_Type = "normal",
  runGSEA = FALSE
)
g2 <- resNorm$compared$VST_corrPlot$corrPlot_disease +
  labs(title = "Normal samples")

resCancer <- correlationAnalyzeR::analyzeGenePairs(
  genesOfInterest = c("BRCA1", "NQO1"),
  Sample_Type = "cancer",
  runGSEA = FALSE
)
g3 <- resCancer$compared$VST_corrPlot$corrPlot_disease +
  labs(title = "Cancer samples")


geneOne <- "BRCA1"
geneTwo <- "NQO1"
titleStr <- "Normal samples"
Rval <- res$compared$VST_corrPlot$Rval
Padj <- res$compared$VST_corrPlot$Padj
resNorm$compared$VST_corrPlot$corrPlot_VST_data %>%
  mutate(Condition = ifelse(Group %in% c(
                                         # "Mammary - Normal",
                                         # "Respiratory - Normal",
                                         # "Thyroid - Normal",
                                         # "Esophagus - Normal",
                                         "Pancreas - Normal",
                                         "Prostate - Normal",
                                         "Skin - Normal",
                                         # "Mammary - Normal",
                                         # "Female Reproductive - Normal",
                                         # "Respiratory - Normal",
                                         # "Kidney - Normal",
                                         "Liver - Normal"
                                         # "Cartilage - Normal",
                                         # "Adispoe - Normal",
                                         # "Muscle - Normal",
                                         # "Brain - Normal",
                                         # "Retina - Normal",
                                         # "Endothelial - Normal"
                                         ), "Correlated",
                            ifelse(Group %in% c(
                              "Prenatal - Normal",
                              "Male Reproductive - Normal",
                              "Stem Like - Normal",
                              "Bone - Normal"
                            ), "Anticorrelated", "Other"
                          ))) %>%
  mutate(Condition = factor(Condition, levels = c(
    "Correlated",
    "Anticorrelated",
    "Other"
  ))) %>%
  arrange(desc(Condition)) %>%
  ggplot2::ggplot(ggplot2::aes_string(x = geneOne,
                                      y = geneTwo,
                                      group="Group",
                                      text = "samples",
                                      color = "Condition")) +
  ggplot2::geom_point(alpha = .8) +
  ggplot2::labs(title = titleStr) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::xlab(paste0(geneOne, " Expression (VST)")) +
  ggplot2::ylab(paste0(geneTwo, " Expression (VST)"))






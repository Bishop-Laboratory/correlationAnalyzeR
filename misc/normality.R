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

# library(correlationAnalyzeR)
# res <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = "SRY",
#                                                whichCompareGroups = "normal",
#                                          crossCompareMode = TRUE)


#
# res <- correlationAnalyzeR::analyzeGenePairs(genesOfInterest = c("BRCA1", "BRCA1"), crossCompareMode = TRUE)
#
# res <- correlationAnalyzeR::geneVsGeneListAnalyze(pairedGenesList = list("BRCA1" = c("ATM", "EZH2", "STAG2")))
#
# res <- correlationAnalyzeR::analyzeGenesetTopology(genesOfInterest = c("ATM", "SETX", "CDKN1A", "CDKN2A", "FUS"))
#
#
#
# correlationAnalyzeR::getAvailableGenes()
#
# genesOfInterest <- c("CDK12", "AURKB", "SFPQ", "NFKB1", "BRCC3", "BRCA2", "PARP1",
#                      "DHX9", "SON", "AURKA", "SETX", "BRCA1", "ATMIN")
# res <- correlationAnalyzeR::analyzeGenesetTopology(genesOfInterest = genesOfInterest,
#                                             Sample_Type = "cancer", returnDataOnly = TRUE,
#                                             Tissue = "brain")
#
# res <- correlationAnalyzeR::getTissueVST(genesOfInterest = c("BRCA1", "ATM"),
#                     Tissues = c("brain", "respiratory"),
#                     Sample_Type = "all",
#                     useBlackList = TRUE)
#
# genesOfInterest <- c("ATM", "SLC7A11")
# res <- correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#                                       GSEA_Type = "simple", returnDataOnly = TRUE,
#                                       Sample_Type = c("normal", "normal"),
#                                       Tissue = c("brain", "brain"))
# genesOfInterest <- c("BRCA1", "BRCA1")
# res <- correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#                                       GSEA_Type = "simple", returnDataOnly = TRUE,
#                                       Sample_Type = c("normal", "cancer"),
#                                       Tissue = c("respiratory", "respiratory"))
# genesOfInterest <- c("NFKB1", "SOX10")
# res <- correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#                                       returnDataOnly = TRUE,
#                                       crossCompareMode = TRUE)
# genesOfInterest <- c("SOX10", "SOX10")
# res <- correlationAnalyzeR::analyzeGenePairs(genesOfInterest = genesOfInterest,
#                                              returnDataOnly = TRUE,
#                                              crossCompareMode = TRUE)
#
#
#
# pairedGenesList <- list("TP53" = c("BRCA1", "CDK12", "PARP1"),
#                         "SON" = c("AURKB", "SFPQ", "DHX9"))
#
# res <- correlationAnalyzeR::geneVsGeneListAnalyze(pairedGenesList = pairedGenesList,
#                                            returnDataOnly = TRUE,
#                                            Sample_Type = "normal",
#                                            Tissue = "brain")
#
#
# genesOfInterest <- c("ATM", "SLC7A11")
# res <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = genesOfInterest,
#                                         returnDataOnly = TRUE,
#                                         GSEA_Type = "simple",
#                                         Sample_Type = c("normal", "cancer"),
#                                         Tissue = c("respiratory", "pancreas"))
#
# genesOfInterest <- c("BRCA1")
# res <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = genesOfInterest,
#                                         GSEA_Type = "simple", returnDataOnly = TRUE,
#                                         crossCompareMode = TRUE,
#                                         whichCompareGroups = "normal")
#
#
# res <- correlationAnalyzeR::getCorrelationData(Sample_Type = "normal",
#                                         Tissue = "kidney",
#                                         geneList = c("ATM", "BRCA1"))
#
#
#
#
# res <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "simple")
# res <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = c("Hallmark", "KEGG"))
#
#
#
# correlationAnalyzeR::getTissueTypes()
#
#
#
# corrDF <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = c("BRCA1"),
#                                                   returnDataOnly = TRUE, runGSEA = FALSE, Sample_Type = "normal")
# ranks <- corrDF$correlations[,1]
# names(ranks) <- rownames(corrDF$correlations)
# TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "simple",
#                                                Species = "hsapiens")
# res <- correlationAnalyzeR::myGSEA(ranks = ranks,
#                             TERM2GENE = TERM2GENE,
#                             plotFile = "GSEA_out", outDir = getwd(),
#                             topPlots = FALSE, returnDataOnly=TRUE, Condition = "GSEA Results")
#
#
#
#
#
#

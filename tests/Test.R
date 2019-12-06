# devtools::install_github("millerh1/correlationAnalyzeR")
# library(correlationAnalyzeR)
#
# genesOfInterest <- "ATM"
# corr <- getCorrelationData(Species = "hsapiens",
#                            Sample_Type = "normal",
#                            Tissue = "blood",
#                            geneList = genesOfInterest)
#
#
# genesOfInterest <- c("ATM", "AKT1", "HIF1A", "PRKAA2", "SLC3A2", "MTOR",
#                      "SLC7A11", "SLC7A5", "G6PD", "HSPB1")
# # Run the function with desired parameters
# res <- analyzeSingleGenes(genesOfInterest = genesOfInterest,
#                           outputPrefix = "tests/singleGeneTesting",
#                           Species = "hsapiens", Tissue = "pancreas",
#                           topPlots = FALSE,
#                           Sample_Type = "cancer")
# corrDF <- res$correlations
# res <- correlationAnalyzeR::analyzeGenePairs(genesOfInterest = c("ATM", "SLC3A2"),
#                                              Tissue = c("pancreas", "pancreas"))
#
#
#
#
# # Check the paired genes
# pairedGeneList <- list("REST" = c("NCOR1", "HSPB1", "BRCA1", "ATM"))
# pairedGenesAnalyzeR(pairedGenesList = pairedGeneList,
#                     Sample_Type = "Normal_Tissues",
#                     outputPrefix = "tests/pairedOut")
#
#
#
# # Check topology analysis methods
# geneList <- c("TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_UP")
# genesOfInterest <- c("ATM", "AKT1", "HIF1A", "PRKAA2", "SLC3A2", "MTOR",
#                      "SLC7A11", "SLC7A5", "G6PD", "HSPB1")
# genesOfInterest <- c("Atm", "Atmin", "Brca1", "Ewsr1", "Slc7a11", "Slc3a2")
# result <- analyzeGenesetTopology(genesOfInterest = genesOfInterest,
#                                  Species = "hsapiens", Sample_Type = "Tumor_Tissues",
#                                  outputPrefix = "tests/topology_test5")
#
#
# genesOfInterest <- "F2"
# analyzeSingleGenes(genesOfInterest = genesOfInterest,
#                    Sample_Type = "Normal_Tissues", outputPrefix = "tests/newTest")
#
#
#
# analyzeSingleGenes(genesOfInterest = genesOfInterest,
#                    outputPrefix = "tests/singleGeneTesting",
#                    Species = "hsapiens",
#                    Sample_Type = "Normal_Tissues")
#
# # load('data/hsapiens_corrSmall_geneNames.rda')
# # mmusculus_corrSmall_geneNames <- geneNames
# # usethis::use_data(mmusculus_corrSmall_geneNames)
#
# #  # Test whether pairedGenesAnalyzeR is working
# # pairedGenesList <- list("SF3B1" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "SF3B2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "SF3B3" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "U2AF2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "U2AF1" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "SRSF2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "AQR" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "XAB2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "WNT5A" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "BCL2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "SOX6" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "FAS" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "EWSR1" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
# #                         "FLI1" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"))
# #
# # outDir <- paste0(getwd(), "/Results/THOC6-TCNER")
# #
# # PathCardsTCNER <- c("DDB1	RPS27A	UBA52	UBB	UBC	CUL4B	CUL4A	RBX1	CCNH	CDK7	ERCC2	ERCC3
# #   GTF2H1	GTF2H2	GTF2H3	GTF2H4	MNAT1	GTF2H5	RPA1	RPA2	RPA3	XPA	ERCC8	EP300
# #   ERCC1	ERCC4	ERCC5	ERCC6	GPS1	HMGN1	PCNA	POLD1	POLD2	POLE	POLE2	POLR2A
# #   POLR2B	POLR2C	POLR2D	POLR2E	POLR2F	POLR2G	POLR2H	POLR2I	POLR2J	POLR2K	POLR2L	RFC1
# #   RFC2	RFC3	RFC4	RFC5	TCEA1	USP7	COPS3	COPS2	AQR	PPIE	POLD3	COPS8
# #   COPS6	COPS5	PRPF19	COPS7A	COPS4	POLK	POLE3	POLE4	XAB2	ISY1	UVSSA	POLD4
# #   COPS7B	ZNF830	MIR1281	PARP1	CETN2	DDB2	LIG1	LIG3	RAD23A	RAD23B	XPC	XRCC1
# #   PARP2	ACTB	ACTL6A	NFRKB	SUMO3	SUMO2	UBE2I	UBE2N	UBE2V2	SUMO1	YY1	PIAS1
# #   RUVBL1	CHD1L	PIAS3	MCRS1	TFPT	INO80	RNF111	INO80D	ACTR5	INO80B	USP45	ACTR8
# #   INO80C	INO80E	MIR6764	ELL")
# # PathCardsTCNER <- unlist(strsplit(PathCardsTCNER, "\n"))
# # PathCardsTCNER <- unlist(strsplit(PathCardsTCNER, "\t"))
# # DDB1 <- c("DDB1", "DCAF1", "DCAF2", "DCAF3", "DCAF4", "DCAF4L1", "DCAF5", "DCAF6",
# #           "DCAF7", "DCAF8", "DCAF9", "DCAF10", "DCAF11", "DCAF12", "DCAF13", "DCAF14",
# #           "DCAF15", "DCAF16", "DDA1", "DDB1", "CRBN")
# # PathCardsTCNERHigh <- PathCardsTCNER[c(1:20)]
# # pairedGenesList <- list("THOC6" = DDB1)
# # pairedGenesList <- list("THOC6" = PathCardsTCNER)
# # PathCardsTCNER <- trimws(PathCardsTCNER)
# #
# # load("Data/geneSets/nrf2.gene.sets.bishop.RData")
# # source("Code/analyzeGenePairs.R")
# #
# # pairedGenesList <- list("ATM" = unique(BISHOP_NRF2_PATHWAY))
# # outDir <- paste0(getwd(), "/Results/ATM_ROS")
# # Res <- pairedGenesAnalyzeR(pairedGenesList = pairedGenesList,
# #                            outDir = outDir, onlyTop = FALSE, topCutoff = .5,
# #                            sigTest = TRUE)
# #
# # ATM_NRF2_MARKERS <- Res$ATM
# # ATM_NRF2_MARKERS <- ATM_NRF2_MARKERS[order(ATM_NRF2_MARKERS)]
# #
# # # Met signature
# # Met1 <- read.table(file = "Data/geneSets/artega_sig.txt")
# # Met2 <- read.table(file = "Data/geneSets/mammaprint_sig_new.txt")
# # Met3 <- read.table(file = "Data/geneSets/werb_49_metastasis_sig.txt")
# # Met <- unique(c(as.character(Met1$V1), as.character(Met2$V1), as.character(Met3$V1)))
# # pairedGenesList <- list("ATM" = Met,
# #                         "SLC3A2" = Met)
# # outDir <- paste0(getwd(), "/Results/ATM_ROS_MET")
# #
# #
# # Res <- pairedGenesAnalyzeR(pairedGenesList = pairedGenesList,
# #                            outDir = outDir, onlyTop = FALSE, topCutoff = .5,
# #                            sigTest = TRUE)
# #
# #
# #
# #
# # Res2 <- pairedGenesAnalyzeR(pairedGenesList = pairedGenesList,
# #                            outDir = outDir, onlyTop = TRUE, topCutoff = .5,
# #                            sigTest = TRUE)
#
#
#

#
#
#
# library(dplyr)
#
#
# Species <- "hsapiens"
# if (Species == "hsapiens") {
#   msigSpec <- "Homo sapiens"
# } else {
#   msigSpec <- "Mus musculus"
# }
# MDF <- msigdbr::msigdbr(species = msigSpec)
# MDF$gs_subcat <- gsub(MDF$gs_subcat, pattern = "CP:", replacement = "", perl = TRUE)
# MDF$gs_cat <- paste0(MDF$gs_cat, ":", MDF$gs_subcat)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = ":$", replacement = "", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C1", replacement = "Cytogenic bands", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C6", replacement = "Oncogenic signatures", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C7", replacement = "Immunological signatures", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C2:", replacement = "", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C5", replacement = "GO", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "H", replacement = "Hallmark", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "CP", replacement = "Canonical pathways", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "CGP", replacement = "Perturbations", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C4:CGN", replacement = "Cancer gene neighborhoods", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C4:CM", replacement = "Cancer modules", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C3:MIR", replacement = "miRNA targets", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C3:TFT", replacement = "TF targets", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "BIOCARTA", replacement = "BioCarta", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "REACTOME", replacement = "Reactome", perl = TRUE)
# MDFHuman <- MDF %>%
#   select(gs_name, gs_cat, human_gene_symbol)
# Species <- "mmusculus"
# if (Species == "hsapiens") {
#   msigSpec <- "Homo sapiens"
# } else {
#   msigSpec <- "Mus musculus"
# }
# MDF <- msigdbr::msigdbr(species = msigSpec)
# MDF$gs_subcat <- gsub(MDF$gs_subcat, pattern = "CP:", replacement = "", perl = TRUE)
# MDF$gs_cat <- paste0(MDF$gs_cat, ":", MDF$gs_subcat)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = ":$", replacement = "", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C1", replacement = "Cytogenic bands", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C6", replacement = "Oncogenic signatures", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C7", replacement = "Immunological signatures", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C2:", replacement = "", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C5", replacement = "GO", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "H", replacement = "Hallmark", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "CP", replacement = "Canonical pathways", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "CGP", replacement = "Perturbations", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C4:CGN", replacement = "Cancer gene neighborhoods", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C4:CM", replacement = "Cancer modules", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C3:MIR", replacement = "miRNA targets", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C3:TFT", replacement = "TF targets", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "BIOCARTA", replacement = "BioCarta", perl = TRUE)
# MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "REACTOME", replacement = "Reactome", perl = TRUE)
# MDFMouse <- MDF %>%
#   select(gs_name, gs_cat, human_gene_symbol, gene_symbol)
#
# MDF <- full_join(MDFHuman, MDFMouse, by = colnames(MDFHuman))
# MDF <- MDF %>% distinct()
# colnames(MDF)[c(3:4)] <- c("human_gene_symbol", "mouse_gene_symbol")
#
#
# # Get TERM2GENE from MDF object
# MDFtoTERM2GENE <- function(MDF, GSEA_Type, species) {
#
#   # # Bug testing
#   # GSEA_Type <- "Basic"
#   # species <- "hsapiens"
#
#   if (species == "hsapiens") {
#     toGrab <- 3
#   } else {
#     toGrab <- 4
#   }
#   print(GSEA_Type)
#   # Filter for pathways of interest
#   optionsNow <- c("Basic", "All", unique(MDF$gs_cat))
#   GSEA_Type <- gsub(GSEA_Type, pattern = "miRNAs",
#                     replacement = "miRNA targets", perl = TRUE)
#   GSEA_Type <- gsub(GSEA_Type, pattern = "Transcription factors",
#                     replacement = "TF targets", perl = TRUE)
#   GSEA_Type <- gsub(GSEA_Type, pattern = "Biological process",
#                     replacement = "GO:BP", perl = TRUE)
#   GSEA_Type <- gsub(GSEA_Type, pattern = "Cellular component",
#                     replacement = "GO:CC", perl = TRUE)
#   GSEA_Type <- gsub(GSEA_Type, pattern = "Molecular function",
#                     replacement = "GO:MF", perl = TRUE)
#   GSEA_Type <- gsub(GSEA_Type, pattern = "Molecular perturbations",
#                     replacement = "Perturbations", perl = TRUE)
#   GSEA_Type <- gsub(GSEA_Type, pattern = "Immuno signatures",
#                     replacement = "Immunological signatures", perl = TRUE)
#   GSEA_Type <- gsub(GSEA_Type, pattern = "Onco signatures",
#                     replacement = "Oncogenic signatures", perl = TRUE)
#   GSEA_Type <- gsub(GSEA_Type, pattern = "Oncogene neighborhoods",
#                     replacement = "Cancer gene neighborhoods", perl = TRUE)
#
#
#   if (! all(GSEA_Type %in% optionsNow)) {
#     stop("\nPlease enter a valid GSEA_Type. Use ?getTERM2GENE to see available options.\n")
#   }
#   categories <- c()
#   if ("Basic" %in% GSEA_Type) {
#     categories <- c(categories, "Hallmark", "Perturbations", "BioCarta",
#                     "GO:BP", "KEGG", "Canonical pathways", "Reactome", "GO:MF", "GO:CC", "PID")
#   }
#   if ("All" %in% GSEA_Type) {
#     categories <- c(categories, optionsNow)
#   }
#   categories <- unique(c(categories, GSEA_Type))
#   colnames(MDF)[toGrab] <- "gene_symbol"
#   TERM2GENE <- MDF %>%
#     filter(.data$gs_cat %in% categories) %>%
#     select(.data$gs_name, .data$gene_symbol) %>%
#     filter(! is.na(.data$gene_symbol)) %>%
#     distinct()
#
#   return(TERM2GENE)
# }
#
# TERM2GENE <- MDFtoTERM2GENE(MDF = MDF, GSEA_Type = "Canonical pathways",
#                             species = "hsapiens")
# pairedRes <- correlationAnalyzeR::analyzeGenePairs(genesOfInterest = c("BRCA1", "BRCA1"),
#                                                    Sample_Type = c("normal", "cancer"),
#                                                    Tissue = c("all", "all"),
#                                                    returnDataOnly = T,
#                                                    TERM2GENE = TERM2GENE,
#                                                    topPlots = F,
#                                                    # nperm = 500, sampler = T,
#                                                    runGSEA = T)
#
#




genesOfInterest <- c("ATM", "AKT1", "HIF1A", "PRKAA2", "SLC3A2", "MTOR",
                     "SLC7A11", "SLC7A5", "G6PD", "HSPB1")
# Run the function with desired parameters
analyzeSingleGenes(genesOfInterest = genesOfInterest,
                                        outputPrefix = "tests/singleGeneTesting",
                                        Species = "hsapiens",
                                        Sample_Type = "Normal_Tissues")


# devtools::install_github("millerh1/correlationAnalyzeR")
library(correlationAnalyzeR)
# Check the paired genes
pairedGeneList <- list("REST" = c("NCOR1", "HSPB1", "BRCA1", "ATM"))
pairedGenesAnalyzeR(pairedGenesList = pairedGeneList,
                    Sample_Type = "Normal_Tissues",
                    outputPrefix = "tests/pairedOut")


genesOfInterest <- "F2"
analyzeSingleGenes(genesOfInterest = genesOfInterest,
                   Sample_Type = "Normal_Tissues", outputPrefix = "tests/newTest")



analyzeSingleGenes(genesOfInterest = genesOfInterest,
                   outputPrefix = "tests/singleGeneTesting",
                   Species = "hsapiens",
                   Sample_Type = "Normal_Tissues")

# load('data/hsapiens_corrSmall_geneNames.rda')
# mmusculus_corrSmall_geneNames <- geneNames
# usethis::use_data(mmusculus_corrSmall_geneNames)

#  # Test whether pairedGenesAnalyzeR is working
# pairedGenesList <- list("SF3B1" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "SF3B2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "SF3B3" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "U2AF2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "U2AF1" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "SRSF2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "AQR" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "XAB2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "WNT5A" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "BCL2" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "SOX6" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "FAS" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "EWSR1" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"),
#                         "FLI1" = c("THOC6", "DNASE1L1", "COG6", "PCYOX1L", "DRP2"))
#
# outDir <- paste0(getwd(), "/Results/THOC6-TCNER")
#
# PathCardsTCNER <- c("DDB1	RPS27A	UBA52	UBB	UBC	CUL4B	CUL4A	RBX1	CCNH	CDK7	ERCC2	ERCC3
#   GTF2H1	GTF2H2	GTF2H3	GTF2H4	MNAT1	GTF2H5	RPA1	RPA2	RPA3	XPA	ERCC8	EP300
#   ERCC1	ERCC4	ERCC5	ERCC6	GPS1	HMGN1	PCNA	POLD1	POLD2	POLE	POLE2	POLR2A
#   POLR2B	POLR2C	POLR2D	POLR2E	POLR2F	POLR2G	POLR2H	POLR2I	POLR2J	POLR2K	POLR2L	RFC1
#   RFC2	RFC3	RFC4	RFC5	TCEA1	USP7	COPS3	COPS2	AQR	PPIE	POLD3	COPS8
#   COPS6	COPS5	PRPF19	COPS7A	COPS4	POLK	POLE3	POLE4	XAB2	ISY1	UVSSA	POLD4
#   COPS7B	ZNF830	MIR1281	PARP1	CETN2	DDB2	LIG1	LIG3	RAD23A	RAD23B	XPC	XRCC1
#   PARP2	ACTB	ACTL6A	NFRKB	SUMO3	SUMO2	UBE2I	UBE2N	UBE2V2	SUMO1	YY1	PIAS1
#   RUVBL1	CHD1L	PIAS3	MCRS1	TFPT	INO80	RNF111	INO80D	ACTR5	INO80B	USP45	ACTR8
#   INO80C	INO80E	MIR6764	ELL")
# PathCardsTCNER <- unlist(strsplit(PathCardsTCNER, "\n"))
# PathCardsTCNER <- unlist(strsplit(PathCardsTCNER, "\t"))
# DDB1 <- c("DDB1", "DCAF1", "DCAF2", "DCAF3", "DCAF4", "DCAF4L1", "DCAF5", "DCAF6",
#           "DCAF7", "DCAF8", "DCAF9", "DCAF10", "DCAF11", "DCAF12", "DCAF13", "DCAF14",
#           "DCAF15", "DCAF16", "DDA1", "DDB1", "CRBN")
# PathCardsTCNERHigh <- PathCardsTCNER[c(1:20)]
# pairedGenesList <- list("THOC6" = DDB1)
# pairedGenesList <- list("THOC6" = PathCardsTCNER)
# PathCardsTCNER <- trimws(PathCardsTCNER)
#
# load("Data/geneSets/nrf2.gene.sets.bishop.RData")
# source("Code/analyzeGenePairs.R")
#
# pairedGenesList <- list("ATM" = unique(BISHOP_NRF2_PATHWAY))
# outDir <- paste0(getwd(), "/Results/ATM_ROS")
# Res <- pairedGenesAnalyzeR(pairedGenesList = pairedGenesList,
#                            outDir = outDir, onlyTop = F, topCutoff = .5,
#                            sigTest = T)
#
# ATM_NRF2_MARKERS <- Res$ATM
# ATM_NRF2_MARKERS <- ATM_NRF2_MARKERS[order(ATM_NRF2_MARKERS)]
#
# # Met signature
# Met1 <- read.table(file = "Data/geneSets/artega_sig.txt")
# Met2 <- read.table(file = "Data/geneSets/mammaprint_sig_new.txt")
# Met3 <- read.table(file = "Data/geneSets/werb_49_metastasis_sig.txt")
# Met <- unique(c(as.character(Met1$V1), as.character(Met2$V1), as.character(Met3$V1)))
# pairedGenesList <- list("ATM" = Met,
#                         "SLC3A2" = Met)
# outDir <- paste0(getwd(), "/Results/ATM_ROS_MET")
#
#
# Res <- pairedGenesAnalyzeR(pairedGenesList = pairedGenesList,
#                            outDir = outDir, onlyTop = F, topCutoff = .5,
#                            sigTest = T)
#
#
#
#
# Res2 <- pairedGenesAnalyzeR(pairedGenesList = pairedGenesList,
#                            outDir = outDir, onlyTop = T, topCutoff = .5,
#                            sigTest = T)




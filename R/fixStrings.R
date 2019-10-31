#' Fix strings (helper function)
#'
#' Convert vector of GSEA (or other) names to publication-ready titles
#'
#' @param StringVec A vector of titles (usually GSEA) to clean for visualizations
#' @return A vector of cleaned string titles in the same order as the input.
#' @export


fixStrings <- function(StringVec) {

  # # Bug testing
  # StringVec <- c("HALLMARK_APOPTOSIS", "GO_MIR21_TARGETS",
  #                "GSE121239_THING_HAPPENED", "CTTGAT_MIR381",
  #                "GTGCAGAG_EZH2", "GSE12309_WHATEVER",
  #                "MEISSNER_NPC_HCP_WITH_H3K4ME2",
  #                "BASSO_CD40_SIGNALING_DN",
  #                "GSE36078_UNTREATED_VS_AD5_INF_IL1R_KO_MOUSE_LUNG_DC_DN",
  #                "KEGG_ACUTE_MYELOID_LEUKEMIA",
  #                "GSE31082_CD4_VS_CD8_SP_THYMOCYTE_UP",
  #                "GSE3920_IFNA_VS_IFNG_TREATED_ENDOTHELIAL_CELL_UP",
  #                "GSE37605_FOXP3_FUSION_GFP_VS_IRES_GFP_TREG_C57BL6_UP",
  #                "GSE41176_UNSTIM_VS_ANTI_IGM_STIM_BCELL_24H_DN",
  #                "MORI_LARGE_PRE_BII_LYMPHOCYTE_UP",
  #                "RYAAAKNNNNNNTTGW_UNKNOWN",
  #                "GGGTGGRR_PAX4_03",
  #                "GGGNNTTTCC_NFKB_Q6_01",
  #                "AAAYWAACM_HFH4_01",
  #                "KANG_DOXORUBICIN_RESISTANCE_UP")

  StringVec <- gsub(StringVec, pattern = "_", replacement = " ")
  StringVec <- tolower(StringVec)
  StringVec <- stringr::str_to_title(StringVec)
  StringVec <- gsub(StringVec, pattern = "Iii", replacement = "III", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ii", replacement = "II", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Of ", replacement = " of ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " To ", replacement = " to ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " In ", replacement = " in ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " With ", replacement = " with ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Without ", replacement = " without ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Upon ", replacement = " upon ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " An ", replacement = " an ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " By ", replacement = " by ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " For ", replacement = " for ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Via ", replacement = " via ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Lof ", replacement = " LOF ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Lof$", replacement = " LOF", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Loh ", replacement = " LOH ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Loh$", replacement = " LOH", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Arms ", replacement = " ARMS ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Erms ", replacement = " ERMS ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nadh ", replacement = " NADH ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nadph ", replacement = " NADPH ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Cll ", replacement = " CLL ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Cml ", replacement = " CML ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Aml ", replacement = " AML ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " All ", replacement = " ALL ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Kim ALL", replacement = "Kim all", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "CMV ALL", replacement = "CMV all", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nhek ", replacement = " NHEK ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ner ", replacement = " NER ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nmda ", replacement = " NMDA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dc ", replacement = " DC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Cd4 ", replacement = " CD4 ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Cd8 ", replacement = " CD8 ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gc ", replacement = " GC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hdl ", replacement = " HDL ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dn$", replacement = " Down", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ldl ", replacement = " LDL ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Tcr ", replacement = " TCR ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Mdc ", replacement = " MDC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Bcr ", replacement = " BCR ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Icp ", replacement = " ICP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hbv ", replacement = " HBV ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Dlbcl", replacement = "DLBCL", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gist ", replacement = " GIST ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gist$", replacement = " GIST", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dp ", replacement = " DP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dn ", replacement = " DN ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " H2o2 ", replacement = " H2O2 ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " With ", replacement = " with ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ntreg", replacement = "nTreg", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Mlr ", replacement = " MLR ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Gfp", replacement = "GFP", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Vs ", replacement = " vs ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " And ", replacement = " and ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Wt ", replacement = " WT ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ros$", replacement = " ROS", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ros ", replacement = " ROS ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ko ", replacement = " KO ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Pdc ", replacement = " PDC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Pdgf", replacement = "PDGF", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Rna ", replacement = " RNA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mrna", replacement = "mRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mirna", replacement = "miRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Sirna", replacement = "siRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Trna", replacement = "tRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ncrna", replacement = "ncRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Snrna", replacement = "snRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Rrna", replacement = "rRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "rna$", replacement = "RNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "rna ", replacement = "RNA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Flii", replacement = "Fli1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hcp ", replacement = " HCP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Tnf ", replacement = " TNF ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Srp ", replacement = " SRP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Utr ", replacement = " UTR ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dna ", replacement = " DNA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Rdna", replacement = "rDNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hiv ", replacement = " HIV ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hiv1 ", replacement = " HIV1 ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "dna", replacement = "DNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Lps ", replacement = " LPS ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gmcsf ", replacement = " GMCSF ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gm Csf ", replacement = " GMCSF ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Balbc", replacement = "BALBc", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Lcmv", replacement = "LCMV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mcmv", replacement = "MCMV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Pcc ", replacement = " PCC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ecm", replacement = "ECM", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "G1s", replacement = "G1S", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " G1 S ", replacement = "G1S", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "G2m", replacement = "G2M", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " G2 M ", replacement = "G2M", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Hcmv", replacement = "HCMV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Pbmc", replacement = "PBMC", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Atp ", replacement = " ATP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Atp$", replacement = " ATP", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Gtp", replacement = "GTP", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Mut ", replacement = " MUT ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Et Al", replacement = "et al", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Cpg", replacement = "CPG", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nkt ", replacement = " NKT ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hsc ", replacement = " HSC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ln ", replacement = " LN ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Cmv", replacement = "CMV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Bm ", replacement = " BM ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Bmdc ", replacement = " BMDC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Esc ", replacement = "ESC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Esc ", replacement = " ESC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mcf10a", replacement = "MCF10A", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Tca", replacement = "TCA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Nkcell", replacement = "NK-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Tcell", replacement = "T-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "T Cell", replacement = "T-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "B Cell", replacement = "B-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Hela", replacement = "HeLa", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Treg", replacement = "T-Reg", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Tconv", replacement = "T-Conv", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Bcell", replacement = "B-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Uv ", replacement = " UV ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Uv$", replacement = " UV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Gse", replacement = "GSE", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Gnf2", replacement = "GNF2", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Gcm", replacement = "GCM", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Morf", replacement = "MORF", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Kegg", replacement = "KEGG", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Pid", replacement = "PID", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Go ", replacement = "GO ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^([GCATSMWNRYK][gcatsmwrnyk]+)( [A-Za-z0-9]+ Q[0-9]$)",
                    replacement = "\\U\\1\\E\\2", perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^([GCATSMWNRYK][gcatsmwrnyk]+)( [A-Za-z0-9]+ Q[0-9] [0-9]+$)",
                    replacement = "\\U\\1\\E\\2", perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^([GCATSMWNRYK][gcatsmwrnyk]+)( [A-Za-z0-9]+ [0-9]+$)",
                    replacement = "\\U\\1\\E\\2", perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^([GCATSMWNRYK][gcatsmwrnyk]+)( Unknown$)",
                    replacement = "\\U\\1\\E\\2", perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mir([0-9]*.*)", replacement = "miR\\1",
                    perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ifn([a-z])", perl = TRUE,
                    replacement = "IFN\\1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Inf([a-z])", perl = TRUE,
                    replacement = "INF\\1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Il([0-9]+)", perl = TRUE,
                    replacement = "IL\\1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "(IL[0-9]+)(r)", perl = TRUE,
                    replacement = "\\1\\U\\2", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Cd([0-9])", perl = TRUE,
                    replacement = "CD\\1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "(Ig)([a-z])", perl = TRUE,
                    replacement = "\\1\\U\\2", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "(H[0-9])(k[0-9]+)", perl = TRUE,
                    replacement = "\\1\\U\\2", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "(B)(i+) ", perl = TRUE,
                    replacement = "\\1\\U\\2 ", ignore.case = FALSE)

  return(StringVec)
}

## correlationAnalyzeR
[![releaseVersion](https://img.shields.io/badge/Version-1.0.0-blue.svg)](https://github.com/millerh1/correlationAnalyzeR)

The correlationAnalyzeR package uses correlation data derived from the [ARCHS4](https://amp.pharm.mssm.edu/archs4/index.html) RNA Seq repository to generate biological insights about a gene or genes of interest.

## Motivation
This project is motivated by a recurring issue that arises during exploratory bioinformatics: *Sometimes little to no information exists about a gene or gene(s) of interest.* 
One way to address this problem is to compute gene expression correlations. These values indicate how genes vary in relation to eachother.

With gene correlation data, it is possible to implement three levels of analysis:
- *Single gene*: Analyses such as [GSEA](http://software.broadinstitute.org/gsea/index.jsp) (Gene Set Enrichment Analysis) which predict
the biological pathways correlated with a gene of interest.
- *Paired gene set*: Statistical approaches to determine if a gene set is significantly correlated with a gene of interest.
- *Gene set toplogy*: Methods to uncover the topology of a gene set by using gene correlation data as the input for tools such as UpSet 
and hierarchical clustering. 

## Contributing [![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/millerh1/correlationAnalyzeR/issues)
Please report issues and feel free to contribute! I hope this project can some day incorporate other sources of data such as ppi networks
and continue to develop as a useful tool for anyone who wants to understand more about their HTS data.

## Installation

First CRAN release coming soon... For now, correlationAnalyzeR can be installed from the development verion:

``` r
## install.packages("devtools")
devtools::install_github("millerh1/correlationAnalyzeR")
```

## Software/tools cited

- Bairoch,A. (2018) The cellosaurus, a cell-line knowledge resource. J. Biomol. Tech.
- Bengtsson,H. (2019) matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors).
- Canty,A. and Ripley,B. (2019) boot: Bootstrap Functions (Originally by Angelo Canty for S).
- Chang,W. et al. (2019) shiny: Web Application Framework for R.
- Dolgalev,I. (2018) msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format.
- Galili et al. (2017) heatmaply: an R package for creating interactive cluster heatmaps for online publishing. Bioinformatics.
- Kassambara,A. (2019) ggpubr: ‘ggplot2’ Based Publication Ready Plots.
- Kolde,R. (2019) pheatmap: Pretty Heatmaps.
- Krijthe,J. and Maaten,L. van der (2018) Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation.
- Lachmann,A. et al. (2018) Massive mining of publicly available RNA-seq data from human and mouse. Nat. Commun., 9.
- Langfelder,P. and Horvath,S. (2008) WGCNA: An R package for weighted correlation network analysis. BMC Bioinformatics.
- Liberzon,A. et al. (2011) Molecular signatures database (MSigDB) 3.0. Bioinformatics.
- Liberzon,A. et al. (2015) The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst.
- Love,M.I. et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol.
- R Core Team (2019) R: A Language and Environment for Statistical Computing.
- Sergushichev,A.A. (2016) An algorithm for fast preranked gene set enrichment analysis using cumulative statistic calculation. bioRxiv.
- Sievert,C. (2018) plotly for R.
- Subramanian,A. et al. (2005) Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U. S. A.
- Wickham,H. (2016) ggplot2: Elegant Graphics for Data Analysis Springer-Verlag New York.
- Xie,Y. et al. (2019) DT: A Wrapper of the JavaScript Library ‘DataTables’.
- Yu,G. et al. (2012) ClusterProfiler: An R package for comparing biological themes among gene clusters. Omi. A J. Integr. Biol.






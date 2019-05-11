## correlationAnalyzeR
[![releaseVersion](https://img.shields.io/badge/Development%20Version-0.9.0-orange.svg)](https://github.com/millerh1/correlationAnalyzeR)

The correlationAnalyzeR package offers biological insights about a gene or genes from correlation data derived from the [ARCHS4](https://amp.pharm.mssm.edu/archs4/index.html) RNA Seq repository.

## Motivation
This project is motivated by a recurring issue that arises during exploratory bioinformatics: *Sometimes little to no information exists about a gene or gene(s) of interest.* 
One way to address this problem is to compute gene expression correlations, the quantitative measurement of how much any particular gene's expression varies with any other gene.

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

First release coming soon... For now, correlationAnalyzeR can be installed from the development verion:

``` r
## install.packages("devtools")
devtools::install_github("millerh1/correlationAnalyzeR")
```


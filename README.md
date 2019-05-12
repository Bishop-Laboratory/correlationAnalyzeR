## correlationAnalyzeR
[![releaseVersion](https://img.shields.io/badge/Development%20Version-0.9.0-orange.svg)](https://github.com/millerh1/correlationAnalyzeR)

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

First release coming soon... For now, correlationAnalyzeR can be installed from the development verion:

``` r
## install.packages("devtools")
devtools::install_github("millerh1/correlationAnalyzeR")
```

## Software/tools cited
[![Rdoc](http://www.rdocumentation.org/badges/version/sva)](http://www.rdocumentation.org/packages/sva)

- Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi: 10.1089/omi.2011.0118.
- Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, Ma’ayan A. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications 9. Article number: 1366 (2018), doi:10.1038/s41467-018-03751-6
- Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE, Storey JD, Zhang Y, Torres LC (2019). sva: Surrogate Variable Analysis. R package version 3.32.0.
- Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, http://ggplot2.org.
- Kassambara A (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. https://rpkgs.datanovia.com/ggpubr/

## License
[![Open Source Love](https://badges.frapsoft.com/os/mit/mit.svg?v=102)](https://opensource.org/licenses/MIT)

MIT License

Copyright (c) 2019 Henry Miller, Bishop Laboratory, UT Health San Antonio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


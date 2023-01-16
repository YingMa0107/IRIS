# IRIS


![IRIS\_pipeline](IRIS_logo.pdf)
We developed a statistical method for Integrative and Reference-Informed Spatial Domain Detection for Spatial Transcriptomics. IRIS is a reference-informed integrative method for detecting spatial domains on multiple tissue slices from spatial transcriptomics with spot-level, single-cell level, or subcellular level resolutions. IRIS builds upon the fact that each spatial domain on the tissue is characterized by a unique composition of cell types and that similar composition is observed for the same spatial domain across different slices of the same tissue 28,29. Consequently, IRIS integrates a reference scRNA-seq data to inform and characterize the cell type composition on the tissue of spatial transcriptomics for accurate and interpretable spatial domain detection. In the process, IRIS accommodates cell type compositional similarity across locations within the same slice and across different slices on the same domain to borrow information both within and between tissue slices for integrative and accurate spatial domain detection. Importantly, IRIS comes with an efficient optimization framework with multiple algebraic innovations for scalable computation and can easily handle multiple spatial transcriptomics datasets with millions of spatial locations and tens of thousands of genes. IRIS is implemented as an open-source R package, freely available at www.xzlab.org/software.html. 

Installation
------------
You can install the released version of IRIS from Github with the following code, for more installation details or solutions that might solve related issues (specifically MacOS system) see the [link](https://yingma0107.github.io/IRIS/documentation/02_installation.html).

## Dependencies 
* R version >= 4.2.2.
* Dependent R packages: Rcpp (>= 1.0.9), RcppArmadillo, SingleCellExperiment, SummarizedExperiment, methods, Matrix, MCMCpack, fields, wrMisc, RANN, stats, ggplot2, grDevices, reshape2


``` r
# install devtools if necessary
install.packages('devtools')

# install the IRIS package
devtools::install_github('YingMa0107/IRIS')

# load package
library(IRIS)

```
The R package has been installed successfully on Operating systems: 
* macOS Catalina 10.15, macOS Monterey 12.4
* Ubuntu 18.04.6 LTS
* Windows 10

# Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible example and also please provide the output of your sessionInfo() in R! 

How to cite `IRIS`
-------------------
Ying Ma, Xiang Zhou. Integrative and Reference-Informed Spatial Domain Detection for Spatial Transcriptomics, 2023. 

How to use `IRIS`
-------------------
Details in [Tutorial](https://yingma0107.github.io/IRIS/)

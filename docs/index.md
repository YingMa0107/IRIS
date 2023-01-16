---
layout: full
homepage: true
disable_anchors: true
description: Integrative and Reference-Informed Spatial Domain Detection for Spatial Transcriptomics 
---
## IRIS Overview
<p align="center">
<img src="IRIS_logo.png" width="300" />
IRIS is a reference-informed integrative method for detecting spatial domains on multiple tissue slices from spatial transcriptomics with spot-level, single-cell level, or subcellular level resolutions. IRIS builds upon the fact that each spatial domain on the tissue is characterized by a unique composition of cell types and that similar composition is observed for the same spatial domain across different slices of the same tissue 28,29. Consequently, IRIS integrates a reference scRNA-seq data to inform and characterize the cell type composition on the tissue of spatial transcriptomics for accurate and interpretable spatial domain detection. In the process, IRIS accommodates cell type compositional similarity across locations within the same slice and across different slices on the same domain to borrow information both within and between tissue slices for integrative and accurate spatial domain detection. Importantly, IRIS comes with an efficient optimization framework with multiple algebraic innovations for scalable computation and can easily handle multiple spatial transcriptomics datasets with millions of spatial locations and tens of thousands of genes. IRIS is implemented as an open-source R package, freely available at [www.xzlab.org/software.html](http://www.xzlab.org/software.html). 

### Example Analysis with IRIS: [here](https://yingma0107.github.io/IRIS/documentation/04_IRIS_Example.html).

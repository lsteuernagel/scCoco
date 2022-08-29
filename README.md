
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scCoco

<!-- badges: start -->
<!-- badges: end -->

Annotation of single cell clusters to brain regions using the common
coordinate framework and ISH data from Allen brain atlas

This package provides wrapper functions around the [cocoframer R
package](https://github.com/AllenInstitute/cocoframer) that allows to
query the [Allen brain atlas API](atlas.brain-map.org).

## Installation

You can install the development version of scCoco like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

<!-- README.md is generated from README.Rmd. Please edit that file -->

# scUtils

<!-- badges: start -->
<!-- badges: end -->

Single-cell utility functions by Lukas Steuernagel

## Installation

Install scUtils using:

``` r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("lsteuernagel/scCoco")
```

Depends on the cocoframer R package.

## Function overview

Load package:

``` r
library(scCoco)
```

List of functions:

``` r
sort(as.character(lsf.str("package:scCoco")))
#>  [1] "apply_DoubletFinder"           "CalculateMultScore"           
#>  [3] "conserved_correlations"        "determine_cluster_resolution" 
#>  [5] "downsample_balanced_iterative" "feature_statistics"           
#>  [7] "find_children"                 "gene_pct_cluster"             
#>  [9] "identify_variable_features"    "infer_sex"                    
#> [11] "match_sample_names"            "rasterize_ggplot"             
#> [13] "Read10xFormat"                 "ReadDGEFormat"                
#> [15] "seurat_recipe"                 "writeList_to_JSON"
```

Check the man pages for details!

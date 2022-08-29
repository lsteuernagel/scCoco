
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

Install scCoco using:

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
#> [1] "aba_ids"                   "coco_expression"          
#> [3] "find_children"             "findRegions_genesets"     
#> [5] "intersect_from_list"       "median_rank"              
#> [7] "region_annotation"         "summariseRegions_genesets"
```

Check the man pages for details!


<!-- README.md is generated from README.Rmd. Please edit that file -->

# ComutR

<!-- badges: start -->
<!-- badges: end -->

Create comut plots in R! Built on top of ComplexHeatmap
[ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
Heavily inspired by the python package comut:
[comut](https://pypi.org/project/comut/)

## Installation

You can install the development version of ComutR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("taprati/ComutR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ComutR)

input_maf <- data.frame(
  Tumor_Sample_Barcode = c("1", "1", "1", "2", "3", "4", "4"),
  Hugo_Symbol = c("A", "B", "C", "C", "A", "A", "B"),
  Variant_Classification = c("Missense_Mutation", "Nonsense_Mutation", "In_Frame_Del", "In_Frame_Del", "Missense_Mutation", "Nonsense_Mutation", "Nonsense_Mutation")
)

comut(data = input_maf)
#> No ids specified. Defaulting to all ids!
```

<img src="man/figures/README-example-1.png" width="100%" />

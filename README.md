
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PhenoR <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/AtilioRausch/PhenoR)](https://github.com/AtilioRausch/PhenoR/issues)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

The goal of `PhenoR` is to infer synteny networks from whole-genome
protein sequence data and analyze them. Anchor pairs from synteny
analyses are treated as an undirected unweighted graph (i.e., a synteny
network), and users can perform:

  - **Synteny detection** using a native implementation of the [MCScanX
    algorithm](https://doi.org/10.1093/nar/gkr1293), a C++ program that
    has been modified and ported to R with Rcpp. This way, users do not
    need to install MCScanX beforehand, because `PhenoR` has its own
    implementation of the same algorithm.
  - **Synteny network inference** by treating anchor pairs as edges of a
    graph;
  - **Network clustering** using the Infomap algorithm;
  - **Phylogenomic profiling**, which consists in identifying which
    species contain which clusters. This analysis can reveal highly
    conserved synteny clusters and taxon-specific ones (e.g., family-
    and order-specific clusters);
  - **Microsynteny-based phylogeny reconstruction** with maximum
    likelihood, which can be achieved by inferring a phylogeny from a
    binary matrix of phylogenomic profiles with IQTREE.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `PhenoR` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("PhenoR")
```

And the development version from
[GitHub](https://github.com/AtilioRausch/PhenoR) with:

``` r
BiocManager::install("AtilioRausch/PhenoR")
```

## Citation

Below is the citation output from using `citation('PhenoR')` in R.
Please run this yourself to check for any updates on how to cite
**PhenoR**.

``` r
print(citation('PhenoR'), bibtex = TRUE)
#> Warning in citation("PhenoR"): could not determine year for 'PhenoR' from
#> package DESCRIPTION file
#> To cite package 'PhenoR' in publications use:
#> 
#>   Rausch A (????). _PhenoR: What the Package Does (One Line, Title
#>   Case)_. R package version 1.0.0.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {PhenoR: What the Package Does (One Line, Title Case)},
#>     author = {Atilio O. Rausch},
#>     note = {R package version 1.0.0},
#>   }
```

Please note that `PhenoR` was only made possible thanks to many other R
and bioinformatics software authors, which are cited either in the
vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `PhenoR` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

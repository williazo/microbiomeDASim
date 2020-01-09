Microbiome Differential Abundance Simulation
================

[![Bioc
build](http://bioconductor.org/shields/build/release/bioc/microbiomeDASim.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/microbiomeDASim/)
[![Bioc
time](http://bioconductor.org/shields/years-in-bioc/microbiomeDASim.svg)](https://bioconductor.org/packages/microbiomeDASim)
[![Bioc
downloads](http://bioconductor.org/shields/downloads/release/microbiomeDASim.svg)](http://bioconductor.org/packages/stats/bioc/microbiomeDASim/)

This package is developed to simulate microbiome data for longitudinal
differential abundance analyses. Microbiome data have a variety of
features that make typical simulation methods inappropriate. For an
in-depth description of the types of problems this simulation package is
designed to solve, plus basics of the functionality please refer to the
F1000 manuscript
[f1000research.20660.1](https://doi.org/10.12688/f1000research.20660.1).

## Installation

To install the `microbiomeDASim` package the latest release version is
available from Bioconductor

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("microbiomeDASim")
```

Alternatively to use the latest development version from Bioconductor
use the following commands

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("microbiomeDASim", version="devel")
```

## Examples

An interactive examples for how to simulate data from a multivariate
normal distribution and fit a trend line using
[metagenomeSeq::fitTimeSeries](https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html)
are available in the `inst/scripts` directory.

This notebook can be run interactively using Google Collab by clicking
the “Open in Colab” marker at the top of the notebook.

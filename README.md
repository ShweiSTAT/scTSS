scTSS for studying transcription start sites using 5’ scRNA-seq data
================
Shiwei Fu, Wei Vivian Li
2024-07-09

## Introduction

scTSS is a package that can leverage the 5’ single-cell RNA-seq data
(both paired-end and single-end data) to study the transcription start
site (TSS) at single-cell resolution. scTSS comprises of two main
stages, TSS prediction and differential TSS usage (DU) test for TSSs
across two or multiple conditions or cell types.

## Installation

`scTSS` is currently available on GitHub. To download and load package
please follow:

``` r
devtools::install_github("https://github.com/ShweiSTAT/scTSS")
library(scTSS)
```

## Usage

For more information, please use this
[link](https://github.com/ShweiSTAT/scTSS/wiki) to see the detailed
tutorial.

## Issues and communications

If you have any issues using this package, please post them
[here](https://github.com/ShweiSTAT/scTSS/issues). Any suggestions and
comments are welcome! For suggestions and comments, please contact
Shiwei Fu (<shiwei.fu@email.ucr.edu>) or Wei Vivian Li (<weil@ucr.edu>).

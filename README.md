scTSS for studying transcription start sites using 5’ scRNA-seq data
================
Shiwei Fu, Wei Vivian Li
2023-02-22

## Introduction

scTSS is a bioinformatics tool that leverages 5’ single-cell RNA
sequencing (scRNA-seq) data to predict transcription start sites (TSSs)
in single cells and identify alternative TSS usage between different
conditions or cell types. scTSS consists of two primary functions: TSS
prediction and differential TSS usage (DU) test.

For the TSS prediction function, scTSS is capable of analyzing both
single-end and paired-end data, yielding robust and accurate TSS
predictions. For the DU test function, scTSS can perform DU tests using
multiple single-cell samples across biological conditions or cell types.

## Installation

`scTSS` is currently available on GitHub. To download and load package
please follow:

``` r
devtools::install_github("https://github.com/ShweiSTAT/scTSS")
library(scTSS)
```

## Usage

For more information, please use this
[link](https://github.com/ShweiSTAT/scTSS/wiki/scTSS-Tutorial) to see
the detailed tutorial.

## Issues and communications

If you have any issues using this package, please post them
[here](https://github.com/ShweiSTAT/scTSS/issues). Any suggestions and
comments are welcome! For suggestions and comments, please contact
Shiwei Fu (<shiwei.fu@email.ucr.edu>) or Wei Vivian Li (<weil@ucr.edu>).

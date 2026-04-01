
## Transfer Learning for Robust Structured Regression with Bi-level Source Detection

<!-- badges: start -->
<!-- badges: end -->

A method to simultaneously tackle data limitation and contamination under the structured regression setting by integrating transfer learning with the robust L2E criterion.

*Jointly developed with [Xiaoqian Liu](https://xiaoqian-liu.github.io/)
*

## Introduction
This repository presents the R package `TransL2E` developed for the manuscript titled "Transfer Learning for Robust Structured Regression with Bi-level Source Detection". It offers transfer learning methods that accounts for contamination in both target and source data while effectively transferring information from source to target. The current package support both sparse regression and group-sparse regression.

## Installation

To install the latest stable version of the package (from GitHub):

``` r
devtools::install_github("haomingsj98/TransL2E")
```

## Implementation

To run TransL2E with sparse regression, use the `TL_L2E_sparse` function with the required input.

To run TransL2E with group-sparse regression, use the `TL_L2E_glasso` function with the required input.

Please see the introductory [demo](https://haomingsj98.github.io/TransL2E/articles/TransL2E_Demo.html) on how to use the TransL2E framework with sparse regression examples.


## Contact

Please reach out to xiaoqian.liu@ucr.edu for any questions. 



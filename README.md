
# Transfer Learning for Robust Structured Regression with Bi-level Source Detection

<!-- badges: start -->
<!-- badges: end -->

Transfer learning has emerged as a powerful tool for mitigating data limitations across a variety of statistical settings. However, a critical practical issue that has been overlooked in this context is data contamination, which may significantly degrade the performance of transfer learning methods. 
To bridge this gap, we propose a novel approach that simultaneously tackle data limitation and contamination under the structured regression setting.  By integrating transfer learning with the robust L2E criterion, we develop the TransL2E method that accounts for contamination in both target and source data while effectively transferring information from source to target. Beyond robust estimation,  TransL2E introduces a data-driven bi-level source detection mechanism, operating at both individual and cohort levels, which possesses multiple advantages over existing source detection approaches.  The proposed TransL2E method is a general and flexible framework that enables transfer learning for robust estimation across a wide spectrum of structural assumptions. %including sparsity, group sparsity, and low-rankness. 
Comprehensive simulation studies and a real data application demonstrate the superior performance of TransL2E in both robust estimation and structure recovery in the presence of data limitation and contamination.

## Installation

To install the latest stable version of the package (from GitHub):

``` r
devtools::install_github("haomingsj98/TranL2E")
```




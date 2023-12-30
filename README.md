
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRCIGVAR

<!-- badges: start -->
<!-- badges: end -->

MRCIGVAR includes 8 classes of econometric models:

    + VAR: vector autoregressive models
    + CIVAR: cointegrated vector autoregressive models which are also known as vector error correction models (VECM)
    + GVAR: global vector autoregressive models
    + CIGVAR: cointegrated global vector autoregressive models
    + MRVAR: multi regime vector autoregressive models
    + MRCIVAR: multi regime cointegrated vector autoregressive models
    + MRGVAR: multi regime global vector autoregressive models
    + MRCIGVAR: multi regime cointegrated global vector autoregressive models

This package contains functions for data generation with flexible model
specifications, parameter estimation, model selection, impulse response
functions, and parameter constraint tests for each model class. Through
sample codes, this introduction demonstrates practical applications of
these functions.

## Installation

You can install the development version of MRCIGVAR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("puchen8229/MRCIGVAR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MRCIGVAR)
## basic example code
```

A multi regime VAR model is defined:

![Y_t = \Big(C_o^{(1)} + \sum\_{j=1}^{p}B\_{j}^{(1)} Y\_{t}+u_t^{(1)}\Big){\bf 1}\_{\[f\_{t-d} \le \tau\]}+\Big(C_o^{(2)} + \sum\_{j=1}^{p}B\_{j}^{(2)} Y\_{t}+u_t^{(2)}\Big){\bf 1}\_{\[f\_{t-d} \> \tau\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_t%20%3D%20%5CBig%28C_o%5E%7B%281%29%7D%20%2B%20%5Csum_%7Bj%3D1%7D%5E%7Bp%7DB_%7Bj%7D%5E%7B%281%29%7D%20Y_%7Bt%7D%2Bu_t%5E%7B%281%29%7D%5CBig%29%7B%5Cbf%201%7D_%7B%5Bf_%7Bt-d%7D%20%5Cle%20%5Ctau%5D%7D%2B%5CBig%28C_o%5E%7B%282%29%7D%20%2B%20%5Csum_%7Bj%3D1%7D%5E%7Bp%7DB_%7Bj%7D%5E%7B%282%29%7D%20Y_%7Bt%7D%2Bu_t%5E%7B%282%29%7D%5CBig%29%7B%5Cbf%201%7D_%7B%5Bf_%7Bt-d%7D%20%3E%20%5Ctau%5D%7D "Y_t = \Big(C_o^{(1)} + \sum_{j=1}^{p}B_{j}^{(1)} Y_{t}+u_t^{(1)}\Big){\bf 1}_{[f_{t-d} \le \tau]}+\Big(C_o^{(2)} + \sum_{j=1}^{p}B_{j}^{(2)} Y_{t}+u_t^{(2)}\Big){\bf 1}_{[f_{t-d} > \tau]}")

where
![u^{(k)}\_t \sim N(0,\Sigma_0^{(k)})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u%5E%7B%28k%29%7D_t%20%5Csim%20N%280%2C%5CSigma_0%5E%7B%28k%29%7D%29 "u^{(k)}_t \sim N(0,\Sigma_0^{(k)})")
for
![k= 1,2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%3D%201%2C2 "k= 1,2").
This model is also called threshold VAR model as the choice of regime
depends on the threshold condition
![f\_{t-d}\le \tau](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_%7Bt-d%7D%5Cle%20%5Ctau "f_{t-d}\le \tau").
Often
![f_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f_t "f_t")
is one component of the endogenous variable
![Y_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_t "Y_t").
In this case, it is called self-exiting threshold VAR model.

The main functions for this class of models are

- MRVARData
- MRVARest
- MRVAR_Select
- irf_MRVAR_NM_CB
- girf_MRVAR_RM_CB

The following sample codes demonstrate how MRVARData is used to generate
data from an MRVAR process.

``` r
p = matrix(c(2,1,0,0),2,2)
res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
max(res_d$Y)
colnames(res_d$Y) = c("R","P")
res_e = MRVARest(res=res_d)
res_e$Summary
```

For more detailed use of MRCIGVAR please read Introduction to
MRCIGVAR.html in the vignettes directory, or click: https://htmlpreview.github.io/?https://github.com/puchen8229/MRCIGVAR/blob/main/vignettes/Introduction to MRCIGVAR.html

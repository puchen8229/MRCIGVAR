
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

## VAR/CIVAR Models

### VAR

A vector autoregressive (VAR) model is defined as follows.

![Y_t = C_0 + B_1Y\_{t-1} + B_2Y\_{t-1} + ... + B_pY\_{t-p} + u_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_t%20%3D%20C_0%20%2B%20B_1Y_%7Bt-1%7D%20%2B%20B_2Y_%7Bt-1%7D%20%2B%20...%20%2B%20B_pY_%7Bt-p%7D%20%2B%20u_t "Y_t = C_0 + B_1Y_{t-1} + B_2Y_{t-1} + ... + B_pY_{t-p} + u_t")

with
![u_t \sim N(0,\Sigma_0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_t%20%5Csim%20N%280%2C%5CSigma_0%29 "u_t \sim N(0,\Sigma_0)").

The main functions for this class of models are

- VARData
- VARest
- irf_VAR_CB

The following codes demonstrate how VARData is used to generate data
from a VAR(p) process with a variety of specifications.

``` r
res_d = VARData(n=2,p=2,T=100,type="const")
res_d = VARData(n=2,p=2,T=100,Co=c(1:2)*0,type="none")
res_d = VARData(n=2,p=2,T=100,Co=c(1:2),  type="const")
res_d = VARData(n=3,p=2,T=200,type="exog1",X=matrix(rnorm(400),200,2))
res_d = VARData(n=3,p=2,T=200,Co=matrix(c(0,0,0,1,2,3,3,2,1),3,3), type="exog0",X=matrix(rnorm(400),200,2))
```

VARData can be used to generate data for given parameters. If parameters
are not provided VARData will generate data by default with randomly
chosen parameters.

The output of VARData is a VAR object containing generated data, used
parameters as well as the residuals. The components of a generated VAR
object correspond to the model variables and parameters as follows.

res_d\$Y:
![(Y_1,Y_2,...,Y_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28Y_1%2CY_2%2C...%2CY_T%29%27 "(Y_1,Y_2,...,Y_T)'")

res_d\$Co:
![C_o](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_o "C_o")

res_d\$B:
![B_1,B_2,...B_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B_1%2CB_2%2C...B_p "B_1,B_2,...B_p")

res_d\$Sigmao:
![\Sigma_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma_0 "\Sigma_0")

res_d\$resid:
![(u_1,u_2,...,u_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28u_1%2Cu_2%2C...%2Cu_T%29%27 "(u_1,u_2,...,u_T)'")

res_d\$p:
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")

The VARest function is used to estimate the parameters for the given
data.

``` r
res_d = VARData(n=2,p=2,T=100,type="const")
res_e = VARest(res=res_d);
IRF_CB = irf_VAR_CB(res=res_e,nstep=20, comb=NA, irf = "gen1", runs = 100, conf = c(0.05, 0.95))
IRF_list <-IRF_graph(IRF_CB)
res_e$Summary
```

The output of VARest is a VAR object with estimated parameters.

For empirical applications it is recommended to follow the steps below:

- Use VARData to generate a VAR object.
- Replace the generated model data with the empirical data.
- Run a model selection procedure to determine the specification of the
  model if necessary.
- Use VARData to re-generate a VAR object with the selected
  specification and replace the model data with the empirical data
- Estimate the selected model.
- Run the analysis.

``` r
library(vars)
#> Warning: package 'vars' was built under R version 4.0.5
#> Loading required package: MASS
#> Loading required package: strucchange
#> Loading required package: zoo
#> 
#> Attaching package: 'zoo'
#> The following objects are masked from 'package:base':
#> 
#>     as.Date, as.Date.numeric
#> Loading required package: sandwich
#> Warning: package 'sandwich' was built under R version 4.0.5
#> Loading required package: urca
#> Loading required package: lmtest
data(Canada)
dim(Canada)
trend = 1:84
# step 1: generating a VAR object
res_d = VARData(n=4,p=2,T=84,Co=matrix(c(0,0,0,0),4,1),type="none")

# step 2: Replacing the model data with the empirical data
res_d$Y = as.matrix(Canada)
res_e = VARest(res=res_d)
res_e$Summary$BIC
res_e$Summary$AIC

# step 3 Model selection based on AIC or BIC

res_d = VARData(n=4,p=2,T=84,Co=matrix(c(1,1,1,1),4,1),type="const")
res_d$Y = as.matrix(Canada)
res_e = VARest(res=res_d)
res_e$Summary$BIC
res_e$Summary$AIC

trend=as.matrix(trend)
colnames(trend) = c("TREND")

res_d = VARData(n=4,p=2,T=84,Co=matrix(c(0,0,0,0,1,1,1,1),4,2),type="exog0",X=trend)
res_d$Y = as.matrix(Canada)
res_e = VARest(res=res_d)
res_e$Summary$BIC
res_e$Summary$AIC

res_d = VARData(n=4,p=2,T=84,Co=matrix(c(1,1,1,1,1,1,1,1),4,2),type="exog1",X=trend)
res_d$Y = as.matrix(Canada)
res_e = VARest(res=res_d)
res_e$Summary$BIC
res_e$Summary$AIC
```

Comparing the AIC values of the 4 different specifications, the model
with constant intercept and no trend is the most appropriate model for
the data. Some useful functions for modelling VARs can be found in the R
package
![vars](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;vars "vars").
In particular, the function VARselect of
![vars](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;vars "vars")
can be used to specify the lag length of a VAR model.

``` r
# step 4: Estimate the chosen model
res_d = VARData(n=4,p=2,T=84,Co=matrix(c(1,1,1,1),4,1),type="const")
res_d$Y = as.matrix(Canada)
res_e = VARest(res=res_d)
res_e$Summary

# step 5: IRF
IRF_CB = irf_VAR_CB(res=res_e,nstep=50, comb=NA, irf = "gen1", runs = 100, conf = c(0.05, 0.95))
IRF_list <-IRF_graph(IRF_CB)
```

In the following example is a VARX model with two exogenous variables.
Using the option irf=“irfX” we obtain impulse response functions to
exogenous shocks. The option Xshk=c(1:2) is used to generate two
one-unit shocks of the exogenous variables.

``` r
res_d = VARData(n=3,p=2,T=200,type="exog1",X=matrix(rnorm(400),200,2))
res_e = VARest(res=res_d);

IRF_CB = irf_VAR_CB(res=res_e,nstep=20, comb=NA, irf = "irfX", runs = 100, conf = c(0.05, 0.95),Xshk=c(1:2))
IRF_list <-IRF_graph(IRF_CB, response=c(1,2,3),Names=c("Y1","Y2","Y3"),IName=c("X1","X2"), impulse=c(1,2),ncol=2)
```

### CIVAR

A cointegrated vector autoregressive model (CIVAR) is defined as
follows.

![Y_t = C_0 + B_1Y\_{t-1} + B_2Y\_{t-1} + ... + B_pY\_{t-p} + u_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_t%20%3D%20C_0%20%2B%20B_1Y_%7Bt-1%7D%20%2B%20B_2Y_%7Bt-1%7D%20%2B%20...%20%2B%20B_pY_%7Bt-p%7D%20%2B%20u_t "Y_t = C_0 + B_1Y_{t-1} + B_2Y_{t-1} + ... + B_pY_{t-p} + u_t")

Because of unit roots and cointegration relations in the VAR process, it
has the following vector error correction form (VECM):

![\Delta Y_t = C_0 + \alpha\beta' Y\_{t-1} + B^\*\_1\Delta Y\_{t-1} + ... + B^\*\_{p-1}\Delta Y\_{t-p+1} + u_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20Y_t%20%3D%20C_0%20%2B%20%5Calpha%5Cbeta%27%20Y_%7Bt-1%7D%20%2B%20B%5E%2A_1%5CDelta%20Y_%7Bt-1%7D%20%2B%20...%20%2B%20B%5E%2A_%7Bp-1%7D%5CDelta%20Y_%7Bt-p%2B1%7D%20%2B%20u_t "\Delta Y_t = C_0 + \alpha\beta' Y_{t-1} + B^*_1\Delta Y_{t-1} + ... + B^*_{p-1}\Delta Y_{t-p+1} + u_t")

where
![\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha "\alpha")
and
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
are
![(n\times ckr)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28n%5Ctimes%20ckr%29 "(n\times ckr)")
matrices and
![u_t \sim N(0,\Sigma_0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u_t%20%5Csim%20N%280%2C%5CSigma_0%29 "u_t \sim N(0,\Sigma_0)").

The main functions for CIVAR models are

- CIVARData
- CIVARest
- irf_CIVAR_CB
- CIVARTest
- CCIVARData
- CCIVARest

The following codes demonstrate how CIVARData can be used to generate
data from a cointegrated VAR process with a variety of specifications.

``` r
res_d = CIVARData(n=2,p=2,T=100,type="const")
res_d = CIVARData(n=3,p=2,T=100,Co=c(1:2)*0,type="none")
res_d = CIVARData(n=4,p=1,T=100,Co=c(1:2)*NA,type="const",crk=2)
STAT(res_d$B)
```

The output of CIVARData is a CIVAR object containing the generated data,
the used parameters as well as the residuals. The components correspond
to the model variables and parameters as follows.

res_d\$Y:
![(Y_1,Y_2,...,Y_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28Y_1%2CY_2%2C...%2CY_T%29%27 "(Y_1,Y_2,...,Y_T)'")

res_d\$Co:
![C_o](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_o "C_o")

res_d\$B:
![B_1,B_2,...B_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B_1%2CB_2%2C...B_p "B_1,B_2,...B_p")

res_d\$Sigmao:
![\Sigma_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma_0 "\Sigma_0")

res_d\$resid:
![(u_1,u_2,...,u_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28u_1%2Cu_2%2C...%2Cu_T%29%27 "(u_1,u_2,...,u_T)'")

res_d\$crk:
![crk](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;crk "crk")

res_d\$p:
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")

CIVARData differs from VARData in that the data generating process
contains unit-roots. The STAT function calculates the roots of
characteristic polynomial in lags or equivalently the eigenvalues of the
companion matrix of the CIVAR model.

The default number of unit-root in CIVARData is one. The option of r_np
matrix can be used to specify the number of unit roots.

CIVARest is used to estimate the parameters of the generated CIVAR data.

``` r
p = 3
n = 4
r_np = matrix(c(1,2,1.5,1.5,2.5,2.5,2,-1.5,2,-4,1.9,-2.1),4,3)
### the r_np matrix shows there is one unit root in the system. Hence the cointegration rank: crk = 3
res_d  = CIVARData(n=4,p=3,T=500,r_np=r_np,Co=matrix(0,n,1) ,type="none",crk=3)
res_e  = CIVARest(res=res_d)
res_e$Summary


B0 = res_d$B
res_d  = CIVARData(n=4,p=3,T=5000,B = B0,Co=matrix(0,n,1) ,type="none",crk=3)
res_e  = CIVARest(res=res_d)
res_e$Summary
sum(abs(res_e$B-res_d$B))
```

sum(abs(res_e\$B-res_d\$B)) can be used to gauge the difference between
the estimated parameters and the true parameters used in the data
generating process. By changing sample size T, we can gauge the
convergence rate of the estimator. res_e\$Summary shows the Johansen
test results. For empirical applications, the following steps are
recommended.

- Use CIVARData to generate a CIVAR object.
- Replace the generated model data with the empirical data.
- Run a model selection and testing procedure to determine the
  specification of the model.
- Use CIVARData to generate a CIVAR object with the selected
  specification and replace the model data with empirical data.
- Estimate the selected model with the replaced data.
- Run the analysis.

``` r

data(Canada)
dim(Canada)
trend = 1:84
# step 1: generating a VAR object
plot(ts(as.matrix(Canada)))

# step 1 level VAR selection

VARselect(Canada)

res_d = CIVARData(n=4,p=2,T=84,type="const",crk=2)
res_d$Y = Canada
plot(ts(res_d$Y))
res_e = CIVARest(res=res_d)
res_e$Summary
### The JH test shows that there is only one cointegration relation

# step 2: Re-specify the DGP and replacing the model data with the empirical data and estimate the model parameters

res_d = CIVARData(n=4,p=2,T=84,type="const",crk=1)
res_d$Y = Canada
res_e = CIVARest(res=res_d)
res_e$Summary

res_e$Summary$BIC
res_e$Summary$AIC

# step 3 Model selection based on AIC or BIC
res_d2 = CIVARData(n=4,p=2,T=84,Co=matrix(c(1,1,1,1),4,1)*0,type="none",crk=1)
res_d2$Y = as.matrix(Canada)
res_e2 = CIVARest(res=res_d2)
res_e2$Summary$BIC
res_e2$Summary$AIC

##### select the specification with the intercept

# step 4 IRF

IRF_CB = irf_CIVAR_CB(res=res_e, nstep=20, comb=NA, irf = "gen1", runs = 20, conf = c(0.05, 0.95))
IRF_list <-IRF_graph(IRF_CB,Names=c("e","prod","rw","U"))
```

### Test of restrictions on ![\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha "\alpha") and ![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")

CIVARTest function calculates maximum likelihood tests on the
restrictions within the cointegration space i.e. restrictions on
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
and
![\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha "\alpha").
The restrictions are formulated as follows.

![\mbox{vec}(\alpha') = G \psi \hspace{2cm} \mbox{vec}(\beta) = H \phi + h](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmbox%7Bvec%7D%28%5Calpha%27%29%20%3D%20G%20%5Cpsi%20%5Chspace%7B2cm%7D%20%5Cmbox%7Bvec%7D%28%5Cbeta%29%20%3D%20H%20%5Cphi%20%2B%20h "\mbox{vec}(\alpha') = G \psi \hspace{2cm} \mbox{vec}(\beta) = H \phi + h")

``` r
res_d  = CIVARData(n=7,p=2,T=207,crk=2,type="const")
colnames(res_d$Y) =  c("w","p","U_l","r","yn","y","fsi")
res_e =  CIVARest(res=res_d)
res_e$Summary

n = 7;  crk = 2
G = diag(n*crk); G0 = G; psi = matrix(1,n*crk,1); psi0=psi;   ### No restrictions on alpha

H = diag(n*crk); H2 = H[,-c(1,8)]; h2 = matrix(0,n*crk,1);
h2[c(1,8),1] = c(1,1);  phi2 = matrix(1,12,1)                 ### normalization

# G0%*%psi0; H2%*%phi2+h2

ABtest = CIVARTest(res=res_d,H=H2,h=h2,phi=phi2,G=G0,Dxflag=0)
ABtest$LR
ABtest$betar
ABtest$alphar
```

The sample codes above do not impose any binding restrictions in
![\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha "\alpha")
and
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta").
The likelihood ratio value is zero.

### Impulse response functions to exogenous shocks

For CIVARX models irf_CIVAR_CB can be used to generate impulse response
functions to exogenous shocks. In the following example there are two
exogenous variables. Using the option irf=“irfX” we obtain impulse
response functions to exogenous shocks.

``` r

T = 100
X=matrix(rnorm(2*T),T,2)
#colnames(X) = c("ex1","ex2")
res_d = CIVARData(n=3,p=2,T=T,Co=matrix(c(2,2,2,1,2,3,3,2,1),3,3), type="exog1",X=X) ;   res_e = CIVARest(res=res_d);res_e$Summary
IRF_CB = irf_CIVAR_CB(res=res_e,nstep=20, comb=NA, irf = "irfX", runs = 100, conf = c(0.05, 0.95),Xshk=c(1:2))
IRF_list <-IRF_graph(IRF_CB, response=c(1,2,3),Names=c("Y1","Y2","Y3"),IName=c("X1","X2"), impulse=c(1,2),ncol=2)

res_e$Co
```

### Conditional cointegration processes

CCIVARData and CCIVARest are functions used to model conditional
cointegrated VAR processes.\[^1\] This conditional cointegrated VAR
process differs from CIVARX in that the conditioning variables are
weakly exogenous I(1) variables that enter the cointegration space,
while the exogenous variables in CIVARX are I(0) variables and hence
they do not enter the cointegration space.

\[^1\] I. Harbo, S. Johansen, B. Nielsen, and A. Rahbek (1998)
Asymptotic Inference on Cointegrating Rank in Partial Systems, Journal
of Business & Economic Statistics, 16(4) (1998), 388-399.

``` r
T = 200
res_d <- CCIVARData(n1=4,n2=3,crk=3,p=3,T=T,type="const")
res_e <- CCIVARest(res=res_d)
res_e$Summary
res_e$tst$beta
```

### Vector error correction models with I(1) and I(0) variables

MIxCIVARData function is used to generate CIVAR data with mixed I(1) and
I(0) variables. It is also used to estimate the mixed CIVAR models. The
following example codes generate a set of data from a mixed CIVAR
process with a total of 9 variables and 5 independent unit roots. The
first 7 columns in Y are I(1) variables and the the last 2 columns are
I(0) variables.

``` r
RR <- MIxCIVARData(n=9,p=2,T=209,r=5,k=2,type="const",Bo=NA,Y=NA,D=NA)
## DGP of I(1) and I(0) mixed VECM
plot(ts(RR$Y))
```

The graph indicates that among the 9 variables there two I(0) variables
and 7 I(1) variables. The mixed CIVAR model implies restrictions on the
cointegration space
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta"):

![\beta' =\left(
  \begin{array}{ccc}
    B^{21} & B^{22} & 0 \\
    0 & 0 & B^{33} \\
  \end{array}
\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta%27%20%3D%5Cleft%28%0A%20%20%5Cbegin%7Barray%7D%7Bccc%7D%0A%20%20%20%20B%5E%7B21%7D%20%26%20B%5E%7B22%7D%20%26%200%20%5C%5C%0A%20%20%20%200%20%26%200%20%26%20B%5E%7B33%7D%20%5C%5C%0A%20%20%5Cend%7Barray%7D%0A%5Cright%29 "\beta' =\left(
  \begin{array}{ccc}
    B^{21} & B^{22} & 0 \\
    0 & 0 & B^{33} \\
  \end{array}
\right)")

![B^{33}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B%5E%7B33%7D "B^{33}")
iS a (2 x 2) submatrix corresponding to the two I(0) variables .
![(B^{21},B^{22})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28B%5E%7B21%7D%2CB%5E%7B22%7D%29 "(B^{21},B^{22})")
is a (2 x 7) cointegration matrix of the 7 I(1) variables. These
restrictions on
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
can be tested using CIVARTest

``` r
#### testing the mixed VECM via testing the restrictions on beta
res_d <- CIVARData(n=9,p=2,T=209,type="const",crk=4)
res_e = CIVARest(res=res_d)
n = 9; crk = 4; k = 2; r = 5
CC  <- c(8,9,17,18)
GG  <- c(19:25,28:34)
G = diag(n*crk); psi=matrix(1,n*crk,1)
### this implies there is no restrictions on the adjustment coefficients alpha
H = diag(n*crk);               H2 = H[,-c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n),CC,GG)]
#### only normalization
h = matrix(0,n*crk,1);         h[c(seq(1,(n-r-k)*n,n),seq((n-r-k+1)*n,(n-r)*n,n)),1] <- 1
phi = matrix(1,ncol(H2),1)
## check consistency of the restrictions    G%*%psi; H2%*%phi + h
#
# G%*%psi; H2%*%phi + h
#
ABtest = CIVARTest(res=res_d,H=H2,h=h,phi=phi,G=G,Dxflag=0)
ABtest$betar
ABtest$LR
1-pchisq(ABtest$LR,14)   ### The Ho of last two are I(0) is rejected

### The null of mixed CIVAR is rejected for pure I(1) data

RR <- MIxCIVARData(n=9,p=2,T=209,r=5,k=2,type="const",Bo=NA,Y=NA,D=NA)
res_d$Y = RR$Y
ABtest = CIVARTest(res=res_d,H=H2,h=h,phi=phi,G=G,Dxflag=0)
ABtest$betar
ABtest$LR
1-pchisq(ABtest$LR,14)   ### The Ho of last two are I(0) can be rejected

RR$GABtest$LR
### The null of mixed CIVAR is not rejected for mixed data generated by MIxCIVARData
```

For CIVAR models, impulse response functions are often decomposed into
shocks with permanent effects and transitory effects. The option
irf=“PTdecomp”can be used to produce such decomposition.

``` r

res_d = CIVARData(n=4,p=2,T=84,type="none",crk=2)
res_e = CIVARest(res=res_d)
res_e$Summary


alpha = t(res_e$tst$estimation$coefficients[1:2,])
G = cbind(alpha_perpf(alpha),alpha)


IRF_CB = irf_CIVAR_CB(res=res_e, nstep=20, comb=NA, irf = "PTdecomp", runs = 20, G=G, conf = c(0.05, 0.95))
IRF_g = IRF_graph(IRF_CB)
```

## MRVAR/MRCIVAR Models

### MRVAR

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

Following sample codes demonstrate how MRVARData is used to generate
data from an MRVAR process.

``` r
p = matrix(c(2,1,0,0),2,2)
res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
max(res_d$Y)
colnames(res_d$Y) = c("R","P")
res_e = MRVARest(res=res_d)
res_e$Summary
```

The output of MRVARData is an MRVAR object containing the generated
data, the used parameters as well as the residuals. The components of an
MRVAR object correspond to the model variables and parameters as
follows.

res_d\$Y:
![(Y_1,Y_2,...,Y_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28Y_1%2CY_2%2C...%2CY_T%29%27 "(Y_1,Y_2,...,Y_T)'")

res_d\$Co:
![C_o^{(1)}, C_o^{(1)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_o%5E%7B%281%29%7D%2C%20C_o%5E%7B%281%29%7D "C_o^{(1)}, C_o^{(1)}")

res_d\$B_0:
![B^{(1)}\_1,B^{(1)}\_2,...,B^{(1)}\_p;B^{(2)}\_1,B^{(2)}\_2,...B^{(2)}\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B%5E%7B%281%29%7D_1%2CB%5E%7B%281%29%7D_2%2C...%2CB%5E%7B%281%29%7D_p%3BB%5E%7B%282%29%7D_1%2CB%5E%7B%282%29%7D_2%2C...B%5E%7B%282%29%7D_p "B^{(1)}_1,B^{(1)}_2,...,B^{(1)}_p;B^{(2)}_1,B^{(2)}_2,...B^{(2)}_p")

res_d\$Sigmao:
![\Sigma^{(1)}\_0,\Sigma^{(2)}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma%5E%7B%281%29%7D_0%2C%5CSigma%5E%7B%282%29%7D_0 "\Sigma^{(1)}_0,\Sigma^{(2)}_0")

res_d\$Uo:
![(u^{(1)}\_1,u^{(1)}\_2,...,u^{(1)}\_T;u^{(2)}\_1,u^{(2)}\_2,...,u^{(2)}\_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28u%5E%7B%281%29%7D_1%2Cu%5E%7B%281%29%7D_2%2C...%2Cu%5E%7B%281%29%7D_T%3Bu%5E%7B%282%29%7D_1%2Cu%5E%7B%282%29%7D_2%2C...%2Cu%5E%7B%282%29%7D_T%29%27 "(u^{(1)}_1,u^{(1)}_2,...,u^{(1)}_T;u^{(2)}_1,u^{(2)}_2,...,u^{(2)}_T)'")

res_d\$TH:
![\tau](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctau "\tau")

res_d\$p:
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")

![p_1= \left(  \begin{array}{cc}  2 & 0 \\  1 & 0 \\  \end{array} \right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_1%3D%20%5Cleft%28%20%20%5Cbegin%7Barray%7D%7Bcc%7D%20%202%20%26%200%20%5C%5C%20%201%20%26%200%20%5C%5C%20%20%5Cend%7Barray%7D%20%5Cright%29 "p_1= \left(  \begin{array}{cc}  2 & 0 \\  1 & 0 \\  \end{array} \right)")
specifies that the lag for the first regime is 2 and there is no
exogenous variables, while for the second regime the lag length is 1 and
there is no exogenous variables too. The function MRVARest is used to
estimate the parameters of the generated data.

``` r
p = matrix(c(2,1,0,0),2,2)
res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="none")
max(res_d$Y)
colnames(res_d$Y) = c("R","P")
res_e = MRVARest(res=res_d)


  
  ### generalized impulse response function with regime migrations
  #RF4 = girf_MRVAR_RM_CB(res_e, shock=c(1,1), R=100, nstep=40, Omega_hist = NA, resid_method = "parametric", conf_level=c(0.05,0.95), N=100)
  #x11()
  #IRF_list <-IRF_graph(RF4,Names =c("R","P") )
  
  ### regime specific impulse response functions without any regime migration
  IRF_CB = irf_MRVAR_NM_CB(res_e,nstep=40,comb=NA,irf="gen",runs=200,conf=c(0.05,0.90))
  #x11()
  IRF_list1 <-IRF_graph(IRF_CB[,,,,1])      
  #x11()
  IRF_list2 <-IRF_graph(IRF_CB[,,,,2])
```

girf_MRVAR_RM_CB is used to produce the generalized impulse response
function with regime migrations.[^1] irf_MRVAR_NM_CB is used to generate
the within-regime impulse response functions [^2]

The following sample codes demonstrate how AIC and BIC criteria can be
used in the model specification.

``` r
Sigma = 1:(4*4*2)
dim(Sigma) = c(4,4,2)
Sigma[,,1] = diag(4)*2
Sigma[,,2] = diag(4)*2
p=matrix(0,2,2)
p[,1] = c(3,2)

res_d = MRVARData(n=4,p=p,T=800,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="none")
max(res_d$Y)
res_e = MRVARest(res=res_d)
res_e$Summary

# providing possible threshold values
TH_v = c(-0.5,0.0,0.5)

# providing the maximum lags for the two regimes
L_v = c(5,5)

Sel = MRVAR_Select(res=res_e,L_V=L_v,TH_V=TH_v)
MRVAR_Select_Summary(Sel)
```

MRVAR_Select calculates AIC and BIC for different lag lengths in
different regimes as well as for different choice of threshold values.
In the example above, the values of information criteria of the MRVAR
AIC and BIC are less than that of the one regime model ORAIC and ORBIC,
which indicates that the two-regime-model is a more proper model for the
data.

Option X in the function MRVARData is used to construct MRVAR models
with exogenous variables. The following examples show how it works.

``` r
#### MRVAR with exogenous variables X is (T x k x S ) array, allowing
#### different exogenous variables in different regimes. k is the number of exogenous variables.

X1=matrix(rnorm(2*300),300,2)
X2=matrix(rnorm(2*300),300,2)

XX = cbind(X1,X2); dim(XX) = c(300,2,2)

p = matrix(c(2,1,2,2),2,2)
res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="exog1",X=XX)
#max(res_d$Y)
res_e = MRVARest(res=res_d)
#res_e$Summary

#### If the number of exogenous variables are different in two regimes, the excess columns 
#### can be filled with zeros.

#X=matrix(rnorm(2*300),300,2)
#XX = cbind(X,X); dim(XX) = c(300,2,2);XX[,2,2]=0

#p = matrix(c(2,1,2,2),2,2)
#res_d = MRVARData(n=2,p=p,T=300,S=2,SESVI=1,type="exog1",X=XX)
#max(res_d$Y)
#res_e = MRVARest(res=res_d)
#res_e$Summary
```

For empirical applications it is recommended to follow the following
steps:

- Use MRVARData to generate an MRVAR object.
- Replace the generated model data with the empirical data.
- Estimate model parameters using MRVARest
- Specify the the threshold value and the lag length using MRVAR_Select
  if necessary.
- Re-estimate the selected model
- Run analysis.

### MRCIVAR

A multi regime cointegrated VAR model is defined as follow.

![Y_t = \Big(C_o^{(1)} + \sum\_{j=1}^{p}B\_{j}^{(1)} Y\_{t-j}+u_t^{(1)}\Big){\bf 1}\_{\[f\_{t-d} \le \tau\]}+\Big(C_o^{(2)} + \sum\_{j=1}^{p}B\_{j}^{(2)} Y\_{t-j}+u_t^{(2)}\Big){\bf 1}\_{\[f\_{t-d} \> \tau\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_t%20%3D%20%5CBig%28C_o%5E%7B%281%29%7D%20%2B%20%5Csum_%7Bj%3D1%7D%5E%7Bp%7DB_%7Bj%7D%5E%7B%281%29%7D%20Y_%7Bt-j%7D%2Bu_t%5E%7B%281%29%7D%5CBig%29%7B%5Cbf%201%7D_%7B%5Bf_%7Bt-d%7D%20%5Cle%20%5Ctau%5D%7D%2B%5CBig%28C_o%5E%7B%282%29%7D%20%2B%20%5Csum_%7Bj%3D1%7D%5E%7Bp%7DB_%7Bj%7D%5E%7B%282%29%7D%20Y_%7Bt-j%7D%2Bu_t%5E%7B%282%29%7D%5CBig%29%7B%5Cbf%201%7D_%7B%5Bf_%7Bt-d%7D%20%3E%20%5Ctau%5D%7D "Y_t = \Big(C_o^{(1)} + \sum_{j=1}^{p}B_{j}^{(1)} Y_{t-j}+u_t^{(1)}\Big){\bf 1}_{[f_{t-d} \le \tau]}+\Big(C_o^{(2)} + \sum_{j=1}^{p}B_{j}^{(2)} Y_{t-j}+u_t^{(2)}\Big){\bf 1}_{[f_{t-d} > \tau]}")

where
![u^{(k)}\_t \sim N(0,\Sigma_0^{(k)})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;u%5E%7B%28k%29%7D_t%20%5Csim%20N%280%2C%5CSigma_0%5E%7B%28k%29%7D%29 "u^{(k)}_t \sim N(0,\Sigma_0^{(k)})")
for
![k= 1,2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%3D%201%2C2 "k= 1,2").
Because of unit roots and cointegration relations the MRCIVAR has a
vector error correction form:

![\Delta Y_t = \Big(C_o^{(1)} + \alpha^{(1)}\beta Y\_{t-1} + \sum\_{j=1}^{p-1}B\_{j}^{\*(1)} \Delta Y\_{t-j}+u_t^{(1)}\Big){\bf 1}\_{\[f\_{t-d} \le \tau\]}+\Big(C_o^{(2)} + \alpha^{(2)}\beta Y\_{t-1} + \sum\_{j=1}^{p-1}B\_{j}^{\*(2)} \Delta Y\_{t-j}+u_t^{(2)}\Big){\bf 1}\_{\[f\_{t-d} \> \tau\]}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20Y_t%20%3D%20%5CBig%28C_o%5E%7B%281%29%7D%20%2B%20%5Calpha%5E%7B%281%29%7D%5Cbeta%20Y_%7Bt-1%7D%20%2B%20%5Csum_%7Bj%3D1%7D%5E%7Bp-1%7DB_%7Bj%7D%5E%7B%2A%281%29%7D%20%5CDelta%20Y_%7Bt-j%7D%2Bu_t%5E%7B%281%29%7D%5CBig%29%7B%5Cbf%201%7D_%7B%5Bf_%7Bt-d%7D%20%5Cle%20%5Ctau%5D%7D%2B%5CBig%28C_o%5E%7B%282%29%7D%20%2B%20%5Calpha%5E%7B%282%29%7D%5Cbeta%20Y_%7Bt-1%7D%20%2B%20%5Csum_%7Bj%3D1%7D%5E%7Bp-1%7DB_%7Bj%7D%5E%7B%2A%282%29%7D%20%5CDelta%20Y_%7Bt-j%7D%2Bu_t%5E%7B%282%29%7D%5CBig%29%7B%5Cbf%201%7D_%7B%5Bf_%7Bt-d%7D%20%3E%20%5Ctau%5D%7D "\Delta Y_t = \Big(C_o^{(1)} + \alpha^{(1)}\beta Y_{t-1} + \sum_{j=1}^{p-1}B_{j}^{*(1)} \Delta Y_{t-j}+u_t^{(1)}\Big){\bf 1}_{[f_{t-d} \le \tau]}+\Big(C_o^{(2)} + \alpha^{(2)}\beta Y_{t-1} + \sum_{j=1}^{p-1}B_{j}^{*(2)} \Delta Y_{t-j}+u_t^{(2)}\Big){\bf 1}_{[f_{t-d} > \tau]}")

where
![\alpha^{(1)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha%5E%7B%281%29%7D "\alpha^{(1)}"),
![\alpha^{(2)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha%5E%7B%282%29%7D "\alpha^{(2)}")
and
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
are
![(n\times ckr)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28n%5Ctimes%20ckr%29 "(n\times ckr)")
matrices. Usually the cointegration relation
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
is assumed to be identical in different regimes while the adjustment
speed
![\alpha^{(k)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha%5E%7B%28k%29%7D "\alpha^{(k)}")
can be regime-dependent.

The main functions for this class of models are

- MRCIVARDatam, MRCIVARestm1
- ABC_MRCIVARestm
- MRCIVAR_Selectm
- irf_MRCIVAR_CB
- girf_MRCIVAR_RM_CB

The following sample codes demonstrate how MRCIVARDatam can be used to
generate data from a two-regime cointegrated VAR process with a variety
of specifications and how MRCIVARestm is used to estimate the
parameters. MRCIVAR_Selectm calculates AIC or BIC criteria for a range
of parameter settings.

``` r
Sigma = 1:(4*4*2)
dim(Sigma) = c(4,4,2)
Sigma[,,1] = diag(4)
Sigma[,,2] = diag(4)
p=matrix(0,2,2)
p[,1] = c(3,2)

res_d = MRCIVARDatam(n=4,p=p,T=261,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=1)
res_e = MRCIVARestm1(res=res_d)
res_e$Summary
```

The output of MRCIVARData is an MRCIVAR object containing generated data
and parameters as well as the residuals. The components of the generated
MRCIVAR object correspond to the model variables and parameters as
follows.

res_d\$Y:
![(Y_1,Y_2,...,Y_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28Y_1%2CY_2%2C...%2CY_T%29%27 "(Y_1,Y_2,...,Y_T)'")

res_d\$Co:
![C_o^{(1)}, C_o^{(2)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_o%5E%7B%281%29%7D%2C%20C_o%5E%7B%282%29%7D "C_o^{(1)}, C_o^{(2)}")

res_d\$B_0:
![B^{(1)}\_1,B^{(1)}\_2,...,B^{(1)}\_p;B^{(2)}\_1,B^{(2)}\_2,...B^{(2)}\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B%5E%7B%281%29%7D_1%2CB%5E%7B%281%29%7D_2%2C...%2CB%5E%7B%281%29%7D_p%3BB%5E%7B%282%29%7D_1%2CB%5E%7B%282%29%7D_2%2C...B%5E%7B%282%29%7D_p "B^{(1)}_1,B^{(1)}_2,...,B^{(1)}_p;B^{(2)}_1,B^{(2)}_2,...B^{(2)}_p")

res_d\$Sigmao:
![\Sigma^{(1)}\_0,\Sigma^{(2)}\_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma%5E%7B%281%29%7D_0%2C%5CSigma%5E%7B%282%29%7D_0 "\Sigma^{(1)}_0,\Sigma^{(2)}_0")

res_d\$Uo:
![(u^{(1)}\_1,u^{(1)}\_2,...,u^{(1)}\_T;u^{(2)}\_1,u^{(2)}\_2,...,u^{(2)}\_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28u%5E%7B%281%29%7D_1%2Cu%5E%7B%281%29%7D_2%2C...%2Cu%5E%7B%281%29%7D_T%3Bu%5E%7B%282%29%7D_1%2Cu%5E%7B%282%29%7D_2%2C...%2Cu%5E%7B%282%29%7D_T%29%27 "(u^{(1)}_1,u^{(1)}_2,...,u^{(1)}_T;u^{(2)}_1,u^{(2)}_2,...,u^{(2)}_T)'")

res_d\$TH:
![\tau](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctau "\tau")

res_d\$p:
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")

res_d\$crk:
![crk](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;crk "crk")

![p= \left(  \begin{array}{cc}  3 & 0 \\  2 & 0 \\  \end{array} \right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3D%20%5Cleft%28%20%20%5Cbegin%7Barray%7D%7Bcc%7D%20%203%20%26%200%20%5C%5C%20%202%20%26%200%20%5C%5C%20%20%5Cend%7Barray%7D%20%5Cright%29 "p= \left(  \begin{array}{cc}  3 & 0 \\  2 & 0 \\  \end{array} \right)")
specifies that the first regime has a lag length of 3 and no exogenous
variables, while the second regime has a lag length of 2 and no
exogenous variables. MRCIVARestm produces an MRCIVAR object containing
estimated results and test statistics.

In the preceding example, r=1, specifies one unit root and thus three
cointegration relations. The regression format of the VECM is contained
in res_e\$Summary: CI11, CI12, and CI13 are three cointegrated error
terms of regime 1 and CI21,CI22,CI23 are the same cointegrated error
terms of regime 2. The vector error correction form has two lags because
regime 1 has three lags in level VAR. dy1.1, dy2.1, dy3.1, dy4.1, dy1.2,
dy2.2, dy3.2, and dy4.2 denote
![\Delta y\_{1,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B1%2Ct-1%7D "\Delta y_{1,t-1}"),
![\Delta y\_{2,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B2%2Ct-1%7D "\Delta y_{2,t-1}"),
![\Delta y\_{3,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B3%2Ct-1%7D "\Delta y_{3,t-1}"),
![\Delta y\_{4,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B4%2Ct-1%7D "\Delta y_{4,t-1}"),
![\Delta y\_{1,t-2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B1%2Ct-2%7D "\Delta y_{1,t-2}"),
![\Delta y\_{2,t-2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B2%2Ct-2%7D "\Delta y_{2,t-2}"),
![\Delta y\_{3,t-2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B3%2Ct-2%7D "\Delta y_{3,t-2}"),
and
![\Delta y\_{4,t-2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B4%2Ct-2%7D "\Delta y_{4,t-2}").
They are the regime 1 regressors. In the second block, dy1.1, dy2.1,
dy3.1, and dy4.1 represent the regime 2 regressors
![\Delta y\_{1,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B1%2Ct-1%7D "\Delta y_{1,t-1}"),
![\Delta y\_{2,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B2%2Ct-1%7D "\Delta y_{2,t-1}"),
![\Delta y\_{3,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B3%2Ct-1%7D "\Delta y_{3,t-1}"),
and
![\Delta y\_{4,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7B4%2Ct-1%7D "\Delta y_{4,t-1}").
Because regime 2 has two lags in level VAR, the error correction form
has only one lag. The result of the Johansen test of the cointegration
rank is also displayed in res_e\$Summary.

B2CIB can be used to calculate the cointegration vectors and the
adjustment vectors from the cointegrated VAR coefficients.

``` r
## calculaion of alpha and beta in the two regimes
B2CIB(res_d$Bo[,,,1])[[2]]
B2CIB(res_d$Bo[,,,2])[[2]]
B2CIB(res_d$Bo[,,,1])[[3]]
B2CIB(res_d$Bo[,,,2])[[3]]

## check the roots of the characteristic polynomials of the two regimes
STAT(res_d$Bo[,,,1])
STAT(res_d$Bo[,,,2])
```

For empirical applications it is recommended to follow the steps:

- Use MRCIVARDatam to generate a MRCIVAR object.
- Replace the generated model data with the empirical data.
- Estimate the model
- Run analysis.

MRCIVAR_Select can be used to assist the model specification based on
the information criteria AIC and BIC.

``` r
Sigma = 1:(4*4*2)
dim(Sigma) = c(4,4,2)
Sigma[,,1] = diag(4)
Sigma[,,2] = diag(4)
p=matrix(0,2,2)
p[,1] = c(3,2)

res_d = MRCIVARDatam(n=4,p=p,T=260,S=2,SESVI=1,TH=0,Sigmao=Sigma,type="const",r=1)
res_e = MRCIVARestm1(res=res_d)


TH_v = c(0,0.1)
L_v = c(6,6)
plot(ts(res_d$Y))

Selm = MRCIVAR_Selectm(res=res_e,L_V=L_v,TH_V=TH_v)
Selm[which.min(Selm[,5]),]
Selm[which.min(Selm[,4]),]
MRVAR_Select_Summary(Selm)
```

The sample codes above demonstrate how model selection criteria can be
used to specify the lag length and threshold values in the two regimes.
The two-regime AIC and BIC values are less than those of the one-regime
models ORAIC and ORBIC. This suggests that the two-regime model is a
better fit for the data. The function MRCIVAR_Select has two arguments:
L_V and L_TH, which specify the maximum lag length of the two regimes
and the threshold values, over which the information criteria AIC and
BIC will be calculated. The MRVAR_Select_Summary function is used to
output the minimum of the information criteria.

## GVAR/CIGVAR Models

### GVAR

A global vector autoregressive (GVAR) model in its country model form
is:

![y\_{i,t} =  c\_{i}+\sum\_{l=1}^{p\_{i}}B\_{i,l}y\_{i,t-l}+\sum\_{l=1}^{q\_{i}}A\_{i,l}y\_{i,t-l}^{\*}+\epsilon\_{i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%20%3D%20%20c_%7Bi%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bp_%7Bi%7D%7DB_%7Bi%2Cl%7Dy_%7Bi%2Ct-l%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bq_%7Bi%7D%7DA_%7Bi%2Cl%7Dy_%7Bi%2Ct-l%7D%5E%7B%2A%7D%2B%5Cepsilon_%7Bi%7D "y_{i,t} =  c_{i}+\sum_{l=1}^{p_{i}}B_{i,l}y_{i,t-l}+\sum_{l=1}^{q_{i}}A_{i,l}y_{i,t-l}^{*}+\epsilon_{i}")

where
![y\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D "y_{i,t}")
is an
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")-vector
of domestic variables with
![y\_{i,t}=(y\_{i,1,t},y\_{i,2,t},...,y\_{i,m,t})'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%3D%28y_%7Bi%2C1%2Ct%7D%2Cy_%7Bi%2C2%2Ct%7D%2C...%2Cy_%7Bi%2Cm%2Ct%7D%29%27 "y_{i,t}=(y_{i,1,t},y_{i,2,t},...,y_{i,m,t})'")
and
![y\_{i,t}^{\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%5E%7B%2A%7D "y_{i,t}^{*}")
is an
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")-vector
of the corresponding foreign variables with
![y\_{i,t}^{\*}=(y\_{i,1,t}^{\*},y\_{i,2,t}^{\*},...,y\_{i,m,t}^{\*})'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%5E%7B%2A%7D%3D%28y_%7Bi%2C1%2Ct%7D%5E%7B%2A%7D%2Cy_%7Bi%2C2%2Ct%7D%5E%7B%2A%7D%2C...%2Cy_%7Bi%2Cm%2Ct%7D%5E%7B%2A%7D%29%27 "y_{i,t}^{*}=(y_{i,1,t}^{*},y_{i,2,t}^{*},...,y_{i,m,t}^{*})'"),
![i=1,2,...,n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%3D1%2C2%2C...%2Cn "i=1,2,...,n").
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
is the number of countries. Stacking all
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
country models together we have GVAR in its global form:

![Y\_{t}=C+\sum\_{l=1}^{P}G\_{l}Y\_{t-l}+\epsilon_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7Bt%7D%3DC%2B%5Csum_%7Bl%3D1%7D%5E%7BP%7DG_%7Bl%7DY_%7Bt-l%7D%2B%5Cepsilon_t "Y_{t}=C+\sum_{l=1}^{P}G_{l}Y_{t-l}+\epsilon_t")

with
![\epsilon_t \sim N(0,\Sigma)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_t%20%5Csim%20N%280%2C%5CSigma%29 "\epsilon_t \sim N(0,\Sigma)"),

![Y\_{t}=(y\_{1,t}',y\_{2,t}',...,y\_{n,t}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7Bt%7D%3D%28y_%7B1%2Ct%7D%27%2Cy_%7B2%2Ct%7D%27%2C...%2Cy_%7Bn%2Ct%7D%27%29%27 "Y_{t}=(y_{1,t}',y_{2,t}',...,y_{n,t}')'")

![\epsilon_t=(\epsilon\_{1,t}',\epsilon\_{2,t}',...,\epsilon\_{n,t}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_t%3D%28%5Cepsilon_%7B1%2Ct%7D%27%2C%5Cepsilon_%7B2%2Ct%7D%27%2C...%2C%5Cepsilon_%7Bn%2Ct%7D%27%29%27 "\epsilon_t=(\epsilon_{1,t}',\epsilon_{2,t}',...,\epsilon_{n,t}')'")

![C=(c\_{1}',c\_{2}',...,c\_{n}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C%3D%28c_%7B1%7D%27%2Cc_%7B2%7D%27%2C...%2Cc_%7Bn%7D%27%29%27 "C=(c_{1}',c_{2}',...,c_{n}')'")

![G\_{l}  =  \left(\begin{array}{cccc}
B\_{1,l} & 0 & \dots & 0\\
0 & B\_{2,l} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & B\_{n,l}
\end{array}\right)+\left(\begin{array}{cccc}
A\_{1,l} & 0 & \dots & 0\\
0 & A\_{2,l} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A\_{n,l}
\end{array}\right)\left(\begin{array}{cccc}
0I_m & w\_{1,2}I_m & \dots & w\_{1,n}I_m\\
w\_{2,1}I_m & 0I_m & \ddots & w\_{2,n}I_m\\
\vdots & \ddots & \ddots & \vdots\\
w\_{n,1}I_m & \dots & w\_{n,n-1}I_m & 0I_m
\end{array}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_%7Bl%7D%20%20%3D%20%20%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0AB_%7B1%2Cl%7D%20%26%200%20%26%20%5Cdots%20%26%200%5C%5C%0A0%20%26%20B_%7B2%2Cl%7D%20%26%20%20%26%20%5Cvdots%5C%5C%0A%5Cvdots%20%26%20%20%26%20%5Cddots%20%26%200%5C%5C%0A0%20%26%20%5Cdots%20%26%200%20%26%20B_%7Bn%2Cl%7D%0A%5Cend%7Barray%7D%5Cright%29%2B%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0AA_%7B1%2Cl%7D%20%26%200%20%26%20%5Cdots%20%26%200%5C%5C%0A0%20%26%20A_%7B2%2Cl%7D%20%26%20%20%26%20%5Cvdots%5C%5C%0A%5Cvdots%20%26%20%20%26%20%5Cddots%20%26%200%5C%5C%0A0%20%26%20%5Cdots%20%26%200%20%26%20A_%7Bn%2Cl%7D%0A%5Cend%7Barray%7D%5Cright%29%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0A0I_m%20%26%20w_%7B1%2C2%7DI_m%20%26%20%5Cdots%20%26%20w_%7B1%2Cn%7DI_m%5C%5C%0Aw_%7B2%2C1%7DI_m%20%26%200I_m%20%26%20%5Cddots%20%26%20w_%7B2%2Cn%7DI_m%5C%5C%0A%5Cvdots%20%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5Cvdots%5C%5C%0Aw_%7Bn%2C1%7DI_m%20%26%20%5Cdots%20%26%20w_%7Bn%2Cn-1%7DI_m%20%26%200I_m%0A%5Cend%7Barray%7D%5Cright%29 "G_{l}  =  \left(\begin{array}{cccc}
B_{1,l} & 0 & \dots & 0\\
0 & B_{2,l} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & B_{n,l}
\end{array}\right)+\left(\begin{array}{cccc}
A_{1,l} & 0 & \dots & 0\\
0 & A_{2,l} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A_{n,l}
\end{array}\right)\left(\begin{array}{cccc}
0I_m & w_{1,2}I_m & \dots & w_{1,n}I_m\\
w_{2,1}I_m & 0I_m & \ddots & w_{2,n}I_m\\
\vdots & \ddots & \ddots & \vdots\\
w_{n,1}I_m & \dots & w_{n,n-1}I_m & 0I_m
\end{array}\right)")

The main functions in the class of GVAR models are:

- GVARData
- GVARest
- GVAR_Select
- irf_GVAR_CB

The following sample codes demonstrate how GVARData is used to generate
data from GVAR processes with a variety of specifications.

``` r
##### A case of m=2 variables, n=5 countries, domestic lag = 2, foreign lag = 1, and no exogenous variables

n = 5; m = 2
p = (1:15)*0; dim(p) = c(5,3)
p[,1] = 2; p[,2]=1;  p[1,2] = 2
res_d = GVARData(m=2,n=5,p=p,T=100,type="const")
max(res_d$Y)
dim(res_d$Y)
res_e = GVARest(res = res_d)
res_e$Summary
```

The output of GVARData is a GVAR object containing generated data and
used parameters as well as the residuals. The components of a generated
MRVAR object correspond to the model variables and parameters as
follows.

res_d\$Y:
![(Y_1,Y_2,...,Y_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28Y_1%2CY_2%2C...%2CY_T%29%27 "(Y_1,Y_2,...,Y_T)'")

res_d\$Uo:
![(\epsilon_1,\epsilon_2,...,\epsilon_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Cepsilon_1%2C%5Cepsilon_2%2C...%2C%5Cepsilon_T%29%27 "(\epsilon_1,\epsilon_2,...,\epsilon_T)'")

res_d\$Co:
![(c_1',c_2',c_3',...,c_n')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28c_1%27%2Cc_2%27%2Cc_3%27%2C...%2Cc_n%27%29%27 "(c_1',c_2',c_3',...,c_n')'")

res_d\$Bo:
![B\_{1,1},...,B\_{1,p},B\_{2,1},...,B\_{2,p},...,B\_{n,1},...B\_{n,p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B_%7B1%2C1%7D%2C...%2CB_%7B1%2Cp%7D%2CB_%7B2%2C1%7D%2C...%2CB_%7B2%2Cp%7D%2C...%2CB_%7Bn%2C1%7D%2C...B_%7Bn%2Cp%7D "B_{1,1},...,B_{1,p},B_{2,1},...,B_{2,p},...,B_{n,1},...B_{n,p}")

res_d\$Ao:
![A\_{1,1},...,A\_{1,q},A\_{2,1},...,A\_{2,q},...,A\_{n,1},...A\_{n,q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A_%7B1%2C1%7D%2C...%2CA_%7B1%2Cq%7D%2CA_%7B2%2C1%7D%2C...%2CA_%7B2%2Cq%7D%2C...%2CA_%7Bn%2C1%7D%2C...A_%7Bn%2Cq%7D "A_{1,1},...,A_{1,q},A_{2,1},...,A_{2,q},...,A_{n,1},...A_{n,q}")

res_d\$C:
![C](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C "C")

res_d\$G:
![G_1,G_2,...,G_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_1%2CG_2%2C...%2CG_p "G_1,G_2,...,G_p")

res_d\$Sigmao:
![\Sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma "\Sigma")

res_d\$p:
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")

res_d\$W:

![W=\left(\begin{array}{cccc}
0 & w\_{1,2} & \dots & w\_{1,n}\\
w\_{2,1} & 0 & \ddots & w\_{2,n}\\
\vdots & \ddots & \ddots & \vdots\\
w\_{N,1} & \dots & w\_{n,n-1} & 0
\end{array}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%3D%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0A0%20%26%20w_%7B1%2C2%7D%20%26%20%5Cdots%20%26%20w_%7B1%2Cn%7D%5C%5C%0Aw_%7B2%2C1%7D%20%26%200%20%26%20%5Cddots%20%26%20w_%7B2%2Cn%7D%5C%5C%0A%5Cvdots%20%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5Cvdots%5C%5C%0Aw_%7BN%2C1%7D%20%26%20%5Cdots%20%26%20w_%7Bn%2Cn-1%7D%20%26%200%0A%5Cend%7Barray%7D%5Cright%29 "W=\left(\begin{array}{cccc}
0 & w_{1,2} & \dots & w_{1,n}\\
w_{2,1} & 0 & \ddots & w_{2,n}\\
\vdots & \ddots & \ddots & \vdots\\
w_{N,1} & \dots & w_{n,n-1} & 0
\end{array}\right)")

![p= \left(  \begin{array}{ccc}  2 & 2 & 0 \\  2 & 1 & 0 \\  2 & 1 & 0 \\  2 & 1 & 0 \\  2 & 1 & 0 \\  \end{array} \right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3D%20%5Cleft%28%20%20%5Cbegin%7Barray%7D%7Bccc%7D%20%202%20%26%202%20%26%200%20%5C%5C%20%202%20%26%201%20%26%200%20%5C%5C%20%202%20%26%201%20%26%200%20%5C%5C%20%202%20%26%201%20%26%200%20%5C%5C%20%202%20%26%201%20%26%200%20%5C%5C%20%20%5Cend%7Barray%7D%20%5Cright%29 "p= \left(  \begin{array}{ccc}  2 & 2 & 0 \\  2 & 1 & 0 \\  2 & 1 & 0 \\  2 & 1 & 0 \\  2 & 1 & 0 \\  \end{array} \right)")
specifies that for country 1 the lag length of the domestic variables is
2 and that of the foreign variables is also 2. For countries 2 to 5, the
lag lengths of domestic variables are 2 and those of the foreign
variables are 1. The output of GVARest, res_e, is a GVAR object
containing estimated results and test statistics. res_e\$Summary gives
the country model regression results and test statistics.

In the country regression output, res_e\$Summary of country 4: Y1.l1,
Y2.l1, Y1.l2, Y2.l2, FY1.1 FY2.1 correspond to the domestic variable and
the foreign variables:
![Y\_{4,1,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7B4%2C1%2Ct-1%7D "Y_{4,1,t-1}"),
![Y\_{4,2,t-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7B4%2C2%2Ct-1%7D "Y_{4,2,t-1}"),
![Y\_{4,1,t-2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7B4%2C1%2Ct-2%7D "Y_{4,1,t-2}"),
![Y\_{4,2,t-2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7B4%2C2%2Ct-2%7D "Y_{4,2,t-2}"),
![Y\_{4,1,t-1}^{\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7B4%2C1%2Ct-1%7D%5E%7B%2A%7D "Y_{4,1,t-1}^{*}"),
![Y\_{4,2,t-1}^{\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7B4%2C2%2Ct-1%7D%5E%7B%2A%7D "Y_{4,2,t-1}^{*}"),respectively.

Option X in GVARData can be used to include exogenous variables in GVAR
models. irf_GVAR_CB calculates impulse response functions of an
estimated GVAR model.

``` r
X1 = matrix(1,200,1)
X2 = matrix(rnorm(200),200,1)
X3 = matrix(rnorm(200),200,1)
X4 = matrix(rnorm(200),200,1)
X  = cbind(X1,X2,X3,X4)
dim(X) = c(200,1,4)
n = 4
p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1;   p[,3]=1; p[2,2]=2;
p

res_d = GVARData(m=2,n=4,p=p,T=200,type="exog0",X=X)
res_e = GVARest(res = res_d)
res_e$Summary

IRF_CB = irf_GVAR_CB(res_e,nstep=10,comb=NA,irf="gen",runs=200,conf=c(0.05,0.95))
dim(IRF_CB)
IRF_g = IRF_graph(IRF_CB,Names=NA,response=c(1,4),impulse=c(1,2,3,4), ncol=4)
```

The package allows a country-wise specification of the lag length.
Information criteria can be applied to select the model country by
country to obtain a proper specification of the lags in a GVAR model.
GVAR_Select calculates AIC and BIC to lags up to the upper bounds given
in L_V for country I.

``` r
#### specification through country-wise model selection
n = 4
p = (1:12)*0; dim(p) = c(4,3);p[,1] = 2; p[,2]=1; p[2:3,2] = 2
res_d = GVARData(m=2,n=4,p=p,T=4000,type="const")

I = 1
L_V = c(4,4)
res_d$p
GVARSelect = GVAR_Select(res=res_d,L_V=c(4,4),I=1)
GVARSelect[which.min(GVARSelect[,3]),]

I = 2
L_V = c(4,4)
res_d$p
GVARSelect = GVAR_Select(res=res_d,L_V=c(4,4),I=2)
GVARSelect[which.min(GVARSelect[,3]),]


I = 3
L_V = c(4,4)
res_d$p
GVARSelect = GVAR_Select(res=res_d,L_V=c(4,4),I=3)
GVARSelect[which.min(GVARSelect[,3]),]
```

For empirical applications it is recommended to follow the steps below:

- Use GVARData to generate a GVAR object.
- Replace the generated model data with the empirical data.
- Run a model selection procedure to determine the specification of the
  model equation by equation.
- Estimate the selected model.
- Run analysis.

### CIGVAR

A cointegrated global vector autoregressive (CIGVAR) model is a GVAR
with unit roots in every country. Each county model is assumed to have
the following vector error correction form:

![\Delta y\_{i,t} =  c\_{i}+\alpha_i\beta_i'y\_{i,t-1}  +\sum\_{l=1}^{p\_{i}-1}B^{\*}\_{i,l}\Delta y\_{i,t-l}+\sum\_{l=1}^{q\_{i}}A\_{i,l}^{\*}\Delta y\_{i,t-l}^{\*}+\epsilon\_{i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7Bi%2Ct%7D%20%3D%20%20c_%7Bi%7D%2B%5Calpha_i%5Cbeta_i%27y_%7Bi%2Ct-1%7D%20%20%2B%5Csum_%7Bl%3D1%7D%5E%7Bp_%7Bi%7D-1%7DB%5E%7B%2A%7D_%7Bi%2Cl%7D%5CDelta%20y_%7Bi%2Ct-l%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bq_%7Bi%7D%7DA_%7Bi%2Cl%7D%5E%7B%2A%7D%5CDelta%20y_%7Bi%2Ct-l%7D%5E%7B%2A%7D%2B%5Cepsilon_%7Bi%7D "\Delta y_{i,t} =  c_{i}+\alpha_i\beta_i'y_{i,t-1}  +\sum_{l=1}^{p_{i}-1}B^{*}_{i,l}\Delta y_{i,t-l}+\sum_{l=1}^{q_{i}}A_{i,l}^{*}\Delta y_{i,t-l}^{*}+\epsilon_{i}")

The country model in level of the variables:

![y\_{i,t} =  c\_{i}+\sum\_{l=1}^{p\_{i}}B\_{i,l}y\_{i,t-l}+\sum\_{l=1}^{q\_{i}}A\_{i,l}y\_{i,t-l}^{\*}+\epsilon\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%20%3D%20%20c_%7Bi%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bp_%7Bi%7D%7DB_%7Bi%2Cl%7Dy_%7Bi%2Ct-l%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bq_%7Bi%7D%7DA_%7Bi%2Cl%7Dy_%7Bi%2Ct-l%7D%5E%7B%2A%7D%2B%5Cepsilon_%7Bi%2Ct%7D "y_{i,t} =  c_{i}+\sum_{l=1}^{p_{i}}B_{i,l}y_{i,t-l}+\sum_{l=1}^{q_{i}}A_{i,l}y_{i,t-l}^{*}+\epsilon_{i,t}")

Stacking all
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
country models together we have the CIGVAR in its global form:

![Y\_{t}=C+\sum\_{l=1}^{P}G\_{l}Y\_{t-l}+\epsilon_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7Bt%7D%3DC%2B%5Csum_%7Bl%3D1%7D%5E%7BP%7DG_%7Bl%7DY_%7Bt-l%7D%2B%5Cepsilon_t "Y_{t}=C+\sum_{l=1}^{P}G_{l}Y_{t-l}+\epsilon_t")

with
![\epsilon_t \sim N(0,\Sigma)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_t%20%5Csim%20N%280%2C%5CSigma%29 "\epsilon_t \sim N(0,\Sigma)"),

![Y\_{t}=(y\_{1,t}',y\_{2,t}',...,y\_{n,t}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7Bt%7D%3D%28y_%7B1%2Ct%7D%27%2Cy_%7B2%2Ct%7D%27%2C...%2Cy_%7Bn%2Ct%7D%27%29%27 "Y_{t}=(y_{1,t}',y_{2,t}',...,y_{n,t}')'")

![\epsilon_t=(\epsilon\_{1,t}',\epsilon\_{2,t}',...,\epsilon\_{n,t}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_t%3D%28%5Cepsilon_%7B1%2Ct%7D%27%2C%5Cepsilon_%7B2%2Ct%7D%27%2C...%2C%5Cepsilon_%7Bn%2Ct%7D%27%29%27 "\epsilon_t=(\epsilon_{1,t}',\epsilon_{2,t}',...,\epsilon_{n,t}')'")

![C=(c\_{1}',c\_{2}',...,c\_{n}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C%3D%28c_%7B1%7D%27%2Cc_%7B2%7D%27%2C...%2Cc_%7Bn%7D%27%29%27 "C=(c_{1}',c_{2}',...,c_{n}')'")

![G\_{l}  =  \left(\begin{array}{cccc}
B\_{1,l} & 0 & \dots & 0\\
0 & B\_{2,l} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & B\_{n,l}
\end{array}\right)+\left(\begin{array}{cccc}
A\_{1,l} & 0 & \dots & 0\\
0 & A\_{2,l} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A\_{n,l}
\end{array}\right)\left(\begin{array}{cccc}
0I_m & w\_{1,2}I_m & \dots & w\_{1,n}I_m\\
w\_{2,1}I_m & 0I_m & \ddots & w\_{2,n}I_m\\
\vdots & \ddots & \ddots & \vdots\\
w\_{N,1}I_m & \dots & w\_{n,n-1}I_m & 0I_m
\end{array}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_%7Bl%7D%20%20%3D%20%20%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0AB_%7B1%2Cl%7D%20%26%200%20%26%20%5Cdots%20%26%200%5C%5C%0A0%20%26%20B_%7B2%2Cl%7D%20%26%20%20%26%20%5Cvdots%5C%5C%0A%5Cvdots%20%26%20%20%26%20%5Cddots%20%26%200%5C%5C%0A0%20%26%20%5Cdots%20%26%200%20%26%20B_%7Bn%2Cl%7D%0A%5Cend%7Barray%7D%5Cright%29%2B%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0AA_%7B1%2Cl%7D%20%26%200%20%26%20%5Cdots%20%26%200%5C%5C%0A0%20%26%20A_%7B2%2Cl%7D%20%26%20%20%26%20%5Cvdots%5C%5C%0A%5Cvdots%20%26%20%20%26%20%5Cddots%20%26%200%5C%5C%0A0%20%26%20%5Cdots%20%26%200%20%26%20A_%7Bn%2Cl%7D%0A%5Cend%7Barray%7D%5Cright%29%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0A0I_m%20%26%20w_%7B1%2C2%7DI_m%20%26%20%5Cdots%20%26%20w_%7B1%2Cn%7DI_m%5C%5C%0Aw_%7B2%2C1%7DI_m%20%26%200I_m%20%26%20%5Cddots%20%26%20w_%7B2%2Cn%7DI_m%5C%5C%0A%5Cvdots%20%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5Cvdots%5C%5C%0Aw_%7BN%2C1%7DI_m%20%26%20%5Cdots%20%26%20w_%7Bn%2Cn-1%7DI_m%20%26%200I_m%0A%5Cend%7Barray%7D%5Cright%29 "G_{l}  =  \left(\begin{array}{cccc}
B_{1,l} & 0 & \dots & 0\\
0 & B_{2,l} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & B_{n,l}
\end{array}\right)+\left(\begin{array}{cccc}
A_{1,l} & 0 & \dots & 0\\
0 & A_{2,l} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A_{n,l}
\end{array}\right)\left(\begin{array}{cccc}
0I_m & w_{1,2}I_m & \dots & w_{1,n}I_m\\
w_{2,1}I_m & 0I_m & \ddots & w_{2,n}I_m\\
\vdots & \ddots & \ddots & \vdots\\
w_{N,1}I_m & \dots & w_{n,n-1}I_m & 0I_m
\end{array}\right)")

The main functions for this class of cointegrated global VAR models are

- CIGVARData
- CIGVARest
- CIGVAR_Select
- irf_CIGVAR_CB

CIGVARData creates a CIGVAR object with provided or randomly generated
parameters and simulated data. The option Ncommfakt allows common
stochastic trends and country specific idiosyncratic stochastic trends.
CIGVARest estimates the parameters of a CIGVAR object.  
The following sample codes demonstrate the use of these five functions.
r_npo array specifies the number of unit roots in each country model.

``` r
n = 5
p = (1:15)*0; dim(p) = c(5,3)
p[,1] = 3; p[,2]=3;
res_d = CIGVARData(m=4,n=5,p=p,T=1000,type="const")
max(abs(res_d$Y))
STAT(res_d$G)
res_e$Summary
## as the default cointegration rank is m-1, we need adjust the specification of the cointegration rank in order to obtain a proper estimation,
res_d$crk = c(3,3,3,3,3)
res_e = CIGVARest(res_d)
res_e$Summary
```

The output of CIGVARData is a CIGVAR object containing generated data
and parameters as well as the residuals. The components of the generated
CIGVARData object correspond to the model variables and parameters as
follows.

res_d\$Y:
![Y = (Y_1,Y_2,...,Y_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%20%3D%20%28Y_1%2CY_2%2C...%2CY_T%29%27 "Y = (Y_1,Y_2,...,Y_T)'"),
![Y_t = (y\_{1,1,t},y\_{1,2,t},...,y\_{1,m,t},y\_{2,1,t},y\_{2,2,t},...y\_{2,m,t},...,y\_{n,m,t})'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_t%20%3D%20%28y_%7B1%2C1%2Ct%7D%2Cy_%7B1%2C2%2Ct%7D%2C...%2Cy_%7B1%2Cm%2Ct%7D%2Cy_%7B2%2C1%2Ct%7D%2Cy_%7B2%2C2%2Ct%7D%2C...y_%7B2%2Cm%2Ct%7D%2C...%2Cy_%7Bn%2Cm%2Ct%7D%29%27 "Y_t = (y_{1,1,t},y_{1,2,t},...,y_{1,m,t},y_{2,1,t},y_{2,2,t},...y_{2,m,t},...,y_{n,m,t})'")

res_d\$Uo:
![U_o = (U_1,U_2,...,U_T)'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;U_o%20%3D%20%28U_1%2CU_2%2C...%2CU_T%29%27 "U_o = (U_1,U_2,...,U_T)'"),
![U_t = (\epsilon\_{1,1,t},\epsilon\_{1,2,t},...,\epsilon\_{1,m,t},\epsilon\_{2,1,t},\epsilon\_{2,2,t},...\epsilon\_{2,m,t},...,\epsilon\_{n,m,t})'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;U_t%20%3D%20%28%5Cepsilon_%7B1%2C1%2Ct%7D%2C%5Cepsilon_%7B1%2C2%2Ct%7D%2C...%2C%5Cepsilon_%7B1%2Cm%2Ct%7D%2C%5Cepsilon_%7B2%2C1%2Ct%7D%2C%5Cepsilon_%7B2%2C2%2Ct%7D%2C...%5Cepsilon_%7B2%2Cm%2Ct%7D%2C...%2C%5Cepsilon_%7Bn%2Cm%2Ct%7D%29%27 "U_t = (\epsilon_{1,1,t},\epsilon_{1,2,t},...,\epsilon_{1,m,t},\epsilon_{2,1,t},\epsilon_{2,2,t},...\epsilon_{2,m,t},...,\epsilon_{n,m,t})'")

res_d\$Co:
![C_o = (c_1',c_2',c_3',...,c_n')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_o%20%3D%20%28c_1%27%2Cc_2%27%2Cc_3%27%2C...%2Cc_n%27%29%27 "C_o = (c_1',c_2',c_3',...,c_n')'")

res_d\$Bo:
![B\_{1,1},...,B\_{1,p},B\_{2,1},...,B\_{2,p},...,B\_{n,1},...B\_{n,p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B_%7B1%2C1%7D%2C...%2CB_%7B1%2Cp%7D%2CB_%7B2%2C1%7D%2C...%2CB_%7B2%2Cp%7D%2C...%2CB_%7Bn%2C1%7D%2C...B_%7Bn%2Cp%7D "B_{1,1},...,B_{1,p},B_{2,1},...,B_{2,p},...,B_{n,1},...B_{n,p}")

res_d\$Ao:
![A\_{1,1},...,A\_{1,q},A\_{2,1},...,A\_{2,q},...,A\_{n,1},...A\_{n,q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A_%7B1%2C1%7D%2C...%2CA_%7B1%2Cq%7D%2CA_%7B2%2C1%7D%2C...%2CA_%7B2%2Cq%7D%2C...%2CA_%7Bn%2C1%7D%2C...A_%7Bn%2Cq%7D "A_{1,1},...,A_{1,q},A_{2,1},...,A_{2,q},...,A_{n,1},...A_{n,q}")

res_d\$C:
![C=(c\_{1}',c\_{2}',...,c\_{n}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C%3D%28c_%7B1%7D%27%2Cc_%7B2%7D%27%2C...%2Cc_%7Bn%7D%27%29%27 "C=(c_{1}',c_{2}',...,c_{n}')'")

res_d\$G:
![G_1,G_2,...,G_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_1%2CG_2%2C...%2CG_p "G_1,G_2,...,G_p")

res_d\$Sigmao:
![\Sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma "\Sigma")

res_d\$p:
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")

res_d\$W:
![W=\left(\begin{array}{cccc} 0 & w\_{1,2} & \dots & w\_{1,N}\\ w\_{2,1} & 0 & \ddots & w\_{2,N}\\ \vdots & \ddots & \ddots & \vdots\\ w\_{N,1} & \dots & w\_{N,N-1} & 0 \end{array}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W%3D%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%200%20%26%20w_%7B1%2C2%7D%20%26%20%5Cdots%20%26%20w_%7B1%2CN%7D%5C%5C%20w_%7B2%2C1%7D%20%26%200%20%26%20%5Cddots%20%26%20w_%7B2%2CN%7D%5C%5C%20%5Cvdots%20%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5Cvdots%5C%5C%20w_%7BN%2C1%7D%20%26%20%5Cdots%20%26%20w_%7BN%2CN-1%7D%20%26%200%20%5Cend%7Barray%7D%5Cright%29 "W=\left(\begin{array}{cccc} 0 & w_{1,2} & \dots & w_{1,N}\\ w_{2,1} & 0 & \ddots & w_{2,N}\\ \vdots & \ddots & \ddots & \vdots\\ w_{N,1} & \dots & w_{N,N-1} & 0 \end{array}\right)")

STAT(res_d\$G) shows how many unit roots there are in the CIGVAR system.
The output of GVARest, res_e, is an estimated object of CIGVAR with
estimated parameters, residuals and a series of test statistics.
res_e\$Summary gives the regression formats of each country’s vector
error correction model. The Johansen test results show that the
cointegration rank is 3 for all 4 countries. Therefore, the
specification of the cointegration rank has to adjust to achieve a
proper estimation.

The option of X in CIGVARData can be used to include exogenous variable
in the CIGVAR model.

``` r

X1 = matrix(rnorm(400),200,2)
X2 = matrix(rnorm(400),200,2)
X3 = matrix(rnorm(400),200,2)
X4 = matrix(rnorm(400),200,2)
X5 = matrix(rnorm(400),200,2)

X  = cbind(X1,X2,X3,X4,X5)
dim(X) = c(200,2,5)

n = 5
p = (1:15)*0; dim(p) = c(5,3)
p[,1] = 3; p[,2]=3; p[,3] = 2

res_d = CIGVARData(m=3,n=5,p=p,T=200,type="exog0",X=X,DFYflag=1,Ncommfakt = 0)
max(res_d$Y)
plot(ts(res_d$Y[,1:10]))
res_e = CIGVARest(res=res_d)
res_e$Summary

res_d = CIGVARData(m=3,n=5,p=p,T=200,type="exog1",X=X,DFYflag=1,Ncommfakt = 0)
max(res_d$Y)
plot(ts(res_d$Y[,1:10]))
res_e = CIGVARest(res=res_d)
res_e$Summary


### model selection



L_V = c(4,4)

res = res_d
#str(res_e)
res=res_e

CIGVARSelecte = CIGVAR_Select(res=res_d,L_V=c(4,4))
CIGVARSelecte[which.min(CIGVARSelecte[,17]),]
CIGVARSelecte[which.min(CIGVARSelecte[,19]),]
pp = CIGVARSelecte[which.min(CIGVARSelecte[,19]),][1:(n*3)]
dim(pp)=c(n,3)
pp
```

The example above shows how CIGVAR_Select is used for the model
selection. The data generating process has a lag length of 3 for every
country’s domestic variables and the foreign variables as well. The
information criteria selection is very close to the true specification.

In global VAR modelling, often global trends play an important role.
CIGVARData function is used to generate data from a CIGVAR process with
global trends.

``` r

n = 5; m = 3
p = (1:15)*0; dim(p) = c(5,3)
p[,1] = 3; p[,2]=3;
r_npo = (1:(3 * 3 * 5))*0; dim(r_npo) = c(3,3,5)
r_npo[,,1] = matrix(c(1,1,3,2,3,3,3,3,3),3,3)   # country 1 contains two unit roots.
r_npo[,,2] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)   # country 2 contains one unit root.
r_npo[,,3] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)   # country 3 contains one unit root.
r_npo[,,4] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)
r_npo[,,5] = matrix(c(1,2,3,2,3,3,3,3,3),3,3)

res_d = CIGVARData(m=3,n=5,p=p,T=500,r_npo=r_npo,type="const",DFYflag=0,Ncommfakt=1)
max(res_d$Y)
STAT(res_d$G)
plot(ts(res_d$Y[,1:10]))
res_e = CIGVARest(res_d)
#res_e$Summary
## adjust the model specification according to the cointegration test results
res_d$crk = c(1,2,2,2,2)
res_e = CIGVARest(res_d)
#res_e$Summary
#### unit roots in the CIGVAR model
STAT(res_e$G)
STAT(res_e$GC)
```

For empirical applications it is recommended to follow the steps below:

- Use CIGVARData to generate a CIGVAR object.
- Replace the generated model data with the empirical data.
- Run a model selection procedure to determine the specification of the
  model if necessary.
- Estimate the selected model.
- Run the analysis.

## MRGVAR/MRCIGVAR Models

### MRGVAR

A multi regime global vector autoregressive model (MRGVAR) is defined as
follows.

![y\_{i,t}  =  c\_{i,S\_{i,t}}+\sum\_{l=1}^{p\_{i,S\_{i,t}}}B\_{i,l,S\_{i,t}}y\_{i,t-l}+\sum\_{l=1}^{q\_{i,S\_{i,t}}}A\_{i,l,S\_{i,t}}y\_{i,t-l}^{\*}+\epsilon\_{i,S\_{i,t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%20%20%3D%20%20c_%7Bi%2CS_%7Bi%2Ct%7D%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bp_%7Bi%2CS_%7Bi%2Ct%7D%7D%7DB_%7Bi%2Cl%2CS_%7Bi%2Ct%7D%7Dy_%7Bi%2Ct-l%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bq_%7Bi%2CS_%7Bi%2Ct%7D%7D%7DA_%7Bi%2Cl%2CS_%7Bi%2Ct%7D%7Dy_%7Bi%2Ct-l%7D%5E%7B%2A%7D%2B%5Cepsilon_%7Bi%2CS_%7Bi%2Ct%7D%7D "y_{i,t}  =  c_{i,S_{i,t}}+\sum_{l=1}^{p_{i,S_{i,t}}}B_{i,l,S_{i,t}}y_{i,t-l}+\sum_{l=1}^{q_{i,S_{i,t}}}A_{i,l,S_{i,t}}y_{i,t-l}^{*}+\epsilon_{i,S_{i,t}}")

where
![y\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D "y_{i,t}")
is an
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")-vector
of domestic variables with
![y\_{i,t}=(y\_{i,1,t},y\_{i,2,t},...,y\_{i,m,t})'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%3D%28y_%7Bi%2C1%2Ct%7D%2Cy_%7Bi%2C2%2Ct%7D%2C...%2Cy_%7Bi%2Cm%2Ct%7D%29%27 "y_{i,t}=(y_{i,1,t},y_{i,2,t},...,y_{i,m,t})'")
and
![y\_{i,t}^{\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%5E%7B%2A%7D "y_{i,t}^{*}")
is an
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")-vector
of the corresponding foreign variables with
![y\_{i,t}^{\*}=(y\_{i,1,t}^{\*},y\_{i,2,t}^{\*},...,y\_{i,m,t}^{\*})'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%5E%7B%2A%7D%3D%28y_%7Bi%2C1%2Ct%7D%5E%7B%2A%7D%2Cy_%7Bi%2C2%2Ct%7D%5E%7B%2A%7D%2C...%2Cy_%7Bi%2Cm%2Ct%7D%5E%7B%2A%7D%29%27 "y_{i,t}^{*}=(y_{i,1,t}^{*},y_{i,2,t}^{*},...,y_{i,m,t}^{*})'"),
![i=1,2,...,n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%3D1%2C2%2C...%2Cn "i=1,2,...,n").
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")
is the number of variables;
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
is the number of countries. Stacking all
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
country models together we have MRGVAR in its global form:

![Y\_{t}=C\_{S\_{t}}+\sum\_{l=1}^{p}G\_{l,S\_{t}}y\_{t-l}+\epsilon\_{S\_{t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7Bt%7D%3DC_%7BS_%7Bt%7D%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bp%7DG_%7Bl%2CS_%7Bt%7D%7Dy_%7Bt-l%7D%2B%5Cepsilon_%7BS_%7Bt%7D%7D "Y_{t}=C_{S_{t}}+\sum_{l=1}^{p}G_{l,S_{t}}y_{t-l}+\epsilon_{S_{t}}")

with
![\epsilon\_{S\_{t}}\sim N(0,\Sigma\_{S\_{t}})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7BS_%7Bt%7D%7D%5Csim%20N%280%2C%5CSigma_%7BS_%7Bt%7D%7D%29 "\epsilon_{S_{t}}\sim N(0,\Sigma_{S_{t}})").
![Y_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_t "Y_t")
is an (mn) vector.
![c\_{S_t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_%7BS_t%7D "c_{S_t}")
is a constant (mn) vector depending on the state of each country at time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t").

![Y\_{t}=(y\_{1,t}',y\_{2,t}',...,y\_{n,t}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7Bt%7D%3D%28y_%7B1%2Ct%7D%27%2Cy_%7B2%2Ct%7D%27%2C...%2Cy_%7Bn%2Ct%7D%27%29%27 "Y_{t}=(y_{1,t}',y_{2,t}',...,y_{n,t}')'")

![\epsilon\_{S\_{t}}=(\epsilon\_{1,S\_{1,t}}',\epsilon\_{2,S\_{2,t}}',...,\epsilon\_{n,S\_{n,t}}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7BS_%7Bt%7D%7D%3D%28%5Cepsilon_%7B1%2CS_%7B1%2Ct%7D%7D%27%2C%5Cepsilon_%7B2%2CS_%7B2%2Ct%7D%7D%27%2C...%2C%5Cepsilon_%7Bn%2CS_%7Bn%2Ct%7D%7D%27%29%27 "\epsilon_{S_{t}}=(\epsilon_{1,S_{1,t}}',\epsilon_{2,S_{2,t}}',...,\epsilon_{n,S_{n,t}}')'")

![C\_{S\_{t}}=(c\_{1,S\_{1,t}}',c\_{2,S\_{2,t}}',...,c\_{n,S\_{n,t}}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_%7BS_%7Bt%7D%7D%3D%28c_%7B1%2CS_%7B1%2Ct%7D%7D%27%2Cc_%7B2%2CS_%7B2%2Ct%7D%7D%27%2C...%2Cc_%7Bn%2CS_%7Bn%2Ct%7D%7D%27%29%27 "C_{S_{t}}=(c_{1,S_{1,t}}',c_{2,S_{2,t}}',...,c_{n,S_{n,t}}')'")

![G\_{l,S\_{t}} =  \left(\begin{array}{cccc}
A\_{1,l,S\_{1,t}} & 0 & \dots & 0\\
0 & A\_{2,l,S\_{2,t}} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A\_{n,l,S\_{N,t}}
\end{array}\right)+\left(\begin{array}{cccc}
A\_{1,l,S\_{1,t}}^{\*} & 0 & \dots & 0\\
0 & A\_{2,l,S\_{2,t}}^{\*} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A\_{n,l,S\_{n,t}}^{\*}
\end{array}\right)\left(\begin{array}{cccc}
0I_m & w\_{1,2}I_m & \dots & w\_{1,n}I_m\\
w\_{2,1}I_m & 0I_m & \ddots & w\_{2,n}I_m\\
\vdots & \ddots & \ddots & \vdots\\
w\_{n,1}I_m & \dots & w\_{n,n-1}I_m & 0I_m
\end{array}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_%7Bl%2CS_%7Bt%7D%7D%20%3D%20%20%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0AA_%7B1%2Cl%2CS_%7B1%2Ct%7D%7D%20%26%200%20%26%20%5Cdots%20%26%200%5C%5C%0A0%20%26%20A_%7B2%2Cl%2CS_%7B2%2Ct%7D%7D%20%26%20%20%26%20%5Cvdots%5C%5C%0A%5Cvdots%20%26%20%20%26%20%5Cddots%20%26%200%5C%5C%0A0%20%26%20%5Cdots%20%26%200%20%26%20A_%7Bn%2Cl%2CS_%7BN%2Ct%7D%7D%0A%5Cend%7Barray%7D%5Cright%29%2B%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0AA_%7B1%2Cl%2CS_%7B1%2Ct%7D%7D%5E%7B%2A%7D%20%26%200%20%26%20%5Cdots%20%26%200%5C%5C%0A0%20%26%20A_%7B2%2Cl%2CS_%7B2%2Ct%7D%7D%5E%7B%2A%7D%20%26%20%20%26%20%5Cvdots%5C%5C%0A%5Cvdots%20%26%20%20%26%20%5Cddots%20%26%200%5C%5C%0A0%20%26%20%5Cdots%20%26%200%20%26%20A_%7Bn%2Cl%2CS_%7Bn%2Ct%7D%7D%5E%7B%2A%7D%0A%5Cend%7Barray%7D%5Cright%29%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0A0I_m%20%26%20w_%7B1%2C2%7DI_m%20%26%20%5Cdots%20%26%20w_%7B1%2Cn%7DI_m%5C%5C%0Aw_%7B2%2C1%7DI_m%20%26%200I_m%20%26%20%5Cddots%20%26%20w_%7B2%2Cn%7DI_m%5C%5C%0A%5Cvdots%20%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5Cvdots%5C%5C%0Aw_%7Bn%2C1%7DI_m%20%26%20%5Cdots%20%26%20w_%7Bn%2Cn-1%7DI_m%20%26%200I_m%0A%5Cend%7Barray%7D%5Cright%29 "G_{l,S_{t}} =  \left(\begin{array}{cccc}
A_{1,l,S_{1,t}} & 0 & \dots & 0\\
0 & A_{2,l,S_{2,t}} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A_{n,l,S_{N,t}}
\end{array}\right)+\left(\begin{array}{cccc}
A_{1,l,S_{1,t}}^{*} & 0 & \dots & 0\\
0 & A_{2,l,S_{2,t}}^{*} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A_{n,l,S_{n,t}}^{*}
\end{array}\right)\left(\begin{array}{cccc}
0I_m & w_{1,2}I_m & \dots & w_{1,n}I_m\\
w_{2,1}I_m & 0I_m & \ddots & w_{2,n}I_m\\
\vdots & \ddots & \ddots & \vdots\\
w_{n,1}I_m & \dots & w_{n,n-1}I_m & 0I_m
\end{array}\right)")

The main functions for this class of models are

The main functions for this class of models are:

- MRGVARData
- MRGVARest
- MRGVAR_Select
- irf_MRGVAR_CB
- girf_MRGVAR_RM_CB

The following sample codes demonstrate how MRGVARData can be used as a
data generating process to generate data from MRGVAR processes with a
variety of specifications.

``` r

## case of n = 10, m = 2, S = 2 
p = rep(1,60); dim(p) = c(10,3,2)
p[1,1,2] = 2; p[2,2,2]=2; p[,3,] = 0
TH = c(1:10)*0; dim(TH) = c(1,10)
res_d <- MRGVARData(m=2,n=10,p=p,TH=TH,T=200,S=2,SESVI=seq(1,20,2),type="const")
max(res_d$Y)

### estimation of the MRGVAR model
colnames(res_d$Y) = c("Pa","Qa","Pb","Qb","Pc","Qc","Pd","Qd","Pe","Qe","Pf","Qf","Pg","Qg","Ph","Qh","Pi","Qi","Pj","Qj")
AName <- c("Pa","Qa","Pb","Qb","Pc","Qc","Pd","Qd","Pe","Qe","Pf","Qf","Pg","Qg","Ph","Qh","Pi","Qi","Pj","Qj")
res_e = MRGVARest(res=res_d)
res_e$Summary
STAT(res_e$Go[,,,1])   # check for the necessary condition of stationarity
STAT(res_e$Go[,,,2])

if (! max(Mod(STAT(res_e$Go[,,,1])),Mod(STAT(res_e$Go[,,,2]))) > 1.001)  {

IRF_CB  = irf_MRGVAR_CB(res=res_e,nstep=10,comb=NA,state=rep(1,10),irf="gen1",runs=20,conf=c(0.05,0.95))
IRF_g = IRF_graph(IRF_CB[[1]], impulse=c(1,3),response=c(2,3,4,6),ncol=2)    #IRF
#IRF_g = IRF_graph(IRF_CB[[2]])   # accumulated IRF

}
```

MRGVARest is used to estimate the parameters of an MRGVAR object with
provided data. irf_MRGVAR_CB and girf_MRGVAR_RM_CB generate regime
specific impulse response functions without regime migration and
generalized impulse response functions with regime migration,
respectively.

### MRCIGVAR

A multi regime cointegrated global vector autoregressive model is an
MRGVAR with unit roots in every country model. Each county model is
assumed to have a vector error correction form:

![\Delta y\_{i,t} =  c\_{i,S\_{i,t}}+\alpha\_{i,S\_{i,t}}\beta_i'y\_{i,t-1}  +\sum\_{l=1}^{p\_{i}-1}B^{\*}\_{i,l,S\_{i,t}}\Delta y\_{i,t-l}+\sum\_{l=1}^{q\_{i}}A^{\*}\_{i,l,S\_{i,t}}\Delta y\_{i,t-l}^{\*}+\epsilon\_{i,S\_{i,t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20y_%7Bi%2Ct%7D%20%3D%20%20c_%7Bi%2CS_%7Bi%2Ct%7D%7D%2B%5Calpha_%7Bi%2CS_%7Bi%2Ct%7D%7D%5Cbeta_i%27y_%7Bi%2Ct-1%7D%20%20%2B%5Csum_%7Bl%3D1%7D%5E%7Bp_%7Bi%7D-1%7DB%5E%7B%2A%7D_%7Bi%2Cl%2CS_%7Bi%2Ct%7D%7D%5CDelta%20y_%7Bi%2Ct-l%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bq_%7Bi%7D%7DA%5E%7B%2A%7D_%7Bi%2Cl%2CS_%7Bi%2Ct%7D%7D%5CDelta%20y_%7Bi%2Ct-l%7D%5E%7B%2A%7D%2B%5Cepsilon_%7Bi%2CS_%7Bi%2Ct%7D%7D "\Delta y_{i,t} =  c_{i,S_{i,t}}+\alpha_{i,S_{i,t}}\beta_i'y_{i,t-1}  +\sum_{l=1}^{p_{i}-1}B^{*}_{i,l,S_{i,t}}\Delta y_{i,t-l}+\sum_{l=1}^{q_{i}}A^{*}_{i,l,S_{i,t}}\Delta y_{i,t-l}^{*}+\epsilon_{i,S_{i,t}}")

where
![y\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D "y_{i,t}")
is an
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")-vector
of domestic variables with
![y\_{i,t}=(y\_{i,1,t},y\_{i,2,t},...,y\_{i,m,t})'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%3D%28y_%7Bi%2C1%2Ct%7D%2Cy_%7Bi%2C2%2Ct%7D%2C...%2Cy_%7Bi%2Cm%2Ct%7D%29%27 "y_{i,t}=(y_{i,1,t},y_{i,2,t},...,y_{i,m,t})'")
and
![y\_{i,t}^{\*}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%5E%7B%2A%7D "y_{i,t}^{*}")
is an
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")-vector
of the corresponding foreign variables with
![y\_{i,t}^{\*}=(y\_{i,1,t}^{\*},y\_{i,2,t}^{\*},...,y\_{i,m,t}^{\*})'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%5E%7B%2A%7D%3D%28y_%7Bi%2C1%2Ct%7D%5E%7B%2A%7D%2Cy_%7Bi%2C2%2Ct%7D%5E%7B%2A%7D%2C...%2Cy_%7Bi%2Cm%2Ct%7D%5E%7B%2A%7D%29%27 "y_{i,t}^{*}=(y_{i,1,t}^{*},y_{i,2,t}^{*},...,y_{i,m,t}^{*})'"),
![i=1,2,...,n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%3D1%2C2%2C...%2Cn "i=1,2,...,n").
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")
is the number of variables;
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
is the number of countries. The country model in level of the variables
is

![y\_{i,t}  =  c\_{i,S\_{i,t}}+\sum\_{l=1}^{p\_{i,S\_{i,t}}}B\_{i,l,S\_{i,t}}y\_{i,t-l}+\sum\_{l=1}^{q\_{i,S\_{i,t}}}A\_{i,l,S\_{i,t}}y\_{i,t-l}^{\*}+\epsilon\_{i,S\_{i,t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D%20%20%3D%20%20c_%7Bi%2CS_%7Bi%2Ct%7D%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bp_%7Bi%2CS_%7Bi%2Ct%7D%7D%7DB_%7Bi%2Cl%2CS_%7Bi%2Ct%7D%7Dy_%7Bi%2Ct-l%7D%2B%5Csum_%7Bl%3D1%7D%5E%7Bq_%7Bi%2CS_%7Bi%2Ct%7D%7D%7DA_%7Bi%2Cl%2CS_%7Bi%2Ct%7D%7Dy_%7Bi%2Ct-l%7D%5E%7B%2A%7D%2B%5Cepsilon_%7Bi%2CS_%7Bi%2Ct%7D%7D "y_{i,t}  =  c_{i,S_{i,t}}+\sum_{l=1}^{p_{i,S_{i,t}}}B_{i,l,S_{i,t}}y_{i,t-l}+\sum_{l=1}^{q_{i,S_{i,t}}}A_{i,l,S_{i,t}}y_{i,t-l}^{*}+\epsilon_{i,S_{i,t}}")

Stacking all
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
country models together we have the MRCIGVAR in its global form:

![Y\_{t}=c\_{S\_{t}}+\sum\_{l=1}^{P}G\_{l,S\_{t}}Y\_{t-l}+\epsilon\_{S\_{t}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7Bt%7D%3Dc_%7BS_%7Bt%7D%7D%2B%5Csum_%7Bl%3D1%7D%5E%7BP%7DG_%7Bl%2CS_%7Bt%7D%7DY_%7Bt-l%7D%2B%5Cepsilon_%7BS_%7Bt%7D%7D "Y_{t}=c_{S_{t}}+\sum_{l=1}^{P}G_{l,S_{t}}Y_{t-l}+\epsilon_{S_{t}}")

with
![\epsilon\_{S\_{t}}\sim N(0,\Sigma\_{S\_{t}})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7BS_%7Bt%7D%7D%5Csim%20N%280%2C%5CSigma_%7BS_%7Bt%7D%7D%29 "\epsilon_{S_{t}}\sim N(0,\Sigma_{S_{t}})")

![Y\_{t}=(y\_{1,t}',y\_{2,t}',...,y\_{n,t}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_%7Bt%7D%3D%28y_%7B1%2Ct%7D%27%2Cy_%7B2%2Ct%7D%27%2C...%2Cy_%7Bn%2Ct%7D%27%29%27 "Y_{t}=(y_{1,t}',y_{2,t}',...,y_{n,t}')'")

![\epsilon\_{S\_{t}}=(\epsilon\_{1,S\_{1,t}}',\epsilon\_{2,S\_{2,t}}',...,\epsilon\_{n,S\_{n,t}}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7BS_%7Bt%7D%7D%3D%28%5Cepsilon_%7B1%2CS_%7B1%2Ct%7D%7D%27%2C%5Cepsilon_%7B2%2CS_%7B2%2Ct%7D%7D%27%2C...%2C%5Cepsilon_%7Bn%2CS_%7Bn%2Ct%7D%7D%27%29%27 "\epsilon_{S_{t}}=(\epsilon_{1,S_{1,t}}',\epsilon_{2,S_{2,t}}',...,\epsilon_{n,S_{n,t}}')'")

![c\_{S\_{t}}=(c\_{1,S\_{1,t}}',c\_{2,S\_{2,t}}',...,c\_{n,S\_{n,t}}')'](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_%7BS_%7Bt%7D%7D%3D%28c_%7B1%2CS_%7B1%2Ct%7D%7D%27%2Cc_%7B2%2CS_%7B2%2Ct%7D%7D%27%2C...%2Cc_%7Bn%2CS_%7Bn%2Ct%7D%7D%27%29%27 "c_{S_{t}}=(c_{1,S_{1,t}}',c_{2,S_{2,t}}',...,c_{n,S_{n,t}}')'")

![G\_{l,S\_{t}} =  \left(\begin{array}{cccc}
B\_{1,l,S\_{1,t}} & 0 & \dots & 0\\
0 & B\_{2,l,S\_{2,t}} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & B\_{n,l,S\_{n,t}}
\end{array}\right)+\left(\begin{array}{cccc}
A\_{1,l,S\_{1,t}}^{\*} & 0 & \dots & 0\\
0 & A\_{2,l,S\_{2,t}}^{\*} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A\_{n,l,S\_{n,t}}^{\*}
\end{array}\right)\left(\begin{array}{cccc}
0I_m & w\_{1,2}I_m & \dots & w\_{1,n}I_m\\
w\_{2,1}I_m & 0I_m & \ddots & w\_{2,n}I_m\\
\vdots & \ddots & \ddots & \vdots\\
w\_{n,1}I_m & \dots & w\_{n,n-1}I_m & 0I_m
\end{array}\right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G_%7Bl%2CS_%7Bt%7D%7D%20%3D%20%20%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0AB_%7B1%2Cl%2CS_%7B1%2Ct%7D%7D%20%26%200%20%26%20%5Cdots%20%26%200%5C%5C%0A0%20%26%20B_%7B2%2Cl%2CS_%7B2%2Ct%7D%7D%20%26%20%20%26%20%5Cvdots%5C%5C%0A%5Cvdots%20%26%20%20%26%20%5Cddots%20%26%200%5C%5C%0A0%20%26%20%5Cdots%20%26%200%20%26%20B_%7Bn%2Cl%2CS_%7Bn%2Ct%7D%7D%0A%5Cend%7Barray%7D%5Cright%29%2B%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0AA_%7B1%2Cl%2CS_%7B1%2Ct%7D%7D%5E%7B%2A%7D%20%26%200%20%26%20%5Cdots%20%26%200%5C%5C%0A0%20%26%20A_%7B2%2Cl%2CS_%7B2%2Ct%7D%7D%5E%7B%2A%7D%20%26%20%20%26%20%5Cvdots%5C%5C%0A%5Cvdots%20%26%20%20%26%20%5Cddots%20%26%200%5C%5C%0A0%20%26%20%5Cdots%20%26%200%20%26%20A_%7Bn%2Cl%2CS_%7Bn%2Ct%7D%7D%5E%7B%2A%7D%0A%5Cend%7Barray%7D%5Cright%29%5Cleft%28%5Cbegin%7Barray%7D%7Bcccc%7D%0A0I_m%20%26%20w_%7B1%2C2%7DI_m%20%26%20%5Cdots%20%26%20w_%7B1%2Cn%7DI_m%5C%5C%0Aw_%7B2%2C1%7DI_m%20%26%200I_m%20%26%20%5Cddots%20%26%20w_%7B2%2Cn%7DI_m%5C%5C%0A%5Cvdots%20%26%20%5Cddots%20%26%20%5Cddots%20%26%20%5Cvdots%5C%5C%0Aw_%7Bn%2C1%7DI_m%20%26%20%5Cdots%20%26%20w_%7Bn%2Cn-1%7DI_m%20%26%200I_m%0A%5Cend%7Barray%7D%5Cright%29 "G_{l,S_{t}} =  \left(\begin{array}{cccc}
B_{1,l,S_{1,t}} & 0 & \dots & 0\\
0 & B_{2,l,S_{2,t}} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & B_{n,l,S_{n,t}}
\end{array}\right)+\left(\begin{array}{cccc}
A_{1,l,S_{1,t}}^{*} & 0 & \dots & 0\\
0 & A_{2,l,S_{2,t}}^{*} &  & \vdots\\
\vdots &  & \ddots & 0\\
0 & \dots & 0 & A_{n,l,S_{n,t}}^{*}
\end{array}\right)\left(\begin{array}{cccc}
0I_m & w_{1,2}I_m & \dots & w_{1,n}I_m\\
w_{2,1}I_m & 0I_m & \ddots & w_{2,n}I_m\\
\vdots & \ddots & \ddots & \vdots\\
w_{n,1}I_m & \dots & w_{n,n-1}I_m & 0I_m
\end{array}\right)")

The main functions for this class of models are

- MRCIGVARData,
- MRCIGVARest
- MRCIGVAR_Select
- irf_MRCIGVAR_CB
- girf_MRCIGVAR_RM_CB

The sample codes below show how to use MRCIGVARData as a data generating
process to generate data from MRCIGVAR processes with a variety of
specifications. With the data provided, MRCIGVARest estimates the
parameters in an MRCIGVAR object. MRCIGVAR_Select can be used to select
an MRCIGVAR model’s model specification.

``` r

m = 3
n = 5
p = c(2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0); dim(p) = c(5,3,2)
p = p[1:n,,]; p[,1,] = 3; p[,2,] = 2

p[2,2,] = 3

TH = c(1:n)*0; dim(TH) = c(1,n)
SESVI=rep(1,4,7,10,13)
r  = rep(1,n)
i = 1
S = 2

T= 200

XX = matrix(rnorm(T*10*2),T*10,2) 
dim(XX) = c(T,2,5,2)

p[,3,]=2


res_d <- MRCIGVARData(m=3,n=5,p=p,TH=TH,T=T,S=2, SESVI=c(1,4,7,10,13),type="exog0",r=rep(1,5),DFYflag=0,Ncommfakt=1,X=XX)    ## m: number of variables, n: number of countries

max(abs(res_d$Y))
plot(ts(res_d$Y[,1:10]))
STAT(res_d$Go[,,,2]);
STAT(res_d$GC[,,,1]);STAT(res_d$GC[,,,2]);
STAT(res_d$Go[,,,1])
max(abs(res_d$Y))                 

Ao=res_d$Ao;Bo=res_d$Bo;Co=res_d$Co;


res_d <- MRCIGVARData(m=3,n=5,p=p,TH=TH,T=T,S=2, SESVI=c(1,4,7,10,13),type="exog0",Ao=Ao,Bo=Bo,Co=Co,r = rep(1,5),DFYflag = 0,Ncommfakt=1,X=XX)    ## m: number of variables, n: number of countries

res_e  = MRCIGVARest(res=res_d)




sum(abs(res_e$Ao-res_d$Ao))        
sum(abs(res_e$Bo-res_d$Bo))
sum(abs(res_e$Co-res_d$Co))    

res_e$Do


res_d <- MRCIGVARData(m=3,n=5,p=p,TH=TH,T=T,S=2, SESVI=c(1,4,7,10,13),type="exog1",Ao=Ao,Bo=Bo,Co=Co,r = rep(1,5),DFYflag = 0,Ncommfakt=1,X=XX)    ## m: number of variables, n: number of countries

res_e  = MRCIGVARest(res=res_d)

sum(abs(res_e$Ao-res_d$Ao))        
sum(abs(res_e$Bo-res_d$Bo))
sum(abs(res_e$Co-res_d$Co))    

res_e$Do

res_e$Summary
```

## Impulse Response functions of VAR, CIVAR, GVAR, CIGVAR, MRVAR, MRCIVAR, MRGVAR and MRCIGVAR

In GVAR, CIGVAR, MRGVAR, and MRCIGVAR conventionally generalized impulse
response functions [^3] are used. The following variations can be
generated:

- country shocks
- global shocks
- regional shocks
- coordinated actions
- global responses
- regional responses

The impulse response function options comb and comb1 allow to calculate
global shocks, regional shocks, coordinated actions, and generate the
system’s global and regional responses.

### Global Shocks

A global shock is defined as the weighted average of shocks of all
countries in the same variable [^4]

Using the weightings which are used in the construction of the foreign
variables in each countries, we can defined a global financial stress
shock in state S:
![v\_{t,S}^{gf}=(w\_{1},0,w\_{2},0,...,w\_{N},0)\epsilon\_{t,S}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_%7Bt%2CS%7D%5E%7Bgf%7D%3D%28w_%7B1%7D%2C0%2Cw_%7B2%7D%2C0%2C...%2Cw_%7BN%7D%2C0%29%5Cepsilon_%7Bt%2CS%7D "v_{t,S}^{gf}=(w_{1},0,w_{2},0,...,w_{N},0)\epsilon_{t,S}")
with
![\sum\_{i=1}^{N}w\_{i}=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csum_%7Bi%3D1%7D%5E%7BN%7Dw_%7Bi%7D%3D1 "\sum_{i=1}^{N}w_{i}=1").
Denoting
![w'=(w\_{1},0,w\_{2},0,...,w\_{N},0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w%27%3D%28w_%7B1%7D%2C0%2Cw_%7B2%7D%2C0%2C...%2Cw_%7BN%7D%2C0%29 "w'=(w_{1},0,w_{2},0,...,w_{N},0)"),
we have
![v\_{t,S}^{gf}=w'\epsilon\_{t,S}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_%7Bt%2CS%7D%5E%7Bgf%7D%3Dw%27%5Cepsilon_%7Bt%2CS%7D "v_{t,S}^{gf}=w'\epsilon_{t,S}")
and the standard error of the global financial shock is  
![\sqrt{w'\Sigma\_{S}w}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csqrt%7Bw%27%5CSigma_%7BS%7Dw%7D "\sqrt{w'\Sigma_{S}w}").

Following Chudik and Pesaran (2018), the generalized impulse response
function to a one standard error global financial shock is given by:

![\mbox{irf}(h,y\|v\_{t,S}^{gf}) = E\big(y\_{t+h}\|v\_{t,S}^{fg}=\sqrt{w'\Sigma\_{S}w},\Omega\_{t-1}\big)-E(y\_{t+h}\|\Omega\_{t-1}) =  \frac{R\_{h}\Sigma\_{S}w}{\sqrt{w'\Sigma\_{S}w}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmbox%7Birf%7D%28h%2Cy%7Cv_%7Bt%2CS%7D%5E%7Bgf%7D%29%20%3D%20E%5Cbig%28y_%7Bt%2Bh%7D%7Cv_%7Bt%2CS%7D%5E%7Bfg%7D%3D%5Csqrt%7Bw%27%5CSigma_%7BS%7Dw%7D%2C%5COmega_%7Bt-1%7D%5Cbig%29-E%28y_%7Bt%2Bh%7D%7C%5COmega_%7Bt-1%7D%29%20%3D%20%20%5Cfrac%7BR_%7Bh%7D%5CSigma_%7BS%7Dw%7D%7B%5Csqrt%7Bw%27%5CSigma_%7BS%7Dw%7D%7D "\mbox{irf}(h,y|v_{t,S}^{gf}) = E\big(y_{t+h}|v_{t,S}^{fg}=\sqrt{w'\Sigma_{S}w},\Omega_{t-1}\big)-E(y_{t+h}|\Omega_{t-1}) =  \frac{R_{h}\Sigma_{S}w}{\sqrt{w'\Sigma_{S}w}}")

where
![R\_{h}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bh%7D "R_{h}")
are obtained recursively as

![R\_{h}=\sum\_{l=1}^{p}G\_{l,S}R\_{h-1}\mbox{ with }R\_{0}=I\mbox{ and }R\_{l}=0\mbox{ for }j\<0,](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bh%7D%3D%5Csum_%7Bl%3D1%7D%5E%7Bp%7DG_%7Bl%2CS%7DR_%7Bh-1%7D%5Cmbox%7B%20with%20%7DR_%7B0%7D%3DI%5Cmbox%7B%20and%20%7DR_%7Bl%7D%3D0%5Cmbox%7B%20for%20%7Dj%3C0%2C "R_{h}=\sum_{l=1}^{p}G_{l,S}R_{h-1}\mbox{ with }R_{0}=I\mbox{ and }R_{l}=0\mbox{ for }j<0,")

![\Omega\_{t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5COmega_%7Bt%7D "\Omega_{t}")
contains all information available up to time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t").
For time period
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
the simultaneous response is

![\mbox{irf}(0,y,v\_{t,S}^{gf})=\frac{\Sigma\_{S}w}{\sqrt{w'\Sigma\_{S}w}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmbox%7Birf%7D%280%2Cy%2Cv_%7Bt%2CS%7D%5E%7Bgf%7D%29%3D%5Cfrac%7B%5CSigma_%7BS%7Dw%7D%7B%5Csqrt%7Bw%27%5CSigma_%7BS%7Dw%7D%7D "\mbox{irf}(0,y,v_{t,S}^{gf})=\frac{\Sigma_{S}w}{\sqrt{w'\Sigma_{S}w}}")

### Regional Shocks

In a similar way, regional shocks can be constructed through appropriate
choices of weightings. For a euro-zone financial stress shock, all the
weightings that correspond to non-euro-zone countries are set to be
zero, and only the weightings that correspond to the financial stress
variables of the euro-zone countries are positive. These positive
weightings have to be re normalized so that they sum up to one. In this
context, a country shock is a specific case of a regional shock, where
the region consists of only one country.

### Coordinated Policy Actions

We define coordinated financial stress shocks in state
![S](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S "S"):
![v\_{t,S}^{cf}=(v\_{1},0,v\_{3},0,...,v\_{2N-1},0)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_%7Bt%2CS%7D%5E%7Bcf%7D%3D%28v_%7B1%7D%2C0%2Cv_%7B3%7D%2C0%2C...%2Cv_%7B2N-1%7D%2C0%29 "v_{t,S}^{cf}=(v_{1},0,v_{3},0,...,v_{2N-1},0)"),
where
![v\_{i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_%7Bi%7D "v_{i}")
is the scale of the financial stress shock
![\epsilon\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bi%2Ct%7D "\epsilon_{i,t}").
Following the idea of generalized impulse response function, we define
the impulse response function for concerted shocks as follows

where
![R\_{h}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bh%7D "R_{h}")
are obtained recursively as

![R\_{h}=\sum\_{l=1}^{p}G\_{l,S}R\_{h-1}\mbox{ with }R\_{0}=I\mbox{ and }R\_{l}=0\mbox{ for }j\<0,](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bh%7D%3D%5Csum_%7Bl%3D1%7D%5E%7Bp%7DG_%7Bl%2CS%7DR_%7Bh-1%7D%5Cmbox%7B%20with%20%7DR_%7B0%7D%3DI%5Cmbox%7B%20and%20%7DR_%7Bl%7D%3D0%5Cmbox%7B%20for%20%7Dj%3C0%2C "R_{h}=\sum_{l=1}^{p}G_{l,S}R_{h-1}\mbox{ with }R_{0}=I\mbox{ and }R_{l}=0\mbox{ for }j<0,")

Defining the impulse response function in this way, the simultaneous
response of
![y\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_%7Bi%2Ct%7D "y_{i,t}")
to the concerted shock
![v\_{t,S}^{cf}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;v_%7Bt%2CS%7D%5E%7Bcf%7D "v_{t,S}^{cf}")
is the sum of conditional expected value of
![\epsilon\_{i,t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bi%2Ct%7D "\epsilon_{i,t}")
given concerted shocks of other countries plus the shock of the own
country.

### Global effects and reginal effects

The global effect of a policy shock can be evaluated as a linear
combination of each country’s response to the shock (see Greenwood,
Nguyen and Shin (2013) for this type of applications.). Using the same
weights which are used to construct the global shocks, we can calculate
the global effects of a policy shock, which can be a global shock or a
coordinated policy shock.

``` r
###   Using data from Chen and Semmler (2018)
###
###   country shocks
###
p        <- p_G
res_d    <- MRGVARData(m=2,n=15,p=p,T=370,S=2,SESVI=c((1:15)*2-1),type="const")
res_d$Y  <- FSIIO
res_d$W  <- Wmat_G
res_d$TH <- t(as.matrix(TH_G))

res_e   <- MRGVARest(res_d)
STAT(res_e$Go[,,,1])
STAT(res_e$Go[,,,2])

IRF_CB = irf_MRGVAR_CB(res=res_e,nstep=20,comb=NA,state=rep(1,15),irf="gen1",runs=20,conf=c(0.05,0.95),NT=1)
# FSI responses to USA FSI shock
IRF_list <- IRF_graph(IRF_CB[[1]],Names=names(res_d$Y),impulse=c(1),response=seq(1,30,2),ncol=5)
IRF_list <- IRF_graph(IRF_CB[[2]],Names=names(res_d$Y),impulse=c(1),response=seq(1,30,2),ncol=5)
# FSI responses to USA IP shock
IRF_list <- IRF_graph(IRF_CB[[1]],Names=names(res_d$Y),impulse=c(2),response=seq(1,30,2),ncol=5)
IRF_list <- IRF_graph(IRF_CB[[2]],Names=names(res_d$Y),impulse=c(2),response=seq(1,30,2),ncol=5)
```

For country shocks different types of the impulse response functions are
listed below

- irf = “gen”: generalized impulse response function with one standard
  deviation shocks
- irf = “gen1”: generalized impulse response function with one unit
  shocks
- irf = “genN1”: generalized impulse response function with negative one
  unit shocks
- irf = “chol”: Cholesky decomposition with one standard deviation
  shocks
- irf = “chol1”: Cholesky decomposition with one unit shocks
- irf = “PTdecom”: permanent and transitory decomposition for
  cointegrated system. In this case option G has to be used.
- comb = NA

The above options are also used in VAR, CIVAR, MRVAR, and MRCIVAR models
to specify the different types of shocks. For country shocks, option
comb is empty or NA.

To generated a global or a regional shock or a concerted policy action
the following type of irf options can be used

- irf = “consert0”: this type is applied with a weighting vector in a
  column of the
  (![n \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n"))
  matrix comb to generate the coordinated policy action. The shocks are
  scaled in standard deviations
- irf = “concert1”: this type is applied with a weighting vector in a
  column of the
  (![n \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n"))
  matrix comb to generate the coordinated policy action. The shocks are
  scaled in units,
- irf = “concertc”: this type is applied with a weighting vector in a
  column of the
  (![n \times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%5Ctimes%20n "n \times n"))
  matrix comb to generate the coordinated policy action.
- irf = “comb”: this is used together with the GVAR weighting matrix to
  create the global, regional or country specific shock. The weights in
  the column of comb must sum up tp one. The shocks are scaled in
  standard deviations.
- irf = “comb1”: this is used together with the GVAR weighting matrix to
  create the global, regional or country specific shock. The weights in
  the column of comb must sum up tp one. The shocks are scaled in
  standard deviations.
- comb = a weighting matrix with column sums of one.

``` r

######## The following weighting vector are taken from the MRGVAR weighting matrix to form global and regional weights.
#######################################################################################################
# c("BEL","AUT","NLD","NOR","SWE","FIN","DNK","CAN","ESP","ITA","GBR","FRA","DEU","JPN","USA")
#### global and regional shocks

p        <- p_G
res_d    <- MRGVARData(m=2,n=15,p=p,T=370,S=2,SESVI=c((1:15)*2-1),type="const")
res_d$Y  <- FSIIO
res_d$W  <- Wmat_G
res_d$TH <- t(as.matrix(TH_G))

res_e   <- MRGVARest(res_d)

SW   = c(16768,     4919, 3730, 2806, 2678, 2149,   1393,   1826, 313, 267, 579, 512, 428,  524 ,853 )
SWEU = c(16768*0,   4919*0, 3730, 2806, 2678*0, 2149,   1393,   1826*0, 313*0, 267, 579*0, 512*0, 428,  524 ,853 )

SWnorth = c(16768*0,4919*0,3730*0,2806*0,2678*0,2149*0,1393*0,1826*0,313*1,267*1,579*1,512*1,428*0,524*0,853*0 )
SWsouth = c(16768*0,4919*0,3730*0,2806*0,2678*0,2149*1,1393*1,1826*0,313*0,267*0,579*0,512*0,428*0,524*0,853*0 )
SWameri = c(16768*1,4919*0,3730*0,2806*0,2678*0,2149*0,1393*0,1826*1,313*0,267*0,579*0,512*0,428*0,524*0,853*0 )

comb_f = SW2comb(SW,15,2,1)
comb_i = SW2comb(SW,15,2,2)
comb_eu_f =  SW2comb(SWEU,15,2,1)
comb_eu_i =  SW2comb(SWEU,15,2,2)
comb_nt_f =  SW2comb(SWnorth,15,2,1)
comb_nt_i =  SW2comb(SWnorth,15,2,2)
comb_st_f =  SW2comb(SWsouth,15,2,1)
comb_st_i =  SW2comb(SWsouth,15,2,2)
comb_am_f =  SW2comb(SWameri,15,2,1)
comb_am_i =  SW2comb(SWameri,15,2,2)

comb_all = comb_f*0
comb_all[,1]  = comb_f[,1]
comb_all[,2]  = comb_i[,1]
comb_all[,3]  = comb_eu_f[,1]
comb_all[,4]  = comb_eu_i[,1]
comb_all[,5]  = comb_nt_f[,1]
comb_all[,6]  = comb_nt_i[,1]
comb_all[,7]  = comb_st_f[,1]
comb_all[,8]  = comb_st_i[,1]
comb_all[,9]  = comb_am_f[,1]
comb_all[,10] = comb_am_i[,1]
colSums(comb_all)
glw = comb_all[,1:2]


###  Global and regional shocks  ##############################################
###  Note that the first two columns of comb_all generate a global FSI shock and an IP shock respectively.
###  The third and the fourth columns generate an EU FSI and an EU IP shock respectively.
###  The fifth and sixth columns the Northern Europe shocks and so on.
###  The seventh and eighth columns the Southern Europe shocks and so on.
###  Then the next two columns are the American FSI and IP 

IRF_CB  = irf_MRGVAR_CB(res=res_e,nstep=20,comb=comb_all,state=rep(2,15),irf="comb1",runs=20,conf=c(0.05,0.95),NT=1)


iNames <- c("Global FSI","Global IP","EU FSI","EU IP","Northern FSI","Northern IP",
                            "Southern FSI","Southern IP", "American FSI","American IP")


# FSI responses to a global one unit FSI shock
IRF_list <- IRF_graph(IRF_CB[[1]],Names=names(res_d$Y),IName=iNames,impulse=c(1),response=seq(1,30,2),ncol=5)
# IP responses to a global one unit FSI  shock
IRF_list <- IRF_graph(IRF_CB[[1]],Names=names(res_d$Y),IName=iNames,impulse=c(1),response=seq(2,30,2),ncol=5)

# FSI responses to an EU one unit FSI shock
IRF_list <- IRF_graph(IRF_CB[[1]],Names=names(res_d$Y),IName=iNames,impulse=c(3),response=seq(1,30,2),ncol=5)
# IP responses to a EU one unit FSI shock
IRF_list <- IRF_graph(IRF_CB[[1]],Names=names(res_d$Y),IName=iNames,impulse=c(3),response=seq(2,30,2),ncol=5)

# FSI responses to a global one unit IP shock
IRF_list <- IRF_graph(IRF_CB[[1]],Names=names(res_d$Y),IName=iNames,impulse=c(2),response=seq(1,30,2),ncol=5)
# IP responses to a global one unit IP  shock
IRF_list <- IRF_graph(IRF_CB[[1]],Names=names(res_d$Y),IName=iNames,impulse=c(2),response=seq(2,30,2),ncol=5)



###############################################################################
#############   country specific negatively scaled shocks

Countryshock = matrix(0,30,30)
for (i in 1:15) Countryshock[2*i-1,2*i-1]=-1.5
for (i in 1:15) Countryshock[2*i,2*i]= 2.5

IRF_CB <- irf_MRGVAR_CB(res=res_e,nstep=20,comb=Countryshock,state=rep(2,15),irf="comb1",runs=20,conf=c(0.05,0.95),NT=1)
#x11()  # FSI responses to USA scaled FSI shock (-1.5)
IRF_list <- IRF_graph(IRF_CB[[1]],impulse=c(1),response=seq(1,30,2),ncol=5)
#x11()  # FSI responses to USA scaled IP shock  (2.5)
IRF_list <- IRF_graph(IRF_CB[[1]],impulse=c(2),response=seq(1,30,2),ncol=5)


###########################################################################
###  concerted policy action IRF
###########################################################################

### concertA matrix specifies different degrees of the concerted actions
concertedA = matrix(0,30,30)
for (i in 1:15 ) for ( j in 1:10 )   concertedA[2*i-1,j] = -1+(j-1)*0.1

### Using irf = "concert0", "concert1" or "concertc" together with comb = concertedA that is a (30 x 30)
### matrix specifying the intensity of the concerted actions, the IRF of
### coordinated policy actions can be generated.


IRF  = irf_MRGVAR(res=res_e,nstep=20,comb=concertedA,state=rep(2,15),irf="concertc")

IRF_CB = irf_MRGVAR_CB(res=res_e,nstep=20,comb=concertedA,state=rep(2,15),irf="concertc",runs=20,conf=c(0.05,0.95),NT=1)

# FSI responses to concerted actions in FSI (1: each country with one unit)
IRF_list <- IRF_graph(IRF_CB[[1]],impulse=c(1),response=seq(1,30,2),ncol=5)
# FSI responses to concerted actions in FSI (2: each country with 0.9 unit)
IRF_list <- IRF_graph(IRF_CB[[1]],impulse=c(2),response=seq(1,30,2),ncol=5)
```

The preceding example shows how we can aggregate different shocks on the
shock side. We can also aggregate on the responses on the response side
using the weighting matrix comb_all to obtain global and regional
responses.

``` r

##########################################################################
### Global and regional responses
##########################################################################

### The same weighting matrix that is used to aggregate global impulses can 
### also be used to aggregate global responses and regional responses. 


Global_Response_CB <- irf_GloabalResponse_CB(IRF_CB,comb_all)


IRF_list <- IRF_graph(Global_Response_CB[[1]],Names=iNames,INames=iNames,impulse=c(1),response=c(1:8),ncol=4)
```

For a more detailed use of MRCIGVAR please read
Introduction-to-MRCIGVAR.html or click the link below:
<https://htmlpreview.github.io/?https://github.com/puchen8229/MRCIGVAR/blob/master/vignettes/Introduction-to-MRCIGVAR.html>

[^1]: G. Koop, M. H. Pesaran, and S. M. Potter (1998) Impulse response
    analysis in nonlinear multivariate models, Journal of Econometrics
    74 (1996), 19-147.

[^2]: M. Ehrmann, M. Ellison, and N. Valla (2003) Regime-dependent
    impulse response functions in a Markov-switching vector
    autoregression model, Economics Letters 78(3) (2003), 295-299.

[^3]: H. Pesaran and Y. Shin （1998）for more details of the generalized
    impulse response function.

[^4]: M. Greenwood-Nimmo, V. Nguyen, and Y. Shin (2013) for more details
    of impulse response functions in GVAR models.

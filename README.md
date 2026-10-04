# sscor-MCD
Code to reproduce the figures of paper "Finite-sample calibration of robust covariance estimation under heavy tails"

> **Note:** This repository contains code and instructions for replicating the results of a paper currently under peer review. Author information is intentionally omitted.

## Abstract
The Minimum Covariance Determinant (MCD) estimator is one of the most widely used high-breakdown estimators of multivariate location and scatter. Its practical success relies on suitable corrections, which remove the asymptotic and small-sample biases introduced by trimming. While explicit consistency factors are available under the multivariate normal distribution and, more recently, under multivariate Student-$t$ models, these corrections only address the asymptotic behaviour of the estimator.
%
In the Gaussian case, \cite{pis+al:02} showed that the MCD also exhibits a non-negligible finite-sample bias, particularly for small sample sizes and relatively large dimensions. Such bias propagates to robust Mahalanobis distances and may affect subsequent procedures, including multivariate outlier detection and robust inference.
%
This paper investigates the finite-sample behaviour of the MCD estimator and of the related robust distances under multivariate Student-$t$ distributions, a potentially relevant scenario in many application fields. Extensive Monte Carlo experiments are carried out over a broad range of sample sizes, dimensions, trimming levels and degrees of freedom. Based on these experiments, empirical correction factors are obtained and approximated by smooth functions of the sample size and trimming proportion, providing a computationally efficient calibration suitable for practical implementations to heavy-tailed data. The proposed corrections are developed and made available both for the scatter bias and the distance quantiles.

## Finite-sample calibration of MCD scatter bias and robust distance quantiles

While asymptotic consistency factors for the Minimum Covariance Determinant (MCD) estimator under multivariate Student-$t$ distributions have recently been established, the estimator exhibits non-negligible finite-sample bias that propagates directly to robust Mahalanobis distances and outlier detection rules. This work provides a comprehensive finite-sample calibration framework for heavy-tailed elliptical models along two complementary directions: (i) an empirical determinant-based correction factor $\hat{\delta}_n$ to remove finite-sample bias from the robust scatter matrix, and (ii) calibrated empirical cutoffs for the distribution of squared robust distances to ensure accurate nominal size in multivariate outlier detection. Using extensive Monte Carlo simulations across varying sample sizes $n$, dimensions $p$, degrees of freedom $\nu$, and trimming proportions $\alpha_0$, the calibration factors are modeled via smooth interpolation functions, enabling fast and accurate calibration for arbitrary sample configurations without requiring additional simulations.

## FSDA and other dependencies

The MCD estimator, which is at the core of the paper, and few other useful
functions rely on the free MATLAB Add-On *FSDA* 
(see https://it.mathworks.com/matlabcentral/fileexchange/72999-fsda-flexible-statistics-data-analysis-toolbox). 
FSDA requires the Statistical and Machine Learning Toolbox. The Parallel 
Processing Toolbox is necessary if the reader needs to replicate the estimates, 
which take lot of time otherwise.

## Results replication: requirements and setup

This repository provides the MATLAB code and the data necessary to replicate the 
main figures and tables from both the main manuscript and the supplementary 
information document. These codes can be used in two modalities: 
- 1) In a standalone MATLAB licensed installation; 
- 2) In the license-free MATLAB Online: https://it.mathworks.com/products/matlab-online.html

### To prepare your local environment 
1. Clone this repository: `git clone https://github.com/UniprJRC/sscor-MCD.git`
2. Navigate to the repository directory: `cd <your path to tdist-MCD>
3. Install FSDA from "Install App" of the standard MATLAB distribution 

### To use MATLAB Online
The free "MATLAB Online" platform can be used up to 20 hours a month, 
which are sufficient to replicate all results except the simulation 
of the KS and AD quantiles. For this, it is sufficient to click on the
link in the table below. Once you are in the MATLAB cloud environment,
go under the "Add Ons section" of the "Home tab" and install FSDA. 

<!---
- R is used version X.X.X or higher
- R packages: [list required packages]
--->

<!---
4. Install required R packages: `Rscript install_packages.R`
--->


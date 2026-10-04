# sscor-MCD
Code to reproduce the figures of paper "Finite-sample calibration of robust covariance estimation under heavy tails"

> **Note:** This repository contains code and instructions for replicating the results of a paper currently under peer review. Author information is intentionally omitted.

## Abstract
The Minimum Covariance Determinant (MCD) estimator is one of the most widely used high-breakdown estimators of multivariate location and scatter. Its practical success relies on suitable corrections, which remove the asymptotic and small-sample biases introduced by trimming. While explicit consistency factors are available under the multivariate normal distribution and, more recently, under multivariate Student-$t$ models, these corrections only address the asymptotic behaviour of the estimator.

In the Gaussian case, \cite{pis+al:02} showed that the MCD also exhibits a non-negligible finite-sample bias, particularly for small sample sizes and relatively large dimensions. Such bias propagates to robust Mahalanobis distances and may affect subsequent procedures, including multivariate outlier detection and robust inference.

This paper investigates the finite-sample behaviour of the MCD estimator and of the related robust distances under multivariate Student-$t$ distributions, a potentially relevant scenario in many application fields. Extensive Monte Carlo experiments are carried out over a broad range of sample sizes, dimensions, trimming levels and degrees of freedom. Based on these experiments, empirical correction factors are obtained and approximated by smooth functions of the sample size and trimming proportion, providing a computationally efficient calibration suitable for practical implementations to heavy-tailed data. The proposed corrections are developed and made available both for the scatter bias and the distance quantiles.

## Practical Impact: Unbiased Scatter Estimation and Accurate Outlier Testing

In finite samples, asymptotic consistency corrections alone are insufficient: they underestimate scatter and render standard outlier tests overly liberal, resulting in substantial rates of false alarms on heavy-tailed data. To solve this, this repository delivers practical finite-sample tools that (i) remove the systematic bias in the MCD scatter matrix via determinant-based correction factors $\hat{\delta}_n$, and (ii) provide calibrated cutoffs for squared robust distances that strictly control the empirical false-positive rate. Through precomputed calibration tables and smooth interpolation functions, practitioners can instantly apply accurate corrections to any combination of sample size $n$, dimension $p$, degrees of freedom $\nu$, and trimming proportion $\alpha_0$ in real-time, eliminating the need for expensive simulation studies.

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



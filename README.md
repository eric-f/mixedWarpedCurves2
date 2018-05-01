R package: mixedWarpedCurves2
================

-   [Installation](#installation)
-   [Model](#model)
-   [Example 1: Curve registration](#example-1-curve-registration)
-   [Example 2: Clustering on phase variation](#example-2-clustering-on-phase-variation)
-   [Reference](#reference)

Installation
============

To install the latest version:

    library(devtools)
    devtools::install_github('eric-f/mixedWarpedCurves2')

Model
=====

We model the functional data as noisy observations of latent curves, ![y\_i](https://latex.codecogs.com/png.latex?y_i "y_i")'s, at time ![t\_{ij}](https://latex.codecogs.com/png.latex?t_%7Bij%7D "t_{ij}"), as ![y\_{ij} = y\_i(t) + \\epsilon\_{ij}.](https://latex.codecogs.com/png.latex?y_%7Bij%7D%20%3D%20y_i%28t%29%20%2B%20%5Cepsilon_%7Bij%7D. "y_{ij} = y_i(t) + \epsilon_{ij}.")

The latent curve, ![y\_i(t) = f\_i \\circ h\_i(t),](https://latex.codecogs.com/png.latex?y_i%28t%29%20%3D%20f_i%20%5Ccirc%20h_i%28t%29%2C "y_i(t) = f_i \circ h_i(t),") is composed of two functions, ![f\_i](https://latex.codecogs.com/png.latex?f_i "f_i") and ![h\_i](https://latex.codecogs.com/png.latex?h_i "h_i"), representing amplitude and phase variations. Specifically, the amplitude function is rank-1: ![f\_i = a\_{i,sh} + a\_{i,sc} f(t)](https://latex.codecogs.com/png.latex?f_i%20%3D%20a_%7Bi%2Csh%7D%20%2B%20a_%7Bi%2Csc%7D%20f%28t%29 "f_i = a_{i,sh} + a_{i,sc} f(t)"), with a common base shape ![f](https://latex.codecogs.com/png.latex?f "f") whereas the warping function ![h\_i:\[0, 1\]\\to\[0, 1\]](https://latex.codecogs.com/png.latex?h_i%3A%5B0%2C%201%5D%5Cto%5B0%2C%201%5D "h_i:[0, 1]\to[0, 1]") is continuous and strictly increasing.

Amplitude function
------------------

The base shape, ![f](https://latex.codecogs.com/png.latex?f "f") is either a known function or modelled as a spline that is to be estimated. The amplitude random effect, ![(a\_{i,sh}, a\_{i,sc})](https://latex.codecogs.com/png.latex?%28a_%7Bi%2Csh%7D%2C%20a_%7Bi%2Csc%7D%29 "(a_{i,sh}, a_{i,sc})"), is Gaussian with mean ![(\\mu\_{sh}, \\mu\_{sc})](https://latex.codecogs.com/png.latex?%28%5Cmu_%7Bsh%7D%2C%20%5Cmu_%7Bsc%7D%29 "(\mu_{sh}, \mu_{sc})") and variance-covariance matrix ![\\Sigma](https://latex.codecogs.com/png.latex?%5CSigma "\Sigma").

`fsim_unimodal()` fit the fixed base shape version where ![f(t) = -4t(1-t)](https://latex.codecogs.com/png.latex?f%28t%29%20%3D%20-4t%281-t%29 "f(t) = -4t(1-t)").

`fsim_mixed_warped_curves()` fit the unknown base shape version of the model. The base shape is modelled using B-spline basis. The order and knot locations of the spline are specified via the options `f_order` and `f_knots` in `control_saem()`.

Warping function
----------------

The warping function is modelled using B-spline basis,

![h\_i(t) = \\sum\_{k=1}^{K\_h} w\_{i,k}B\_k(t),](https://latex.codecogs.com/png.latex?h_i%28t%29%20%3D%20%5Csum_%7Bk%3D1%7D%5E%7BK_h%7D%20w_%7Bi%2Ck%7DB_k%28t%29%2C "h_i(t) = \sum_{k=1}^{K_h} w_{i,k}B_k(t),")

 where ![0\\equiv w\_{i,1} \\leq \\cdots \\leq w\_{i,K\_h} \\equiv=1](https://latex.codecogs.com/png.latex?0%5Cequiv%20w_%7Bi%2C1%7D%20%5Cleq%20%5Ccdots%20%5Cleq%20w_%7Bi%2CK_h%7D%20%5Cequiv%3D1 "0\equiv w_{i,1} \leq \cdots \leq w_{i,K_h} \equiv=1"). The first difference of the random coefficients, ![(w\_{i,2}, w\_{i,3}-w\_{i,2}, \\ldots, 1-w\_{i,K\_h-1}),](https://latex.codecogs.com/png.latex?%28w_%7Bi%2C2%7D%2C%20w_%7Bi%2C3%7D-w_%7Bi%2C2%7D%2C%20%5Cldots%2C%201-w_%7Bi%2CK_h-1%7D%29%2C "(w_{i,2}, w_{i,3}-w_{i,2}, \ldots, 1-w_{i,K_h-1}),") follows a Dirichlet distribution whose mean is fixed such that ![E(h\_i(t))=t](https://latex.codecogs.com/png.latex?E%28h_i%28t%29%29%3Dt "E(h_i(t))=t").

The order and knot locations of the spline are specified via the options `h_order` and `h_knots` in `control_saem()`.

Example 1: Curve registration
=============================

![](README_files/figure-markdown_github-ascii_identifiers/example-1-data-1.png)![](README_files/figure-markdown_github-ascii_identifiers/example-1-data-2.png)

Estimate model parameters and warping functions by SAEM
-------------------------------------------------------

``` r
out <- fsim_mixed_warped_curves(
  y = sim_data$y,
  obs_time = sim_data$x,
  curve_id = sim_data$id,
  saem_control = control_saem(
    n_saem_iter = 1000,
    n_saem_burn = 100,
    n_mcmc_burn = 5,
    h_knots = seq(0, 1, length=5),
    f_knots = seq(0, 1, length=7)
  ))
```

    ## [1] "Randomizing initial cluster labels..."
    ## Import parameters...
    ## Generate Cholesky centering matrix...
    ## Import data...
    ## Initialize basis evaluation matrices...
    ## Entering SAEM loop...
    ## 0.0%...5.0%...Initialize clustering with user inputs...
    ## cluster_size
    ##    50.0000
    ## 
    ## p_clusters
    ##    1.0000
    ## 
    ## 10.0%...15.0%...20.0%...25.0%...30.0%...35.0%...40.0%...45.0%...50.0%...55.0%...60.0%...65.0%...70.0%...75.0%...80.0%...85.0%...90.0%...95.0%...(Done)

Fitted curves and registered curves
-----------------------------------

![](README_files/figure-markdown_github-ascii_identifiers/example-1-plots-1.png)![](README_files/figure-markdown_github-ascii_identifiers/example-1-plots-2.png)

![](README_files/figure-markdown_github-ascii_identifiers/example-1-tracks-1.png)

Example 2: Clustering on phase variation
========================================

![](README_files/figure-markdown_github-ascii_identifiers/example-2-data-1.png)![](README_files/figure-markdown_github-ascii_identifiers/example-2-data-2.png)

Fit mixture of warping function by SAEM
---------------------------------------

### Known unimodal shape

``` r
clust_out <- fsim_unimodal(
  y = mix_data$y,
  obs_time = mix_data$x,
  curve_id = mix_data$id,
  n_clust = 4,
  saem_control = control_saem(
    n_saem_iter = 2000,
    n_saem_burn = 200,
    n_mcmc_burn = 5,
    h_knots = seq(0, 1, length=5)
  ))
```

    ## [1] "Randomizing initial cluster labels..."
    ## Entering SAEM loop...
    ## 0.0%...5.0%...Initialize clustering with user inputs...
    ## cluster_size
    ##    43.0000
    ##    53.0000
    ##    48.0000
    ##    56.0000
    ## 
    ## p_clusters
    ##    0.2150
    ##    0.2650
    ##    0.2400
    ##    0.2800
    ## 
    ## 10.0%...15.0%...20.0%...25.0%...30.0%...35.0%...40.0%...45.0%...50.0%...55.0%...60.0%...65.0%...70.0%...75.0%...80.0%...85.0%...90.0%...95.0%...(Done)

### Unknown base shape

``` r
flex_clust_out <- fsim_mixed_warped_curves(
  y = mix_data$y,
  obs_time = mix_data$x,
  curve_id = mix_data$id,
  n_clust = 4,
  saem_control = control_saem(
    n_saem_iter = 2000,
    n_saem_burn = 200,
    n_mcmc_burn = 5,
    h_knots = seq(0, 1, length=7),
    f_knots = seq(0, 1, length=3)
  ))
```

    ## [1] "Randomizing initial cluster labels..."
    ## Import parameters...
    ## Generate Cholesky centering matrix...
    ## Import data...
    ## Initialize basis evaluation matrices...
    ## Entering SAEM loop...
    ## 0.0%...5.0%...Initialize clustering with user inputs...
    ## cluster_size
    ##    59.0000
    ##    50.0000
    ##    52.0000
    ##    39.0000
    ## 
    ## p_clusters
    ##    0.2950
    ##    0.2500
    ##    0.2600
    ##    0.1950
    ## 
    ## 10.0%...15.0%...20.0%...25.0%...30.0%...35.0%...40.0%...45.0%...50.0%...55.0%...60.0%...65.0%...70.0%...75.0%...80.0%...85.0%...90.0%...95.0%...(Done)

Fitted curves and clustering on phase variation
-----------------------------------------------

### Known unimodal shape

![](README_files/figure-markdown_github-ascii_identifiers/example-2-unimodal-plots-1.png)![](README_files/figure-markdown_github-ascii_identifiers/example-2-unimodal-plots-2.png)![](README_files/figure-markdown_github-ascii_identifiers/example-2-unimodal-plots-3.png)

### Unknown base shape

![](README_files/figure-markdown_github-ascii_identifiers/example-2-flexShape-plots-1.png)![](README_files/figure-markdown_github-ascii_identifiers/example-2-flexShape-plots-2.png)![](README_files/figure-markdown_github-ascii_identifiers/example-2-flexShape-plots-3.png)

Reference
=========

Fu, E. and Heckman, N. (2017). Model-based curve registration via stochastic approximation EM algorithm. <https://arxiv.org/abs/1712.07265>

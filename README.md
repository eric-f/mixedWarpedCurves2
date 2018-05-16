mixedWarpedCurves2
================

[![Travis-CI Build Status](https://travis-ci.org/eric-f/mixedWarpedCurves2.svg?branch=master)](https://travis-ci.org/eric-f/mixedWarpedCurves2) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/eric-f/mixedWarpedCurves2?branch=master&svg=true)](https://ci.appveyor.com/project/eric-f/mixedWarpedCurves2)

Dependencies
============

-   GNU Scientific Library (GSL)
-   R packages:
    -   splines
    -   gtools
    -   Rcpp
    -   RcppArmadillo
    -   RcppGSL
    -   RcppProgress
    -   BH

On Windows, the R package will be built using the GSL at <https://github.com/rwinlib/gsl>. On other platforms, you will need to have the library installed manually before this R package can be compiled.

On OSX, you can install GSL via [Homebrew](https://brew.sh). If you do not already have Homebrew, go to <https://brew.sh> for installation instructions. After that, run the following command in your terminal

    brew install gsl

On Debian or Ubuntu, install from terminal by running

    sudo apt-get install libgsl0-dev

Installation
============

Once all the dependencies have been installed, run the following lines in R to install the latest version of this package:

    library(devtools)
    devtools::install_github('eric-f/mixedWarpedCurves2')

Overview
========

The main functions in this package to fit the warped curves model via Stochastic approximation EM algorithms are `fsim()` and `fsim_unimodal()`. Control arguments for the model and algorithm can be specificy via `control_saem()`.

Warped curves model
===================

We model the functional data as noisy observations of latent curves, ![y\_i](https://latex.codecogs.com/png.latex?y_i "y_i")'s, at time ![t\_{ij}](https://latex.codecogs.com/png.latex?t_%7Bij%7D "t_{ij}"), as ![y\_{ij} = y\_i(t\_{ij}) + \\epsilon\_{ij}](https://latex.codecogs.com/png.latex?y_%7Bij%7D%20%3D%20y_i%28t_%7Bij%7D%29%20%2B%20%5Cepsilon_%7Bij%7D "y_{ij} = y_i(t_{ij}) + \epsilon_{ij}").

The latent curve, ![y\_i(t) = f\_i \\circ h\_i(t),](https://latex.codecogs.com/png.latex?y_i%28t%29%20%3D%20f_i%20%5Ccirc%20h_i%28t%29%2C "y_i(t) = f_i \circ h_i(t),") is composed of two functions, ![f\_i](https://latex.codecogs.com/png.latex?f_i "f_i") and ![h\_i](https://latex.codecogs.com/png.latex?h_i "h_i"), representing amplitude and phase variations respectively. Specifically, the amplitude function has the form: ![f\_i = a\_{i,sh} + a\_{i,sc} f(t)](https://latex.codecogs.com/png.latex?f_i%20%3D%20a_%7Bi%2Csh%7D%20%2B%20a_%7Bi%2Csc%7D%20f%28t%29 "f_i = a_{i,sh} + a_{i,sc} f(t)"), where ![f](https://latex.codecogs.com/png.latex?f "f") is a common base shape whereas the warping function is continuous, strictly increasing and satisfies ![h\_i(0)\\equiv 0](https://latex.codecogs.com/png.latex?h_i%280%29%5Cequiv%200 "h_i(0)\equiv 0") and ![h\_i(1)\\equiv 1](https://latex.codecogs.com/png.latex?h_i%281%29%5Cequiv%201 "h_i(1)\equiv 1").

Amplitude component
-------------------

The base shape, ![f](https://latex.codecogs.com/png.latex?f "f"), is either a known function or modelled as a spline that is to be estimated. The amplitude random effect, ![(a\_{i,sh}, a\_{i,sc})](https://latex.codecogs.com/png.latex?%28a_%7Bi%2Csh%7D%2C%20a_%7Bi%2Csc%7D%29 "(a_{i,sh}, a_{i,sc})"), is Gaussian with mean ![(\\mu\_{sh}, \\mu\_{sc})](https://latex.codecogs.com/png.latex?%28%5Cmu_%7Bsh%7D%2C%20%5Cmu_%7Bsc%7D%29 "(\mu_{sh}, \mu_{sc})") and variance-covariance matrix ![\\Sigma](https://latex.codecogs.com/png.latex?%5CSigma "\Sigma").

The fixed base shape version is implemented in `fsim_unimodal()` where ![f(t) = -4t(1-t)](https://latex.codecogs.com/png.latex?f%28t%29%20%3D%20-4t%281-t%29 "f(t) = -4t(1-t)").

The unknown base shape version is implemented `fsim_mixed_warped_curves()` where the base shape is modelled using a B-spline. The order and knot locations of the B-spline are specified by the arguments `f_order` and `f_knots` via `control_saem()`.

Phase component
---------------

The warping functions are modelled using B-splines,

![h\_i(t) = \\sum\_{k=1}^{K\_h} w\_{i,k}B\_k(t),](https://latex.codecogs.com/png.latex?h_i%28t%29%20%3D%20%5Csum_%7Bk%3D1%7D%5E%7BK_h%7D%20w_%7Bi%2Ck%7DB_k%28t%29%2C "h_i(t) = \sum_{k=1}^{K_h} w_{i,k}B_k(t),")

 where ![0\\equiv w\_{i,1} \\leq \\cdots \\leq w\_{i,K\_h} \\equiv 1](https://latex.codecogs.com/png.latex?0%5Cequiv%20w_%7Bi%2C1%7D%20%5Cleq%20%5Ccdots%20%5Cleq%20w_%7Bi%2CK_h%7D%20%5Cequiv%201 "0\equiv w_{i,1} \leq \cdots \leq w_{i,K_h} \equiv 1"). The first difference of the random coefficients, ![(w\_{i,2}, w\_{i,3}-w\_{i,2}, \\ldots, 1-w\_{i,K\_h-1}),](https://latex.codecogs.com/png.latex?%28w_%7Bi%2C2%7D%2C%20w_%7Bi%2C3%7D-w_%7Bi%2C2%7D%2C%20%5Cldots%2C%201-w_%7Bi%2CK_h-1%7D%29%2C "(w_{i,2}, w_{i,3}-w_{i,2}, \ldots, 1-w_{i,K_h-1}),") follows a Dirichlet distribution whose mean is fixed such that ![E(h\_i(t))=t](https://latex.codecogs.com/png.latex?E%28h_i%28t%29%29%3Dt "E(h_i(t))=t").

The order and knot locations of the B-splines are specified by the arguments `h_order` and `h_knots` via `control_saem()`.

Usage
=====

This package comes with a simulated dataset, `data_ex1`, of 50 curves observed on 101 equally spaced time points.

``` r
library(mixedWarpedCurves2)
library(ggplot2)
library(plyr)
library(dplyr)

## Plot base shape
data_ex1 %>%
  dplyr::filter(id=="1") %>%
  ggplot() + 
  geom_line(aes(x=warped_x, y=(y0 - sh)/sc)) + 
  labs(x="Time, t", y=expression(paste("Base shape, ", f(t))))

## Plot observed curves
ggplot(data_ex1) + 
  geom_line(aes(x=x, y=y, col=id), size=0.5, show.legend = F) + 
  labs(x=expression(paste("Time, ", t[ij])), 
       y=expression(paste("Observed curve, ", y[ij])))

## Plot underlying amplitude functions
ggplot(data_ex1) + 
  geom_line(aes(x=warped_x, y=y0, col=id), show.legend = F) + 
  labs(x="Time, t", y=expression(paste("Amplitude function, ", f[i](t))))

## Plot underlying warping functions
ggplot(data_ex1) + 
  geom_line(aes(x=x, y=warped_x, col=id), show.legend = F) + 
  labs(x="Time, t", 
       y=expression(paste("Warping function, ", h[i](t))))
```

![](man/figures/README-example-1-data-1.png)![](man/figures/README-example-1-data-2.png)![](man/figures/README-example-1-data-3.png)![](man/figures/README-example-1-data-4.png) The figures above show the base shape (top left), the observe curves (top right), the amplitude functions, ![f\_i](https://latex.codecogs.com/png.latex?f_i "f_i"), (bottom left) and the warping functions, ![h\_i](https://latex.codecogs.com/png.latex?h_i "h_i"), (bottom right).

Call `fsim()` to fit the model by maximum likelihood via a stochastic approximation expectation-maximization (SAEM) algorithm:

``` r
out <- fsim(
  y = data_ex1$y,
  obs_time = data_ex1$x,
  curve_id = data_ex1$id,
  saem_control = control_saem(
    seed = 0,
    n_saem_iter = 1000,
    n_saem_burn = 100,
    n_mcmc_burn = 5,
    h_knots = seq(0, 1, length=5),
    f_knots = seq(0, 1, length=7)
  ))
```

    ## Entering SAEM loop...
    ## 0.0%...5.0%...10.0%...15.0%...20.0%...25.0%...30.0%...35.0%...40.0%...45.0%...50.0%...55.0%...60.0%...65.0%...70.0%...75.0%...80.0%...85.0%...90.0%...95.0%...(Done)

The returned objects has the following components:

``` r
names(out)
```

    ## [1] "pars"             "curves"           "aux"             
    ## [4] "pars_track"       "se_info"          "y_scaling_factor"

See documentation in R for details of the output:

    help(fsim)

The fitted curves and the warping functions can be extract by:

``` r
out_df <- ldply(out$curves, function(crv){
  data.frame(
    id = as.character(c(crv$curve_id)),
    x = c(crv$x),
    y = c(crv$y),
    warped_x = c(crv$warped_x),
    fitted_y = c(crv$fitted_y)
  )
})
```

The figures below shows the fitted curves (left) and the registered curves (right).

``` r
## Plot fitted curves
ggplot(out_df) + 
  geom_line(aes(x=x, y=y, group=id), size=0.3) + 
  geom_line(aes(x=x, y=fitted_y, col=id), 
            size=0.3, show.legend = F) + 
  labs(x="Time, t", 
       y=expression(paste("Fitted curve, ", hat(y)[i](t))))

## Plot registered curves
ggplot(out_df) + 
  geom_line(aes(x=warped_x, y=y, col=id),
            show.legend = F) + 
  labs(x=expression(paste("Predicted warped time, ", hat(h)[i](t))), 
       y=expression(y[i](t)))
```

![](man/figures/README-example-1-plots-1.png)![](man/figures/README-example-1-plots-2.png)

The sequences of estimated parameters are returned in the `pars_track` component. The figure below shows the estimated error variances as a simple convergence diagnostics:

``` r
sigma2_track <- out$pars_track$sigma2_track
qplot(x=seq(along=sigma2_track),
      y=sigma2_track, 
      geom="line") + 
  labs(x="Iterations", 
       y=expression(paste("Estimated ", sigma[epsilon]^2))) + 
  coord_cartesian(ylim=c(sigma2_track[50], min(sigma2_track)))
```

![](man/figures/README-example-1-tracks-1.png)

Reference
=========

Fu, E. and Heckman, N. (2017). Model-based curve registration via stochastic approximation EM algorithm. <https://arxiv.org/abs/1712.07265>

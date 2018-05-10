# rm(list=ls())
require(fda)
require(gtools)

#' Make base shape by standarding a symmetric beta-density
#'
#' Internal function to be called by sim_warping_mixture()
#' @param x numeric vector between 0 and 1 where the function is evaluated
#' @param a non-negative parameters of the Beta distribution
#' @import stats
base_shape <- function(x, a=2){
  -dbeta(x, a, a) / dbeta(0.5, a, a)
}

#' Make base shape from a beta-density
#'
#' Internal function to be called by sim_warping_mixture()
#' @param x numeric vector between 0 and 1 where the function is evaluated
#' @param a non-negative parameters of the Beta distribution
#' @param b non-negative parameters of the Beta distribution
#' @import stats
beta_shape <- function(x, a=2, b=3){
  -dbeta(x, a, b) / dbeta((a-1)/(a+b-2), a, b)
}


#' Simulate curves with mixture of warping
#'
#' This function simulate observations from the mixture of warping model
#' with a unimodal base shape. The base shape is negative of a standardized beta density
#' function with minimum equals to -1
#' @param n integer, number of curves
#' @param p numeric vector, mixture proportions
#' @param kappa G x K nuermic matrix, means of the Dirichlet distribution
#' for the G mixture component, stacked in rows
#' @param tau numeric, common concentration parameter of the Dirichlet components
#' @param ni numeric, number of points per curve
#' @param mu_sh numeric, mean of amplitude shifting effect
#' @param mu_sc numeric, mean of amplitude scaling effect
#' @param sd_sh numeric, standard deviation of amplitude shifting
#' @param sd_sc numeric, standard deviation of amplitude scaling
#' @param sd_err numeric, standard deviation of the residual term
#' @param shape_a numeric, non-negative parameters of the Beta distribution
#' @param shape_b numeric, non-negative parameters of the Beta distribution
#' @import stats splines gtools
#' @examples
#' \dontrun{
#' require(ggplot2)
#' require(fda)
#' kappa0 <- c(1,2,2,1)
#' kappa0 <- kappa0 / sum(kappa0)
#' kappa1 <- c(0,4,4,0)
#' kappa1 <- kappa1 / sum(kappa1)
#' kappa2 <- c(4,1,1,4)
#' kappa2 <- kappa2 / sum(kappa2)
#' kappa3 <- c(0,0,1,1)
#' kappa3 <- kappa3 / sum(kappa3)
#' kappa4 <- c(1,1,0,0)
#' kappa4 <- kappa4 / sum(kappa4)
#'
#' x0 <- seq(0, 1, length=1001)
#' h0_basis <- create.bspline.basis(c(0, 1), 5)
#' h0_basis_mat <- eval.basis(h0_basis, x0)
#' h1 <- h0_basis_mat %*% cumsum(c(0, kappa0))
#' plot(h1~x0, type="l")
#'
#' dat0 <- sim_warping_mixture(100, rep(1/3, 3),
#'                     rbind(kappa1,
#'                           kappa2,
#'                           kappa3),
#'                     ni = 201,
#'                     tau = 30,
#'                     mu_sh = -25, mu_sc = 500,
#'                     sd_sh = 10, sd_sc=50, sd_err = 5)
#'
#' ggplot(dat0) +
#'   geom_line(aes(x=x, y=y, group=id, col=clust))
#' ggplot(dat0) +
#'   geom_line(aes(x=x, y=warped_x, group=id, col=clust))
#' }
#' @export
sim_warping_mixture <- function(n, p, kappa,
                                tau=1, ni=1001,
                                mu_sh = -15, mu_sc = 500,
                                sd_sh = 5, sd_sc=50, sd_err = 5,
                                shape_a=2, shape_b=2){
  x0 <- seq(0, 1, length=ni)
  n_clust <- length(p)
  clust_lbl <- sample(n_clust, n, prob = p, replace = TRUE)
  n_basis <- ncol(kappa) + 1

  w0 <- kappa[clust_lbl,] * tau
  w <- apply(w0, 1, function(k){rdirichlet(1, k)})
  w <- apply(rbind(0, w), 2, cumsum)

  h_basis_mat <- splines::bs(x = x0,
                             degree = 3,
                             df = n_basis,
                             Boundary.knots = c(0, 1),
                             intercept = TRUE)
  ht <- h_basis_mat %*% w
  y0 <- beta_shape(ht, shape_a, shape_b)

  err <- matrix(rnorm(n*ni, 0, sd_err), ni, n)

  a_sh <- rnorm(n, mu_sh, sd_sh)
  a_sc <- rnorm(n, mu_sc, sd_sc)
  y = t(a_sh + a_sc * t(y0)) + err

  dat <- data.frame(
    id = as.factor(rep(1:n, each=ni)),
    clust = as.factor(rep(clust_lbl, each=ni)),
    sh = rep(c(a_sh), each=ni),
    sc = rep(c(a_sc), each=ni),
    x = rep(x0, n),
    warped_x = c(ht),
    y0 = c(y0),
    y = c(y))
  attr(dat, "w") <- w
  return(dat)
}





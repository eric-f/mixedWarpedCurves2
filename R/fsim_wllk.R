#' Function to compute the mixture of warping likelihood
#'
#' This function computes the likelihood of the warping component
#' for the mixture of warping function. The predicted warping
#' coefficients by stochastic approximation is treated as the data;
#' the likelihood of the mixture of Dirichlet distributions is
#' returned
#' @param fsim_obj : object returned by the mixture of warping model with fixed unimodel template
#' @export
fsim_wllk <- function(fsim_obj){

  kappa_hat <- fsim_obj$pars$kappa_clusters
  p_hat <- fsim_obj$pars$p_clusters
  dw_mat <- apply(sapply(fsim_obj$curves, function(crv){crv$sapprox_w}), 2, diff)
  wllks <- mixture_of_dirichlet_logLik(dw_mat, kappa_hat, p_hat)
  wllk = sum(wllks)
  n <- length(fsim_obj$curves)
  npars <- length(c(kappa_hat)) + length(c(p_hat)) - 1
  out <- data.frame(wllk = wllk,
           wbic = npars * log(n) - 2 * wllk,
           waic = npars * 2 - 2 * wllk)
  cat("LLK, AIC and BIC for fsim_unimodal fit\n")
  return(out)
}

# out <- readRDS("~/org/project/Simulations/module/sim-20180220_unimodal_mixture_brillinger_start_multi_K_3-200x1000/output/unimodal_mixture_brillinger_start_multi_K_3-200x1000_set_000_K=3.rds")
# fsim_obj <- out$saemOut
# fsim_wllk(fsim_sfobj)


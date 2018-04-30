#' Function to compute the mixture of warping likelihood
#'
#' This function computes the likelihood of the warping component
#' for the mixture of warping function. The predicted warping
#' coefficients by stochastic approximation is treated as the data;
#' the likelihood of the mixture of Dirichlet distributions is
#' returned
#' @param md_obj : object returned by mixture_of_dirichlet()
#' @export
mixture_of_dirichlet_wllk <- function(md_obj){

  kappa_hat <- md_obj$alpha_hat
  p_hat <- md_obj$p_hat
  wllk = md_obj$llk
  n <- nrow(md_obj$post_p)
  npars <- length(c(kappa_hat)) + length(c(p_hat)) - 1
  out <- data.frame(wllk = wllk,
           wbic = npars * log(n) - 2 * wllk,
           waic = npars * 2 - 2 * wllk)
  cat("LLK, AIC and BIC for mixture_of_dirichlet fit\n")
  return(out)
}

# out <- readRDS("~/org/project/Simulations/module/sim-20180220_unimodal_mixture_brillinger_start_multi_K_3-200x1000/output/unimodal_mixture_brillinger_start_multi_K_3-200x1000_set_000_K=3.rds")
# fsim_obj <- out$saemOut
# dw_mat <- apply(sapply(fsim_obj$curves, function(crv){crv$sapprox_w}), 2, diff)
# md_obj <- mixture_of_dirichlet(dw_mat, 3)
# out2$llk
# mixture_of_dirichlet_wllk(md_obj)

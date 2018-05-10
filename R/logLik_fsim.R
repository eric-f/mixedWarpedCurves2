#' Observed data log-likelihood of the registration model
#'
#' Compute log-likelihood of the fitted model by sampling based approximation
#' @param saemObj fitted registration model
#' @param method method for computing the log-likelihood
#' @param n_sim number of Monte Carlo sample
#' @param seed seed for random number generator
logLik_fsim <- function(saemObj, method=c("joint_mc", "joint_is", "marginal_is", "mcmc"), n_sim=10, seed=0){
  # kappa (<-> identity)
  tmp_y <- tmp_x <- seq(0, 1, length=1000)
  h_knots <- sort(c(rep(range(saemObj$aux$h_break_points), saemObj$aux$h_order-1), saemObj$aux$h_break_points))
  bhx <- splineDesign(h_knots, tmp_x, saemObj$aux$h_order, tmp_x)
  warping_ols <- lm(tmp_y ~ bhx - 1)
  saemObj$pars$kappa_id <- diff(unname(warping_ols$coefficients))
  # number of cluster
  saemObj$pars$num_clusters = ncol(saemObj$pars$kappa_clusters)
  saemObj$pars$prop_tau = mean(colSums(saemObj$pars$kappa_clusters))
  saemObj$aux$prop_sigma = rev(saemObj$aux$proposal_sigma_history)[1]
  sampling_method = match(method, c("joint_mc", "joint_is", "marginal_is", "mcmc"))[1]
  if(is.na(sampling_method)) sampling_method = 1
  out <- logLik_unimodal(saemObj$pars, saemObj$curves, saemObj$aux, sampling_method, n_sim, seed)
  return(out)
}

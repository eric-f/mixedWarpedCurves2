#' Function to compute the log-likelihood, AIC and BIC for the
#' mixture of dirichlet distribution
#'
#' @param dw_mat k x n matrix of k-dimensional vector of compositional data
#' @param kappa_hat k x M matrix of shape parameters for the M clusters
#' @param p_hat M x 1 vector of mixing proportion
mixture_of_dirichlet_logLik <- function(dw_mat, kappa_hat, p_hat){
  apply(dw_mat, 2, function(dw){
    log(sum(p_hat * apply(kappa_hat, 2, function(kappa){ddirichlet(dw, kappa)})))
  })
}

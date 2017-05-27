#' Auxiliary for Specifying the control parameters for the stochastic approximation
#'
#' This function help specify the control parameters for the stochastic approximation.
#' @param n_saem_iter
#' @param n_saem_burn
#' @param n_mcmc_burn
#' @param n_core
#' @param saem_step_seq_pow
#' @param prop_sigma
#' @param need_centering
#' @param accept_rate_lb
#' @param accept_rate_ub
#' @param accept_rate_window
#' @param f_order
#' @param h_order
#' @param f_knots
#' @param h_knots
#' @param f_n_knots
#' @param h_n_knots
#' @export
control_saem <- function(n_saem_iter=1,
                         n_saem_burn=1,
                         n_mcmc_burn=1,
                         n_core=1,
                         saem_step_seq_pow=1,
                         prop_sigma = 1e-2,
                         need_centering=FALSE,
                         accept_rate_lb=0.17,
                         accept_rate_ub=0.33,
                         accept_rate_window=5,
                         f_order=4,
                         h_order=4,
                         f_knots=NULL,
                         h_knots=NULL,
                         f_n_knots=NULL,
                         h_n_knots=NULL){
  n_saem_iter = max(1, n_saem_iter)
  n_saem_burn = max(1, n_saem_burn)
  n_mcmc_burn = max(1, n_mcmc_burn)
  saem_step_seq_pow = max(0.5, min(saem_step_seq_pow, 1))
  prop_sigma = max(1e-16, prop_sigma)

  if(is.null(f_knots)){
    if(!is.null(f_n_knots)){
      if(f_n_knots <= 0)
        f_knots = c(0, 0.5, 1)
      else
        f_knots = seq(0, 1, length=f_n_knots)
    }
    else
      f_knots = c(0, 0.5, 1)
  }

  if(is.null(h_knots)){
    if(!is.null(h_n_knots)){
      if(h_n_knots <= 0)
        h_knots = c(0, 0.5, 1)
      else
        h_knots = seq(0, 1, length=h_n_knots)
    }
    else{
      h_knots = c(0, 0.5, 1)
    }
  }

  if(any(range(f_knots) != c(0, 1)))
    stop("f_knots needs to be within 0 and 1")
  if(any(range(h_knots) != c(0, 1)))
    stop("h_knots needs to be within 0 and 1")

  list(n_saem_iter = n_saem_iter,
       n_saem_burn = n_saem_burn,
       n_mcmc_burn = n_mcmc_burn,
       n_core = n_core,
       saem_step_seq_pow = saem_step_seq_pow,
       prop_sigma = prop_sigma,
       need_centering = need_centering,
       accept_rate_lb = accept_rate_lb,
       accept_rate_ub = accept_rate_ub,
       accept_rate_window = accept_rate_window,
       f_order = f_order,
       h_order = h_order,
       f_knots = f_knots,
       h_knots = h_knots)
}

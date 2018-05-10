#' Auxiliary for Specifying the control parameters for the stochastic approximation
#'
#' This function help specify the control parameters for the stochastic approximation.
#' @param n_saem_iter Number of EM iterations (in addition to SAEM burn)
#' @param n_saem_burn Number of SAEM burn-in steps
#' @param n_mcmc_burn Number of MCMC steps
#' @param n_core      Number of cores for parallel computing
#' @param saem_step_seq_pow  Power of the harmonic series of SAEM step sizes
#' @param prop_sigma         Initial value for the variance of the MH proposal
#' @param need_centering     Logical variable of whether the sampled amplitude effects should be centered to have sample mean equals to mu
#' @param accept_rate_lb     Lower bound of the desired MH acceptance rate
#' @param accept_rate_ub     Upper bound of the desired MH acceptance rate
#' @param accept_rate_window Number of E-steps between successive calibration of MH proposal
#' @param f_order    Order of the spline for the base shape
#' @param h_order    Order of the splines for the warping functions
#' @param f_knots    Knot locations of the base shape spline (including boundary knots)
#' @param h_knots    Knot locations of the warping function splines (including boundary knots)
#' @param f_n_knots  Number of equally spaced knots, from 0 to 1, to use if f_knots is not specified
#' @param h_n_knots  Number of equally spaced knots, from 0 to 1, to use if h_knots is not specified
#' @param ind_amp    TRUE if the amplitude shifting and scaling effects are assumed independent, FALSE otherwise
#' @param seed       Random seed for MCMC sampler in C++
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
                         h_n_knots=NULL,
                         ind_amp=FALSE,
                         seed=NULL){
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

  if(is.null(seed)){
    seed <- round(runif(1, 1, 2^8))
  }
  else{
    if(!is.integer(seed))
      seed <- round(runif(1, 1, 2^8))
  }

  list(seed = seed,
       n_saem_iter = n_saem_iter,
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
       h_knots = h_knots,
       ind_amp = ind_amp)
}

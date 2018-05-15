#' Model-based curve registration with unknown base shape
#'
#' This function fits the model \deqn{Y_i(t) = a_{i,sh} + a_{i,sc} f \circ h_i(y) + error}
#' by maximum likelihood via a stochastic approximation (SA) EM algorithm. In the model,
#' $f$ is a B-spline representing a common shape whereas \eqn{h_i:[0, 1] \to [0, 1]}
#' is a monotone B-spline representing a random time transformation, referred to as a
#' warping function or registration function. The vector of coefficients of $h_i$ follows a
#' Dirichlet distributions.
#' @param y vector of observed curves at times obs_time
#' @param obs_time vector of the observation times
#' @param curve_id vector of curve IDs
#' @param n_clust number of clusters. \emph{(Deprecated!! Please use fsim_mixed_warped_curves() for simultaneous registration and clustering)}
#' @param saem_control list of values to control MCMC and stochastic approximation
#' @param trace if TRUE, tracing information of the estimated parameters are printed
#' @return \describe{
#'   \item{pars}{List of estimated or fixed model parameters
#'     \describe{
#'       \item{mu0}{Fixed mean vector for the Gaussian amplitude effects}
#'       \item{kappa0}{Fixed Dirichlet mean so that the mean warping function is the identity function}
#'       \item{alpha}{Estimated B-spline basis coefficient for the common base shape}
#'       \item{sigma2}{Estimated variance of the error term}
#'       \item{big_sigma}{Estimated variance-covariance matrix of the Gaussian amplitude effect}
#'       \item{p_clusters}{Estimated mixing proportion. Identical to 1 in the no clustering case}
#'       \item{tau_clusters}{Overall concentration of the Dirichlet random effect distribution}
#'       \item{kappa_clusters}{Estimated means of the Dirichlet components. Identical to kappa0 in the no clustering case}
#'     }
#'   }
#'   \item{curves}{List of curves with stochastic approximation to sufficient statistics, each curves has the following components
#'     \describe{
#'       \item{curve_id}{Curve ID. (Caution: this might be different from the inputted curve id, if the original id's is not a sequence from 1 to n.}
#'       \item{x}{Inputted observation time}
#'       \item{y}{Inputted observed curves}
#'       \item{warped_x}{Estimated warped time}
#'       \item{fitted_y}{Fitted curve}
#'       \item{sapprox_a}{Stochastic approximation to the conditional expectation of amplitude effects given data}
#'       \item{sapprox_w}{Stochastic approximation to the conditional expectation of warping coefficients given data}
#'       \item{sapprox_warped_f_basis_mat}{Stochastic approximation to sufficient statistics for SAEM}
#'       \item{sapprox_aug_warped_f_basis_m}{Stochastic approximation to sufficient statistics for SAEM}
#'       \item{sapprox_hat_mat}{Stochastic approximation to sufficient statistics for SAEM}
#'       \item{sapprox_sigma_a}{Stochastic approximation to sufficient statistics for SAEM}
#'       \item{sapprox_log_dw}{Stochastic approximation to sufficient statistics for SAEM}
#'       \item{sapprox_cluster_membership}{(Deprecated!!) Stochastic approximation to predictive probabilities of cluster membership}
#'     }
#'   }
#'   \item{aux}{List of auxiliary information and intermediate variables for MCMC-SAEM}
#'   \item{pars_track}{Sequence of estimated parameters for convergence diagnostics}
#'   \item{se_info}{(Experimental!!) Stochastic approximation to Fishers information matrix}
#'   \item{y_scaling_factor}{Maximum absolute value of the observed curve}
#' }
#' @references Fu, E. and Heckman, N. (2017). Model-based curve registration via stochastic approximation EM algorithm. https://arxiv.org/abs/1712.07265
#' @useDynLib mixedWarpedCurves2
#' @importFrom splines splineDesign
#' @export
fsim <- function(y,
                 obs_time,
                 curve_id,
                 n_clust=1,
                 saem_control = control_saem(),
                 trace=FALSE){
  ## --------------------------------------------------------------------------
  ## Scale reponses and pack data into data.frame
  ## --------------------------------------------------------------------------
  if(min(obs_time)!=0 | max(obs_time)!=1){
    stop("observation time needs to be within 0 and 1")
  }
  y_scaling_factor <- max(abs(y))
  data <- data.frame(y = y / y_scaling_factor,
                     x = obs_time,
                     id = curve_id)

  ## --------------------------------------------------------------------------
  ## Initialize model parameters ----------------------------------------------
  # --------------------------------------------------------------------------
  pars <- NULL
  # mu (fixed)
  pars$mu <- c(0, 1)
  # kappa (fixed)
  tmp_y <- tmp_x <- seq(0,1,length=1000)
  h_knots <- sort(c(rep(range(saem_control$h_knots), saem_control$h_order-1), saem_control$h_knots))
  bhx <- splineDesign(h_knots, tmp_x, saem_control$h_order, rep(0, 1000))
  warping_ols <- lm(tmp_y~bhx-1, data=data)
  pars$kappa <- diff(unname(warping_ols$coefficients))
  # alpha and sigma2
  f_knots <- sort(c(rep(range(saem_control$f_knots), saem_control$f_order-1), saem_control$f_knots))
  bfx <- splineDesign(f_knots, data$x, saem_control$f_order, rep(0, length(data$x)))
  shape_ols <- lm(y~bfx-1, data=data)
  pars$alpha <- unname(shape_ols$coefficients)
  pars$sigma2 <- var(shape_ols$residuals)
  pars$num_clusters <- n_clust
  # big_sigma and tau to be initialized in cpp...

  ## --------------------------------------------------------------------------
  ## Initialize common auxiliary objects --------------------------------------
  ## --------------------------------------------------------------------------
  saem_control$n_total = length(curve_id)
  saem_control$n_curve = length(unique(curve_id))

  ## --------------------------------------------------------------------------
  ## Convert data to list of curves -------------------------------------------
  ## --------------------------------------------------------------------------
  data_lst <- split(data, data$id)

  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## Iterate SAEM -------------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  out <- saem_fit(data_lst, pars, saem_control, y_scaling_factor, trace)
  out$y_scaling_factor = y_scaling_factor

  return(out)
}

#' Model-based curve registration and clustering with unknown base shape
#'
#' This function fits the model \deqn{Y_i(t) = a_{i,sh} + a_{i,sc} f \circ h_i(y) + error}
#' by maximum likelihood via a stochastic approximation EM algorithm. In the model,
#' $f$ is a B-spline representing a common shape whereas \eqn{h_i:[0, 1] \to [0, 1]}
#' is a monotone B-spline representing a random time transformation, referred to as a
#' warping function or registration function. The vector of (warping) coefficients of $h_i$
#' follows a Dirichlet distributions. The function also supports a finite mixture of
#' Dirichlet distributions for the warping coefficient for simultaneous curve clustering and
#' registration.
#' @param y vector of observed curves
#' @param obs_time vector of the observation times
#' @param curve_id vector of curve IDs
#' @param init_clust vector of inital clutering label with length equals to the number of curves
#' @param n_clust integer, number of clusters (set to 1 if not clustering)
#' @param saem_control a list of values to control the MCMC and stochastic approximation. See control_saem().
#' @param trace if TRUE, tracing information of the estimated parameters are printed
#' @return \describe{
#'   \item{pars}{List of estimated or fixed model parameters
#'     \describe{
#'       \item{alpha}{Estimated B-spline basis coefficient for the common base shape}
#'       \item{mu_a}{Fixed mean vector for the Gaussian amplitude effects}
#'       \item{sigma2}{Estimated error variance}
#'       \item{sigma2_a}{Estimated variance-covariance matrix of the Gaussian amplitude effect}
#'       \item{p_clusters}{Estimated mixing proportion}
#'       \item{kappa_id}{Fixed Dirichlet mean vector of the first component where the mean warping function is the identity function}
#'       \item{tau1}{Overall concentration of the Dirichlet distribution for the first component}
#'       \item{kappa_clusters}{Estimated Dirichlet concentration parameters}
#'     }
#'   }
#'   \item{curves}{List of curves with stochastic approximation to sufficient statistics, each curves has the following components
#'     \describe{
#'       \item{curve_id}{Curve ID. (Caution: this might be different from the inputted curve id, if the original id's is not a sequence from 1 to n.}
#'       \item{x}{Inputted observation time}
#'       \item{y}{Inputted observed curves}
#'       \item{y}{Inputted or random initial cluster label}
#'       \item{warped_x}{Estimated warped time}
#'       \item{fitted_y}{Fitted curve}
#'       \item{sapprox_residual_sum_of_squares}{Stochastic approximation to residual sum of squares}
#'       \item{sapprox_a}{Stochastic approximation to the conditional expectation of amplitude effects given data}
#'       \item{sapprox_w}{Stochastic approximation to the conditional expectation of warping coefficients given data}
#'       \item{sapprox_log_dw}{Stochastic approximation to sufficient statistics for SAEM}
#'       \item{sapprox_cluster_membership}{Stochastic approximation to predictive probabilities of cluster membership}
#'     }
#'   }
#'   \item{aux}{List of auxiliary information and intermediate variables for MCMC-SAEM}
#'   \item{pars_track}{Sequence of estimated parameters for convergence diagnostics}
#'   \item{se_info}{Not currently implemented}
#'   \item{y_scaling_factor}{Maximum absolute value of the observed curve}
#' }
#' @seealso See https://github.com/eric-f/mixedWarpedCurves2 for examples
#' @references Fu, E. and Heckman, N. (2017). Model-based curve registration via stochastic approximation EM algorithm. https://arxiv.org/abs/1712.07265
#' @useDynLib mixedWarpedCurves2
#' @importFrom splines splineDesign
#' @export
fsim_mixed_warped_curves <- function(y,
                                     obs_time,
                                     curve_id,
                                     init_clust=NULL,
                                     n_clust=1,
                                     saem_control = control_saem(),
                                     trace=FALSE){
  ## --------------------------------------------------------------------------
  ## Scale reponses and pack data into a data.frame
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
  ## --------------------------------------------------------------------------
  pars <- NULL
  # mu (fixed)
  pars$mu <- c(0, 1)
  # kappa (<-> identity)
  tmp_y <- tmp_x <- seq(0, 1, length=1000)
  h_knots <- sort(c(rep(range(saem_control$h_knots), saem_control$h_order-1), saem_control$h_knots))
  bhx <- splineDesign(h_knots, tmp_x, saem_control$h_order, rep(0, 1000))
  warping_ols <- lm(tmp_y ~ bhx - 1, data=data)
  pars$kappa <- diff(unname(warping_ols$coefficients))
  pars$kappa <- pars$kappa / sum(pars$kappa)
  # f and sigma2
  f_knots <- sort(c(rep(range(saem_control$f_knots), saem_control$f_order-1), saem_control$f_knots))
  bfx <- splineDesign(f_knots, data$x, saem_control$f_order, rep(0, length(data$x)))
  shape_ols <- lm(y ~ bfx - 1, data=data)
  pars$alpha <- unname(shape_ols$coefficients)
  pars$sigma2 <- var(shape_ols$residuals)
  # big_sigma to be initialized in C++...
  # tau_1 and (tau_2, kappa_2), ..., (tau_M, kappa_M)
  # to be initialized in C++...
  pars$num_clusters <- n_clust

  ## --------------------------------------------------------------------------
  ## Initialize common auxiliary objects --------------------------------------
  ## --------------------------------------------------------------------------
  saem_control$n_total = length(curve_id)
  saem_control$n_curve = length(unique(curve_id))

  ## --------------------------------------------------------------------------
  ## Convert data to list of curves -------------------------------------------
  ## --------------------------------------------------------------------------
  data_lst <- split(data, data$id)
  n_curve <- length(data_lst)
  if(is.null(init_clust)){
    print("Randomizing initial cluster labels...")
    init_clust <- sample(n_clust, n_curve, replace=TRUE)
  }
  else{
    if(length(init_clust) != n_curve |
       any(is.na(init_clust)) |
       min(init_clust, na.rm = T) < 1 |
       max(init_clust, na.rm = T) > n_clust){
      print("Invalid initial cluster labels. Will use random clustering configuration as starting value...")
      init_clust <- sample(n_clust, n_curve, replace=TRUE)
    }
  }
  for(i in seq(along=data_lst)){
    data_lst[[i]]$init_clust = init_clust[i]
  }

  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## Run SAEM in C++ ----------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  out <- saem_fit_mixed_warped_curves(data_lst, pars, saem_control,
                                      y_scaling_factor, trace)
  out$y_scaling_factor = y_scaling_factor

  return(out)
}

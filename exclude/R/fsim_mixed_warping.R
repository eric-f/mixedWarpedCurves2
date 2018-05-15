#' Fitting shape-invariant model with flexible time transformations that follows
#' a finite mixture of dirichlet distribution
#' (obsolete implementation, see fsim_mixed_warped_curves)
#'
#' This function fits the model \deqn{Y_i(t) = a_{i,sh} + a_{i,sc} f \circ h_i(y)}
#' by maximum likelihood via a stochastic approximation EM algorithm. In the model,
#' $f$ is a B-spline representing a common shape whereas \eqn{h_i:[0, 1] \to [0, 1]}
#' is a monotone B-spline representing a random time transformation
#' @param y a vector of observed curves
#' @param obs_time a vector of the observation times
#' @param curve_id a vector of curve IDs
#' @param init_clust vector of inital clutering label with length equals to the number of curves
#' @param n_clust integer, number of clusters
#' @param saem_control a list of values to control the MCMC and stochastic approximation. See control_saem().
#' @param trace if TRUE, tracing information of the estimated parameters are printed
#' @importFrom Rcpp evalCpp
#' @noRd
#' @useDynLib mixedWarpedCurves2
fsim_mixed_warping <- function(y,
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
  if(is.null(init_clust) |
     length(init_clust) != n_curve |
     min(init_clust) < 1 |
     max(init_clust) > n_clust){
    print("Randomizing initial cluster labels...")
    init_clust <- sample(n_clust, n_curve, replace=TRUE)
  }
  for(i in seq(along=data_lst)){
    data_lst[[i]]$init_clust = init_clust[i]
  }

  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## Run SAEM in C++ ----------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  out <- saem_fit_mixed_warping(data_lst, pars, saem_control,
                                y_scaling_factor, trace)
  out$y_scaling_factor = y_scaling_factor

  return(out)
}

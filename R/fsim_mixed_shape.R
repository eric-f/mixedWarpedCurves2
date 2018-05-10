#' Fitting mixture of shapes model with amplitude shifting and scaling by SAEM
#'
#' This function fits the model \deqn{Y_i(t) = a_{i,sh} + a_{i,sc} f_{m_i}(t)}
#' by maximum likelihood via a stochastic approximation EM algorithm. In the model,
#' $f_m$ is a B-spline representing a common shape for cluster $m$. This model assume
#' no warpings.
#' @param y a vector of observed curves
#' @param obs_time a vector of the observation times
#' @param curve_id a vector of curve IDs
#' @param init_clust vector of inital clutering label with length equals to the number of curves
#' @param n_clust integer, number of clusters
#' @param saem_control a list of values to control the MCMC and stochastic approximation. See control_saem().
#' @param trace if TRUE, tracing information of the estimated parameters are printed
#' @importFrom splines splineDesign
#' @importFrom Rcpp evalCpp
#' @useDynLib mixedWarpedCurves2
fsim_mixed_shape <- function(y,
                             obs_time,
                             curve_id,
                             init_clust=NULL,
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
                     id = curve_id,
                     init_clust = init_clust)

  ## --------------------------------------------------------------------------
  ## Initialize model parameters ----------------------------------------------
  # --------------------------------------------------------------------------
  pars <- NULL
  # mu (fixed)
  pars$mu <- c(0, 1)
  # alpha and sigma2
  f_knots <- sort(c(rep(range(saem_control$f_knots), saem_control$f_order-1), saem_control$f_knots))
  pars$alpha <- sapply(split(data, data$init_clust), function(sub_data){
    bfx <- splineDesign(f_knots, sub_data$x, saem_control$f_order, rep(0, length(sub_data$x)))
    shape_ols <- lm(y~bfx-1, data=sub_data)
    return(unname(shape_ols$coefficients))
  })
  bfx <- splineDesign(f_knots, data$x, saem_control$f_order, rep(0, length(data$x)))
  shape_ols <- lm(y~bfx-1, data=data)
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
  print("Iterate SAEM")
  out <- saem_fit_mixed_shape(data_lst, pars, saem_control, y_scaling_factor, trace)
  out$y_scaling_factor = y_scaling_factor

  return(out)
}

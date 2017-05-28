#' Fitting shape-invariant model with flexible time transformations
#'
#' This function fits the model \deqn{Y_i(t) = a_{i,sh} + a_{i,sc} f \circ h_i(y)}
#' by maximum likelihood via a stochastic approximation EM algorithm. In the model,
#' $f$ is a B-spline representing a common shape whereas \eqn{h_i:[0, 1] \to [0, 1]}
#' is a monotone B-spline representing a random time transformation
#' @param y a vector of observed curves
#' @param obs_time a vector of the observation times
#' @param curve_id a vector of curve IDs
#' @param init_pars a list of values to initialize parameters
#' @param pars.control a list of values to control which parameters will be treated as known
#' @param basis.control a list of values to specify the B-splines for the common shape and the warping functions
#' @param sa.control a list of values to control the MCMC and stochastic approximation
#' @param track_pars a logical variable for storing the whole sequence of estimates parameters.
#' @useDynLib mixedWarpedCurves
#' @export
fsim <- function(y,
                 obs_time,
                 curve_id,
                 saem_control = control_saem()){
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
  ## --------------------------------------------------------------------------
  pars <- NULL
  # mu (fixed)
  pars$mu <- c(0, 1)
  # kappa (fixed)
  tmp_y <- tmp_x <- seq(0,1,length=1000)
  h_knots <- sort(c(rep(range(saem_control$h_knots), saem_control$h_order-1), saem_control$h_knots))
  bhx <- splineDesign(h_knots, tmp_x, saem_control$h_order, 0)
  warping_ols <- lm(tmp_y~bhx-1, data=data)
  pars$kappa <- diff(unname(warping_ols$coefficients))
  # alpha and sigma2
  f_knots <- sort(c(rep(range(saem_control$f_knots), saem_control$f_order-1), saem_control$f_knots))
  bfx <- splineDesign(f_knots, data$x, saem_control$f_order, 0)
  shape_ols <- lm(y~bfx-1, data=data)
  pars$alpha <- unname(shape_ols$coefficients)
  pars$sigma2 <- var(shape_ols$residuals)
  # big_sigma and tau to be initialized in cpp...

  ## --------------------------------------------------------------------------
  ## Initialize common auxiliary objects --------------------------------------
  ## --------------------------------------------------------------------------
  cat("  Initialize common auxiliary objects ... ")
  saem_control$n_total = length(curve_id)
  saem_control$n_curve = length(unique(curve_id))
  cat("(done)\n")

  ## --------------------------------------------------------------------------
  ## Convert data to list of curves -------------------------------------------
  ## --------------------------------------------------------------------------
  cat("  Converting data to list of curves ... ")
  data_lst <- split(data, data$id)
  cat("(done)\n")

  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## Iterate SAEM -------------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  cat("Stochastic approximation EM algorithm ... \n")
  out <- saem_fit(data_lst, pars, saem_control, y_scaling_factor)
  out$y_scaling_factor = y_scaling_factor


  return(out)
}

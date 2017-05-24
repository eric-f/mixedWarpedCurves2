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
                 init_pars = NULL,
                 pars.control = control_model_param(),
                 basis.control = control_bspline_basis(),
                 sa.control = control_sa(),
                 track_pars = FALSE){
  ## --------------------------------------------------------------------------
  ## Scale reponses and pack data into data.frame
  ## --------------------------------------------------------------------------
  mean_y <- mean(y)
  sd_y <- sd(y)
  data <- data.frame(y = (y - mean_y)  / sd_y,
                     x = obs_time,
                     id = curve_id)

  ## --------------------------------------------------------------------------
  ## Create B-spline basis objects --------------------------------------------
  ##  - using create.bspline.basis(fda) and project.basis(fda) ----------------
  ## --------------------------------------------------------------------------
  cat("Generating f_basis from basis.control  ... \n")
  f_basis <-
    create.bspline.basis(
      rangeval = basis.control$f.rangeval,
      norder = basis.control$f.degree + 1,
      nbasis = basis.control$f.degree + 1 + basis.control$f.nknots)
  cat("Generating h_basis from basis.control  ... \n")
  h_basis <-
    create.bspline.basis(
      rangeval = basis.control$h.rangeval,
      norder = basis.control$h.degree + 1,
      nbasis = basis.control$h.degree + 1 + basis.control$h.nknots)
  cat("Generating default init_pars ... \n")
  if (is.null(init_pars) || is.null(init_pars$mu)) {
    init_pars$mu <- c(0, 1)
  }
  if (is.null(init_pars) || is.null(init_pars$Sigma)) {
    init_pars$Sigma <- 1e-2 * diag(2)
  }
  if (is.null(init_pars) || is.null(init_pars$kappa)){
    tmp_y <- tmp_x <- seq(0,1,length=1000)
    init_pars$kappa <- diff(c(project.basis(tmp_y, tmp_x, h_basis)))
  }
  if (is.null(init_pars) || is.null(init_pars$tau)){
    init_pars$tau <- 1e3
  }
  if (is.null(init_pars) || is.null(init_pars$alp)) {
    init_pars$alp <- project.basis(data$y, data$x, basisobj = f_basis)
  }


  ## --------------------------------------------------------------------------
  ## Initialize model parameters ----------------------------------------------
  ## --------------------------------------------------------------------------
  cat("Initializing...\n")
  pars <- init_pars  ## need to check if all parameters are provided and valid


  ## --------------------------------------------------------------------------
  ## Initialize common auxiliary objects --------------------------------------
  ## pkg func invoked: solve_exact, generate_L_matrix
  ## extern func invoked: fd
  ## --------------------------------------------------------------------------
  cat("  Initialize common auxiliary objects ... ")
  pars$aux <- list(
    n_total = nrow(data),
    n_curve = length(unique(data$id)),
    f_basis = f_basis,
    f_order = basis.control$f.degree + 1,
    f_full_knots = c(rep(f_basis$rangeval[1], basis.control$f.degree+1),
                     f_basis$params,
                     rep(f_basis$rangeval[2], basis.control$f.degree+1)),
    f_break_points = c(f_basis$rangeval[1], f_basis$params, f_basis$rangeval[2]),
    h_basis = h_basis,
    Sigma_inv = solve_exact(pars$Sigma),
    D = length(pars$kappa),
    L = generate_L_matrix(length(pars$kappa))
  )
  cat("(done)\n")


  ## --------------------------------------------------------------------------
  ## Convert data to list of curves -------------------------------------------
  ## entern func invoked: dlply, fd, predict.fd
  ## --------------------------------------------------------------------------
  cat("  Converting data to list of curves ... ")
  data_obj_lst <-
    dlply(data, .(id),
          function(curve, f_basis, h_basis){
            ## Generate B-spline basis matrix for the warping functions
            h_fd <- fd(diag(h_basis$nbasis), h_basis)
            h_basis_mat <- predict(h_fd, curve$x)
            return(list(data = curve, h_basis_mat = h_basis_mat))
          },
          f_basis=f_basis, h_basis=h_basis)
  cat("(done)\n")


  ## --------------------------------------------------------------------------
  ## Initialize Markov Chain of the random effect -----------------------------
  ## pkg func invoked: get_warped_time, get_warped_shape
  ## entern func invoked: llply
  ## --------------------------------------------------------------------------
  cat("  Initialize Markov chain of random effects ... ")
  data_obj_lst <-
    llply(data_obj_lst,
          function(data_obj, pars){
            curr_warped_time = get_warped_time(cumsum(c(0, pars$kappa)), data_obj$h_basis_mat)
            curr_warped_f = get_warped_shape(curr_warped_time, pars)
            re_mc_state = c(pars$mu, cumsum(c(0, pars$kappa)))
            return(c(data_obj,
                     list(re_mc_state = re_mc_state,
                          curr_warped_time = curr_warped_time,
                          curr_warped_f = curr_warped_f)))
          },
          pars = pars)
  cat("(done)\n")


  ## --------------------------------------------------------------------------
  ## Initialize sufficient statistics -----------------------------------------
  ## pkg func invoked: initialize_stat
  ## entern func invoked: llply
  ## --------------------------------------------------------------------------
  cat("  Initialize sufficient statistics ... ")
  data_obj_lst <-
    llply(data_obj_lst, .fun = initialize_stat,
          pars = pars,
          sa.control = sa.control,
          .parallel = sa.control$nCore>1)
  cat("(done)\n")


  ## --------------------------------------------------------------------------
  ## Initialize basis coef. for amp. func. and var. of error term -------------
  ## pkg func invoked: maximization_step, update_DG_aux
  ## --------------------------------------------------------------------------
  cat("  Initialize basis coef. for amp. func. and var. of error term ... ")
  out <- maximization_step(pars = pars,
                           data_obj_lst = data_obj_lst,
                           pars.control = control_model_param(
                             fixed.alp = FALSE, fixed.sigma2 = FALSE,
                             fixed.mu = TRUE, fixed.Sigma = TRUE,
                             fixed.kappa = TRUE, fixed.tau = TRUE))
  pars <- out$pars
  pars <- aux_m_step(pars)
  cat("(done)\n")


  ## --------------------------------------------------------------------------
  ## Initialize stochastic approximates for SE --------------------------------
  ## pkg func invoked: initialize_complete_data_gradient_and_hessian
  ## --------------------------------------------------------------------------
  infoMat <- initialize_complete_data_gradient_and_hessian(pars=pars, pars.control=pars.control)


  ## --------------------------------------------------------------------------
  ## Step size schedule -------------------------------------------------------
  ## --------------------------------------------------------------------------
  gks <- c(rep(1, sa.control$nBurnSAEM),
           1 / (1:sa.control$nIter)) ^ sa.control$alphaSAEM


  ## --------------------------------------------------------------------------
  ## Allocation memeory for tracking objects ----------------------------------
  ## --------------------------------------------------------------------------
  if (track_pars){
    pars_track <- list(
      alp = array(NA, c(sa.control$nIter + sa.control$nBurnSAEM, length(init_pars$alp))),
      sigma2 = array(NA, c(sa.control$nIter + sa.control$nBurnSAEM)),
      Sigma = array(NA, c(sa.control$nIter + sa.control$nBurnSAEM, 2, 2)),
      tau = array(NA, c(sa.control$nIter + sa.control$nBurnSAEM))
    )
    sa_track <- array(NA, c((sa.control$nIter + sa.control$nBurnSAEM),
                            pars$aux$n_curve,
                            length(data_obj_lst[[1]]$re_s_approx)))
    mcmc_track <- array(NA, c((sa.control$nIter + sa.control$nBurnSAEM),
                              pars$aux$n_curve,
                              length(data_obj_lst[[1]]$re_mc_state)))
    score_track <- array(NA, c((sa.control$nIter + sa.control$nBurnSAEM),
                               length(pars$alp) +
                                 sum(upper.tri(pars$Sigma, TRUE)) +
                                 length(pars$tau) +
                                 length(pars$sigma)))
  }
  else{
    pars_track <- sa_track <- mcmc_track <- score_track <- NULL
  }
  obs_llk_update_track <- numeric(sa.control$nIter + sa.control$nBurnSAEM)
  obs_llk_track <-numeric(sa.control$nIter + sa.control$nBurnSAEM)




  ## --------------------------------------------------------------------------
  ## Initialize acceptance rate trackers
  ## --------------------------------------------------------------------------
  sa.control$accept_rate <- numeric(sa.control$nIter + sa.control$nBurnSAEM)
  sa.control$past_accept_rates <- numeric(sa.control$n_accept_rates)
  sa.control$prop_sigma_track <- numeric(2*sa.control$nBurnSAEM)
  sa.control$ar_counter <- 1

  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## Iterate SAEM -------------------------------------------------------------
  ## --------------------------------------------------------------------------
  ## --------------------------------------------------------------------------
  cat("Stochastic approximation EM algorithm ... \n")

  # pars: mu, Sigma, kappa, tau, alp, aux, sigma2
  saem_rcpp(curve_list=data_obj_lst,
            pars=pars,
            sa_control=sa.control,
            pars_control=pars.control,
            step_size=gks,
            infoMat=infoMat,
            obs_llk_track=obs_llk_track,
            pars_track=pars_track)

  cat("(done)\n")



  ## --------------------------------------------------------------------------
  ## Combine to get approximated Fisher information matrix --------------------
  ## --------------------------------------------------------------------------
  infoMat$I <- - (infoMat$SA_H - infoMat$SA_C + infoMat$SA_G %*% t(infoMat$SA_G))

  ## --------------------------------------------------------------------------
  ## Return
  ## --------------------------------------------------------------------------
  return(list(pars = pars,
              data_obj_lst = data_obj_lst,
              data = data,
              init_pars = init_pars,
              infoMat = infoMat,
              ## --------------------------------------------------------------
              pars_track = pars_track,
              sa_track = sa_track,
              mcmc_track = mcmc_track,
              score_track = score_track,
              prop_sigma_track = sa.control$prop_sigma_track,
              ## --------------------------------------------------------------
              obs_llk_track = obs_llk_track,
              obs_llk_update_track = obs_llk_update_track,
              ## --------------------------------------------------------------
              input = list(init_pars = init_pars,
                           f_basis = f_basis,
                           h_basis = h_basis,
                           pars.control = pars.control,
                           basis.control = basis.control,
                           sa.control = sa.control,
                           gks = gks)
  ))
}

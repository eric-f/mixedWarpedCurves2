#' Function to compute the mixture of warping likelihood
#'
#' This function computes the likelihood of the warping component
#' for the mixture of warping function. The predicted warping
#' coefficients by stochastic approximation is treated as the data;
#' the likelihood of the mixture of Dirichlet distributions is
#' returned
#' @param y : n x p matrix of the n curves
#' @param bs97_obj : object returned by mixed_shape that fits Brillinger's model
#' @export
mixed_shape_wllk <- function(y, bs97_obj){
  n <- nrow(bs97_obj$p_jk)
  npars <- length(c(bs97_obj$a_k)) +
    length(c(bs97_obj$s2)) +
    length(c(bs97_obj$p_k))
  a_k <- bs97_obj$a_k * bs97_obj$y_max
  s2 <- bs97_obj$s2 * bs97_obj$y_max^2
  p_k <- bs97_obj$p_k
  mllk <- sum(apply(y, 2, function(crv){
    llk_k <- apply(a_k, 2, function(a){-sum((crv-a)^2/2/s2)})
    max_llk_k <- max(llk_k)
    ## Wrong implementation!!
    ## Need to have the varaince term!!!
    llk <- log(sum(exp(llk_k - max_llk_k) * p_k)) + max_llk_k +
      -sum(log(s2))/2
  }))
  out <- data.frame(mllk = mllk,
                    bic = npars * log(n) - 2 * mllk,
                    aic = npars * 2 - 2 * mllk)
  cat("LLK, AIC and BIC for Brillinger's model fit\n")
  return(out)
}


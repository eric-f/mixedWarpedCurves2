#' Fit the mixture model of Brillinger and Stewart (1997)
#'
#' This function fits the mixture of shape model using the EM algorithm
#' described in Brillinger and Stewart (1997). This is a multi-try
#' version of fit_mixed_shape()
#' @param y matrix, numeric matrix containing depth. Each column represents one dive.
#' @param x numeric, common time grid at which the depths are recorded
#' @param nCat integer, number of dive categories
#' @param nIter integer, number of EM iterations
#' @param nTry integer, number of multi-random start
#' @references Brillinger, D. R. and Stewart, B. S. (1997). Elephant seal movements: dive types and their sequences. In Modelling Longitudinal and Spatially Correlated Data (pp. 275–288). Springer.
#' @export
mixed_shape <- function(y, x=NULL, nCat=9,
                        nIter=100, nTry=10){
  iTry <- 1
  aux <- numeric(nTry)
  out <- fit_mixed_shape(y, x, nCat, nIter, plot.est=FALSE)
  aux[iTry] <- out$aux
  for(l in seq(along=numeric(nTry-1))){
    iTry <- iTry + 1
    fit <- fit_mixed_shape(y, x, nCat, nIter, plot.est=FALSE)
    aux[iTry] <- fit$aux
    # if(fit$mlik > out$mlik)
    if(fit$aux > out$aux)
      out <- fit
  }
  out$aux_ntry <- aux
  return(out)
}

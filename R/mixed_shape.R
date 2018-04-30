#' Fit the mixture model of Brillinger and Stewart (1997)
#'
#' This function fits the mixture of shape model using the EM algorithm
#' described in Brillinger and Stewart (1997).
#' @param y matrix, numeric matrix containing depth. Each column represents one dive.
#' @param x numeric, common time grid at which the depths are recorded
#' @param nCat integer, number of dive categories
#' @param nIter integer, number of EM iterations
#' @param nTry integer, number of multi-random start
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

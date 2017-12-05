devtools::use_package("gtools")

#' Fit the mixture model of Brillinger and Stewart (1997)
#'
#' This function fits the mixture of shape model using the EM algorithm
#' described in Brillinger and Stewart (1997).
#' @param y matrix, numeric matrix containing depth. Each column represents one dive.
#' @param x numeric, common time grid at which the depths are recorded
#' @param nCat integer, number of dive categories
fit_mixed_shape <- function(y, x=NULL, nCat=9, nIter=100, plot.est=FALSE){
  nx <- nrow(y)
  nid <- ncol(y)
  a_k <- matrix(NA, nrow=nx, ncol=nCat)
  p_k <- numeric(nCat)
  mlik_track <- rep(0, nIter)
  # initialize p_jk
  p_jk <- rdirichlet(nid, alpha = rep(1,nCat))
  for(i in 1:nIter){
    a_k <- y %*% p_jk %*% diag(1/colSums(p_jk))
    s2 <- rep(0, nx)
    for(j in 1:nid)
      s2 <- s2 + c((y[,j] - a_k)^2 %*% p_jk[j,]) / nid
    p_k <- colMeans(p_jk)
    for(j in 1:nid){
      z <- -colSums((y[,j] - a_k)^2/s2)/2
      z <- z - max(z) + 400
      tmp_p <- p_k * exp(z)
      tmp_p <- tmp_p / sum(tmp_p)
      p_jk[j,] <- tmp_p
    }
    if (plot.est)
      matplot(a_k, x=x, type="l")
    for(k in 1:nid){
      mlik_k <- log(sum(exp(colSums(-(y[,k] - a_k)^2 / s2 / 2)) * p_k)) -
        log(2*pi)/2 - sum(log(s2))/2
      mlik_track[i] <- mlik_track[i] + mlik_k / nid
    }
  }
  mlik <- mlik_track[nIter]
  nPars <- length(a_k) + length(s2) + length(p_k)
  BIC <- -2*mlik + nPars*log(length(y))/nid
  AIC <- -2*mlik + 2*nPars/nid
  return(list(a_k=a_k, s2=s2, p_k=p_k, p_jk=p_jk, mlik=mlik,
              AIC=AIC, BIC=BIC, mlik_track=mlik_track))
}




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
                        nIter=30, nTry=10){
  out <- fit_mixed_shape(y, x, nCat, nIter, plot.est=FALSE)
  for(l in seq(along=numeric(nTry-1))){
    fit <- fit_mixed_shape(y, x, nCat, nIter, plot.est=FALSE)
    if(fit$mlik > out$mlik)
      out <- fit
  }
  return(out)
}

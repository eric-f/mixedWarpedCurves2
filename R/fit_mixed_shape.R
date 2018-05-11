#' Fit the mixture model of Brillinger and Stewart (1997)
#'
#' This function fits the mixture of shape model using the EM algorithm
#' described in Brillinger and Stewart (1997).
#' @param y matrix, numeric matrix containing depth. Each column represents one dive.
#' @param x numeric, common time grid at which the depths are recorded
#' @param nCat integer, number of dive categories
#' @param nIter integer, number of maximum EM iterations
#' @param plot.est if TRUE, plot estimated shapes
#' @importFrom graphics matplot
fit_mixed_shape <- function(y, x=NULL, nCat=9, nIter=100, plot.est=FALSE){
  y_max <- max(abs(y))
  y = y / y_max
  nx <- nrow(y)
  nid <- ncol(y)
  a_k <- matrix(NA, nrow=nx, ncol=nCat)
  p_k <- numeric(nCat)
  mlik_track <- rep(0, nIter)
  aux_track <- numeric(nIter)
  # Initialize p_jk
  p_jk <- matrix(runif(nCat*nid), nid, nCat)
  p_jk <- p_jk / rowSums(p_jk)
  for(i in 1:nIter){
    ## Fitted shapes
    a_k <- y %*% p_jk %*% diag(1/colSums(p_jk))
    ## Estimate common variance
    s2 <- numeric(nx)
    for(j in 1:nid)
      s2 <- s2 + c((y[,j] - a_k)^2 %*% p_jk[j,]) / nid
    ## Update p_jk
    p_k <- colMeans(p_jk)
    for(j in 1:nid){
      ## Sum of squares
      z <- -colSums((y[,j] - a_k)^2/s2)/2
      z <- z - max(z)
      ## Joint density
      tmp_p <- p_k * exp(z)
      ## Posterior density
      p_jk[j,] <- tmp_p / sum(tmp_p)
    }
    if (plot.est)
      matplot(a_k, x=x, type="l")
    # ## Track marginal likelihood
    # for(k in 1:nid){
    #   mlik_k <- log(sum(exp(colSums(-(y[,k] - a_k)^2 / s2 / 2)) * p_k)) -
    #     log(2*pi)/2 - sum(log(s2))/2
    #   mlik_track[i] <- mlik_track[i] + mlik_k / nid
    # }
    ## Track complete data log-likelihood
    for(k in 1:nid){
      aux_tmp <- p_jk[k,] %*% (log(p_k)/nx-colMeans((y[,k] - a_k)^2/2/s2))
      aux_track[i] <- aux_track[i] + (aux_tmp - mean(log(2*pi*s2)/2))
    }
    if(i>1){
      if((aux_track[i] - aux_track[i-1])/aux_track[i-1] < 1e-6)
        break
    }
  }
  aux <- aux_track[i]
  # mlik <- mlik_track[nIter]
  nPars <- length(a_k) + length(s2) + length(p_k)
  # BIC <- -2*mlik + nPars*log(length(y))/nid
  # AIC <- -2*mlik + 2*nPars/nid
  return(list(a_k=a_k, s2=s2, p_k=p_k, p_jk=p_jk,
              # mlik=mlik,
              # AIC=AIC, BIC=BIC, mlik_track=mlik_track
              aux=aux, aux_track=aux_track, y_max=y_max))
}

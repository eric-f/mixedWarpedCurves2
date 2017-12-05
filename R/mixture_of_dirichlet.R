#' Fit finite mixture of Dirichlet by EM algorithm
#'
#' @param w data
#' @param nclust number of mixture component
#' @param init_clust vector of initial clustering label
#' @param maxit maximum number of iterations for EM algorithm
#' @param nstart number of random start if init_clust is not specified
#' @export
mixture_of_dirichlet <- function(w, nclust, init_clust=NULL, maxit=20, nstart=1){
  if(!is.null(init_clust)){
    return(em_mixture_of_dirichlet(w, init_clust, nclust, maxit))
  }
  else{
    out <- NULL
    for(idx in 1:nstart){
      rand_clust <- sample(1:nclust, ncol(w), replace = T)
      candidate <- em_mixture_of_dirichlet(w, rand_clust, nclust, maxit)
      if(!is.null(out)){
        if(out$llk < candidate$llk)
          # print(out$llk)
          out <- candidate
      }
      else{
        out <- candidate
      }
    }
  }
  return(out)
}

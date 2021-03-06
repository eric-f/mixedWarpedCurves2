% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{em_mixture_of_dirichlet}
\alias{em_mixture_of_dirichlet}
\title{Internal function for fitting mixture of Dirichlet by EM}
\usage{
em_mixture_of_dirichlet(data, init_clust, dim_m, maxit)
}
\arguments{
\item{data}{matrix of data point on simplex}

\item{init_clust}{initial cluster label}

\item{dim_m}{number of mixture components}

\item{maxit}{number of EM iterations}
}
\description{
Internal function for fitting mixture of Dirichlet by EM
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//util.h
#ifndef UTIL_H
#define UTIL_H

arma::vec nllk_dirichlet(double log_tau, arma::vec log_dw, int n_curve, arma::vec kappa);
double compute_llk_dw(arma::vec dw, arma::vec tau_kappa);

#endif

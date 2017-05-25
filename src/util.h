// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//util.h
#ifndef UTIL_H
#define UTIL_H

arma::vec nllk_dirichlet_rcpp (double log_tau, arma::vec log_dw, int n_curve, arma::vec kappa);
void squish(arma::vec *x, double left_bound, double right_bound);
double compute_llk_dw(arma::vec dw, arma::vec tau_kappa);

#endif

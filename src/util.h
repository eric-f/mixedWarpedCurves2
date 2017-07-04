// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//util.h
#ifndef UTIL_H
#define UTIL_H

arma::vec nllk_dirichlet(double log_tau, arma::vec log_dw, int n_curve, arma::vec kappa);
double compute_llk_dw(arma::vec dw, arma::vec tau_kappa);
double newton_step_dirichlet_fixed_mean (double log_tau, arma::vec kappa, arma::vec mean_log_dw);
arma::vec newton_step_dirichlet_free_mean (arma::vec log_tau_kappa, arma::vec mean_log_dw);
#endif

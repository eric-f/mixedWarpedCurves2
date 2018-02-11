// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include "class_warping_function.h"
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>

// Constructor
Warping_Function::Warping_Function(Rcpp::List data_, int id_){
  w = Rcpp::as<arma::vec>(data_["w"]);
  dim_w = w.size();
  dw = arma::diff(w);
  log_dw = log(dw);
}

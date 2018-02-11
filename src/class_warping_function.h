// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>

// class_warping_function.h
#ifndef CLASS_WARPING_FUNCTION_H
#define CLASS_WARPING_FUNCTION_H

class Warping_Function {
public:
  int dim_w;
  arma::vec w;
  arma::vec dw;
  arma::vec log_dw;

  Warping_Function(Rcpp::List data_, int id_);
};

#endif

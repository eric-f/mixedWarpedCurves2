// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>

class Warping_Function;

// class_warping_function.h
#ifndef CLASS_WARPING_MODEL_H
#define CLASS_WARPING_MODEL_H

class Warping_Model {
public:
  int dim_kappa;
  int n;
  arma::vec alpha;
  arma::vec kappa;
  double tau;

  Warping_Model(int _dim_kappa, int _n);
  void Update_Suff_Stat(std::vector<Warping_Function*>* mydata);
  double Do_One_Newton_Update(int max_inner_iter);
  void Find_MLE(int max_outer_iter, int max_inner_iter, double tol_thres);
  void Update_LogLikelihood();
  Rcpp::List Return_Estimates();

private:
  double log_gamma_sum_alpha;
  arma::vec log_gamma_alpha;
  double log_likelihood;
  arma::vec q;
  arma::vec g;
  double z;
  double a;
  double b;
  arma::vec update_step;
  arma::vec new_alpha;
  arma::vec mean_log_dw;
  arma::vec mean_dw;
  arma::vec mean_dw2;
  int outer_iter;
  int inner_iter;
  double tmp_tol;
  double tol;
};

#endif

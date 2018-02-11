// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include "class_warping_model.h"
#include "class_warping_function.h"
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>

// Contructor
Warping_Model::Warping_Model(int _dim_kappa, int _n): dim_kappa(_dim_kappa), n(_n) {
  // Control and counter for Newton_Raphson
  inner_iter = 0;
  outer_iter = 0;
  tmp_tol = 0.0;
  tol = 0.0;

  // Allocate space for estimates and sufficient statistics
  kappa = arma::zeros(dim_kappa);
  tau = 0.0;
  alpha = tau * kappa;
  mean_log_dw = arma::zeros(dim_kappa);
  mean_dw = arma::zeros(dim_kappa);
  mean_dw2 = arma::zeros(dim_kappa);
  // Initialize working variables
  log_gamma_sum_alpha = 0.0;
  log_gamma_alpha = arma::zeros(dim_kappa);
  log_likelihood = 0.0;
  q = arma::zeros(dim_kappa);
  g = arma::zeros(dim_kappa);
  a = 0.0;
  b = 0.0;
  z = 0.0;
  update_step = arma::zeros(dim_kappa);
  new_alpha = arma::zeros(dim_kappa);
}


void Warping_Model::Update_Suff_Stat(std::vector<Warping_Function*>* mydata) {
  n = mydata->size();
  for(std::vector<Warping_Function*>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    mean_log_dw += (*it)->log_dw / n;
    mean_dw += (*it)->dw / n;
    mean_dw2 += arma::square((*it)->dw) / n;
  }
}


double Warping_Model::Do_One_Newton_Update(int max_inner_iter) {
  // Compute gradient, hessian and update step
  a = boost::math::digamma(arma::sum(alpha));
  z = boost::math::trigamma(arma::sum(alpha));
  for(int i = 0; i < dim_kappa; ++i){
    q(i) = - boost::math::trigamma(alpha(i));
    g(i) = a - boost::math::digamma(alpha(i)) + mean_log_dw(i);
  }
  b = arma::sum(g / q) / (1/z+arma::sum(1/q));
  update_step = (g - b * arma::ones(dim_kappa)) / q;

  // Step halving if new estimate is out of bound
  inner_iter = 0;
  do{
    new_alpha = alpha - update_step;
    update_step /= 2;
    ++inner_iter;
  } while (arma::any(new_alpha <= 0) && inner_iter < max_inner_iter);
  tmp_tol = arma::max(arma::abs(new_alpha - alpha));

  // Check successful update
  if(arma::all(new_alpha > 0)){
    alpha = new_alpha;
    tau = arma::sum(alpha);
    kappa = alpha / tau;
  }
  else{
    Rcpp::Rcout << "Some alpha's are still negative even after maximum attempts of step-halfings." << std::endl;
    Rcpp::Rcout << "update_step: " << update_step.t() << std::endl;
    Rcpp::Rcout << "mean_log_dw: " << mean_log_dw.t() << std::endl;
    Rcpp::Rcout << "q: " << q.t() << std::endl;
    Rcpp::Rcout << "g: " << g.t() << std::endl;
    Rcpp::Rcout << "a: " << a << std::endl;
    Rcpp::Rcout << "z: " << z << std::endl;
    tmp_tol = -99; // Terminate Newton-Raphson
    // throw;
  }
  return tmp_tol;
}


void Warping_Model::Find_MLE(int max_outer_iter, int max_inner_iter, double tol_thres){
  // Initialize estimates
  kappa = mean_dw;
  tau = (mean_dw(0) - mean_dw2(0)) / (mean_dw2(0) - pow(mean_dw(0), 2.0));
  alpha = kappa * tau;

  // Newton-Raphson
  outer_iter = 0;
  do{
    tol = Do_One_Newton_Update(max_inner_iter);
    ++outer_iter;
  } while (tol > tol_thres & outer_iter < max_outer_iter);

  // Check convergence
  if(tol <= tol_thres && tol >= 0){
    Rcpp::Rcout << "Converged, tol = " << tol << " iter = " << outer_iter << std::endl;
  }
  else if(tol < 0){
    Rcpp::Rcout << "Not converged, maximum step halving reached" << std::endl;
  }
  else{
    Rcpp::Rcout << "Not converged, tol = " << tol << " iter = " << outer_iter << std::endl;
  }
}


void Warping_Model::Update_LogLikelihood(){
  log_gamma_sum_alpha = boost::math::lgamma(arma::sum(alpha));
  for(int i = 0; i < dim_kappa; ++i){
    log_gamma_alpha(i) = boost::math::lgamma(alpha(i));
  }
  log_likelihood = n * log_gamma_sum_alpha -
    n*sum(log_gamma_alpha) +
    arma::as_scalar((alpha - arma::ones(alpha.size())).t() * mean_log_dw);
}


Rcpp::List Warping_Model::Return_Estimates() {
  // Update loglikelihood;
  Update_LogLikelihood();

  // Compute gradient and hessian
  arma::mat hess(dim_kappa, dim_kappa, arma::fill::zeros);
  a = boost::math::digamma(arma::sum(alpha));
  z = boost::math::trigamma(arma::sum(alpha));
  for(int i = 0; i < dim_kappa; ++i){
    q(i) = -boost::math::trigamma(alpha(i));
    g(i) = a - boost::math::digamma(alpha(i)) + mean_log_dw(i);
  }
  g *= n;
  q *= n;
  z *= n;
  hess = z * arma::ones(dim_kappa, dim_kappa) + arma::diagmat(q);

  // Return as R-List
  return Rcpp::List::create(
    Rcpp::Named("alpha", Rcpp::wrap(alpha)),
    Rcpp::Named("loglik", Rcpp::wrap(log_likelihood)),
    Rcpp::Named("grad", Rcpp::wrap(g)),
    Rcpp::Named("hessian", Rcpp::wrap(hess)),
    Rcpp::Named("kappa", Rcpp::wrap(kappa)),
    Rcpp::Named("tau", Rcpp::wrap(tau)),
    Rcpp::Named("n", Rcpp::wrap(n)),
    Rcpp::Named("mean_log_dw", Rcpp::wrap(mean_log_dw))
  );
}

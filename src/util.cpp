// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>

// // Function to compute minus log likelihood and derivatives
// // for Dirichlet distribution with fixed mean (old version, no mixture)
// arma::vec nllk_dirichlet (double log_tau,
//                           arma::vec log_dw,
//                           int n_curve,
//                           arma::vec kappa) {
//   double tau = std::exp(log_tau);
//   arma::vec tau_kappa = tau * kappa;
//
//   // log gamma tau, log gamma tau kappa
//   double lgamma_tau = boost::math::lgamma(tau);
//   arma::vec lgamma_tau_kappa(kappa.size(), arma::fill::zeros);
//   for (int i=0; i<kappa.size(); i++) {
//     lgamma_tau_kappa(i) = boost::math::lgamma(tau_kappa(i));
//   }
//
//   // digamma tau, digamma tau kappa
//   double digamma_tau = boost::math::digamma(tau);
//   arma::vec digamma_tau_kappa(kappa.size(), arma::fill::zeros);
//   for (int i=0; i<kappa.size(); i++) {
//     digamma_tau_kappa(i) = boost::math::digamma(tau_kappa(i));
//   }
//
//   // trigamma tau, trigamma tau kappa
//   double trigamma_tau = boost::math::trigamma(tau);
//   arma::vec trigamma_tau_kappa(kappa.size(), arma::fill::zeros);
//   for (int i=0; i<kappa.size(); i++) {
//     trigamma_tau_kappa[i] = boost::math::trigamma(tau_kappa(i));
//   }
//
//   // log-likelihood
//   double drv0 =  arma::as_scalar((tau_kappa - arma::ones(tau_kappa.size())).t() * log_dw);
//   drv0 -= n_curve * (sum(lgamma_tau_kappa) -  lgamma_tau);
//
//   // first derivative (gradient)
//   double drv1 = arma::as_scalar(tau_kappa.t() *  log_dw);
//   drv1 -= n_curve * (sum(tau_kappa % digamma_tau_kappa) - tau * digamma_tau);
//
//   // second derivative (hessian)
//   double drv2 = arma::as_scalar(tau_kappa.t() *  log_dw);
//   drv2 -= n_curve * (sum(tau_kappa % digamma_tau_kappa) +
//     sum(pow(tau_kappa, 2.0) % trigamma_tau_kappa) -
//     tau * digamma_tau -
//     std::pow(tau, 2.0) * trigamma_tau);
//
//   // return
//   arma::vec nllk(3);
//   nllk(0) = -drv0;
//   nllk(1) = -drv1;
//   nllk(2) = -drv2;
//
//   return nllk;
// }



// Function to compute Dirichlet log-likelihood
double compute_llk_dw(arma::vec dw, arma::vec tau_kappa) {
  int dim_dw;
  bool any_degenerate(false);
  double log_d_terms;
  double log_w_terms = 0;
  double density;

  dim_dw = dw.size();
  // log reciprocal of generalized beta function
  log_d_terms = - boost::math::lgamma(arma::sum(tau_kappa));
  for(int i=0; i<dim_dw; i++){
    log_d_terms += boost::math::lgamma(tau_kappa(i));
  };
  // check degenerate case
  for(int i=0; i<dim_dw; i++){
    if(dw(i)==0 && tau_kappa(i)==1){
      any_degenerate = true;
    }
  }
  // density kernel
  if(any_degenerate){
    log_w_terms = R_NegInf;
  } else{
    for(int i=0; i<dim_dw; i++){
      log_w_terms += (tau_kappa(i) - 1) * log(dw(i));
    }
  }
  density = log_w_terms - log_d_terms;
  return density;
}



// Function to compute minus log likelihood and derivatives
// for Dirichlet distribution with fixed mean (new version, for mixture)
// Parameter: 
double newton_step_dirichlet_fixed_mean (double log_tau,
                                         arma::vec kappa,
                                         arma::vec mean_log_dw) {
  double gradient;
  double hessian;
  double tau = exp(log_tau);
  arma::vec tau_kappa = tau * kappa;
  int dim_kappa = kappa.size();
  // Digamma
  double digamma_tau = boost::math::digamma(tau);
  arma::vec digamma_tau_kappa(dim_kappa, arma::fill::zeros);
  for (int i = 0; i < dim_kappa; ++i) {
    digamma_tau_kappa(i) = boost::math::digamma(tau_kappa(i));
  }
  // Trigamma
  double trigamma_tau = boost::math::trigamma(tau);
  arma::vec trigamma_tau_kappa(dim_kappa, arma::fill::zeros);
  for (int i = 0; i < dim_kappa; ++i) {
    trigamma_tau_kappa(i) = boost::math::trigamma(tau_kappa(i));
  }
  // Gradient
  gradient = tau * digamma_tau -
    arma::as_scalar(tau_kappa.t() * digamma_tau_kappa) +
    arma::as_scalar(tau_kappa.t() * mean_log_dw);
  // Hessian
  hessian = gradient +
    pow(tau, 2) * trigamma_tau -
    arma::as_scalar(arma::square(tau_kappa).t() * trigamma_tau_kappa);

  return(gradient / hessian);
}


// Function to compute minus log likelihood and derivatives
// for Dirichlet distribution with fixed mean (new version, for mixture)
// [OBSOLETE] check alternative implementation with step-halving 
//            and no log-transformation
arma::vec newton_step_dirichlet_free_mean (arma::vec log_tau_kappa,
                                           arma::vec mean_log_dw) {
  arma::vec gradient;
  arma::mat hessian;
  int dim_kappa = log_tau_kappa.size();
  arma::vec tau_kappa = exp(log_tau_kappa);
  double tau = arma::as_scalar(arma::sum(tau_kappa));
  // Digamma
  double digamma_tau = boost::math::digamma(tau);
  arma::vec digamma_tau_kappa(dim_kappa, arma::fill::zeros);
  for (int i = 0; i < dim_kappa; i++) {
    digamma_tau_kappa(i) = boost::math::digamma(tau_kappa(i));
  }
  // Trigamma
  double trigamma_tau = boost::math::trigamma(tau);
  arma::vec trigamma_tau_kappa(dim_kappa, arma::fill::zeros);
  for (int i = 0; i < dim_kappa; i++) {
    trigamma_tau_kappa(i) = boost::math::trigamma(tau_kappa(i));
  }
  // Gradient
  gradient = tau_kappa % (digamma_tau - digamma_tau_kappa + mean_log_dw);
  // Hessian
  hessian = tau_kappa * tau_kappa.t() * trigamma_tau +
    arma::diagmat(gradient - arma::square(tau_kappa) % trigamma_tau_kappa);

  return(arma::solve(hessian, gradient));
}

















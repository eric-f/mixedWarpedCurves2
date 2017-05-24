// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>

class Curve;

// class_pars.h
#ifndef CLASS_PARS_H
#define CLASS_PARS_H

class Pars {
public:
  // SAEM iterations
  int saem_counter;
  arma::vec saem_step_sizes;

  // Model Parameters
  arma::vec alpha;     // k_f x 1
  double sigma2;
  arma::vec mu;        // 2 x 1
  arma::mat big_sigma; // 2 x 2
  arma::mat big_sigma_inverse; // 2 x 2
  arma::vec kappa;     // (dim_w - 1) x 1
  double tau;

  // Dimension
  int dim_a;
  int dim_w; // dim_w - 1 = dimension of kappa;
  int dim_alpha;

  // Auxiliary variables
  int n_total;                 // total number of data points
  int n_curve;                 // number of curves
  int f_order;                 // order of the base curve spline
  arma::vec f_full_knots;      // (repeated) knot locations of the base curve
  RcppGSL::Vector f_break_points;    // boundary and internal knot locations of the base curve
  double f_left_bound;
  double f_right_bound;
  arma::mat chol_centering_mat; // (dim_w - 1) x (dim_w - 2), Cholesky decomposition of the I - 1/(dim_w - 1) * J matrix.
  arma::mat identity_cor_mat; // dim_z x dim_z identity matrix (dim_z = dim_w - 2)

  // SA-MCMC settings
  double sa_step_size_mod;
  int n_iterations;
  int n_burn_saem;
  int calibrate_period;
  int n_burn_mcmc;
  int n_core;
  double proposal_sigma;
  bool need_centering;
  double mh_accept_rate_lb;
  double mh_accept_rate_ub;
  int mh_accept_rates_table_ncol;
  arma::mat mh_accept_rate_table;
  int mh_accept_rate_table_counter;
  int mh_accept_rate_history_counter_ncol;
  arma::vec mh_accept_rate_history;
  int mh_accept_rate_history_counter;

  // Temporary variables for centering step
  arma::mat current_a_mat;

  // Constructor
  Pars(Rcpp::List pars_list,
       Rcpp::List control_list,
       RcppGSL::Vector break_points);

  // Method
  void post_simulation_housekeeping();
  void update_parameter_estimates(std::vector<Curve>* mydata);
  void advance_iteration_counter();
  void print_estimates(int interval);
  Rcpp::List return_list();
};

#endif

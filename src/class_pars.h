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
  // Model Parameters
  arma::vec alpha;     // k_f x 1
  double sigma2;
  arma::vec mu0;        // 2 x 1
  arma::mat big_sigma; // 2 x 2
  arma::mat big_sigma_inverse; // 2 x 2
  arma::vec kappa0;     // (dim_w - 1) x 1
  arma::mat kappa_clusters; // (dim_w - 1) x num_clusters
  arma::vec p_clusters;     // num_clusters x 1
  arma::vec tau_clusters;

  // Dimension
  int dim_a;
  int dim_w; // dim_w - 1 = dimension of kappa;
  int dim_alpha;
  int num_clusters;

  // Auxiliary variables
  int n_total;                    // total number of data points
  int n_curve;                    // number of curves
  int f_order;                    // order of the base curve spline
  RcppGSL::Vector f_break_points; // boundary and internal knot locations of the base curve spline
  double f_left_bound;
  double f_right_bound;
  int h_order;                    // order of the warping function splines
  RcppGSL::Vector h_break_points; // boundary and internal knot locations of the warping function splines
  double h_left_bound;
  double h_right_bound;
  arma::mat chol_centering_mat;   // (dim_w - 1) x (dim_w - 2), Cholesky decomposition of the I - 1/(dim_w - 1) * J matrix.
  arma::mat identity_cor_mat;     // dim_z x dim_z identity matrix (dim_z = dim_w - 2)

  // Diagonal var-cov matrix for amplitude effect
  bool diag_big_sigma;

  // SA-MCMC settings
  int n_core;
  int n_iterations;
  int saem_counter;
  int n_burn_saem;
  int calibrate_period;
  // SA
  double sa_step_size_mod;
  arma::vec saem_step_sizes;
  // MCMC
  int n_burn_mcmc;
  double proposal_sigma;
  bool need_centering;
  double mh_accept_rate_lb;
  double mh_accept_rate_ub;
  int mh_accept_rates_table_ncol;
  arma::mat mh_accept_rate_table;
  int mh_accept_rate_table_counter;
  arma::vec mh_accept_rate_history;
  arma::vec proposal_sigma_history;

  // Temporary variables for centering step
  arma::mat current_a_mat;
  arma::ivec current_m_vec;

  // Tracker
  arma::mat alpha_track;
  arma::vec sigma2_track;
  arma::mat big_sigma_track;
  arma::mat tau_clusters_track;
  arma::cube kappa_clusters_track;
  arma::imat sampled_m_track;

  // Stochastic approximation of Fisher information
  int num_pars;
  arma::mat sapprox_H;
  arma::vec sapprox_G;
  arma::mat sapprox_C;
  arma::mat current_H;
  arma::vec current_G;

  // Constructor
  Pars(Rcpp::List pars_list,
       Rcpp::List control_list,
       RcppGSL::Vector f_break_points_r,
       RcppGSL::Vector h_break_points_r);

  // Method
  void generate_chol_centering_mat();
  void track_mh_acceptance_and_calibrate_proposal();
  void update_parameter_estimates(std::vector<Curve>* mydata);
  void update_fisher_information_approx(std::vector<Curve>* mydata);
  void advance_iteration_counter();
  void track_estimates();
  void print_estimates(int interval);
  Rcpp::List return_pars();
  Rcpp::List return_pars(double y_scaling_factor);
  Rcpp::List return_aux();
  Rcpp::List return_pars_tracker();
  Rcpp::List return_pars_tracker(double y_scaling_factor);
  Rcpp::List return_fisher_pieces(double y_scaling_factor);
};

#endif

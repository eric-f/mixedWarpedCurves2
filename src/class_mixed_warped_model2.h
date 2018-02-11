// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>

class Mixed_Warped_Curve;

// class_pars.h
#ifndef CLASS_MIXED_WARPED_MODEL_H
#define CLASS_MIXED_WARPED_MODEL_H

class Mixed_Warped_Model {
public:

  // Model Parameters to be imported from R
  arma::vec alpha;             // dim_f x 1
  double sigma2;               // 1 x 1
  arma::vec mu_a;              // 2 x 1
  arma::vec kappa_id;          // (dim_w - 1) x 1
  int n_cluster;

  // Dimension
  int dim_f;
  int dim_a;
  int dim_w;
  int dim_kappa; // dim_kappa = dim_w - 1 = dimension of kappa;
  double dim_z; // dim_z = dim_w - 2;

  // Model Parameters to be initialize in C
  arma::vec alpha_aug;         // dim_f + 1 x 1
  arma::vec sigma2_a;          // 2 x 1
  arma::mat sigma2_a_mat;      // 2 x 2
  arma::mat sigma2_a_inv;      // 2 x 2
  arma::vec p_clusters;        // n_cluster x 1
  double tau1;                 // 1 x 1 [tau of primary component]
  arma::mat kappa_clusters;    // (dim_w - 1) x n_cluster

  // Auxiliary variables
  int n_total;                    // total number of data points
  int n_curve;                    // number of curves
  int f_order;                    // order of base shape splines
  int h_order;                    // order of the warping function splines
  double h_left_bound;
  double h_right_bound;
  arma::mat chol_centering_mat;   // (dim_w - 1) x (dim_w - 2), Cholesky decomposition of the I - 1/(dim_w - 1) * J matrix.
  arma::mat identity_cor_mat;     // dim_z x dim_z identity matrix (dim_z = dim_w - 2)
  // GSL objects...
  RcppGSL::Vector f_break_points; // boundary and internal knot locations of base shape splines
  RcppGSL::Vector h_break_points; // boundary and internal knot locations of the warping function splines

  // SA-MCMC Control parameters to be imported from R
  int n_core;
  int n_burn_saem;
  int n_burn_mcmc;
  int n_iterations;
  bool need_centering;
  double sa_step_size_mod;

  // Stochastic approximation step sizes
  arma::vec saem_step_sizes;

  // Variable for Tuning proposal sigma over SAEM burning step
  double proposal_sigma;
  int calibrate_period;
  double mh_accept_rate_lb;
  double mh_accept_rate_ub;
  arma::mat mh_accept_rate_table;
  arma::vec mh_accept_rate_history;
  arma::vec proposal_sigma_history;
  int mh_accept_rates_table_ncol;
  int mh_accept_rate_table_counter;

  // Temporary variables for tracking MH acceptance rate
  double tmp_accept_rate;

  // Sufficient Statistics for M-step
  arma::mat SS_mean_hat_mat;
  arma::mat SS_mean_sigma_a;
  arma::mat SS_mean_log_dw;
  arma::vec SS_mean_cluster_sizes;

  // Variables for keeping track of cluster membership
  arma::ivec current_m_vec;
  // Temporary variables for centering step
  arma::mat current_a_mat;

  // Counters
  int saem_counter;
  int generic_idx;
  int cluster_idx;
  int newton_idx;
  int newton_inner_idx;
  int idx_d;
  int idx_r;
  std::vector<Mixed_Warped_Curve>::iterator it;

  // Variables for Newton-Raphson steps
  // For full Dirichlet
  arma::vec newton_q;
  arma::vec newton_g;
  double newton_z;
  double newton_b;
  arma::vec newton_update_step_kappa;
  arma::vec newton_new_kappa;
  // For Dirichlet with fixed mean
  double newton_sum_digamma_k1;
  double newton_sum_trigamma_k1;
  double newton_grad1;
  double newton_hess1;
  double newton_update_step_tau1;
  double newton_new_tau1;
  // Control parameters
  int newton_max_iter;
  int newton_max_inner_iter;
  double newton_abs_tol;


  // Tracker
  arma::mat alpha_track;
  arma::vec sigma2_track;
  arma::mat sigma2_a_track;
  arma::imat sampled_m_track;
  arma::cube kappa_clusters_track;

  // Stochastic approximation of Fisher information
  int num_pars;
  arma::mat sapprox_H;
  arma::vec sapprox_G;
  arma::mat sapprox_C;
  arma::mat current_H;
  arma::vec current_G;

  // Constructor
  Mixed_Warped_Model(Rcpp::List pars_list,
                     Rcpp::List control_list,
                     RcppGSL::Vector f_break_points_r,
                     RcppGSL::Vector h_break_points_r);

  // Method
  void generate_chol_centering_mat();
  void track_mh_acceptance_and_calibrate_proposal();
  void initialize_clustering_with_user_inputs(std::vector<Mixed_Warped_Curve>* mydata);
  void gather_sufficient_statistics(std::vector<Mixed_Warped_Curve>* mydata);
  void update_estimates();  
  void update_fisher_information_approx(std::vector<Mixed_Warped_Curve>* mydata);
  void advance_iteration_counter();
  void track_estimates();
  void print_estimates(int interval);
  Rcpp::List return_pars(double y_scaling_factor);
  Rcpp::List return_aux();
  Rcpp::List return_pars_tracker(double y_scaling_factor);
  Rcpp::List return_fisher_pieces(double y_scaling_factor);

private:
  void update_observation_model();
  void update_amplitude_model();
  void update_mixture_proportion();
  void update_primary_warping_component();
  void synchronize_components();
  void update_other_warping_components();

};

#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>

class Unimodal_Curve;

// class_pars.h
#ifndef CLASS_UNIMODAL_MODEL_H
#define CLASS_UNIMODAL_MODEL_H

class Unimodal_Model {
public:
  // Model Parameters
  double sigma2;               // 1 x 1
  arma::vec mu_a;              // 2 x 1
  arma::vec sigma2_a;          // 2 x 1
  arma::mat sigma2_a_mat;      // 2 x 2
  arma::mat sigma2_a_inv;      // 2 x 2
  arma::vec p_clusters;        // n_cluster x 1
  arma::vec cluster_sizes;     // n_cluster x 1
  arma::mat kappa_clusters;    // (dim_w - 1) x n_cluster


  // Dimension
  int dim_a;
  int dim_w;
  int dim_kappa; // dim_kappa = dim_w - 1 = dimension of kappa;
  double dim_z; // dim_z = dim_w - 2;
  int n_cluster;

  // Auxiliary variables
  int n_total;                    // total number of data points
  int n_curve;                    // number of curves
  int h_order;                    // order of the warping function splines
  RcppGSL::Vector h_break_points; // boundary and internal knot locations of the warping function splines
  double h_left_bound;
  double h_right_bound;
  arma::mat chol_centering_mat;   // (dim_w - 1) x (dim_w - 2), Cholesky decomposition of the I - 1/(dim_w - 1) * J matrix.
  arma::mat identity_cor_mat;     // dim_z x dim_z identity matrix (dim_z = dim_w - 2)
  arma::vec kappa_id;             // (dim_w - 1) x 1

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

  // Temporary variables
  arma::ivec current_m_vec;
  double tmp_accept_rate;
  arma::vec minka_q;
  arma::vec minka_g;
  double minka_z;
  double minka_b;
  arma::vec minka_update_step;
  int cluster_idx;
  int generic_idx;
  int newton_idx;
  arma::mat tmp_sum_log_dw;

  // Kmeans & Newton-Raphson control parameters
  double tmp_min_km_ss;
  double tmp_new_km_ss;
  int newton_max_iter;
  double newton_abs_tol;

  // Tracker
  arma::vec sigma2_track;
  arma::mat mu_a_track;
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
  Unimodal_Model(Rcpp::List pars_list,
                 Rcpp::List control_list,
                 RcppGSL::Vector h_break_points_r);

  // Method
  void generate_chol_centering_mat();
  void track_mh_acceptance_and_calibrate_proposal();
  void update_estimates_data_mod(std::vector<Unimodal_Curve>* mydata);
  void update_estimates_amp_mod(std::vector<Unimodal_Curve>* mydata);
  void update_estimates_single_warp_mod(std::vector<Unimodal_Curve>* mydata);
  void initialize_clustering_with_user_inputs(std::vector<Unimodal_Curve>* mydata);
  void initialize_clustering_with_kmeans(std::vector<Unimodal_Curve>* mydata);
  void update_estimates_mixture_warp_mod(std::vector<Unimodal_Curve>* mydata);
  void update_fisher_information_approx(std::vector<Unimodal_Curve>* mydata);
  void advance_iteration_counter();
  void track_estimates();
  void print_estimates(int interval);
  Rcpp::List return_pars(double y_scaling_factor);
  Rcpp::List return_aux();
  Rcpp::List return_pars_tracker(double y_scaling_factor);
  Rcpp::List return_fisher_pieces(double y_scaling_factor);
};

#endif

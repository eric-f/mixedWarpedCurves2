// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>

class Mixed_Warped_Model;

// class_curve.h
#ifndef CLASS_MIXED_WARPED_CURVE_H
#define CLASS_MIXED_WARPED_CURVE_H

class Mixed_Warped_Curve {
public:

  // id for internal use
  int curve_id;

  // Pointer to a common Pars class objects
  Mixed_Warped_Model *common_pars;

  // Raw Data
  arma::vec y;
  arma::vec x;
  int init_clust;

  // Dimensions
  int n_i;   // number of points per curve
  int dim_f; // number of basis coefficients for the base shape
  int dim_a; // number of amplitude effects (default = 2)
  int dim_w; // number of basis coefficients for the warping function
  int dim_z; // dim_w - 2
  int n_cluster;

  // Basis Evaluation Matrix - warping function
  arma::mat h_basis_mat;



  // MCMC working variables
  int current_m;                             // cluster membership
  arma::vec current_a;                       // dim_a x 1
  arma::vec current_dw;                      // (dim_w - 1) x 1
  arma::vec current_w;                       // dim_w x 1

  // Linking current_w to random walk in R^dim_z
  arma::vec current_z;                       // (dim_w - 2) x 1

  // Fitted object based on current warping
  arma::vec current_warped_x;                // n_i x 1
  arma::mat current_warped_f_basis_mat;      // n_i x dim_f;

  arma::vec proposed_w;                      // dim_w x 1
  arma::vec proposed_dw;                     // (dim_w - 1) x 1
  arma::vec proposed_z;                      // (dim_w - 2) x 1
  arma::vec proposed_warped_x;               // n_i x 1
  arma::mat proposed_warped_f_basis_mat;     // n_i x dim_f;

  // Sufficient statistics based on the current MCMC draw
  arma::mat current_log_dw;                  // (dim_w - 1) x n_cluster
  arma::vec current_cluster_membership;      // n_cluster x 1
  arma::mat current_aug_warped_f_basis_mat;  // n_i x (dim_f + 1)
  arma::mat current_hat_mat;                 // (dim_f + 1) x (dim_f + 1)
  arma::mat current_cov_a;                   // dim_a x dim_a

  // Stochastic approximated sufficient statistics
  arma::vec sapprox_a;                       // dim_a x 1
  arma::vec sapprox_w;                       // dim_w x 1
  arma::mat sapprox_log_dw;                  // (dim_w - 1) x n_cluster
  arma::vec sapprox_cluster_membership;      // n_cluster x 1
  arma::mat sapprox_aug_warped_f_basis_mat;  // n_i x (dim_f + 1)          [Psi]
  arma::mat sapprox_hat_mat;                 // (dim_f + 1) x (dim_f + 1)  [C_mat]
  arma::mat sapprox_cov_a;                   // dim_a x dim_a              [Sigma_a]  

  // Constructor
  Mixed_Warped_Curve(Rcpp::List data, Mixed_Warped_Model* pars, int id, int seed);

  // Methods
  void initialize_h_basis_mat();
  void initialize_current_f_basis_mat();
  void do_simulation_step();
  void center_current_a();
  void update_sufficient_statistics_approximates();
  Rcpp::List return_list();
  Rcpp::List return_list(double y_scaling_factor);
  void print_random_number();

private:
  // random number generator
  gsl_rng *rng_gen;

  // Temporary variables
  // ... for initialization of current_z
  arma::vec tmp_log_current_dw;
  arma::vec tmp_log_centered_dw;
  // ... for initialize_h_basis_mat()
  gsl_vector *tmp_b_vec;
  gsl_bspline_workspace *tmp_bw;
  // ... for draw_new_m()
  double tmp_u;
  arma::vec pred_prob_clusters;
  int cluster_idx;
  // ... for draw_new_a()
  arma::mat tmp_f_mat;
  arma::vec tmp_mu_post;
  arma::mat tmp_sigma_post;
  //... for propose_new_w()
  arma::vec tmp_dw;
  // ... for compute_log_mh_ratio()
  arma::vec proposed_warped_f;
  arma::vec proposed_fitted_y;
  arma::vec current_warped_f;
  arma::vec current_fitted_y;
  double proposed_residual_sum_of_squares;
  double current_residual_sum_of_squares;
  double proposed_minus_current_llk_data;
  double current_llk_w;
  double proposed_llk_w;
  double log_jacobian_term;
  double mh_log_accept_prob;
  // ... for mh_accept_reject()
  double mh_randu;
  // ... for center_current_a()
  arma::mat centering_mat;
  arma::vec mean_current_a;
  // ... for update_sufficient_statistics_approximates()
  double current_step_size;
  arma::mat tmp_half_hat_mat;
  // ... for do_simulation_step()
  int mcmc_iter;

  // Internal functions for the MH-within-Gibbs sampler
  void propose_new_w();
  void compute_proposed_warping_and_f_basis_mat();
  void compute_log_mh_ratio();
  void mh_accept_reject();
  void draw_new_a();
  void draw_new_m();
};

#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>

class Mixed_Warping_Model;

// class_curve.h
#ifndef CLASS_MIXED_WARPING_CURVE_H
#define CLASS_MIXED_WARPING_CURVE_H

class Mixed_Warping_Curve {
public:

  // Pointer to a common Pars class objects
  Mixed_Warping_Model *common_pars;

  // id for internal use
  int curve_id;

  // Raw Data
  arma::vec y;
  arma::vec x;
  int init_clust;

  // Dimensions
  int n_i;  // number of points per curve
  int dim_f; // number of basisi coefficients for the shape function
  int dim_w; // number of basis coefficients for the warping function
  int dim_a; // number of amplitude effects (default = 2)
  int dim_z; // dim_w - 2
  int n_cluster;

  // Basis Evaluation Matrix - warping function
  arma::mat h_basis_mat;                // n_i x dim_w

  // Metropolis-Hasting
  arma::mat current_a;                    // dim_a x n_cluster
  arma::mat current_dw;                   // (dim_w - 1) x n_cluster
  arma::mat current_w;                    // dim_w x n_cluster
  arma::mat current_log_dw;               // (dim_w - 1) x n_cluster
  arma::mat current_z;                    // dim_z x n_cluster
  arma::mat current_warped_x;             // n_i x n_cluster
  arma::cube current_warped_f_basis_mat;  // n_i x dim_f

  arma::vec proposed_dw;                  // (dim_w - 1) x 1
  arma::vec proposed_w;                   // dim_w x 1
  arma::vec proposed_z;                   // (dim_w - 2) x 1
  arma::vec proposed_warped_x;            // n_i x 1
  arma::mat proposed_warped_f_basis_mat;  // n_i x dim_f

  // Sufficient Statistics
  arma::vec sapprox_residual_sum_of_squares; // n_cluster
  arma::mat sapprox_a;                       // dim_a x n_cluster
  arma::mat sapprox_sq_a;                    // dim_a x n_cluster
  arma::mat sapprox_w;                       // dim_w x n_cluster
  arma::mat sapprox_dw;                      // (dim_w - 1) x n_cluster
  arma::mat sapprox_log_dw;                  // (dim_w - 1) x n_cluster
  arma::cube sapprox_warped_f_basis_mat;     // n_i x dim_f x n_cluster
  arma::mat sapprox_warped_f;                // n_i x n_cluster
  arma::mat sapprox_fitted_y;                // n_i x n_cluster

  arma::vec current_yy;                      // n_cluster x 1
  arma::mat current_Xy;                      // dim_f x n_cluster
  arma::cube current_XX;                     // dim_f x dim_f x n_cluster

  arma::vec sapprox_yy;                      // n_cluster x 1
  arma::mat sapprox_Xy;                      // dim_f x n_cluster
  arma::cube sapprox_XX;                     // dim_f x dim_f x n_cluster

  arma::vec pred_prob_clusters;              // n_cluster x 1

  // Constructor
  Mixed_Warping_Curve(Rcpp::List data, Mixed_Warping_Model* pars, int id, int seed);

  // Methods
  void initialize_h_basis_mat();
  void do_simulation_step();
  void update_sufficient_statistics_approximates();
  Rcpp::List return_list(double y_scaling_factor);

private:
  // Random number generator
  gsl_rng * rng_gen;

  // Temporary variables
  // ... for update_sufficient_statistics_approximates()
  double current_step_size;

  // Internal functions for the MH-within-Gibbs sampler
  void propose_new_w(int cluster_idx);
  arma::vec tmp_dw;
  // --------------------------------------------------
  void compute_proposed_warping_and_f_basis_mat();
  gsl_vector *tmp_b_vec;
  gsl_bspline_workspace *tmp_bw;
  // --------------------------------------------------
  void compute_log_mh_ratio(int cluster_idx);
  arma::vec proposed_warped_f;                   // n_i x 1
  arma::vec proposed_fitted_y;                   // n_i x 1
  arma::mat current_warped_f;                    // n_i x n_cluster
  arma::mat current_fitted_y;                    // n_i x n_cluster
  double proposed_residual_sum_of_squares;
  arma::vec current_residual_sum_of_squares;     // n_cluster x 1
  double proposed_minus_current_llk_data;
  arma::vec current_llk_w;                       // n_cluster x 1
  double proposed_llk_w;
  double log_jacobian_term;
  // --------------------------------------------------
  void mh_accept_reject(int cluster_idx);
  double mh_randu;
  double tmp_mh_log_accept_prob;
  // --------------------------------------------------
  void draw_new_a(int cluster_idx);
  arma::mat tmp_f_mat;
  arma::vec tmp_mu_post;
  arma::mat tmp_sigma_post;
  // --------------------------------------------------
  void center_current_a();
  arma::mat centering_mat;
  arma::vec mean_current_a;
  // --------------------------------------------------
  void update_post_prob();
  arma::vec sapprox_llk_y; // n_cluster x 1
  arma::vec sapprox_llk_a; // n_cluster x 1
  arma::vec sapprox_llk_w; // n_cluster x 1
};

#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>

class Unimodal_Model;

// class_curve.h
#ifndef CLASS_UNIMODAL_CURVE_H
#define CLASS_UNIMODAL_CURVE_H

class Unimodal_Curve {
public:

  // Pointer to a common Pars class objects
  Unimodal_Model *common_pars;

  // id for internal use
  int curve_id;

  // Raw Data
  arma::vec y;
  arma::vec x;
  int init_clust;

  // Basis Evaluation Matrix - warping function
  arma::mat h_basis_mat;

  // Dimensions
  int n_i;  // number of points per curve
  int dim_a; // number of amplitude effects (default = 2)
  int dim_w; // number of basis coefficients for the warping function
  int dim_z; // dim_w - 2
  int n_cluster;

  // Metropolis-Hasting
  int current_m;                        // cluster membership
  arma::vec current_a;                  // dim_a x 1
  arma::mat current_log_dw;             // (dim_w - 1) x n_cluster
  arma::vec current_cluster_membership; // n_cluster x 1

  arma::vec current_w;                  // dim_w x 1
  arma::vec current_dw;                 // (dim_w - 1) x 1
  arma::vec current_z;                  // (dim_w - 2) x 1
  arma::vec current_warped_x;           // n_i x 1

  arma::vec proposed_w;                 // dim_w x 1
  arma::vec proposed_dw;                // (dim_w - 1) x 1
  arma::vec proposed_z;                 // (dim_w - 2) x 1
  arma::vec proposed_warped_x;          // n_i x 1


  // Sufficient Statistics
  double sapprox_residual_sum_of_squares; // 1 x 1
  arma::vec sapprox_a;                       // dim_a x 1
  arma::mat sapprox_sq_a;                    // dim_a x 1
  arma::mat sapprox_log_dw;                  // (dim_w - 1) x n_cluster
  arma::vec sapprox_cluster_membership;      // n_cluster x 1

  // S-approx for prediction
  arma::vec sapprox_w;                       // dim_w x 1
  arma::mat sapprox_warped_f;                // n_i x 1
  arma::mat sapprox_fitted_y;                // n_i x 1


  // Constructor
  Unimodal_Curve(Rcpp::List data, Unimodal_Model* pars, int id, int seed);

  // Methods
  void initialize_h_basis_mat();
  void do_simulation_step();
  void update_sufficient_statistics_approximates();
  Rcpp::List return_list();
  Rcpp::List return_list(double y_scaling_factor);
  void print_random_number();

private:
  // random number generator
  gsl_rng * rng_gen;

  // Temporary variables
  // ... for draw_new_m()
  double tmp_u;
  arma::vec tmp_pred_prob_clusters;
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
  // ... for mh_accept_reject()
  double mh_randu;
  // ... for center_current_a()
  arma::mat centering_mat;
  arma::vec mean_current_a;
  // ... for update_sufficient_statistics_approximates()
  double current_step_size;
  arma::mat tmp_half_hat_mat;
  // ...for mh_accept_reject()
  double tmp_mh_log_accept_prob;

  // Internal functions for the MH-within-Gibbs sampler
  void propose_new_w();
  void compute_proposed_warping();
  void compute_log_mh_ratio();
  void mh_accept_reject();
  void draw_new_a();
  void draw_new_m();
};

#endif

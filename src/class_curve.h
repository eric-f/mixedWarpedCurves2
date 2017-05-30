// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <random>

class Pars;

// class_curve.h
#ifndef CLASS_CURVE_H
#define CLASS_CURVE_H

class Curve {
public:

  // Pointer to a common Pars class objects
  Pars *common_pars;

  // id for internal use
  int curve_id;

  // Raw Data
  arma::vec y;
  arma::vec x;

  // Basis Evaluation Matrix - warping function
  arma::mat h_basis_mat;

  // Dimensions
  int n_i;  // number of points per curve
  int dim_a; // number of amplitude effects (default = 2)
  int dim_w; // number of basis coefficients for the warping function
  int dim_z; // dim_w - 2
  int dim_alpha; // number of basis coefficients for the base curve

  // Sufficient Statistics
  arma::vec current_a; // dim_a x 1

  arma::vec current_w; // dim_w x 1
  arma::vec current_dw; // dim_w x 1
  arma::vec current_z; // dim_w x 1
  arma::vec current_warped_x;  // n_i x 1
  arma::mat current_warped_f_basis_mat; // n_i x dim_alpha

  arma::vec proposed_w; // dim_w x 1
  arma::vec proposed_dw; // dim_w x 1
  arma::vec proposed_z; // dim_w x 1
  arma::vec proposed_warped_x; // n_i x 1
  arma::mat proposed_warped_f_basis_mat; // n_i x dim_alpha

  arma::vec sapprox_a; // dim_a x 1
  arma::vec sapprox_w; // dim_w x 1
  arma::mat sapprox_warped_f_basis_mat; // n_i x dim_alpha

  arma::mat sapprox_aug_warped_f_basis_mat; // n_i x (dim_alpha + 1)  [Psi]
  arma::mat sapprox_hat_mat;   // (dim_alpha + 1) x (dim_alpha + 1)   [C_mat]
  arma::mat sapprox_sigma_a; // dim_a x dim_a                         [Sigma_a]
  arma::vec sapprox_log_dw;  // (dim_w - 1) x 1

  arma::mat current_aug_warped_f_basis_mat;   // n_i x (dim_alpha + 1)
  arma::mat current_hat_mat;   // (dim_alpha + 1) x (dim_alpha + 1)
  arma::mat current_sigma_a; // dim_a x dim_a
  arma::vec current_log_dw;  // (dim_w - 1) x 1

  // Constructor
  Curve(Rcpp::List data, Pars* pars, int id, int seed);

  // Methods
  void initialize_h_basis_mat();
  void initialize_current_f_basis_mat();
  void draw_new_a();
  void do_simulation_step();
  void center_current_a();
  void update_sufficient_statistics_approximates();
  Rcpp::List return_list();
  Rcpp::List return_list(double y_scaling_factor);

private:
  // random number generator
  std::mt19937 gen;
  std::uniform_real_distribution<double> dist;
  // Temporary variables
  gsl_vector *tmp_b_vec;
  gsl_bspline_workspace *tmp_bw;
  // ... for draw_new_a()
  arma::mat tmp_f_mat;
  arma::vec tmp_mu_post;
  arma::mat tmp_sigma_post;
  //... for propose_new_w()
  arma::vec tmp_dw;
  // ... for compute_log_mh_ratio()
  arma::vec proposed_warped_f;
  arma::vec proposed_warped_y;
  arma::vec current_warped_f;
  arma::vec current_warped_y;
  arma::vec proposed_residual_sum_of_squares;
  arma::vec current_residual_sum_of_squares;
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

  void propose_new_w();
  void compute_proposed_warping_and_f_basis_mat();
  void compute_log_mh_ratio();
  void mh_accept_reject();
};

#endif

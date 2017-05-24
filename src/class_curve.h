// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

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
  Curve(Rcpp::List curve_obj, Pars* pars, int id);

  // Methods
  void initialize_current_f_basis_mat();
  void do_simulation_step();
  void center_current_a();
  void update_sufficient_statistics_approximates();
  Rcpp::List return_list();

private:
  void propose_new_w();
  void compute_proposed_warping_and_f_basis_mat();
  double compute_log_mh_ratio();
  void mh_accept_reject();
  void draw_new_a();
};

#endif

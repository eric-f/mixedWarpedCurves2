// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
// #include <random>
// #include "ctime"

class Mixture_A_Pars;

// class_curve.h
#ifndef CLASS_MIXTURE_A_CURVE_H
#define CLASS_MIXTURE_A_CURVE_H

class Mixture_A_Curve {
public:

  // Pointer to a common Pars class objects
  Mixture_A_Pars *common_pars;

  // id for internal use
  int curve_id;

  // Raw Data
  arma::vec y;
  arma::vec x;

  // Basis Evaluation Matrix - warping function
  arma::mat f_basis_mat;
  arma::mat aug_f_basis_mat;

  // Dimensions
  int n_i;          // number of points per curve
  int dim_a;        // number of amplitude effects (const. = 2)
  int dim_alpha;    // number of basis coefficients for the base curve
  int num_clusters;

  // MCMC state
  int initial_m;
  int current_m;       // cluster membership
  arma::mat current_a; // amplitude effect conditioning on cluster

  // Sufficient Statistics
  arma::vec current_cluster_membership; // num_clusters x 1
  arma::vec sapprox_cluster_membership; // num_clusters x 1
  double SS_sigma2_sh;
  double SS_sigma2_sc;
  arma::cube SS_XtX;
  arma::mat SS_XtY;
  double SS_sigma2;
  arma::vec SS_post_prob;

  // Constructor
  Mixture_A_Curve(Rcpp::List data, Mixture_A_Pars* pars, int id, int seed);

  // Methods
  void initialize_f_basis_mat();
  void do_MH_simulation_step();
  void do_Gibbs_expectation_step();
  void do_expectation_step();
  Rcpp::List return_list(double y_scaling_factor);

private:
  // random number generator
  gsl_rng * rng_gen;

  // Temporary variables
  gsl_vector *tmp_b_vec;
  gsl_bspline_workspace *tmp_bw;
  // ... for update_pieces_for_m_step()
  double current_step_size;
  arma::mat BtB;
  double YtY;
  double Yt1;
  arma::vec Ft1;
  arma::vec FtY;
  arma::vec FtF;
  // ... for compute_a_posterior()
  arma::mat fitted_f;                      // n_i x num_clusters
  arma::mat fitted_1f_slice;               // n_i x 2
  arma::mat fitted_1f_star_slice;          // n_i x 2
  arma::mat post_amp_mu;                 // num_clusters x 2
  arma::mat post_amp_sigma_slice;        // 2 x 2
  arma::cube post_amp_sigma;             // 2 x 2 x num_clusters
  arma::vec post_amp_shsh;               //  num_clusters x 2
  arma::vec post_amp_shsc;               //  num_clusters x 2
  arma::vec post_amp_scsc;               //  num_clusters x 2
  arma::mat post_y_mu;                   // n_i x num_clusters
  arma::cube post_y_sigma_inv;           // n_i x n_i x num_clusters
  arma::vec post_y_sigma_log_det_scaled; // num_clusters
  arma::vec post_y_sigma_log_det_sign;   // num_clusters
  arma::vec post_y_llk;                  // num_clusters x 1
  // ... for draw_new_m()

  // Internal functions for the MH-within-Gibbs sampler
  void compute_a_posterior();
  void draw_new_a();
  void draw_new_m();
  void compute_conditional_post_prob();
  void compute_marginal_post_prob();
  void update_sapprox_amp_effect();
  void update_sapprox_cluster_membership();
  void update_pieces_for_m_step();

};

#endif

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>

class Mixed_Shape_Model;

// class_curve.h
#ifndef CLASS_MIXED_SHAPE_CURVE_H
#define CLASS_MIXED_SHAPE_CURVE_H

class Mixed_Shape_Curve {
public:

  // Pointer to a common Pars class objects
  Mixed_Shape_Model *common_pars;

  // id for internal use
  int curve_id;

  // Raw Data
  arma::vec y;
  arma::vec x;
  double y_scl;

  // Basis Evaluation Matrix - warping function
  arma::mat f_basis_mat;

  // Dimensions
  int n_i;          // number of points per curve
  int dim_a;        // number of amplitude effects (const. = 2)
  int dim_alpha;    // number of basis coefficients for the base curve
  int num_clusters;

  // MCMC state
  int initial_m;
  int current_m;                        // cluster membership
  arma::vec current_cluster_membership; // num_clusters x 1
  arma::vec sapprox_cluster_membership; // num_clusters x 1
  arma::mat current_a;                  // amplitude effect conditioning on cluster

  // Sufficient Statistics
  double SS_sigma2_sh;
  double SS_sigma2_sc;
  arma::cube SS_XtX;
  arma::mat SS_XtY;
  double SS_sigma2;
  arma::vec SS_post_prob;

  // Constructor
  Mixed_Shape_Curve(Rcpp::List data,
                    Mixed_Shape_Model* pars,
                    int id,
                    double y_scl_r,
                    int seed);

  // Methods
  void initialize_f_basis_mat();
  void do_MH_simulation_step();
  Rcpp::List return_list();

private:
  // random number generator
  gsl_rng * rng_gen;

  // Temporary variables
  // ... B-spline workspace
  gsl_vector *tmp_b_vec;
  gsl_bspline_workspace *tmp_bw;
  // ... Fixed terms in the likelihood function
  double YtY;
  double Yt1;
  arma::mat BtB;
  arma::vec Ft1;
  arma::vec FtY;
  arma::vec FtF;
  // ... Fitted base shapes
  arma::mat fitted_f;                      // n_i x num_clusters
  arma::mat fitted_1f_slice;               // n_i x 2
  arma::mat fitted_1f_star_slice;          // n_i x 2
  // ... Posterior moments and likelihood
  arma::mat post_amp_mu;                 // num_clusters x 2
  arma::mat post_amp_sigma_slice;        // 2 x 2
  arma::cube post_amp_sigma;             // 2 x 2 x num_clusters
  arma::vec post_amp_shsh;               //  num_clusters x 2
  arma::vec post_amp_shsc;               //  num_clusters x 2
  arma::vec post_amp_scsc;               //  num_clusters x 2
  arma::vec post_y_llk;                  // num_clusters x 1
  // ... Stochastic approximation step size
  double current_step_size;
  // ... counter
  int clust_idx;

  // Internal functions for the MH-within-Gibbs sampler
  void compute_a_posterior();
  void draw_new_a();
  void draw_new_m();
  void compute_conditional_post_prob();
  void update_sapprox_cluster_membership();
  void update_pieces_for_m_step();

};

#endif

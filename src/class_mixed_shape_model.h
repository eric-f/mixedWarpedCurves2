// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>

class Mixed_Shape_Curve;

// class_pars.h
#ifndef CLASS_MIXED_SHAPE_MODEL_H
#define CLASS_MIXED_SHAPE_MODEL_H

class Mixed_Shape_Model {
public:

  // Counter
  int saem_counter;
  int clust_idx;

  // Model Parameters
  arma::mat alpha;             // k_f x num_clusters
  double sigma2;
  arma::vec mu0;               // 2 x 1
  arma::mat amp_sigma;         // 2 x 2 (Diagonal)
  arma::mat amp_sigma_inverse; // 2 x 2 (Diagonal)
  arma::vec p_clusters;        // num_clusters x 1

  // Dimension
  int dim_a;        // 2
  int dim_alpha;    // k_f
  int num_clusters; // M

  // Auxiliary variables
  int n_total;                    // total number of data points
  int n_curve;                    // number of curves
  int f_order;                    // order of the base curve spline
  double y_scl;

  // B-spline
  RcppGSL::Vector f_break_points; // boundary and internal knot locations of the base curve spline
  double f_left_bound;
  double f_right_bound;

  // SA-MCMC settings
  int n_core;
  int n_iterations;
  int n_burn_saem;
  int calibrate_period;
  // SA
  double sa_step_size_mod;
  arma::vec saem_step_sizes;
  // MCMC
  int n_burn_mcmc;
  bool need_centering;
  arma::ivec current_m_vec;
  // MCMC - Temporary variables for centering step
  arma::mat current_a_mat;
  arma::cube XtX;
  arma::mat XtY;

  // Tracker
  arma::cube alpha_track;
  arma::vec sigma2_track;
  arma::mat amp_sigma_track;
  arma::mat p_clusters_track;
  arma::imat sampled_m_track;

  // Stochastic approximation of Fisher information
  int num_pars;
  arma::mat sapprox_H;
  arma::vec sapprox_G;
  arma::mat sapprox_C;
  arma::mat current_H;
  arma::vec current_G;

  // Constructor
  Mixed_Shape_Model(Rcpp::List pars_list,
                    Rcpp::List control_list,
                    double y_scl_r,
                    RcppGSL::Vector f_break_points_r);

  // Method
  void update_parameter_estimates(std::vector<Mixed_Shape_Curve*>* mydata);
  void advance_iteration_counter();
  void track_estimates();
  void print_estimates(int interval);
  Rcpp::List return_pars();
  Rcpp::List return_aux();
  Rcpp::List return_pars_tracker();
};

#endif

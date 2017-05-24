// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// #include <vector>
// #include "class_curve.h"
//
// // class_model.h
// #ifndef CLASS_MODEL_H
// #define CLASS_MODEL_H
//
// class Model {
// public:
//   // Model Parameters
//   arma::vec alpha;     // k_f x 1
//   double sigma2;
//   arma::vec mu;        // 2 x 1
//   arma::mat big_sigma; // 2 x 2
//   arma::vec kappa;     // (dim_w - 1) x 1
//   double tau;
//
//   // Auxiliary variables
//   int n_total;                 // total number of data points
//   int n_curve;                 // number of curves
//   int f_order;                 // order of the base curve spline
//   arma::vec f_full_knots;      // (repeated) knot locations of the base curve
//   arma::vec f_break_points;    // internal knot locations of the base curve
//   arma::mat big_sigma_inverse; // 2 x 2
//   int D;                       // (dim_w - 1) - dimension of kappa
//   arma::mat L; // D x (D-1), Matrix for transforming a D-dimensional Gaussian proposal to the D-simplex
//
//   // Curve objects
//   std::vector<Curve> curves;
//
//   // SA-MCMC settings
//   double sa_step_size_mod;
//   int n_iterations;
//   int n_burn_saem;
//   int calibrate_period;
//   int n_burn_mcmc;
//   int n_core;
//   double proposal_sigma;
//   bool need_centering;
//   double mh_accept_rate_lb;
//   double mh_accept_rate_ub;
//   int mh_accept_rates_window;
//   arma::vec mh_accept_rate_rolling_avg;
//   arma::vec mh_accept_rate_buffer;
//   int mh_accept_rate_counter;
//
//   // Constructor
//   Model(Rcpp::List pars,
//         Rcpp::List curves_list,
//         Rcpp::List control_list);
//
//   // Methods
//   void do_simulation_step(); // update current_a and current_w
//   void do_centering_step(); // center current_a
//   void do_stochastic_approximation_step(); // update and gather sufficient statistics
//   void do_update_likelihood_and_information_matrix();
//   void do_maximization_step(); // update_model_parameters
// };
//
// #endif

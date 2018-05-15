// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include "class_mixed_shape_model.h"
#include "class_mixed_shape_curve.h"

// Constructor
Mixed_Shape_Model::Mixed_Shape_Model(Rcpp::List pars_list,
                                     Rcpp::List control_list,
                                     double y_scl_r,
                                     RcppGSL::Vector f_break_points_r) : y_scl(y_scl_r), f_break_points(f_break_points_r){

  // Iteration counter
  saem_counter = 0;
  clust_idx = 0;

  // Fixed parameters
  mu0 = Rcpp::as<arma::vec>(pars_list["mu"]);

  // Shape parameter and variance of error term
  alpha = Rcpp::as<arma::mat>(pars_list["alpha"]);
  sigma2 = Rcpp::as<double>(pars_list["sigma2"]);

  // Dimension
  dim_a = mu0.size();
  dim_alpha = alpha.n_rows;
  num_clusters = Rcpp::as<int>(pars_list["num_clusters"]);

  // Random effect parameters
  amp_sigma = arma::eye(dim_a, dim_a) * 0.01;
  amp_sigma_inverse = arma::eye(dim_a, dim_a) * 100;
  p_clusters = arma::ones(num_clusters) / num_clusters;

  // Auxiliary variables
  n_total = Rcpp::as<int>(control_list["n_total"]);
  n_curve = Rcpp::as<int>(control_list["n_curve"]);
  f_order = Rcpp::as<int>(control_list["f_order"]);
  f_left_bound = gsl_vector_min(f_break_points);
  f_right_bound = gsl_vector_max(f_break_points);

  // SA-MCMC Control parameters
  n_burn_saem = Rcpp::as<int>(control_list["n_saem_burn"]);
  n_iterations = Rcpp::as<double>(control_list["n_saem_iter"]) + 2 * n_burn_saem;
  n_burn_mcmc = Rcpp::as<int>(control_list["n_mcmc_burn"]);
  n_core = Rcpp::as<int>(control_list["n_core"]);
  sa_step_size_mod = Rcpp::as<double>(control_list["saem_step_seq_pow"]);

  // Generate step sizes
  // Searching stage
  saem_step_sizes = arma::ones(n_iterations);
  // Zone in stage
  saem_step_sizes(arma::span(n_burn_saem, 2 * n_burn_saem - 1)) =
    arma::pow(arma::linspace(1, n_burn_saem, n_burn_saem), -sa_step_size_mod);
  // Averaging stage
  saem_step_sizes(arma::span(2 * n_burn_saem, n_iterations - 1)) =
    arma::pow(arma::linspace(1, n_iterations - 2 * n_burn_saem, n_iterations - 2 * n_burn_saem), -sa_step_size_mod);

  // Temporary variables for centering step
  current_a_mat = arma::zeros(dim_a, n_curve);
  current_m_vec = arma::zeros<arma::ivec>(n_curve);
  XtX = arma::zeros(dim_alpha, dim_alpha, num_clusters);
  XtY = arma::zeros(dim_alpha, num_clusters);

  // Trackers
  alpha_track = arma::zeros(dim_alpha, num_clusters, n_iterations);
  sigma2_track = arma::zeros(n_iterations);
  amp_sigma_track = arma::zeros(dim_a, n_iterations);
  p_clusters_track = arma::zeros(num_clusters, n_iterations);
  sampled_m_track = arma::zeros<arma::imat>(n_curve, n_iterations);

  // Stochastic approximation of Fisher information
  num_pars = dim_alpha * num_clusters + 1 + dim_a + (num_clusters - 1);
  sapprox_H = arma::zeros(num_pars, num_pars);
  sapprox_C = arma::zeros(num_pars, num_pars);
  sapprox_G = arma::zeros(num_pars);
  current_H = arma::zeros(num_pars, num_pars);
  current_G = arma::zeros(num_pars);
}





// Gather stochastic approximates of the sufficient statistics and perform an M-step
// Depends on:
// Changes:
void Mixed_Shape_Model::update_parameter_estimates(std::vector<Mixed_Shape_Curve*>* mydata){
  // Rcpp::Rcout << "update_parameter_estimates" << std::endl;
  // Reset parameters
  XtX.zeros();
  XtY.zeros();
  p_clusters.zeros();
  amp_sigma.zeros();
  sigma2 = 0.0;

  // Gather sufficient statistics
  for(std::vector<Mixed_Shape_Curve*>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    XtX += (*it)->SS_XtX;
    XtY += (*it)->SS_XtY;
    p_clusters += (*it)->SS_post_prob / n_curve;
    amp_sigma(0,0) += (*it)->SS_sigma2_sh;
    amp_sigma(1,1) += (*it)->SS_sigma2_sc;
    // Rcpp::Rcout << "sigma2: " << sigma2 << std::endl;
    // Rcpp::Rcout << "(*it)->SS_sigma2: " << (*it)->SS_sigma2 << std::endl;
    sigma2 += (*it)->SS_sigma2;
  }

  // Update alpha (by cluster)
  for(clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    alpha.col(clust_idx) = arma::solve(XtX.slice(clust_idx), XtY.col(clust_idx));
  }

  return;
}



// Increase saem iteration counter
// Depends on: Nil
// Changes: saem_counter
void Mixed_Shape_Model::advance_iteration_counter(){
  ++saem_counter;
  return;
}



// Store current estiamtes for tracking
// Depends on: saem_counter, alpha, sigma2, big_sigma, tau
// Changes: alpha_track, sigma2_track, big_sigma_track, tau_track
void Mixed_Shape_Model::track_estimates(){
  // Track estimates
  alpha_track.slice(saem_counter) = alpha;
  sigma2_track(saem_counter) = sigma2;
  amp_sigma_track.col(saem_counter) = amp_sigma.diag();
  p_clusters_track.col(saem_counter) = p_clusters;
  sampled_m_track.col(saem_counter) = current_m_vec;
}



// Print estimates for monitoring
// Depends on: saem_counter, proposal_sigma, big_sigma, alpha, sigma2, tau
//             mh_accept_rate_history
// Changes: Nil
void Mixed_Shape_Model::print_estimates(int interval){
  if(saem_counter % interval == 0 & saem_counter < n_iterations){
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "=================================" << std::endl;
    Rcpp::Rcout << "Iteration: " << std::endl << saem_counter << std::endl;
    Rcpp::Rcout << "big_sigma: " << std::endl << amp_sigma << std::endl;
    Rcpp::Rcout << "alpha: " << std::endl << alpha << std::endl;
    Rcpp::Rcout << "sigma2: " << sigma2 << std::endl;
    Rcpp::Rcout << "p_clusters: " << std::endl << p_clusters.t() << std::endl;
    Rcpp::Rcout << "membership state: " << std::endl << current_m_vec.t() << std::endl;
    Rcpp::Rcout << "=================================" << std::endl;
    Rcpp::Rcout << std::endl;
  }
  return;
}



// Return estimated parameters as R list
Rcpp::List Mixed_Shape_Model::return_pars(){
  // Scaling matrices
  arma::mat a_scaling_mat = arma::eye(2, 2);
  a_scaling_mat(0, 0) = y_scl;
  return Rcpp::List::create(
    Rcpp::Named("mu0", Rcpp::wrap(mu0)),
    Rcpp::Named("alpha", Rcpp::wrap(alpha * y_scl)),
    Rcpp::Named("sigma2", Rcpp::wrap(sigma2 * std::pow(y_scl, 2))),
    Rcpp::Named("amp_sigma", Rcpp::wrap(a_scaling_mat * amp_sigma * a_scaling_mat)),
    Rcpp::Named("p_clusters", Rcpp::wrap(p_clusters))
  );
};



// Return auxiliary information as R list
// Note: Rcpp::List::create seems to limit the number of items in the list
Rcpp::List Mixed_Shape_Model::return_aux(){
  return Rcpp::List::create(
    Rcpp::Named("n_total", Rcpp::wrap(n_total)),
    Rcpp::Named("n_curve", Rcpp::wrap(n_curve)),
    Rcpp::Named("f_order", Rcpp::wrap(f_order)),
    Rcpp::Named("f_break_points", Rcpp::wrap(f_break_points)),
    Rcpp::Named("n_burn_saem", Rcpp::wrap(n_burn_saem)),
    Rcpp::Named("n_iterations", Rcpp::wrap(n_iterations)),
    Rcpp::Named("n_burn_mcmc", Rcpp::wrap(n_burn_mcmc)),
    Rcpp::Named("n_core", Rcpp::wrap(n_core)),
    Rcpp::Named("saem_step_sizes", Rcpp::wrap(saem_step_sizes))
  );
};



// Return sequence of estimated parameters as R list
Rcpp::List Mixed_Shape_Model::return_pars_tracker(){
  // Scaling matrices
  arma::mat scaling_mat = arma::eye(2, 2);
  scaling_mat(0, 0) = pow(y_scl, 2);
  scaling_mat(1, 1) = 1;
  return Rcpp::List::create(
    Rcpp::Named("alpha_track", Rcpp::wrap(alpha_track * y_scl)),
    Rcpp::Named("sigma2_track", Rcpp::wrap(sigma2_track * std::pow(y_scl, 2))),
    Rcpp::Named("amp_sigma_track", Rcpp::wrap(scaling_mat * amp_sigma_track)),
    Rcpp::Named("p_clusters_track", Rcpp::wrap(p_clusters_track)),
    Rcpp::Named("sampled_m_track", Rcpp::wrap(sampled_m_track))
  );
};



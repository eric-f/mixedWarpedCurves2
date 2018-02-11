// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>
#include "class_mixed_warped_model2.h"
#include "class_mixed_warped_curve2.h"
#include "util.h"

// Constructor
Mixed_Warped_Model::Mixed_Warped_Model(Rcpp::List pars_list,
                                       Rcpp::List control_list,
                                       RcppGSL::Vector f_break_points_r,
                                       RcppGSL::Vector h_break_points_r) :
  f_break_points(f_break_points_r), h_break_points(h_break_points_r){

  // Iteration counter
  saem_counter = 0;

  // Model Parameters to be imported from R
  // Shape parameter and variance of error term
  alpha = Rcpp::as<arma::vec>(pars_list["alpha"]);
  sigma2 = Rcpp::as<double>(pars_list["sigma2"]);
  // Random effect parameters
  mu_a = Rcpp::as<arma::vec>(pars_list["mu"]);
  kappa_id = Rcpp::as<arma::vec>(pars_list["kappa"]);
  n_cluster = Rcpp::as<int>(pars_list["num_clusters"]);

  // Dimension
  dim_f = alpha.size();
  dim_a = mu_a.size();
  dim_w = kappa_id.size() + 1;
  dim_kappa = dim_w - 1;
  dim_z = dim_w - 2;

  // Model Parameters to be initialize in C
  alpha_aug = arma::ones(dim_f + 1);
  sigma2_a = arma::ones(dim_a) * 0.01;
  sigma2_a_mat = arma::eye(dim_a, dim_a) * 0.01;
  sigma2_a_inv = arma::eye(dim_a, dim_a) / 0.01;
  p_clusters = arma::ones(n_cluster) / n_cluster;
  tau1 = 100.0;
  kappa_clusters = arma::zeros(dim_kappa, n_cluster);
  for(generic_idx = 0; generic_idx < n_cluster; ++generic_idx){
    kappa_clusters.col(generic_idx) = kappa_id * tau1;
  }

  // Auxiliary variables
  n_total = Rcpp::as<int>(control_list["n_total"]);
  n_curve = Rcpp::as<int>(control_list["n_curve"]);
  f_order = Rcpp::as<int>(control_list["f_order"]);
  h_order = Rcpp::as<int>(control_list["h_order"]);
  h_left_bound = gsl_vector_min(h_break_points);
  h_right_bound = gsl_vector_max(h_break_points);
  chol_centering_mat = arma::zeros(dim_kappa, dim_z);
  identity_cor_mat = arma::eye(dim_z, dim_z);

  // SA-MCMC Control parameters to be imported from R
  n_core = Rcpp::as<int>(control_list["n_core"]);
  n_burn_saem = Rcpp::as<int>(control_list["n_saem_burn"]);
  n_burn_mcmc = Rcpp::as<int>(control_list["n_mcmc_burn"]);
  n_iterations = Rcpp::as<double>(control_list["n_saem_iter"]) + 2 * n_burn_saem;
  need_centering = Rcpp::as<bool>(control_list["need_centering"]);
  sa_step_size_mod = Rcpp::as<double>(control_list["saem_step_seq_pow"]);

  // Stochastic approximation step sizes
  // Searching stage
  saem_step_sizes = arma::ones(n_iterations);
  // Zone in stage
  saem_step_sizes(arma::span(n_burn_saem, 2 * n_burn_saem - 1)) =
    arma::pow(arma::linspace(1, n_burn_saem, n_burn_saem), -sa_step_size_mod);
  // Averaging stage
  saem_step_sizes(arma::span(2 * n_burn_saem, n_iterations - 1)) =
    arma::pow(arma::linspace(1, n_iterations - 2 * n_burn_saem, n_iterations - 2 * n_burn_saem), -sa_step_size_mod);

  // Variable for Tuning proposal sigma over SAEM burning step
  proposal_sigma = Rcpp::as<double>(control_list["prop_sigma"]);
  calibrate_period = std::min(2 * n_burn_saem, n_iterations);
  mh_accept_rate_lb = Rcpp::as<double>(control_list["accept_rate_lb"]);
  mh_accept_rate_ub = Rcpp::as<double>(control_list["accept_rate_ub"]);
  mh_accept_rate_table = arma::zeros(n_curve, Rcpp::as<int>(control_list["accept_rate_window"]));
  mh_accept_rate_history = arma::zeros(n_iterations);
  proposal_sigma_history = arma::zeros(n_iterations);
  mh_accept_rates_table_ncol = 0;
  mh_accept_rate_table_counter = 0;

  // Temporary variables for tracking MH acceptance rate
  tmp_accept_rate = 0.0;

  // Sufficient Statistics for M-step
  SS_mean_hat_mat = arma::zeros(dim_f + 1, dim_f + 1);
  SS_mean_sigma_a = arma::zeros(dim_a, dim_a);
  SS_mean_log_dw = arma::zeros(dim_kappa, n_cluster);
  SS_mean_cluster_sizes = arma::zeros(n_cluster);

  // Variables for keeping track of cluster membership
  current_m_vec = arma::zeros<arma::ivec>(n_curve);
  // Temporary variables for centering step
  current_a_mat = arma::zeros(dim_a, n_curve);

  // Counters
  saem_counter = 0;
  generic_idx = 0;
  cluster_idx = 0;
  newton_idx = 0;
  newton_inner_idx = 0;

  // Variables for Newton-Raphson steps
  // For full Dirichlet
  newton_q = arma::zeros(dim_kappa);
  newton_g = arma::zeros(dim_kappa);
  newton_z = 0.0;
  newton_b = 0.0;
  newton_update_step_kappa = arma::zeros(dim_kappa);
  newton_new_kappa = arma::zeros(dim_kappa);
  // For Dirichlet with fixed mean
  newton_sum_digamma_k1 = 0.0;
  newton_sum_trigamma_k1 = 0.0;
  newton_grad1 = 0.0;
  newton_hess1 = 0.0;
  newton_update_step_tau1 = 0.0;
  newton_new_tau1 = 0.0;
  // Control parameters
  newton_max_iter = 1000;
  newton_max_inner_iter = 1000;
  newton_abs_tol = 1.0e-8;

  // Trackers
  alpha_track = arma::zeros(dim_f, n_iterations);
  sigma2_track = arma::zeros(n_iterations);
  sigma2_a_track = arma::zeros(dim_a, n_iterations);
  sampled_m_track = arma::zeros<arma::imat>(n_curve, n_iterations);
  kappa_clusters_track = arma::zeros(dim_kappa, n_cluster, n_iterations);

  // Stochastic approximation of Fisher information
  num_pars =  dim_f + 1 + dim_a * 2 + n_cluster * (dim_kappa);
  sapprox_H = arma::zeros(num_pars, num_pars);
  sapprox_C = arma::zeros(num_pars, num_pars);
  sapprox_G = arma::zeros(num_pars);
  current_H = arma::zeros(num_pars, num_pars);
  current_G = arma::zeros(num_pars);

}


// Generate the cholesky decomposition of the I - 1 / (dim_kappa) J matrix with a dimension of (dim_kappa) x (dim_kappa)
// Depends on: dim_kappa
// Changes: chol_centering_mat
void Mixed_Warped_Model::generate_chol_centering_mat(){
  if(dim_w < 2)
    throw("Dimension of the warping function must be bigger than 2");
  for(idx_d = 0; idx_d < dim_z; ++idx_d){
    for(idx_r = idx_d; idx_r < dim_z + 1; ++idx_r){
      if(idx_d == idx_r){
        chol_centering_mat(idx_r, idx_d) = std::sqrt((dim_z - idx_d) / (dim_z - idx_d + 1));
      }
      else{
        chol_centering_mat(idx_r, idx_d) = -std::sqrt(1 / (dim_z - idx_d) / (dim_z - idx_d + 1));
      }
    }
  }
  return;
}


// Keep track of the acceptance rate of the metropolis-hastings sampler and
// adapt the scale parameter of the proposal distribution
// Depends on: saem_counter, mh_accept_rate_table_counter, calibrate_period,
//             proposal_sigma, mh_accept_rate_table, n_burn_mcmc
// Changes: mh_accept_rate_history, proposal_sigma_history,
//          mh_accept_rate_table, proposal_sigma
// Counter increment: mh_accept_rate_table_counter
void Mixed_Warped_Model::track_mh_acceptance_and_calibrate_proposal(){
  // Record acceptance rate
  mh_accept_rate_history(saem_counter) =
    arma::mean(mh_accept_rate_table.col(mh_accept_rate_table_counter)) / n_burn_mcmc;
  proposal_sigma_history(saem_counter) = proposal_sigma;

  // Advance counters
  mh_accept_rate_table_counter = (mh_accept_rate_table_counter + 1) % mh_accept_rate_table.n_cols;

  if(mh_accept_rate_table_counter == 0){
    // Calibrate if still adapting
    if((saem_counter < calibrate_period)){
      tmp_accept_rate = accu(mh_accept_rate_table) / mh_accept_rate_table.n_cols / mh_accept_rate_table.n_rows / n_burn_mcmc;
      if(tmp_accept_rate < mh_accept_rate_lb){
        proposal_sigma /= 2;
      }
      if(tmp_accept_rate > mh_accept_rate_ub){
        proposal_sigma *= 2;
      }
    }
    // Reset table
    mh_accept_rate_table.zeros(mh_accept_rate_table.n_rows, mh_accept_rate_table.n_cols);
  }
  return;
}


void Mixed_Warped_Model::initialize_clustering_with_user_inputs(std::vector<Mixed_Warped_Curve>* mydata){
  Rcpp::Rcout << "Initialize clustering with user inputs..." << std::endl;
  // Re-initialize cluster membership
  cluster_idx = 0;
  for(it = mydata->begin(); it != mydata->end(); ++it){
    it->current_m = it->init_clust;
    it->current_cluster_membership.zeros();
    it->current_cluster_membership(it->current_m) = 1;
    ++cluster_idx;
  }

  // Initialize p_clusters
  p_clusters.zeros();
  for(it = mydata->begin(); it != mydata->end(); ++it){
    ++p_clusters(it->current_m);
  }
  Rcpp::Rcout << "cluster_size" << std::endl << p_clusters << std::endl;
  p_clusters /=  n_curve;
  Rcpp::Rcout << "p_clusters" << std::endl << p_clusters << std::endl;

  return;
}


void Mixed_Warped_Model::gather_sufficient_statistics(std::vector<Mixed_Warped_Curve>* mydata){
  // Reset variables
  SS_mean_hat_mat.zeros();
  SS_mean_sigma_a.zeros();
  SS_mean_log_dw.zeros();
  SS_mean_cluster_sizes.zeros();
  // Gather sufficient statistics
  for(it = mydata->begin(); it != mydata->end(); ++it){
    SS_mean_hat_mat += it->sapprox_hat_mat / n_curve; // Divide by n_curve only to avoid overflow.
    SS_mean_sigma_a += it->sapprox_cov_a / n_curve;
    SS_mean_log_dw += it->sapprox_log_dw / n_curve;
    SS_mean_cluster_sizes += it->sapprox_cluster_membership / n_curve;
  }
  return;
}


// Gather stochastic approximates of the sufficient statistics and perform an M-step - alpha, sigma2
// Depends on:
// Changes:
void Mixed_Warped_Model::update_observation_model(){
  // Rcpp::Rcout << "Update alpha... " << std::endl;
  // Update alpha
  alpha = arma::solve(
    SS_mean_hat_mat(arma::span(1, dim_f), arma::span(1, dim_f)),
    SS_mean_hat_mat(arma::span(1, dim_f), arma::span(0, 0)));
  alpha_aug(arma::span(1, dim_f)) = -alpha;
  // Rcpp::Rcout << "Update sigma2... " << std::endl;
  // Update sigma2
  // Rcpp::Rcout << "dim_f: " << dim_f << std::endl;
  // Rcpp::Rcout << "alpha_aug.size(): " << alpha_aug.size() << std::endl;
  // Rcpp::Rcout << "SS_mean_hat_mat.size(): " << SS_mean_hat_mat.size() << std::endl;
  sigma2 = arma::as_scalar(alpha_aug.t() * SS_mean_hat_mat * alpha_aug) * n_curve / n_total;
  return;
}


// Gather stochastic approximates of the sufficient statistics and perform an M-step - sigma2_a
// Depends on:
// Changes:
void Mixed_Warped_Model::update_amplitude_model(){
  // Update sigma2_a and matrix forms
  sigma2_a_mat = arma::diagmat(SS_mean_sigma_a);
  sigma2_a = sigma2_a_mat.diag();
  sigma2_a_inv = arma::inv(sigma2_a_mat);
  return;
}


// Gather stochastic approximates of the sufficient statistics and perform an M-step -
// Depends on:
// Changes:
void Mixed_Warped_Model::update_primary_warping_component(){
  // Update kappa_clusters.col(0) (by Newton-Raphson)
  if(SS_mean_cluster_sizes(0) > 0){
    for (newton_idx = 0; newton_idx < newton_max_iter; ++newton_idx) {
      // reset temporary variables
      newton_sum_digamma_k1 = 0.0;
      newton_sum_trigamma_k1 = 0.0;
      // compute hessian and gradient
      for(generic_idx = 0; generic_idx < dim_kappa; ++generic_idx){
        newton_sum_digamma_k1 += boost::math::digamma(kappa_id(generic_idx) * tau1) * kappa_id(generic_idx);
        newton_sum_trigamma_k1 += boost::math::trigamma(kappa_id(generic_idx) * tau1) * std::pow(kappa_id(generic_idx), 2);
      }
      // Rcpp::Rcout << "kappa_id: " << kappa_id << std::endl;
      // Rcpp::Rcout << "SS_mean_log_dw.col(0): " << SS_mean_log_dw.col(0) << std::endl;
      // Rcpp::Rcout << "SS_mean_cluster_sizes(0): " << SS_mean_cluster_sizes(0) << std::endl;
      newton_grad1 = arma::as_scalar(kappa_id.t() * SS_mean_log_dw.col(0)) - SS_mean_cluster_sizes(0) * (newton_sum_digamma_k1 - boost::math::digamma(tau1));
      newton_hess1 = - newton_sum_trigamma_k1 + boost::math::trigamma(tau1);
      
      // step size
      newton_update_step_tau1 = newton_grad1 / newton_hess1;
      
      // step halving if out-of-bound
      for(newton_inner_idx = 0; newton_inner_idx < newton_max_inner_iter; ++newton_inner_idx){
        newton_new_tau1 = tau1 - newton_update_step_tau1;
        if(newton_new_tau1 > 0){
          break;
        }
        else{
          newton_update_step_tau1 /= 2.0;
        }
      }

      // check validity and update
      if(newton_new_tau1 > 0){
        tau1 = newton_new_tau1;
        kappa_clusters.col(0) = kappa_id * tau1;
      }
      else{
        Rcpp::Rcout << "Maximum step-halving reached. Estimated tau1 still negative." << std::endl;
        throw;
      }

      // check convergence
      if(std::abs(newton_update_step_tau1) < newton_abs_tol) break;
    }
  }
  else{
      // Rcpp::Rcout << "Empty primary cluster" << std::endl;
  }
  return;
}


// Set other warping components as the same as the primary component
// Depends on:
// Changes:
void Mixed_Warped_Model::synchronize_components(){
  // Copy estimates to other clusters
  for(cluster_idx = 1; cluster_idx < n_cluster; ++cluster_idx){
    kappa_clusters.col(cluster_idx) = kappa_clusters.col(0);
  }
  return;
}


void Mixed_Warped_Model::update_mixture_proportion(){
  // Update p_clusters
  p_clusters = SS_mean_cluster_sizes;
  return;
}



// Gather stochastic approximates of the sufficient statistics and perform an M-step -
// Depends on:
// Changes:
void Mixed_Warped_Model::update_other_warping_components(){
  // Update kappa_clusters (by Newton-Raphson)
  for(cluster_idx = 1; cluster_idx < n_cluster; ++cluster_idx){
    if(SS_mean_cluster_sizes(cluster_idx) > 0){
      for(newton_idx = 0; newton_idx < newton_max_iter; ++newton_idx) {
        // compute hessian and gradient
        newton_z = SS_mean_cluster_sizes(cluster_idx) * boost::math::trigamma(arma::sum(kappa_clusters.col(cluster_idx)));
        newton_g = SS_mean_log_dw.col(cluster_idx) + SS_mean_cluster_sizes(cluster_idx) * boost::math::digamma(arma::sum(kappa_clusters.col(cluster_idx)));
        for(generic_idx = 0; generic_idx < dim_kappa; ++generic_idx) {
          newton_q(generic_idx) = -SS_mean_cluster_sizes(cluster_idx) * boost::math::trigamma(kappa_clusters(generic_idx, cluster_idx));
          newton_g(generic_idx) -= SS_mean_cluster_sizes(cluster_idx) * boost::math::digamma(kappa_clusters(generic_idx, cluster_idx));
        }
        newton_b = arma::sum(newton_g / newton_q) / (1/newton_z + arma::sum(1 / newton_q));

        // step size
        newton_update_step_kappa = (newton_g - arma::ones(dim_kappa) * newton_b) / newton_q;

        // step halving if out-of-bound
        for(newton_inner_idx = 0; newton_inner_idx < newton_max_inner_iter; ++newton_inner_idx){
          newton_new_kappa = kappa_clusters.col(cluster_idx) - newton_update_step_kappa;
          if(arma::all(newton_new_kappa > 0)){
            break;
          }
          else{
            newton_update_step_kappa /= 2;
          }
        }

        // Check validity and update
        if(arma::all(newton_new_kappa > 0)){
          kappa_clusters.col(cluster_idx) = newton_new_kappa;
        }
        else{
          Rcpp::Rcout << "Maximum step-halving reached. Some estimated kappa still negative." << std::endl;
          throw;
        }

        // Check convergence
        if(arma::max(arma::abs(newton_update_step_kappa)) < newton_abs_tol) {
          break;
        }
      }
    }
    else{
      // Rcpp::Rcout << "Empty cluster " << cluster_idx << std::endl;
    }
  }
  return;
}


void Mixed_Warped_Model::update_estimates(){
  // Rcpp::Rcout << "Calling update_observation_model..." << std::endl;
  update_observation_model();
  // Rcpp::Rcout << "Calling update_amplitude_model..." << std::endl;
  update_amplitude_model();
  // Rcpp::Rcout << "Calling update_mixture_proportion..." << std::endl;
  update_mixture_proportion();
  // Rcpp::Rcout << "Calling update_primary_warping_component..." << std::endl;
  update_primary_warping_component();
  if (saem_counter < n_burn_saem) {
    // Rcpp::Rcout << "Calling synchronize_components..." << std::endl;
    synchronize_components();
  }
  else {
    // Rcpp::Rcout << "Calling update_other_warping_components..." << std::endl;
    update_other_warping_components();
  }
  return;
}

// Increase saem iteration counter
// Depends on: Nil
// Changes: saem_counter
void Mixed_Warped_Model::advance_iteration_counter(){
  ++saem_counter;
  return;
}


// Store current estiamtes for tracking
// Depends on: saem_counter, alpha, sigma2, sigma2_a
// Changes: sigma2_track, mu_a_track, sigma2_a_track, sampled_m_track, kappa_clusters_track
void Mixed_Warped_Model::track_estimates(){
  // Track estimates
  alpha_track.col(saem_counter) = alpha;
  sigma2_track(saem_counter) = sigma2;
  sigma2_a_track.col(saem_counter) = sigma2_a;
  sampled_m_track.col(saem_counter) = current_m_vec;
  kappa_clusters_track.slice(saem_counter) = kappa_clusters;
}


// Print estimates for monitoring
// Depends on: saem_counter, proposal_sigma, big_sigma, alpha, sigma2
//             mh_accept_rate_history
// Changes: Nil
void Mixed_Warped_Model::print_estimates(int interval){
  if(saem_counter % interval == 0 & saem_counter < n_iterations){
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "=================================" << std::endl;
    Rcpp::Rcout << "Iteration: " << std::endl << saem_counter << std::endl;
    Rcpp::Rcout << "Acceptance rate: " << mh_accept_rate_history(saem_counter) << std::endl;
    Rcpp::Rcout << "Proposal sigma: " << proposal_sigma << std::endl;
    Rcpp::Rcout << "sigma2: " << sigma2 << std::endl;
    Rcpp::Rcout << "mu_a: " << std::endl << mu_a << std::endl;
    Rcpp::Rcout << "sigma2_a: " << std::endl << sigma2_a << std::endl;
    Rcpp::Rcout << "kappa_clusters: " << std::endl << kappa_clusters.t() << std::endl;
    Rcpp::Rcout << "p_clusters: " << std::endl << p_clusters.t() << std::endl;
    Rcpp::Rcout << "membership state: " << std::endl << current_m_vec.t() << std::endl;
    Rcpp::Rcout << "=================================" << std::endl;
    Rcpp::Rcout << std::endl;
  }
  return;
}


// Return estimated parameters as R list
Rcpp::List Mixed_Warped_Model::return_pars(double y_scaling_factor){
  alpha *= y_scaling_factor;
  sigma2 *= std::pow(y_scaling_factor, 2);
  sigma2_a(0) *= std::pow(y_scaling_factor, 2);
  return Rcpp::List::create(
    Rcpp::Named("alpha", Rcpp::wrap(alpha)),
    Rcpp::Named("sigma2", Rcpp::wrap(sigma2)),
    Rcpp::Named("mu_a", Rcpp::wrap(mu_a)),
    Rcpp::Named("sigma2_a", Rcpp::wrap(sigma2_a)),
    Rcpp::Named("p_clusters", Rcpp::wrap(p_clusters)),
    Rcpp::Named("tau1", Rcpp::wrap(tau1)),
    Rcpp::Named("kappa_id", Rcpp::wrap(kappa_id)),
    Rcpp::Named("kappa_clusters", Rcpp::wrap(kappa_clusters))
  );
};


// Return auxiliary information as R list
// Note: Rcpp::List::create seems to limit the number of items in the list
Rcpp::List Mixed_Warped_Model::return_aux(){
  return Rcpp::List::create(
    Rcpp::Named("n_total", Rcpp::wrap(n_total)),
    Rcpp::Named("n_curve", Rcpp::wrap(n_curve)),
    Rcpp::Named("h_order", Rcpp::wrap(h_order)),
    Rcpp::Named("h_break_points", Rcpp::wrap(h_break_points)),
    Rcpp::Named("f_order", Rcpp::wrap(f_order)),
    Rcpp::Named("f_break_points", Rcpp::wrap(f_break_points)),
    Rcpp::Named("chol_centering_mat", Rcpp::wrap(chol_centering_mat)),
    Rcpp::Named("identity_cor_mat", Rcpp::wrap(identity_cor_mat)),
    Rcpp::Named("n_burn_saem", Rcpp::wrap(n_burn_saem)),
    Rcpp::Named("n_iterations", Rcpp::wrap(n_iterations)),
    Rcpp::Named("mh_accept_rate_history", Rcpp::wrap(mh_accept_rate_history)),
    Rcpp::Named("n_burn_mcmc", Rcpp::wrap(n_burn_mcmc)),
    Rcpp::Named("n_core", Rcpp::wrap(n_core)),
    Rcpp::Named("need_centering", Rcpp::wrap(need_centering)),
    Rcpp::Named("mh_accept_rate_lb", Rcpp::wrap(mh_accept_rate_lb)),
    Rcpp::Named("mh_accept_rate_ub", Rcpp::wrap(mh_accept_rate_ub)),
    Rcpp::Named("calibrate_period", Rcpp::wrap(calibrate_period)),
    Rcpp::Named("proposal_sigma_history", Rcpp::wrap(proposal_sigma_history)),
    Rcpp::Named("kappa_id", Rcpp::wrap(kappa_id)),
    Rcpp::Named("saem_step_sizes", Rcpp::wrap(saem_step_sizes))
  );
};


// Return sequence of estimated parameters as R list
Rcpp::List Mixed_Warped_Model::return_pars_tracker(double y_scaling_factor){
  alpha_track *= y_scaling_factor;
  sigma2_track *= std::pow(y_scaling_factor, 2);
  sigma2_a_track.row(0) *= std::pow(y_scaling_factor, 2);
  return Rcpp::List::create(
    Rcpp::Named("alpha_track", Rcpp::wrap(alpha_track)),
    Rcpp::Named("sigma2_track", Rcpp::wrap(sigma2_track)),
    Rcpp::Named("sigma2_a_track", Rcpp::wrap(sigma2_a_track)),
    Rcpp::Named("sampled_m_track", Rcpp::wrap(sampled_m_track)),
    Rcpp::Named("kappa_clusters_track", Rcpp::wrap(kappa_clusters_track))
  );
};


// Return sequence of estimated parameters as R list
Rcpp::List Mixed_Warped_Model::return_fisher_pieces(double y_scaling_factor){
  return Rcpp::List();
};

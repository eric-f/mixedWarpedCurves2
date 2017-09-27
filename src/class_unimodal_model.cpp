// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>
#include "class_unimodal_model.h"
#include "class_unimodal_curve.h"
#include "util.h"

// Constructor
Unimodal_Model::Unimodal_Model(Rcpp::List pars_list,
                               Rcpp::List control_list,
                               RcppGSL::Vector h_break_points_r) : h_break_points(h_break_points_r){

  // Iteration counter
  saem_counter = 0;

  // Random effect parameters
  mu_a = Rcpp::as<arma::vec>(pars_list["mu"]);
  kappa_id = Rcpp::as<arma::vec>(pars_list["kappa"]);

  // Shape parameter and variance of error term
  sigma2 = Rcpp::as<double>(pars_list["sigma2"]);

  // Dimension
  dim_a = mu_a.size();
  dim_w = kappa_id.size() + 1;
  dim_kappa = dim_w - 1;
  dim_z = dim_w - 2.0;
  n_cluster = Rcpp::as<int>(pars_list["num_clusters"]);

  // Random effect parameters
  sigma2_a = arma::vec(dim_a) * 0.01;
  sigma2_a_mat = arma::eye(dim_a, dim_a) * 0.01;
  sigma2_a_inv = arma::eye(dim_a, dim_a) / 0.01;
  p_clusters = arma::ones(n_cluster) / n_cluster;
  cluster_sizes = arma::zeros(n_cluster);
  kappa_clusters = arma::zeros(dim_kappa, n_cluster);
  for(generic_idx = 0; generic_idx < n_cluster; ++generic_idx){
    kappa_clusters.col(generic_idx) = kappa_id * 100;
  }

  // Auxiliary variables
  n_total = Rcpp::as<int>(control_list["n_total"]);
  n_curve = Rcpp::as<int>(control_list["n_curve"]);
  h_order = Rcpp::as<int>(control_list["h_order"]);
  h_left_bound = gsl_vector_min(h_break_points);
  h_right_bound = gsl_vector_max(h_break_points);
  chol_centering_mat = arma::zeros(dim_kappa, dim_z);
  identity_cor_mat = arma::eye(dim_z, dim_z);

  // SA-MCMC Control parameters
  n_burn_saem = Rcpp::as<int>(control_list["n_saem_burn"]);
  n_iterations = Rcpp::as<double>(control_list["n_saem_iter"]) + 2 * n_burn_saem;
  n_burn_mcmc = Rcpp::as<int>(control_list["n_mcmc_burn"]);
  n_core = Rcpp::as<int>(control_list["n_core"]);
  need_centering = Rcpp::as<bool>(control_list["need_centering"]);
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

  // Variable for Tuning proposal sigma over SAEM burning step
  proposal_sigma = Rcpp::as<double>(control_list["prop_sigma"]);
  calibrate_period = std::min(2 * n_burn_saem, n_iterations);
  mh_accept_rate_lb = Rcpp::as<double>(control_list["accept_rate_lb"]);
  mh_accept_rate_ub = Rcpp::as<double>(control_list["accept_rate_ub"]);
  mh_accept_rate_table = arma::zeros(n_curve, Rcpp::as<int>(control_list["accept_rate_window"]));
  mh_accept_rate_table_counter = 0;
  mh_accept_rate_history = arma::zeros(n_iterations);
  proposal_sigma_history = arma::zeros(n_iterations);

  // Temporary variables for keeping track of cluster membership
  current_m_vec = arma::zeros<arma::ivec>(n_curve);

  // Temporary variables for tracking MH acceptance rate
  tmp_accept_rate = 0.0;

  // Temporary variables for Dirichlet maximization
  tmp_sum_log_dw = arma::zeros(dim_kappa, n_cluster);
  minka_q = arma::zeros(dim_kappa);
  minka_g = arma::zeros(dim_kappa);
  minka_z = 0.0;
  minka_b = 0.0;
  minka_update_step = arma::zeros(dim_kappa);

  // Trackers
  sigma2_track = arma::zeros(n_iterations);
  mu_a_track = arma::zeros(dim_a, n_iterations);
  sigma2_a_track = arma::zeros(dim_a, n_iterations);
  sampled_m_track = arma::zeros<arma::imat>(n_curve, n_iterations);
  kappa_clusters_track = arma::zeros(dim_kappa, n_cluster, n_iterations);

  // Stochastic approximation of Fisher information
  num_pars =  1 + dim_a * 2 + n_cluster * (dim_kappa);
  sapprox_H = arma::zeros(num_pars, num_pars);
  sapprox_C = arma::zeros(num_pars, num_pars);
  sapprox_G = arma::zeros(num_pars);
  current_H = arma::zeros(num_pars, num_pars);
  current_G = arma::zeros(num_pars);
}



// Generate the cholesky decomposition of the I - 1 / (dim_kappa) J matrix with a dimension of (dim_kappa) x (dim_kappa)
// Depends on: dim_kappa
// Changes: chol_centering_mat
void Unimodal_Model::generate_chol_centering_mat(){
  if(dim_w < 2)
    throw("Dimension of the warping function must be bigger than 2");
  for(int d = 0; d < dim_z; ++d){
    for(int idx_r = d; idx_r < dim_z + 1; ++idx_r){
      if(d == idx_r){
        chol_centering_mat(idx_r, d) = std::sqrt((dim_z - d) / (dim_z - d + 1));
      }
      else{
        chol_centering_mat(idx_r, d) = -std::sqrt(1 / (dim_z - d) / (dim_z - d + 1));
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
void Unimodal_Model::track_mh_acceptance_and_calibrate_proposal(){
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




// Gather stochastic approximates of the sufficient statistics and perform an M-step
// Depends on: mydata, n_curve, n_total, kappa
// Changes: sigma2_a, sigma2_a_mat, sigma2_a_inv, sigma2

void Unimodal_Model::update_estimates_data_mod(std::vector<Unimodal_Curve>* mydata){
  // Reset estimates
  sigma2 = 0.0;
  // Gather sufficient statistics
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    sigma2 += it->sapprox_residual_sum_of_squares / n_total;
  }

}

void Unimodal_Model::update_estimates_amp_mod(std::vector<Unimodal_Curve>* mydata){
  // Reset estimates
  mu_a.zeros();
  sigma2_a.zeros();

  // Gather sufficient statistics
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    mu_a += it->sapprox_a / n_curve;
    sigma2_a += it->sapprox_sq_a / n_curve;
  }

  // Update sigma2_a and matrix forms
  sigma2_a -= arma::square(mu_a);
  sigma2_a_mat = arma::diagmat(sigma2_a);
  sigma2_a_inv = arma::inv(sigma2_a_mat);
}

void Unimodal_Model::update_estimates_single_warp_mod(std::vector<Unimodal_Curve>* mydata){
  // Reset estimates
  cluster_sizes.zeros();

  // Reset temporary variable
  tmp_sum_log_dw.zeros();

  // Gather sufficient statistics
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    cluster_sizes += it->sapprox_cluster_membership;
    tmp_sum_log_dw += it->sapprox_log_dw;
  }

  // Update kappa_clusters (by Newton-Raphson)
  // ... for when all curves are in group 0 ...
  newton_max_iter = 1000;
  newton_abs_tol = 1.0e-8;

  cluster_idx = 0;

  for(newton_idx = 0; newton_idx < newton_max_iter; ++newton_idx) {
    minka_z = cluster_sizes(cluster_idx) * boost::math::trigamma(arma::sum(kappa_clusters.col(cluster_idx)));
    minka_g = tmp_sum_log_dw.col(cluster_idx) + cluster_sizes(cluster_idx) * boost::math::digamma(arma::sum(kappa_clusters.col(cluster_idx)));
    for(generic_idx = 0; generic_idx < dim_kappa; ++generic_idx) {
      minka_q(generic_idx) = -cluster_sizes(cluster_idx) * boost::math::trigamma(kappa_clusters(generic_idx, cluster_idx));
      minka_g(generic_idx) -= cluster_sizes(cluster_idx) * boost::math::digamma(kappa_clusters(generic_idx, cluster_idx));
    }
    minka_b = arma::sum(minka_g / minka_q) / (1/minka_z + arma::sum(1 / minka_q));
    minka_update_step = (minka_g - arma::ones(dim_kappa) * minka_b) / minka_q;
    kappa_clusters.col(cluster_idx) -= minka_update_step;
    if(arma::max(arma::abs(minka_update_step)) < newton_abs_tol) {
      break;
    }
  }

  // Copy estimates to other clusters
  ++cluster_idx;
  while(cluster_idx < n_cluster) {
    kappa_clusters.col(cluster_idx) = kappa_clusters.col(0);
    ++cluster_idx;
  }
}

void Unimodal_Model::initialize_clustering_with_user_inputs(std::vector<Unimodal_Curve>* mydata){
  Rcpp::Rcout << "Initialize clustering with user inputs..." << std::endl;
  // Re-initialize cluster membership
  cluster_idx = 0;
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    it->current_m = it->init_clust;
    it->current_cluster_membership.zeros();
    it->current_cluster_membership(it->current_m) = 1;
    ++cluster_idx;
  }

  // Print kmeans cluster membership
  arma::vec cluster_size(n_cluster, arma::fill::zeros);
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    ++cluster_size(it->current_m);
  }
  Rcpp::Rcout << "cluster_size" << std::endl << cluster_size << std::endl;
  p_clusters = cluster_size / n_curve;
  Rcpp::Rcout << "p_clusters" << std::endl << p_clusters << std::endl;

  return;
}

void Unimodal_Model::initialize_clustering_with_kmeans(std::vector<Unimodal_Curve>* mydata){

  Rcpp::Rcout << "Initialize clustering with K-means..." << std::endl;

  // Temporary variables
  arma::mat tmp_sapprox_w_mat(dim_w, n_curve, arma::fill::zeros);
  arma::mat tmp_centroid_w;
  arma::vec tmp_cluster_membership(n_curve, arma::fill::zeros);
  arma::mat tmp_km_xx_mat(n_curve, n_cluster, arma::fill::zeros);
  arma::mat tmp_km_xm_mat(n_curve, n_cluster, arma::fill::zeros);
  arma::mat tmp_km_mm_mat(n_curve, n_cluster, arma::fill::zeros);
  arma::mat tmp_km_dist_mat(n_curve, n_cluster, arma::fill::zeros);
  arma::mat tmp_best_km_dist_mat(n_curve, n_cluster, arma::fill::zeros);
  arma::uvec tmp_best_km_min_dist_idx(n_curve, arma::fill::zeros);

  // Extract warpings
  cluster_idx = 0;
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    tmp_sapprox_w_mat.col(cluster_idx) = it->sapprox_w;
    ++cluster_idx;
  }
  // Rcpp::Rcout << "tmp_sapprox_w_mat:" << std::endl;
  // Rcpp::Rcout << tmp_sapprox_w_mat << std::endl;

  // Kmeans (random start 30 times)
  tmp_min_km_ss = arma::accu(arma::square(tmp_sapprox_w_mat)); // initialize SS at some large value
  tmp_km_xx_mat = arma::square(tmp_sapprox_w_mat).t() * arma::ones(dim_w, n_cluster);
  // Rcpp::Rcout << "tmp_km_xx_mat:" << std::endl;
  // Rcpp::Rcout << tmp_km_xx_mat << std::endl;
  for (int km_iter = 0; km_iter < 30; ++km_iter) {
    rand();
    bool status = arma::kmeans(tmp_centroid_w, tmp_sapprox_w_mat, n_cluster, arma::random_subset, 50, false);
    if (status == false) {
      Rcpp::Rcout << "Kmeans clustering failed" << std::endl;
    }
    else {
      tmp_km_xm_mat = tmp_sapprox_w_mat.t() * tmp_centroid_w;
      // Rcpp::Rcout << "tmp_km_xm_mat:" << std::endl;
      // Rcpp::Rcout << tmp_km_xm_mat << std::endl;
      tmp_km_mm_mat = arma::ones(n_curve, dim_w) * arma::square(tmp_centroid_w);
      // Rcpp::Rcout << "tmp_km_mm_mat:" << std::endl;
      // Rcpp::Rcout << tmp_km_mm_mat << std::endl;
      tmp_km_dist_mat = tmp_km_xx_mat - 2 * tmp_km_xm_mat + tmp_km_mm_mat;
      tmp_new_km_ss = arma::accu(arma::min(tmp_km_dist_mat, 1));
      if (tmp_new_km_ss < tmp_min_km_ss) {
        Rcpp::Rcout << "Found a better kmeans result" << std::endl;
        Rcpp::Rcout << "tmp_min_km_ss:" << std::endl;
        Rcpp::Rcout << tmp_min_km_ss << std::endl;
        Rcpp::Rcout << "tmp_new_km_ss:" << std::endl;
        Rcpp::Rcout << tmp_new_km_ss << std::endl;
        tmp_best_km_dist_mat = tmp_km_dist_mat;
        tmp_min_km_ss = tmp_new_km_ss;
      }
    }
  }

  // Re-initialize cluster membership
  tmp_best_km_min_dist_idx = arma::index_min(tmp_best_km_dist_mat, 1);
  cluster_idx = 0;
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    it->init_clust = tmp_best_km_min_dist_idx(cluster_idx);
    it->current_m = tmp_best_km_min_dist_idx(cluster_idx);
    it->current_cluster_membership.zeros();
    it->current_cluster_membership(tmp_best_km_min_dist_idx(cluster_idx)) = 1;
    ++cluster_idx;
  }

  // Print kmeans cluster membership
  arma::vec cluster_size(n_cluster, arma::fill::zeros);
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    ++cluster_size(it->current_m);
  }
  Rcpp::Rcout << "cluster_size" << std::endl << cluster_size << std::endl;
  p_clusters = cluster_size / n_curve;
  Rcpp::Rcout << "p_clusters" << std::endl << p_clusters << std::endl;

  return;
}



void Unimodal_Model::update_estimates_mixture_warp_mod(std::vector<Unimodal_Curve>* mydata){
  // Reset estimates
  cluster_sizes.zeros();

  // Reset temporary variable
  tmp_sum_log_dw.zeros();

  // Gather sufficient statistics
  for(std::vector<Unimodal_Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    cluster_sizes += it->sapprox_cluster_membership;
    tmp_sum_log_dw += it->sapprox_log_dw;
  }

  // Update p_clusters
  p_clusters = cluster_sizes / n_curve;

  // Update kappa_clusters (by Newton-Raphson)
  newton_max_iter = 1000;
  newton_abs_tol = 1.0e-8;

  for(cluster_idx = 0; cluster_idx < n_cluster; ++cluster_idx) {
    if(cluster_sizes(cluster_idx) > 0){
      for(newton_idx = 0; newton_idx < newton_max_iter; ++newton_idx) {
        minka_z = cluster_sizes(cluster_idx) * boost::math::trigamma(arma::sum(kappa_clusters.col(cluster_idx)));
        minka_g = tmp_sum_log_dw.col(cluster_idx) + cluster_sizes(cluster_idx) * boost::math::digamma(arma::sum(kappa_clusters.col(cluster_idx)));
        for(generic_idx = 0; generic_idx < dim_kappa; ++generic_idx) {
          minka_q(generic_idx) = -cluster_sizes(cluster_idx) * boost::math::trigamma(kappa_clusters(generic_idx, cluster_idx));
          minka_g(generic_idx) -= cluster_sizes(cluster_idx) * boost::math::digamma(kappa_clusters(generic_idx, cluster_idx));
        }
        minka_b = arma::sum(minka_g / minka_q) / (1/minka_z + arma::sum(1 / minka_q));
        minka_update_step = (minka_g - arma::ones(dim_kappa) * minka_b) / minka_q;
        kappa_clusters.col(cluster_idx) -= minka_update_step;
        if(arma::max(arma::abs(minka_update_step)) < newton_abs_tol) {
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




// Update the stoastic approximation of the fisher information
// Depends on:
// Changes:
void Unimodal_Model::update_fisher_information_approx(std::vector<Unimodal_Curve>* mydata){
  return;
}



// Increase saem iteration counter
// Depends on: Nil
// Changes: saem_counter
void Unimodal_Model::advance_iteration_counter(){
  ++saem_counter;
  return;
}



// Store current estiamtes for tracking
// Depends on: saem_counter, alpha, sigma2, sigma2_a
// Changes: sigma2_track, mu_a_track, sigma2_a_track, sampled_m_track, kappa_clusters_track
void Unimodal_Model::track_estimates(){
  // Track estimates
  sigma2_track(saem_counter) = sigma2;
  mu_a_track.col(saem_counter) = mu_a;
  sigma2_a_track.col(saem_counter) = sigma2_a;
  sampled_m_track.col(saem_counter) = current_m_vec;
  kappa_clusters_track.slice(saem_counter) = kappa_clusters;
}



// Print estimates for monitoring
// Depends on: saem_counter, proposal_sigma, big_sigma, alpha, sigma2
//             mh_accept_rate_history
// Changes: Nil
void Unimodal_Model::print_estimates(int interval){
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
Rcpp::List Unimodal_Model::return_pars(double y_scaling_factor){
  return Rcpp::List::create(
    Rcpp::Named("sigma2", Rcpp::wrap(std::pow(y_scaling_factor, 2) * sigma2)),
    Rcpp::Named("mu_a", Rcpp::wrap(y_scaling_factor * mu_a)),
    Rcpp::Named("sigma2_a", Rcpp::wrap(std::pow(y_scaling_factor, 2) * sigma2_a)),
    Rcpp::Named("p_clusters", Rcpp::wrap(p_clusters)),
    Rcpp::Named("kappa_clusters", Rcpp::wrap(kappa_clusters))
  );
};



// Return auxiliary information as R list
// Note: Rcpp::List::create seems to limit the number of items in the list
Rcpp::List Unimodal_Model::return_aux(){
  return Rcpp::List::create(
    Rcpp::Named("n_total", Rcpp::wrap(n_total)),
    Rcpp::Named("n_curve", Rcpp::wrap(n_curve)),
    Rcpp::Named("h_order", Rcpp::wrap(h_order)),
    Rcpp::Named("h_break_points", Rcpp::wrap(h_break_points)),
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
Rcpp::List Unimodal_Model::return_pars_tracker(double y_scaling_factor){
  return Rcpp::List::create(
    Rcpp::Named("sigma2_track", Rcpp::wrap(sigma2_track * std::pow(y_scaling_factor, 2))),
    Rcpp::Named("mu_a_track", Rcpp::wrap(mu_a_track * y_scaling_factor)),
    Rcpp::Named("sigma2_a_track", Rcpp::wrap(sigma2_a_track * std::pow(y_scaling_factor, 2))),
    Rcpp::Named("sampled_m_track", Rcpp::wrap(sampled_m_track)),
    Rcpp::Named("kappa_clusters_track", Rcpp::wrap(kappa_clusters_track))
  );
};



// Return sequence of estimated parameters as R list
Rcpp::List Unimodal_Model::return_fisher_pieces(double y_scaling_factor){
  return Rcpp::List();
};

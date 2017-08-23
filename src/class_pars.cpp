// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>
#include "class_pars.h"
#include "class_curve.h"
#include "util.h"

// Constructor
Pars::Pars(Rcpp::List pars_list,
           Rcpp::List control_list,
           RcppGSL::Vector f_break_points_r,
           RcppGSL::Vector h_break_points_r) : f_break_points(f_break_points_r),
                                               h_break_points(h_break_points_r){

  // Iteration counter
  saem_counter = 0;

  // Fixed parameters
  mu0 = Rcpp::as<arma::vec>(pars_list["mu"]);
  kappa0 = Rcpp::as<arma::vec>(pars_list["kappa"]);

  // Shape parameter and variance of error term
  alpha = Rcpp::as<arma::vec>(pars_list["alpha"]);
  sigma2 = Rcpp::as<double>(pars_list["sigma2"]);

  // Dimension
  dim_a = mu0.size();
  dim_w = kappa0.size() + 1;
  dim_alpha = alpha.size();
  num_clusters = Rcpp::as<int>(pars_list["num_clusters"]);

  // Random effect parameters
  big_sigma = arma::eye(dim_a, dim_a) * 0.01;
  big_sigma_inverse = arma::eye(dim_a, dim_a) * 100;
  tau_clusters = arma::ones(num_clusters) * 100;
  p_clusters = arma::ones(num_clusters) / num_clusters;
  kappa_clusters = arma::zeros(dim_w - 1, num_clusters);
  for(int i = 0; i < num_clusters; ++i){
    kappa_clusters.col(i) = kappa0;
  }

  // Auxiliary variables
  n_total = Rcpp::as<int>(control_list["n_total"]);
  n_curve = Rcpp::as<int>(control_list["n_curve"]);
  f_order = Rcpp::as<int>(control_list["f_order"]);
  h_order = Rcpp::as<int>(control_list["h_order"]);
  f_left_bound = gsl_vector_min(f_break_points);
  f_right_bound = gsl_vector_max(f_break_points);
  h_left_bound = gsl_vector_min(h_break_points);
  h_right_bound = gsl_vector_max(h_break_points);
  chol_centering_mat = arma::zeros(dim_w - 1, dim_w - 2);
  identity_cor_mat = arma::eye(dim_w - 2, dim_w - 2);

  // Diagonal var-cov matrix for amplitude effect
  diag_big_sigma = Rcpp::as<bool>(control_list["ind_amp"]);

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

  // saem_step_sizes = arma::linspace(1, n_iterations, n_iterations) - n_burn_saem;
  // saem_step_sizes = arma::clamp(saem_step_sizes, 1, arma::datum::inf);
  // saem_step_sizes = arma::pow(saem_step_sizes, -sa_step_size_mod);

  // Variable for Tuning proposal sigma over SAEM burning step
  proposal_sigma = Rcpp::as<double>(control_list["prop_sigma"]);
  calibrate_period = std::min(2 * n_burn_saem, n_iterations);
  mh_accept_rate_lb = Rcpp::as<double>(control_list["accept_rate_lb"]);
  mh_accept_rate_ub = Rcpp::as<double>(control_list["accept_rate_ub"]);
  mh_accept_rate_table = arma::zeros(n_curve, Rcpp::as<int>(control_list["accept_rate_window"]));
  mh_accept_rate_table_counter = 0;
  mh_accept_rate_history = arma::zeros(n_iterations);
  proposal_sigma_history = arma::zeros(n_iterations);

  // Temporary variables for centering step
  current_a_mat = arma::zeros(dim_a, n_curve);
  current_m_vec = arma::zeros<arma::ivec>(n_curve);

  // Trackers
  alpha_track = arma::zeros(dim_alpha, n_iterations);
  sigma2_track = arma::zeros(n_iterations);
  big_sigma_track = arma::zeros(dim_a * dim_a, n_iterations);
  tau_clusters_track = arma::zeros(num_clusters, n_iterations);
  kappa_clusters_track = arma::zeros(dim_w - 1, num_clusters, n_iterations);
  sampled_m_track = arma::zeros<arma::imat>(n_curve, n_iterations);

  // Stochastic approximation of Fisher information
  num_pars = dim_alpha + 5 + (num_clusters - 1) * (dim_w - 1);
  if(diag_big_sigma)
    num_pars -= 1;
  sapprox_H = arma::zeros(num_pars, num_pars);
  sapprox_C = arma::zeros(num_pars, num_pars);
  sapprox_G = arma::zeros(num_pars);
  current_H = arma::zeros(num_pars, num_pars);
  current_G = arma::zeros(num_pars);
}



// Generate the cholesky decomposition of the I - 1/(dim_w - 1) J matrix with a dimension of (dim_w-1)x(dim_w-1)
// Depends on: dim_w
// Changes: chol_centering_mat
void Pars::generate_chol_centering_mat(){
  double dim_z = dim_w - 2.0;
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
void Pars::track_mh_acceptance_and_calibrate_proposal(){
  double tmp_accept_rate;
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
// Changes: big_sigma, big_sigma_inverse, alpha, sigma2, tau
void Pars::update_parameter_estimates(std::vector<Curve>* mydata){
  // Temporary variables
  arma::mat tmp_mean_sigma_a(dim_a, dim_a, arma::fill::zeros);
  arma::mat tmp_scaled_hat_mat(dim_alpha + 1, dim_alpha + 1, arma::fill::zeros);
  arma::mat tmp_mean_hat_By(dim_alpha, 1, arma::fill::zeros);
  arma::mat tmp_mean_hat_BB(dim_alpha, dim_alpha, arma::fill::zeros);
  arma::vec tmp_alpha_aug(dim_alpha + 1, arma::fill::ones);
  arma::mat tmp_mean_log_dw(dim_w-1, num_clusters, arma::fill::zeros);
  arma::vec tmp_num_in_clusters(num_clusters, arma::fill::zeros);

  int newton_max_iter = 1000;
  double newton_update_step;
  arma::vec newton_update_step_vec;
  double tmp_log_tau;
  arma::vec tmp_log_tau_kappa;
  arma::vec tmp_tau_kappa;

  // Gather sufficient statistics
  for(std::vector<Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    tmp_mean_sigma_a += it->sapprox_sigma_a / n_curve;
    tmp_scaled_hat_mat += it->sapprox_hat_mat / n_curve; // Divide by n_curve only to avoid overflow.
    tmp_mean_log_dw += it->sapprox_log_dw;
    tmp_num_in_clusters += it->sapprox_cluster_membership;
  }
  for(int i = 0; i < tmp_mean_log_dw.n_cols; ++i){
    tmp_mean_log_dw.col(i) = tmp_mean_log_dw.col(i) /  tmp_num_in_clusters(i);
  }

  // Update big_sigma
  if(diag_big_sigma){
    big_sigma = arma::diagmat(tmp_mean_sigma_a);
    big_sigma_inverse = arma::inv(arma::diagmat(big_sigma));
  }
  else{
    big_sigma = tmp_mean_sigma_a;
    big_sigma_inverse = arma::inv_sympd(big_sigma);
  }

  // Update alpha
  tmp_mean_hat_By = tmp_scaled_hat_mat(arma::span(1, dim_alpha), arma::span(0, 0));
  tmp_mean_hat_BB = tmp_scaled_hat_mat(arma::span(1, dim_alpha), arma::span(1, dim_alpha));
  alpha = arma::solve(tmp_mean_hat_BB, tmp_mean_hat_By);

  // Update sigma2
  tmp_alpha_aug(arma::span(1, dim_alpha)) = -alpha;
  sigma2 = arma::as_scalar(tmp_alpha_aug.t() * tmp_scaled_hat_mat * tmp_alpha_aug) * n_curve / n_total;

  // Update tau_1 with fixed kappa (by Newton-Raphson)
  if(tmp_num_in_clusters(0) > 0){
    tmp_log_tau = std::log(tau_clusters(0));
    for (int i = 0; i < newton_max_iter; ++i) {
      newton_update_step = newton_step_dirichlet_fixed_mean(tmp_log_tau, kappa_clusters.col(0), tmp_mean_log_dw.col(0));
      tmp_log_tau -= newton_update_step;
      if(std::abs(newton_update_step) < 1e-6) break; // convergence by step-size
    }
    tau_clusters(0) = std::exp(tmp_log_tau);
    if(tau_clusters(0) != tau_clusters(0)) throw("estimate of tau blown up...");
  }

  // Update tau_m with free kappa if num_cluster > 1
  if(num_clusters > 1){
    for(int cluster_idx = 1; cluster_idx < num_clusters; ++cluster_idx){
      if(tmp_num_in_clusters(cluster_idx) > 0){
        tmp_log_tau_kappa = log(tau_clusters(cluster_idx) * kappa_clusters.col(cluster_idx));
        for (int i = 0; i < newton_max_iter; ++i) {
          newton_update_step_vec = newton_step_dirichlet_free_mean(tmp_log_tau_kappa, tmp_mean_log_dw.col(cluster_idx));
          tmp_log_tau_kappa -= newton_update_step_vec;
          if(max(abs(newton_update_step_vec)) < 1e-6) break; // convergence by step-size
        }
        tmp_tau_kappa = exp(tmp_log_tau_kappa);
        tau_clusters(cluster_idx) = sum(tmp_tau_kappa);
        kappa_clusters.col(cluster_idx) = tmp_tau_kappa / tau_clusters(cluster_idx);
      }
      else{
        // Rcpp::Rcout << "Warning: Empty cluster..." << std::endl;
      }
    }
  }

  // Update p_clusters
  p_clusters = tmp_num_in_clusters / n_curve;

  return;
}



// Update the stoastic approximation of the fisher information
// Depends on:
// Changes: SA_H, SA_G, SA_C
void Pars::update_fisher_information_approx(std::vector<Curve>* mydata){
  double current_step_size = saem_step_sizes(saem_counter);
  // Temporary variables
  arma::mat tmp_sum_sigma_a(dim_a, dim_a, arma::fill::zeros);
  arma::mat tmp_sum_hat_mat(dim_alpha + 1, dim_alpha + 1, arma::fill::zeros);
  arma::mat tmp_sum_hat_yy(1, 1, arma::fill::zeros);
  arma::mat tmp_sum_hat_By(dim_alpha, 1, arma::fill::zeros);
  arma::mat tmp_sum_hat_BB(dim_alpha, dim_alpha, arma::fill::zeros);
  arma::mat tmp_sum_log_dw(dim_w-1, num_clusters, arma::fill::zeros);
  arma::vec tmp_num_in_clusters(num_clusters, arma::fill::zeros);
  double tmp_rss = 0;
  arma::mat D2;
  if(diag_big_sigma){
    D2 = arma::zeros(4, 2);
    D2(0, 0) = 1;
    D2(3, 1) = 1;
  }
  else{
    D2 = arma::zeros(4, 3);
    D2(0, 0) = 1;
    D2(1, 1) = 1;
    D2(2, 1) = 1;
    D2(3, 2) = 1;
  }

  if (num_clusters == 1) {
    current_H.zeros();
    current_G.zeros();

    // Gather sufficient statistics
    // (Duplicating code chunk in update_parameter_estimates)
    for(std::vector<Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
      tmp_sum_sigma_a += it->current_sigma_a;
      tmp_sum_hat_mat += it->current_hat_mat;
      tmp_sum_log_dw += it->current_log_dw.col(0);
    }
    tmp_sum_hat_yy = tmp_sum_hat_mat(arma::span(0, 0), arma::span(0, 0));
    tmp_sum_hat_By = tmp_sum_hat_mat(arma::span(1, dim_alpha), arma::span(0, 0));
    tmp_sum_hat_BB = tmp_sum_hat_mat(arma::span(1, dim_alpha), arma::span(1, dim_alpha));

    // Prep digamma and trigamma
    arma::vec tau_kappa0 = tau_clusters(0) * kappa0;
    arma::vec digamma_tau_kappa0(dim_w - 1, arma::fill::zeros);
    arma::vec trigamma_tau_kappa0(dim_w - 1, arma::fill::zeros);
    for (int i = 0; i < dim_w - 1; ++i) {
      digamma_tau_kappa0(i) = boost::math::digamma(tau_kappa0(i));
      trigamma_tau_kappa0(i) = boost::math::trigamma(tau_kappa0(i));
    }

    // G1: Gradient wrt alpha
    current_G(arma::span(0, dim_alpha - 1)) = (tmp_sum_hat_By - tmp_sum_hat_BB * alpha) / sigma2;
    // G2: Gradient wrt sigma2
    tmp_rss = arma::as_scalar(tmp_sum_hat_yy - 2 * tmp_sum_hat_By.t() * alpha + alpha.t() * tmp_sum_hat_BB * alpha);
    current_G(dim_alpha) =
      (tmp_rss / pow(sigma2, 2) - n_total / sigma2) / 2;
    // G3: Gradient wrt vech(big_sigma)
    current_G(arma::span(dim_alpha + 1, num_pars - 2)) =
      D2.t() * arma::vectorise(big_sigma_inverse * tmp_sum_sigma_a * big_sigma_inverse - n_curve * big_sigma_inverse) / 2;
    // G4: Gradient wrt tau0
    current_G(num_pars - 1) =
      arma::as_scalar(kappa0.t() * tmp_sum_log_dw -
      n_curve * (kappa0.t() * digamma_tau_kappa0 - boost::math::digamma(tau_clusters(0))));

    // H1: Hessian alpha alpha
    current_H(arma::span(0, dim_alpha - 1), arma::span(0, dim_alpha - 1)) =
      - tmp_sum_hat_BB / sigma2;
    // H2: Hessian alpha sigma2
    current_H(arma::span(0, dim_alpha - 1), arma::span(dim_alpha, dim_alpha)) =
      - (tmp_sum_hat_By - tmp_sum_hat_BB * alpha) / pow(sigma2, 2);
    current_H(arma::span(dim_alpha, dim_alpha), arma::span(0, dim_alpha - 1)) =
      current_H(arma::span(0, dim_alpha - 1), arma::span(dim_alpha, dim_alpha)).t();
    // H3: Hessian sigma2 sigma2
    current_H(dim_alpha, dim_alpha) =
      (n_total / 2 / pow(sigma2, 2) - tmp_rss / pow(sigma2, 3));
    // H4: Hessian big_sigma big_sigma
    current_H(arma::span(dim_alpha + 1, num_pars - 2), arma::span(dim_alpha + 1, num_pars - 2)) =
      - n_curve / 2 * D2.t() *
      kron(big_sigma_inverse,
           big_sigma_inverse * (2 * tmp_sum_sigma_a / n_curve - big_sigma) *
             big_sigma_inverse) * D2;
    // H5: Hessian tau0 tau0
    current_H(num_pars - 1, num_pars - 1) =
      - n_curve * arma::as_scalar(square(kappa0).t() * trigamma_tau_kappa0 -
      boost::math::trigamma(tau_clusters(0)));

    sapprox_H = (1 - current_step_size) * sapprox_H + current_step_size * current_H;
    sapprox_G = (1 - current_step_size) * sapprox_G + current_step_size * current_G;
    sapprox_C = (1 - current_step_size) * sapprox_C + current_step_size * current_G * current_G.t();
  }
  else{
    // not implemented for num_cluster > 1
  }
  return;
}



// Increase saem iteration counter
// Depends on: Nil
// Changes: saem_counter
void Pars::advance_iteration_counter(){
  ++saem_counter;
  return;
}



// Store current estiamtes for tracking
// Depends on: saem_counter, alpha, sigma2, big_sigma, tau
// Changes: alpha_track, sigma2_track, big_sigma_track, tau_track
void Pars::track_estimates(){
  // Track estimates
  alpha_track.col(saem_counter) = alpha;
  sigma2_track(saem_counter) = sigma2;
  big_sigma_track.col(saem_counter) = arma::vectorise(big_sigma);
  tau_clusters_track.col(saem_counter) = tau_clusters;
  kappa_clusters_track.slice(saem_counter) = kappa_clusters;
  sampled_m_track.col(saem_counter) = current_m_vec;
}



// Print estimates for monitoring
// Depends on: saem_counter, proposal_sigma, big_sigma, alpha, sigma2, tau
//             mh_accept_rate_history
// Changes: Nil
void Pars::print_estimates(int interval){
  if(saem_counter % interval == 0 & saem_counter < n_iterations){
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "=================================" << std::endl;
    Rcpp::Rcout << "Iteration: " << std::endl << saem_counter << std::endl;
    Rcpp::Rcout << "Acceptance rate: " << mh_accept_rate_history(saem_counter) << std::endl;
    Rcpp::Rcout << "Proposal sigma: " << proposal_sigma << std::endl;
    Rcpp::Rcout << "big_sigma: " << std::endl << big_sigma << std::endl;
    Rcpp::Rcout << "alpha: " << std::endl << alpha.t() << std::endl;
    Rcpp::Rcout << "sigma2: " << sigma2 << std::endl;
    Rcpp::Rcout << "kappa_clusters: " << std::endl << kappa_clusters.t() << std::endl;
    Rcpp::Rcout << "tau_clusters: " << std::endl << tau_clusters.t() << std::endl;
    Rcpp::Rcout << "p_clusters: " << std::endl << p_clusters.t() << std::endl;
    Rcpp::Rcout << "membership state: " << std::endl << current_m_vec.t() << std::endl;
    Rcpp::Rcout << "=================================" << std::endl;
    Rcpp::Rcout << std::endl;
  }
  return;
}



// Return estimated parameters as R list
Rcpp::List Pars::return_pars(){
  return Rcpp::List::create(
    Rcpp::Named("mu0", Rcpp::wrap(mu0)),
    Rcpp::Named("kappa0", Rcpp::wrap(kappa0)),
    Rcpp::Named("alpha", Rcpp::wrap(alpha)),
    Rcpp::Named("sigma2", Rcpp::wrap(sigma2)),
    Rcpp::Named("big_sigma", Rcpp::wrap(big_sigma)),
    Rcpp::Named("tau_clusters", Rcpp::wrap(tau_clusters)),
    Rcpp::Named("p_clusters", Rcpp::wrap(p_clusters)),
    Rcpp::Named("kappa_clusters", Rcpp::wrap(kappa_clusters))
  );
};



// Return estimated parameters as R list
Rcpp::List Pars::return_pars(double y_scaling_factor){
  // Scaling matrices
  arma::mat a_scaling_mat = arma::eye(2, 2);
  a_scaling_mat(0, 0) = y_scaling_factor;
  return Rcpp::List::create(
    Rcpp::Named("mu0", Rcpp::wrap(mu0)),
    Rcpp::Named("kappa0", Rcpp::wrap(kappa0)),
    Rcpp::Named("alpha", Rcpp::wrap(alpha * y_scaling_factor)),
    Rcpp::Named("sigma2", Rcpp::wrap(sigma2 * std::pow(y_scaling_factor, 2))),
    Rcpp::Named("big_sigma", Rcpp::wrap(a_scaling_mat * big_sigma * a_scaling_mat)),
    Rcpp::Named("p_clusters", Rcpp::wrap(p_clusters)),
    Rcpp::Named("tau_clusters", Rcpp::wrap(tau_clusters)),
    Rcpp::Named("kappa_clusters", Rcpp::wrap(kappa_clusters))
  );
};



// Return auxiliary information as R list
// Note: Rcpp::List::create seems to limit the number of items in the list
Rcpp::List Pars::return_aux(){
  return Rcpp::List::create(
    Rcpp::Named("n_total", Rcpp::wrap(n_total)),
    Rcpp::Named("n_curve", Rcpp::wrap(n_curve)),
    Rcpp::Named("f_order", Rcpp::wrap(f_order)),
    Rcpp::Named("h_order", Rcpp::wrap(h_order)),
    Rcpp::Named("f_break_points", Rcpp::wrap(f_break_points)),
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
    Rcpp::Named("saem_step_sizes", Rcpp::wrap(saem_step_sizes))
  );
};



// Return sequence of estimated parameters as R list
Rcpp::List Pars::return_pars_tracker(){
  return Rcpp::List::create(
    Rcpp::Named("alpha_track", Rcpp::wrap(alpha_track)),
    Rcpp::Named("sigma2_track", Rcpp::wrap(sigma2_track)),
    Rcpp::Named("big_sigma_track", Rcpp::wrap(big_sigma_track)),
    Rcpp::Named("tau_clusters_track", Rcpp::wrap(tau_clusters_track)),
    Rcpp::Named("kappa_clusters_track", Rcpp::wrap(kappa_clusters_track)),
    Rcpp::Named("sampled_m_track", Rcpp::wrap(sampled_m_track))
  );
};



// Return sequence of estimated parameters as R list
Rcpp::List Pars::return_pars_tracker(double y_scaling_factor){
  // Scaling matrices
  arma::mat scaling_mat = arma::eye(4, 4);
  scaling_mat(0, 0) = pow(y_scaling_factor, 2);
  scaling_mat(1, 1) = y_scaling_factor;
  scaling_mat(2, 2) = 1;
  return Rcpp::List::create(
    Rcpp::Named("alpha_track", Rcpp::wrap(alpha_track * y_scaling_factor)),
    Rcpp::Named("sigma2_track", Rcpp::wrap(sigma2_track * std::pow(y_scaling_factor, 2))),
    Rcpp::Named("big_sigma_track", Rcpp::wrap(scaling_mat * big_sigma_track)),
    Rcpp::Named("tau_clusters_track", Rcpp::wrap(tau_clusters_track)),
    Rcpp::Named("kappa_clusters_track", Rcpp::wrap(kappa_clusters_track)),
    Rcpp::Named("sampled_m_track", Rcpp::wrap(sampled_m_track))
  );
};



// Return sequence of estimated parameters as R list
Rcpp::List Pars::return_fisher_pieces(double y_scaling_factor){
  arma::mat scaling_mat(num_pars, num_pars, arma::fill::zeros);
  arma::mat varcov(num_pars, num_pars);

  // Construct rescaling matrix
  for(int idx = 0; idx < dim_alpha; ++idx){
    scaling_mat(idx, idx) = y_scaling_factor;
  }
  scaling_mat(dim_alpha, dim_alpha) = pow(y_scaling_factor, 2); // sigma2
  scaling_mat(dim_alpha + 1, dim_alpha + 1) = pow(y_scaling_factor, 2); // big_sigma
  scaling_mat(dim_alpha + 2, dim_alpha + 2) = y_scaling_factor;
  scaling_mat(num_pars - 2, num_pars - 2) = 1;
  scaling_mat(num_pars - 1, num_pars - 1) = 1; // tau_0

  // Compute variance covariance matrix
  varcov = scaling_mat * arma::inv_sympd(- sapprox_H + sapprox_C - sapprox_G * sapprox_G.t()) * scaling_mat;
  return Rcpp::List::create(
    Rcpp::Named("SA_H", Rcpp::wrap(sapprox_H)),
    Rcpp::Named("SA_C", Rcpp::wrap(sapprox_C)),
    Rcpp::Named("SA_G", Rcpp::wrap(sapprox_G)),
    Rcpp::Named("scaling_mat", Rcpp::wrap(scaling_mat)),
    Rcpp::Named("varcov", Rcpp::wrap(varcov))
  );
};

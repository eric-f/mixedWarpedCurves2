// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
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
  mu = Rcpp::as<arma::vec>(pars_list["mu"]);
  kappa = Rcpp::as<arma::vec>(pars_list["kappa"]);
  alpha = Rcpp::as<arma::vec>(pars_list["alpha"]);

  // Dimension
  dim_a = mu.size();
  dim_w = kappa.size() + 1;
  dim_alpha = alpha.size();

  // Model parameters
  sigma2 = Rcpp::as<double>(pars_list["sigma2"]);
  // big_sigma = Rcpp::as<arma::mat>(pars_list["Sigma"]);
  big_sigma = arma::eye(dim_a, dim_a) * 0.01;
  big_sigma_inverse = arma::eye(dim_a, dim_a) * 100;
  // tau = Rcpp::as<double>(pars_list["tau"]);
  tau = 1000;

  // Auxiliary variables
  n_total = Rcpp::as<int>(control_list["n_total"]);
  n_curve = Rcpp::as<int>(control_list["n_curve"]);
  f_order = Rcpp::as<int>(control_list["f_order"]);
  h_order = Rcpp::as<int>(control_list["h_order"]);
  // f_break_points = Rcpp::as<RcppGSL::Vector>(control_list["f_break_points"]); // Need to explicitly initiate above
  // h_break_points = Rcpp::as<RcppGSL::Vector>(control_list["h_break_points"]); // Need to explicitly initiate above
  f_left_bound = gsl_vector_min(f_break_points);
  f_right_bound = gsl_vector_max(f_break_points);
  h_left_bound = gsl_vector_min(h_break_points);
  h_right_bound = gsl_vector_max(h_break_points);
  // chol_centering_mat = Rcpp::as<arma::mat>(control_list["L"]); // chol_centering_mat
  chol_centering_mat = arma::zeros(dim_w - 1, dim_w - 2);
  identity_cor_mat = arma::eye(dim_w - 2, dim_w - 2);

  // SA-MCMC Control parameters
  n_burn_saem = Rcpp::as<int>(control_list["n_saem_burn"]);
  n_iterations = Rcpp::as<double>(control_list["n_saem_iter"]) + n_burn_saem;
  n_burn_mcmc = Rcpp::as<int>(control_list["n_mcmc_burn"]);
  n_core = Rcpp::as<int>(control_list["n_core"]);
  need_centering = Rcpp::as<bool>(control_list["need_centering"]);
  sa_step_size_mod = Rcpp::as<double>(control_list["saem_step_seq_pow"]);

  // Generate step sizes
  saem_step_sizes = clamp(arma::linspace(1, n_iterations, n_iterations) - n_burn_saem, 1, arma::datum::inf);
  saem_step_sizes = pow(saem_step_sizes, -sa_step_size_mod);

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

  // Trackers
  alpha_track = arma::zeros(dim_alpha, n_iterations);
  sigma2_track = arma::zeros(n_iterations);
  big_sigma_track = arma::zeros(dim_a * dim_a, n_iterations);
  tau_track = arma::zeros(n_iterations);
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
        chol_centering_mat(idx_r, d) = std::sqrt((dim_z-d)/(dim_z-d+1));
      }
      else{
        chol_centering_mat(idx_r, d) = -std::sqrt(1/(dim_z-d)/(dim_z-d+1));
      }
    }
  }
  return;
}



// Keep track of the acceptance rate of the metropolis-hastings sampler and
// adapt the scale parameter of the proposal distribution
// Depends on: saem_counter, n_iterations,
//             mh_accept_rate_history,
// Changes: ...
// Counter increment: mh_accept_rate_table_counter,
//                    saem_counter
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
// Depends on:
// Changes:
void Pars::update_parameter_estimates(std::vector<Curve>* mydata){
  // Temporary variables
  arma::mat tmp_mean_sigma_a(dim_a, dim_a, arma::fill::zeros);
  arma::mat tmp_scaled_hat_mat(dim_alpha + 1, dim_alpha + 1, arma::fill::zeros);
  arma::vec tmp_sum_log_dw(dim_w-1, arma::fill::zeros);

  arma::mat tmp_mean_hat_By(dim_alpha, 1, arma::fill::zeros);
  arma::mat tmp_mean_hat_BB(dim_alpha, dim_alpha, arma::fill::zeros);
  arma::vec tmp_alpha_aug(dim_alpha + 1, arma::fill::ones);

  int newton_max_iter = 1000;
  double newton_update_step;
  arma::vec newton_workspace;
  double tmp_log_tau;

  // Gather sufficient statistics
  for(std::vector<Curve>::iterator it = mydata->begin(); it != mydata->end(); ++it){
    tmp_mean_sigma_a += it->sapprox_sigma_a / n_curve;
    tmp_scaled_hat_mat += it->sapprox_hat_mat / n_curve; // Divide by n_curve only to avoid overflow.
    tmp_sum_log_dw += it->sapprox_log_dw;
  }

  // Update big_sigma
  big_sigma = tmp_mean_sigma_a;
  big_sigma_inverse = arma::inv_sympd(big_sigma);

  // Update alpha
  tmp_mean_hat_By = tmp_scaled_hat_mat(arma::span(1, dim_alpha), arma::span(0, 0));
  tmp_mean_hat_BB = tmp_scaled_hat_mat(arma::span(1, dim_alpha), arma::span(1, dim_alpha));
  alpha = arma::solve(tmp_mean_hat_BB, tmp_mean_hat_By);

  // Update sigma2
  tmp_alpha_aug(arma::span(1, dim_alpha)) = -alpha;
  sigma2 = arma::as_scalar(tmp_alpha_aug.t() * tmp_scaled_hat_mat * tmp_alpha_aug) * n_curve / n_total;

  // Update tau (by Newton-Raphson)
  tmp_log_tau = std::log(tau);
  for (int i = 0; i < newton_max_iter; i++) {
    newton_workspace = nllk_dirichlet(tmp_log_tau, tmp_sum_log_dw, n_curve, kappa);
    newton_update_step = newton_workspace(1) / newton_workspace(2);
    tmp_log_tau -= newton_update_step;
    if(std::abs(newton_update_step) < 1e-6) break; // convergence by step-size
  }
  tau = std::exp(tmp_log_tau);
  if(tau != tau) throw("estimate of tau blown up...");
  return;
}



// Increase saem iteration counter
// Depends on: Nil
// Changes: saem_counter
void Pars::advance_iteration_counter(){
  ++saem_counter;
  return;
}



void Pars::track_estimates(){
  // Track estimates
  alpha_track.col(saem_counter) = alpha;
  sigma2_track(saem_counter) = sigma2;
  big_sigma_track.col(saem_counter) = arma::vectorise(big_sigma);
  tau_track(saem_counter) = tau;
}



// Print estimates for monitoring
// Depends on: saem_counter, proposal_sigma, big_sigma, alpha, sigma2, tau
//             mh_accept_rate_history
// Changes: Nil
void Pars::print_estimates(int interval){
  if(saem_counter % interval == 0){
    Rcpp::Rcout << "=================================" << std::endl;
    Rcpp::Rcout << "Iteration: " << std::endl << saem_counter << std::endl;
    Rcpp::Rcout << "Acceptance rate: " << mh_accept_rate_history(saem_counter) << std::endl;
    Rcpp::Rcout << "Proposal sigma: " << proposal_sigma << std::endl;
    Rcpp::Rcout << "big_sigma: " << arma::vectorise(big_sigma).t() << std::endl;
    Rcpp::Rcout << "alpha: " << alpha.t() << std::endl;
    Rcpp::Rcout << "sigma2: " << sigma2 << std::endl;
    Rcpp::Rcout << "new tau: " << tau << std::endl;
    Rcpp::Rcout << "=================================" << std::endl;
  }
  return;
}



// Return estimated parameters as R list
Rcpp::List Pars::return_pars(){
  return Rcpp::List::create(
    Rcpp::Named("mu", Rcpp::wrap(mu)),
    Rcpp::Named("kappa", Rcpp::wrap(kappa)),
    Rcpp::Named("alpha", Rcpp::wrap(alpha)),
    Rcpp::Named("sigma2", Rcpp::wrap(sigma2)),
    Rcpp::Named("big_sigma", Rcpp::wrap(big_sigma)),
    Rcpp::Named("tau", Rcpp::wrap(tau))
  );
};



// Return auxiliary information as R list
Rcpp::List Pars::return_aux(){
  return Rcpp::List::create(
    Rcpp::Named("n_total", Rcpp::wrap(n_total)),
    Rcpp::Named("n_curve", Rcpp::wrap(n_curve)),
    // Rcpp::Named("dim_a", Rcpp::wrap(dim_a)),
    // Rcpp::Named("dim_w", Rcpp::wrap(dim_w)),
    // Rcpp::Named("dim_alpha", Rcpp::wrap(dim_alpha)),
    Rcpp::Named("f_order", Rcpp::wrap(f_order)),
    Rcpp::Named("h_order", Rcpp::wrap(h_order)),
    Rcpp::Named("f_break_points", Rcpp::wrap(f_break_points)),
    Rcpp::Named("h_break_points", Rcpp::wrap(h_break_points)),
    Rcpp::Named("chol_centering_mat", Rcpp::wrap(chol_centering_mat)),
    Rcpp::Named("identity_cor_mat", Rcpp::wrap(identity_cor_mat)),
    Rcpp::Named("n_burn_saem", Rcpp::wrap(n_burn_saem)),
    Rcpp::Named("n_iterations", Rcpp::wrap(n_iterations)),
    Rcpp::Named("n_burn_mcmc", Rcpp::wrap(n_burn_mcmc)),
    Rcpp::Named("n_core", Rcpp::wrap(n_core)),
    Rcpp::Named("need_centering", Rcpp::wrap(need_centering)),
    Rcpp::Named("sa_step_size_mod", Rcpp::wrap(sa_step_size_mod)),
    Rcpp::Named("mh_accept_rate_lb", Rcpp::wrap(mh_accept_rate_lb)),
    Rcpp::Named("mh_accept_rate_ub", Rcpp::wrap(mh_accept_rate_ub)),
    // Rcpp::Named("mh_accept_rate_table", Rcpp::wrap(mh_accept_rate_table)),
    Rcpp::Named("calibrate_period", Rcpp::wrap(calibrate_period)),
    // Rcpp::Named("proposal_sigma", Rcpp::wrap(proposal_sigma)),
    Rcpp::Named("mh_accept_rate_history", Rcpp::wrap(mh_accept_rate_history)),
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
    Rcpp::Named("tau_track", Rcpp::wrap(tau_track))
  );
};

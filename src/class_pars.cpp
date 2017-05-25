// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include "class_pars.h"
#include "class_curve.h"
#include "util.h"

// Constructor
Pars::Pars(Rcpp::List pars_list,
           Rcpp::List control_list,
           RcppGSL::Vector break_points) : f_break_points(break_points){

  // Extract R objects
  Rcpp::List aux = pars_list["aux"];

  // Iteration counter
  saem_counter = 0;

  // Model Parameters
  alpha = Rcpp::as<arma::vec>(pars_list["alp"]);
  sigma2 = Rcpp::as<double>(pars_list["sigma2"]);
  mu = Rcpp::as<arma::vec>(pars_list["mu"]);
  big_sigma = Rcpp::as<arma::mat>(pars_list["Sigma"]);
  kappa = Rcpp::as<arma::vec>(pars_list["kappa"]);
  tau = Rcpp::as<double>(pars_list["tau"]);

  // Dimension
  dim_a = mu.size();
  dim_w = kappa.size() + 1;
  dim_alpha = alpha.size();

  // Auxiliary variables
  n_total = Rcpp::as<int>(aux["n_total"]);
  n_curve = Rcpp::as<int>(aux["n_curve"]);
  f_order = Rcpp::as<int>(aux["f_order"]);
  f_full_knots = Rcpp::as<arma::vec>(aux["f_full_knots"]);
  // f_break_points = Rcpp::as<RcppGSL::Vector>(aux["f_break_points"]); // Need to explicitly initiate above
  f_left_bound = gsl_vector_min(f_break_points);
  f_right_bound = gsl_vector_max(f_break_points);
  big_sigma_inverse = Rcpp::as<arma::mat>(aux["Sigma_inv"]);
  chol_centering_mat = Rcpp::as<arma::mat>(aux["L"]); // chol_centering_mat
  identity_cor_mat = arma::eye(dim_w - 2, dim_w - 2);

  // SA-MCMC Control parameters
  n_burn_saem = Rcpp::as<int>(control_list["nBurnSAEM"]);
  n_iterations = Rcpp::as<double>(control_list["nIter"]) + n_burn_saem;
  n_burn_mcmc = Rcpp::as<int>(control_list["nBurnMCMC"]);
  n_core = Rcpp::as<int>(control_list["nCore"]);
  need_centering = Rcpp::as<bool>(control_list["centering"]);
  sa_step_size_mod = Rcpp::as<double>(control_list["alphaSAEM"]);

  // Generate step sizes
  saem_step_sizes = clamp(arma::linspace(1, n_iterations, n_iterations) - n_burn_saem, 1, arma::datum::inf);
  saem_step_sizes = pow(saem_step_sizes, -sa_step_size_mod);

  // Variable for Tuning proposal sigma over SAEM burning step
  proposal_sigma = Rcpp::as<double>(control_list["prop_sigma"]);
  calibrate_period = std::min(2 * n_burn_saem, n_iterations);
  mh_accept_rate_lb = Rcpp::as<double>(control_list["accept_rate_lb"]);
  mh_accept_rate_ub = Rcpp::as<double>(control_list["accept_rate_ub"]);
  mh_accept_rate_table = arma::zeros(n_curve, Rcpp::as<int>(control_list["n_accept_rates"]));
  mh_accept_rate_table_counter = 0;
  mh_accept_rate_history = arma::zeros(n_iterations);
  mh_accept_rate_history_counter = 0;

  // Temporary variables for centering step
  current_a_mat = arma::zeros(dim_a, n_curve);
}



// Keep track of the acceptance rate of the metropolis-hastings sampler and
// adapt the scale parameter of the proposal distribution
// Depends on: mh_accept_rate_history_counter, n_iterations,
//             mh_accept_rate_history,
// Changes: ...
// Counter increment: mh_accept_rate_history_counter,
//                    mh_accept_rate_table_counter,
//                    saem_counter
void Pars::post_simulation_housekeeping(){
  double tmp_accept_rate;
  // Record acceptance rate
  if(mh_accept_rate_history_counter < n_iterations){
    mh_accept_rate_history(mh_accept_rate_history_counter) =
      arma::mean(mh_accept_rate_table.col(mh_accept_rate_table_counter)) / n_curve / n_burn_mcmc;

    // Table is full, calibrate and reset tallies (TO-DO: print and check table)
    if((mh_accept_rate_table_counter % mh_accept_rate_table.n_cols) == (mh_accept_rate_table.n_cols - 1)){
      // Calibrate if still adapting
      if(saem_counter < calibrate_period){
        tmp_accept_rate = accu(mh_accept_rate_table) / mh_accept_rate_table.n_cols / mh_accept_rate_table.n_rows / n_burn_mcmc;
        if(tmp_accept_rate < mh_accept_rate_lb){
          proposal_sigma /= 2;
        }
        if(tmp_accept_rate > mh_accept_rate_ub){
          proposal_sigma *= 2;
        }
      }
      // Reset tallies
      mh_accept_rate_table.zeros();
    }
    // Advance counters
    ++mh_accept_rate_history_counter;
    mh_accept_rate_table_counter = (mh_accept_rate_table_counter + 1) % mh_accept_rate_table.n_cols;
    return;
  }
  return;
}



// Gather stochastic approximates of the sufficient statistics and perform an M-step
// Depends on:
// Changes:
void Pars::update_parameter_estimates(std::vector<Curve>* mydata){
  // Temporary variables
  arma::mat tmp_mean_sigma_a(dim_a, dim_a, arma::fill::zeros);
  arma::mat tmp_mean_hat_mat(dim_alpha + 1, dim_alpha + 1, arma::fill::zeros);
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
    tmp_mean_hat_mat += it->sapprox_hat_mat / n_curve;
    tmp_sum_log_dw += it->sapprox_log_dw;
  }

  // Update big_sigma
  big_sigma = tmp_mean_sigma_a;
  big_sigma_inverse = arma::inv_sympd(big_sigma);

  // Update alpha
  tmp_mean_hat_By = tmp_mean_hat_mat(arma::span(1, dim_alpha), arma::span(0, 0));
  tmp_mean_hat_BB = tmp_mean_hat_mat(arma::span(1, dim_alpha), arma::span(1, dim_alpha));
  alpha = arma::solve(tmp_mean_hat_BB, tmp_mean_hat_By);

  // Update sigma2
  tmp_alpha_aug(arma::span(1, dim_alpha)) = alpha;
  sigma2 = arma::as_scalar(tmp_alpha_aug.t() * tmp_mean_hat_mat * tmp_alpha_aug * n_curve / n_total);

  // Update tau (by Newton-Raphson)
  tmp_log_tau = std::log(tau);
  for (int i = 0; i < newton_max_iter; i++) {
    newton_workspace = nllk_dirichlet_rcpp(tmp_log_tau, tmp_sum_log_dw, n_curve, kappa);
    newton_update_step = newton_workspace(1) / newton_workspace(2);
    tmp_log_tau -= newton_update_step;
    if(std::abs(newton_update_step) < 1e-6) break; // convergence by step-size
  }
  tau = std::exp(tmp_log_tau);
  if(tau != tau) throw("estimate of tau blown up..."); // convergence by step-size

}



// Increase saem iteration counter
// Depends on: Nil
// Changes: saem_counter
void Pars::advance_iteration_counter(){
  ++saem_counter;
}



// Print estimates for monitoring
// Depends on: saem_counter, proposal_sigma, big_sigma, alpha, sigma2, tau
//             mh_accept_rate_history, mh_accept_rate_history_counter
// Changes: Nil
void Pars::print_estimates(int interval){
  if(saem_counter % interval == 0){
    Rcpp::Rcout << "=================================" << std::endl;
    Rcpp::Rcout << "Iteration: " << std::endl << saem_counter << std::endl;
    Rcpp::Rcout << "Acceptance rate: " << mh_accept_rate_history(mh_accept_rate_history_counter-1) << std::endl;
    Rcpp::Rcout << "Proposal sigma: " << proposal_sigma << std::endl;
    Rcpp::Rcout << "big_sigma: " << arma::vectorise(big_sigma).t() << std::endl;
    Rcpp::Rcout << "alpha: " << alpha.t() << std::endl;
    Rcpp::Rcout << "sigma2: " << sigma2 << std::endl;
    Rcpp::Rcout << "new tau: " << tau << std::endl;
    Rcpp::Rcout << "=================================" << std::endl;
  }
}



// Return as R list
Rcpp::List Pars::return_list(){
  return Rcpp::List::create(
    Rcpp::Named("alpha", Rcpp::wrap(alpha)),
    Rcpp::Named("sigma2", Rcpp::wrap(sigma2)),
    Rcpp::Named("big_sigma", Rcpp::wrap(big_sigma)),
    Rcpp::Named("tau", Rcpp::wrap(tau))
  );
};

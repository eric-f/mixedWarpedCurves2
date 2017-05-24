// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// #include "class_model.h"
//
// // Constructor
// Model::Model(Rcpp::List pars,
//              Rcpp::List curve_list,
//              Rcpp::List control_list){
//
//   // Extract R objects
//   Rcpp::List aux = pars["aux"];
//
//   Rcpp::Rcout << "Importing model parameters" << std::endl;
//
//   // Model Parameters
//   alpha = Rcpp::as<arma::vec>(pars["alp"]);
//   sigma2 = Rcpp::as<double>(pars["sigma2"]);
//   mu = Rcpp::as<arma::vec>(pars["mu"]);
//   big_sigma = Rcpp::as<arma::mat>(pars["Sigma"]);
//   kappa = Rcpp::as<arma::vec>(pars["kappa"]);
//   tau = Rcpp::as<double>(pars["tau"]);
//
//   Rcpp::Rcout << "Importing auxiliary variables" << std::endl;
//
//   // Auxiliary variables
//   n_total = Rcpp::as<int>(aux["n_total"]);
//   n_curve = Rcpp::as<int>(aux["n_curve"]);
//   f_order = Rcpp::as<int>(aux["f_order"]);
//   f_full_knots = Rcpp::as<arma::vec>(aux["f_full_knots"]);
//   f_break_points = Rcpp::as<arma::vec>(aux["f_break_points"]);
//   big_sigma_inverse = Rcpp::as<arma::mat>(aux["Sigma_inv"]);
//   D = Rcpp::as<int>(aux["D"]);
//   L = Rcpp::as<arma::mat>(aux["L"]);
//
//   Rcpp::Rcout << "Importing control variables" << std::endl;
//
//   // SA-MCMC Control parameters
//   sa_step_size_mod = Rcpp::as<double>(control_list["alphaSAEM"]);
//   n_burn_saem = Rcpp::as<int>(control_list["nBurnSAEM"]);
//   calibrate_period = 2 * n_burn_saem;
//   n_iterations = Rcpp::as<double>(control_list["nIter"]);
//   n_burn_mcmc = Rcpp::as<int>(control_list["nBurnMCMC"]);
//   n_core = Rcpp::as<int>(control_list["nCore"]);
//   proposal_sigma = Rcpp::as<double>(control_list["prop_sigma"]);
//   need_centering = Rcpp::as<bool>(control_list["centering"]);
//
//   // Variable for Tuning proposal sigma over SAEM burning step
//   mh_accept_rate_lb = Rcpp::as<double>(control_list["accept_rate_lb"]);
//   mh_accept_rate_ub = Rcpp::as<double>(control_list["accept_rate_ub"]);
//   mh_accept_rates_window = Rcpp::as<int>(control_list["n_accept_rates"]);
//   mh_accept_rate_rolling_avg = Rcpp::as<arma::vec>(control_list["accept_rate"]);
//   mh_accept_rate_buffer = Rcpp::as<arma::vec>(control_list["past_accept_rates"]);
//   mh_accept_rate_counter = 0;
//
//   Rcpp::Rcout << "Importing curve objects" << std::endl;
//
//   // Curve objects
//   for(Rcpp::List::iterator it = curve_list.begin();
//       it != curve_list.end(); ++it) {
//     curves.push_back(Curve(*it, alpha.size(), L, proposal_sigma));
//   }
// }

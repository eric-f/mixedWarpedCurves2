// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include "class_curve.h"
#include "class_pars.h"


// [[Rcpp::export]]
void test_oopc1(Rcpp::List curve_list,
                Rcpp::List pars_list,
                Rcpp::List control_list) {


  // Special handling of f_break_points as an RcppGSL::Vector object to be propergated to the classes.
  RcppGSL::Vector f_break_points = Rcpp::as< RcppGSL::vector<double> >(Rcpp::as<Rcpp::List>(pars_list["aux"])["f_break_points"]);
  RcppGSL::Vector h_break_points = Rcpp::as< RcppGSL::vector<double> >(Rcpp::as<Rcpp::List>(pars_list["aux"])["h_break_points"]);

  Pars pars(pars_list, control_list, f_break_points, h_break_points);
  pars.generate_chol_centering_mat();

  std::vector<Curve> data;
  int id = 0;
  for(Rcpp::List::iterator it = curve_list.begin(); it != curve_list.end(); ++it) {
    data.push_back(Curve(*it, &pars, id));
    ++id;
  }

  Rcpp::Rcout << "Initialize current_f_basis_mat " << std::endl;
  for(std::vector<Curve>::iterator it = data.begin(); it != data.end(); ++it){
    for(int j = 0; j < 5; ++j){
      it->initialize_current_f_basis_mat();
    }
  }

  // for(std::vector<Curve>::iterator it = data.begin();
  //     it != data.end(); ++it) {
  //   Rcpp::Rcout << it->current_a << std::endl;
  // }

  Rcpp::Rcout << "Check SAEM step sizes..." << std::endl;
  Rcpp::Rcout << pars.saem_step_sizes;
  //
  // Rcpp::Rcout << "Check separation of R and Cpp data..." << std::endl;
  // Rcpp::Rcout << "data[0].current_a: " << std::endl << data[0].current_a << std::endl;
  // Rcpp::Rcout << "Set  data[0].current_a = arma::zeros(2)..." << std::endl;
  // data[0].current_a = arma::zeros(2);
  // Rcpp::Rcout << "data[0].current_a: " << std::endl << data[0].current_a << std::endl;
  // Rcpp::Rcout << "curve_list[0]['re_mc_state']: " << std::endl << Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(curve_list[0])["re_mc_state"]) << std::endl;
  //
  //
  // Rcpp::Rcout << "Testing pointer to common pars..." << std::endl;
  // Rcpp::Rcout << "Curve 1 proposal sigma: " << std::endl << data[0].common_pars->proposal_sigma << std::endl;
  // Rcpp::Rcout << "Curve 2 proposal sigma: " << std::endl << data[1].common_pars->proposal_sigma << std::endl;
  // Rcpp::Rcout << "Common proposal sigma: " << std::endl << pars.proposal_sigma << std::endl;
  // Rcpp::Rcout << "Change pars.proposal_sigma = 1000" << std::endl << std::endl;
  // pars.proposal_sigma = 1000;
  // Rcpp::Rcout << "Curve 1 proposal sigma: " << std::endl << data[0].common_pars->proposal_sigma << std::endl;
  // Rcpp::Rcout << "Curve 2 proposal sigma: " << std::endl << data[1].common_pars->proposal_sigma << std::endl;
  // Rcpp::Rcout << "Common proposal sigma: " << std::endl << pars.proposal_sigma << std::endl;
  // Rcpp::Rcout << "Change data[0].common_pars->proposal_sigma = 0.2" << std::endl << std::endl;
  // data[0].common_pars->proposal_sigma = 0.2;
  // Rcpp::Rcout << "Curve 1 proposal sigma: " << std::endl << data[0].common_pars->proposal_sigma << std::endl;
  // Rcpp::Rcout << "Curve 2 proposal sigma: " << std::endl << data[1].common_pars->proposal_sigma << std::endl;
  // Rcpp::Rcout << "Common proposal sigma: " << std::endl << pars.proposal_sigma << std::endl;

  // Rcpp::Rcout << "Testing propose_new_w()" << std::endl;
  // Rcpp::Rcout << "dim_w: " << data[1].dim_w << std::endl;
  // Rcpp::Rcout << "chol_centering_mat: " << std::endl << data[1].common_pars->chol_centering_mat << std::endl;
  // Rcpp::Rcout << "identity_cor_mat: " << std::endl << data[1].common_pars->identity_cor_mat << std::endl;
  // Rcpp::Rcout << "propose_w: " << std::endl << data[0].proposed_w << std::endl;
  // data[0].propose_new_w();
  // Rcpp::Rcout << "propose_w: " << std::endl << data[0].proposed_w << std::endl;
  //
  // Rcpp::Rcout << "Testing draw_new_a()" << std::endl;
  // Rcpp::Rcout << "current_a: " << std::endl << data[0].current_a << std::endl;
  // data[0].draw_new_a();
  // Rcpp::Rcout << "current_a: " << std::endl << data[0].current_a << std::endl;
  //
  // Rcpp::Rcout << "Testing B-spline component" << std::endl;
  // Rcpp::Rcout << pars.f_break_points << std::endl;
  // Rcpp::Rcout << data[0].common_pars->f_break_points << std::endl;
  // Rcpp::Rcout << "data[0].dim_alpha: " << std::endl;
  // Rcpp::Rcout << data[0].dim_alpha << std::endl;
  // Rcpp::Rcout << "data[0].current_warped_f_basis_mat: " << std::endl;
  // Rcpp::Rcout << data[0].current_warped_f_basis_mat(arma::span(0,3), arma::span(0,1)) << std::endl;
  // Rcpp::Rcout << "data[0].proposed_warped_f_basis_mat: " << std::endl;
  // Rcpp::Rcout << data[0].proposed_warped_f_basis_mat(arma::span(0,3), arma::span(0,1)) << std::endl;
  // Rcpp::Rcout << "Update basis evaluation matrix...: " << std::endl;
  // data[0].compute_proposed_warping_and_f_basis_mat();
  // Rcpp::Rcout << "data[0].warped_f_basis_mat: " << std::endl;
  // Rcpp::Rcout << data[0].proposed_warped_f_basis_mat(arma::span(0,3), arma::span(0,3)) << std::endl;
  // Rcpp::Rcout << "Propose new w...: " << std::endl;
  // data[0].propose_new_w();
  // data[0].compute_proposed_warping_and_f_basis_mat();
  // Rcpp::Rcout << "data[0].warped_f_basis_mat: " << std::endl;
  // Rcpp::Rcout << data[0].proposed_warped_f_basis_mat(arma::span(0,3), arma::span(0,3)) << std::endl;
  // Rcpp::Rcout << "Propose new a...: " << std::endl;
  // data[0].draw_new_a();
  // data[0].compute_proposed_warping_and_f_basis_mat();
  // Rcpp::Rcout << "data[0].warped_f_basis_mat: " << std::endl;
  // Rcpp::Rcout << data[0].proposed_warped_f_basis_mat(arma::span(0,3), arma::span(0,3)) << std::endl;

  Rcpp::Rcout << "Test SAEM" << std::endl;
  for(int sim_idx = 0; sim_idx < 100; ++sim_idx){
    for(std::vector<Curve>::iterator it = data.begin(); it != data.end(); ++it){
      // Rcpp::Rcout << "Simulation" << std::endl;
      it->do_simulation_step();
      // Rcpp::Rcout << "Centering" << std::endl;
      it->center_current_a();
      Rcpp::Rcout << "Curve id: " << it->curve_id << std::endl;
      Rcpp::Rcout << "Centered a: " << std::endl << it->current_a << std::endl;
      // Rcpp::Rcout << "Suff. Stat." << std::endl;
      it->update_sufficient_statistics_approximates();
    }
    // Rcpp::Rcout << "Acceptance rate" << std::endl;
    pars.track_mh_acceptance_and_calibrate_proposal();
    // Rcpp::Rcout << "Maximization" << std::endl;
    pars.update_parameter_estimates(&data);
    // Rcpp::Rcout << "Counter increment" << std::endl;
    pars.advance_iteration_counter();
  }

  Rcpp::Rcout << "alpha: " << std::endl << pars.alpha << std::endl;
  Rcpp::Rcout << "sigma2: " << std::endl << pars.sigma2 << std::endl;
  Rcpp::Rcout << "big_sigma: " << std::endl << pars.big_sigma << std::endl;
  Rcpp::Rcout << "tau: " << std::endl << pars.tau << std::endl;
}





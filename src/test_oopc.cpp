// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include "class_curve.h"
#include "class_pars.h"


// [[Rcpp::export]]
void test_oopc(Rcpp::List curve_list,
               Rcpp::List pars_list,
               Rcpp::List control_list) {

  // Special handling of f_break_points as an RcppGSL::Vector object to be propergated to the classes.
  RcppGSL::Vector f_break_points = Rcpp::as< RcppGSL::vector<double> >(Rcpp::as<Rcpp::List>(pars_list["aux"])["f_break_points"]);
  Pars pars(pars_list, control_list, f_break_points);
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

  Rcpp::Rcout << "Test SAEM" << std::endl;
  for(int sim_idx = 0; sim_idx < pars.n_iterations; ++sim_idx){
    for(std::vector<Curve>::iterator it = data.begin(); it != data.end(); ++it){
      it->do_simulation_step();
    }
    for(std::vector<Curve>::iterator it = data.begin(); it != data.end(); ++it){
      it->center_current_a();
      it->update_sufficient_statistics_approximates();
    }
    pars.post_simulation_housekeeping();
    pars.update_parameter_estimates(&data);
    pars.advance_iteration_counter();
  }

  Rcpp::Rcout << "alpha: " << std::endl << pars.alpha << std::endl;
  Rcpp::Rcout << "sigma2: " << std::endl << pars.sigma2 << std::endl;
  Rcpp::Rcout << "big_sigma: " << std::endl << pars.big_sigma << std::endl;
  Rcpp::Rcout << "tau: " << std::endl << pars.tau << std::endl;
}





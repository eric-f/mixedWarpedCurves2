// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "class_curve.h"
#include "class_pars.h"


// [[Rcpp::export]]
Rcpp::List test_oopc(Rcpp::List curve_list,
               Rcpp::List pars_list,
               Rcpp::List control_list) {

  // Special handling of f_break_points as an RcppGSL::Vector object to be propergated to the classes.
  RcppGSL::Vector f_break_points = Rcpp::as< RcppGSL::vector<double> >(Rcpp::as<Rcpp::List>(pars_list["aux"])["f_break_points"]);

  // Import parameters
  Pars pars(pars_list, control_list, f_break_points);

  // Import data
  std::vector<Curve> data;
  int id = 0;
  for(Rcpp::List::iterator it = curve_list.begin(); it != curve_list.end(); ++it) {
    data.push_back(Curve(*it, &pars, id));
    ++id;
  }

  // Initialize warped base shape basis evaluation matrices
  for(std::vector<Curve>::iterator it = data.begin(); it != data.end(); ++it){
    for(int j = 0; j < 5; ++j){
      it->initialize_current_f_basis_mat();
    }
  }

  // Run SAEM
  int tick = std::max(1, (int) pars.n_iterations / 20);
  Progress p(0, false);
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
    // pars.print_estimates(10);
    if (pars.saem_counter % tick == 0) {
      double progress = (double) pars.saem_counter / pars.n_iterations * 100;
      Rprintf("%3.1f%%...", progress);
    }
    if (Progress::check_abort())
      return R_NilValue;
  }

  // Return as R object
  Rcpp::List curves(pars.n_curve);
  Rcpp::List fit;

  for(int i = 0; i < pars.n_curve; ++i){
    curves[i] = data[i].return_list();
  }

  fit = pars.return_list();

  return(Rcpp::List::create(
    Rcpp::Named("fit", fit),
    Rcpp::Named("curves", curves)));

}





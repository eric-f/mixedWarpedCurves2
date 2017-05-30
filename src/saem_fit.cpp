// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include "class_curve.h"
#include "class_pars.h"
#include "ctime"

//' Internal function for fitting the model by SAEM
//'
//' @param curve_list List of curve objects, each should have a component called data, with subcomponents, x and y,
//'        for the time and value of the functional data
//' @param pars_list List with components: mu, kappa, tau and sigma2
//' @param control_list List generated by control_saem()
//' @param y_scaling_factor Numeric value for scaling back the estimates and predictions to the original scale of y
// [[Rcpp::export]]
Rcpp::List saem_fit(Rcpp::List curve_list,
                    Rcpp::List pars_list,
                    Rcpp::List control_list,
                    double y_scaling_factor) {

  int seed0 = std::time(0);

  // Special handling of f_break_points and h_break_points as
  // an RcppGSL::Vector object to be propergated to the classes
  RcppGSL::Vector f_break_points = Rcpp::as< RcppGSL::vector<double> >(control_list["f_knots"]);
  RcppGSL::Vector h_break_points = Rcpp::as< RcppGSL::vector<double> >(control_list["h_knots"]);

  // Import parameters
  Pars* pars = new Pars(pars_list, control_list, f_break_points, h_break_points);
  pars->generate_chol_centering_mat();

  // Import data
  std::vector<Curve>* data = new std::vector<Curve>;
  int id = 0;
  for(Rcpp::List::iterator it = curve_list.begin(); it != curve_list.end(); ++it) {
    data->push_back(Curve(*it, pars, id, seed0));
    ++id;
    ++seed0;
  }

  // Initialize basis evaluation matrices
  for(std::vector<Curve>::iterator it = data->begin(); it != data->end(); ++it){
      it->initialize_h_basis_mat();
      it->initialize_current_f_basis_mat();
  }

  // Setup multi-threading
#ifdef _OPENMP
  omp_set_num_threads(pars->n_core);
  REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif


  // Setup progress tracking
  int tick = std::max(1, (int) pars->n_iterations / 20);
  Progress p(0, false);
  double progress;
  int saem_idx = 0;
  int curve_idx = 0;

  // SAEM loop
  for(saem_idx = 0; saem_idx < pars->n_iterations; ++saem_idx){
    // Simulation step
#pragma omp parallel for
    for(curve_idx = 0; curve_idx < data->size(); ++curve_idx){
      data->at(curve_idx).do_simulation_step();
    }
    // Centering and stochastic approximation step
    for(curve_idx = 0; curve_idx < data->size(); ++curve_idx){
      data->at(curve_idx).center_current_a();
      data->at(curve_idx).update_sufficient_statistics_approximates();
    }
    // MH Calibration
    pars->track_mh_acceptance_and_calibrate_proposal();
    // Maximization step
    pars->update_parameter_estimates(data);
    pars->track_estimates();
    pars->advance_iteration_counter();
    if (pars->saem_counter % tick == 0) {
      progress = (double) pars->saem_counter / pars->n_iterations * 100;
      Rprintf("%3.1f%%...", progress);
      // pars.print_estimates(10);
    }
    if (Progress::check_abort())
      saem_idx = data->size(); // Bump index to exit the loop
  }
  Rcpp::Rcout << "(Done)" << std::endl;


  // Wrap up as R objects
  Rcpp::List curves(pars->n_curve);
  Rcpp::List fit;
  Rcpp::List aux;
  Rcpp::List pars_track;
  for(int i = 0; i < pars->n_curve; ++i){
    curves[i] = data->at(i).return_list(y_scaling_factor);
  }
  fit = pars->return_pars(y_scaling_factor);
  aux = pars->return_aux();
  pars_track = pars->return_pars_tracker(y_scaling_factor);

  // Free memory
  free(data);
  free(pars);

  // Return to R
  return(Rcpp::List::create(
      Rcpp::Named("fit", fit),
      Rcpp::Named("curves", curves),
      Rcpp::Named("aux", aux),
      Rcpp::Named("pars_track", pars_track)));

}





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
#include "class_unimodal_curve.h"
#include "class_unimodal_model.h"

//' Internal function for fitting the model by SAEM
//'
//' @param curve_list List of curve objects, each should have a component called data, with subcomponents, x and y,
//'        for the time and value of the functional data
//' @param pars_list List with components: mu, kappa and sigma2
//' @param control_list List generated by control_saem()
//' @param y_scaling_factor Numeric value for scaling back the estimates and predictions to the original scale of y
//' @param trace Logical value, if TRUE tracing information of the estimated parameters are printed
//' @noRd
// [[Rcpp::export]]
Rcpp::List saem_fit_unimodal(Rcpp::List curve_list,
                             Rcpp::List pars_list,
                             Rcpp::List control_list,
                             double y_scaling_factor,
                             bool trace) {

  int seed = Rcpp::as<int>(control_list["seed"]);

  // Special handling of h_break_points as
  // an RcppGSL::Vector object to be propergated to the classes
  RcppGSL::Vector h_break_points = Rcpp::as< RcppGSL::vector<double> >(control_list["h_knots"]);

  // Import parameters
  Unimodal_Model* pars = new Unimodal_Model(pars_list, control_list, h_break_points);
  pars->generate_chol_centering_mat();

  // Import data
  std::vector<Unimodal_Curve>* data = new std::vector<Unimodal_Curve>;
  int id = 0;
  for(Rcpp::List::iterator it = curve_list.begin(); it != curve_list.end(); ++it) {
    data->push_back(Unimodal_Curve(*it, pars, id, seed));
    ++id;
  }

  // Initialize basis evaluation matrices
  for(std::vector<Unimodal_Curve>::iterator it = data->begin(); it != data->end(); ++it){
    it->initialize_h_basis_mat();
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

  Rcpp::Rcout << "Entering SAEM loop..." << std::endl;
  // SAEM loop
  for(saem_idx = 0; saem_idx < pars->n_iterations; ++saem_idx){
    // Simulation step
#pragma omp parallel for
    for(curve_idx = 0; curve_idx < data->size(); ++curve_idx){
      data->at(curve_idx).do_simulation_step();
    }
    // Initialize clustering configuration by K-mean
    if(pars->saem_counter == pars->n_burn_saem){
      // pars->initialize_clustering_with_kmeans(data);
      pars->initialize_clustering_with_user_inputs(data);
    }
    // Stochastic approximation step
#pragma omp parallel for
    for(curve_idx = 0; curve_idx < data->size(); ++curve_idx){
      data->at(curve_idx).update_sufficient_statistics_approximates();
    }
    // MH Calibration
    pars->track_mh_acceptance_and_calibrate_proposal();
    // Maximization step
    pars->update_estimates_data_mod(data);
    pars->update_estimates_amp_mod(data);
    if (pars->saem_counter < pars->n_burn_saem) {
      pars->update_estimates_single_warp_mod(data);
    }
    else {
      pars->update_estimates_mixture_warp_mod(data);
    }
    // Stochastic approximation to Fisher information matrix
    pars->track_estimates();
    // Progress report
    if (pars->saem_counter % tick == 0) {
      progress = (double) pars->saem_counter / pars->n_iterations * 100;
      Rprintf("%3.1f%%...", progress);
      if(trace){
        pars->print_estimates(1);
      }
    }
    pars->advance_iteration_counter();
    if (Progress::check_abort())
      saem_idx = pars->n_iterations; // Bump index to exit the loop
  }
  Rcpp::Rcout << "(Done)" << std::endl;


  // Wrap up as R objects
  Rcpp::List curves(pars->n_curve);
  Rcpp::List est_pars;
  Rcpp::List aux;
  Rcpp::List pars_track;
  Rcpp::List se_info;
  for(int i = 0; i < pars->n_curve; ++i){
    curves[i] = data->at(i).return_list(y_scaling_factor);
  }
  est_pars = pars->return_pars(y_scaling_factor);
  aux = pars->return_aux();
  pars_track = pars->return_pars_tracker(y_scaling_factor);
  // se_info = pars->return_fisher_pieces(y_scaling_factor);

  // Free memory
  free(data);
  free(pars);

  // Return to R
  return(Rcpp::List::create(
      Rcpp::Named("pars", est_pars),
      Rcpp::Named("curves", curves),
      Rcpp::Named("aux", aux),
      Rcpp::Named("pars_track", pars_track),
      Rcpp::Named("se_info", se_info)));

}





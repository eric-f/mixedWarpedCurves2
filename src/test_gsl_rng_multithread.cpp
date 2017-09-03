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

//' Internal function for fitting the model by SAEM
//'
//' @export
//[[Rcpp::export]]
void test_gsl_rng_multithread(Rcpp::List curve_list,
                              Rcpp::List pars_list,
                              Rcpp::List control_list,
                              double y_scaling_factor) {

  int seed = Rcpp::as<int>(control_list["seed"]);

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
    data->push_back(Curve(*it, pars, id, seed));
    ++id;
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

  // Draw random number
#pragma omp parallel for
  for(int curve_idx = 0; curve_idx < data->size(); ++curve_idx){
    data->at(curve_idx).print_random_number();
  }

}

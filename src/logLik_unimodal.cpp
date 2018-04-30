// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include "class_unimodal_lite.h"

// [[Rcpp::export]]
Rcpp::List logLik_unimodal(Rcpp::List pars_list,
                           Rcpp::List curve_list,
                           Rcpp::List aux_list,
                           int mc_mode,
                           int n_mc,
                           int seed) {
  // method: 1=joint_mc, 2=joint_is, 3=marginal_is
  // Special handling of h_break_points as
  // an RcppGSL::Vector object to be propergated to the classes
  RcppGSL::Vector h_break_points = Rcpp::as< RcppGSL::vector<double> >(aux_list["h_break_points"]);

  // Import parameters
  Unimodal_Model_Lite* pars = new Unimodal_Model_Lite(pars_list, aux_list, h_break_points);
  pars->generate_chol_centering_mat();

  // Import data
  std::vector<Unimodal_Curve_Lite>* data = new std::vector<Unimodal_Curve_Lite>;
  int id = 0;
  for(Rcpp::List::iterator it = curve_list.begin(); it != curve_list.end(); ++it) {
    data->push_back(Unimodal_Curve_Lite(*it, pars, mc_mode, n_mc, id, seed));
    // Initialize basis evaluation matrices
    data->at(id).initialize_h_basis_mat();
    ++id;
  }
  // Direct Monte Carlo
  if(mc_mode==1){
    for(int curve_idx = 0; curve_idx < data->size(); ++curve_idx){
      data->at(curve_idx).approx_logLik_full_mc();
    }
  }
  // Importance sampling
  if(mc_mode==2){
    for(int curve_idx = 0; curve_idx < data->size(); ++curve_idx){
      data->at(curve_idx).approx_logLik_full_is();
    }
  }
  // MCMC
  if(mc_mode==4){
    for(int curve_idx = 0; curve_idx < data->size(); ++curve_idx){
      data->at(curve_idx).approx_logLik_full_mcmc();
    }
  }

  // Wrap up as R objects
  Rcpp::List curves(pars->n_curve);
  Rcpp::List est_pars;
  Rcpp::List aux;
  for(int i = 0; i < pars->n_curve; ++i){
    curves[i] = data->at(i).return_obj();
  }
  est_pars = pars->return_pars();
  aux = pars->return_aux();

  // Free memory
  free(data);
  free(pars);

  return(Rcpp::List::create(
      Rcpp::Named("pars", est_pars),
      Rcpp::Named("curves", curves),
      Rcpp::Named("aux", aux)));

}


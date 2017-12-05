// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include "class_warping_function.h"
#include "class_warping_model.h"
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>

//' Internal function for fitting the model by SAEM
//'
//' @param w_list R List - list of compositional data points
//' @export
// [[Rcpp::export]]
Rcpp::List dirichlet_mle(Rcpp::List w_list) {

  int id;
  Rcpp::List estimates;

  // Import data
  std::vector<Warping_Function*>* data = new std::vector<Warping_Function*>;
  id = 0;
  for(Rcpp::List::iterator it = w_list.begin(); it != w_list.end(); ++it) {
    data->push_back(new Warping_Function(*it, id));
    ++id;
  }

  // Initialize Model
  Warping_Model* model = new Warping_Model(data->at(0)->dim_w - 1,
                                           data->size());

  // Update Sufficient Statistics
  Rcpp::Rcout << "Update Sufficient Statistics" << std::endl;
  model->Update_Suff_Stat(data);

  // Maximum Likelihood Estimation
  Rcpp::Rcout << "MLE by Newton Raphson" << std::endl;
  model->Find_MLE(100, 100, 1e-6); // args: max_outer_iter, max_inner_iter, tol_thres

  // Return Estimates
  estimates = model->Return_Estimates();

  // Deconstruct
  for(std::vector<Warping_Function*>::iterator it = data->begin(); it != data->end(); ++it){
    delete *it;
  }
  delete model;

  // Return Estimates to R
  return estimates;
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include "class_pars.h"
#include "class_curve.h"

void squish(arma::vec *x, double left_bound, double right_bound) {
  int n = x->size();

  for (int i = 0; i < n; ++i) {
    if (x->at(i) < left_bound) {
      x->at(i) = left_bound;
    } else if (x->at(i) > right_bound) {
      x->at(i) = right_bound;
    }
  }
  return;
}



double compute_llk_dw(arma::vec dw, arma::vec tau_kappa) {
  int dim_dw;
  bool any_degenerate(false);
  double log_d_terms;
  double log_w_terms = 0;
  double density;

  dim_dw = dw.size();
  // log reciprocal of generalized beta function
  log_d_terms = - lgamma(sum(tau_kappa));
  for(int i=0; i<dim_dw; i++){
    log_d_terms += lgamma(tau_kappa(i));
  };
  // check degenerate case
  for(int i=0; i<dim_dw; i++){
    if(dw[i]==0 && tau_kappa(i)==1){
      any_degenerate = true;
    }
  }
  // density kernel
  if(any_degenerate){
    log_w_terms = R_NegInf;
  } else{
    for(int i=0; i<dim_dw; i++){
      log_w_terms += (tau_kappa(i) - 1) * log(dw(i));
    }
  }
  density = log_w_terms - log_d_terms;
  return density;
}



// Constructor
Curve::Curve(Rcpp::List curve_obj, Pars* pars, int id) : curve_id(id){

  // Temporary variables
  arma::vec tmp_log_current_dw;
  arma::vec tmp_log_centered_dw;


  // Extract components from R list
  Rcpp::List data = curve_obj["data"];
  arma::vec current_aw = Rcpp::as<arma::vec>(curve_obj["re_mc_state"]);


  // Rcpp::Rcout << "Pointer to common_pars...";

  // Point to common pars
  common_pars = pars;

  // Rcpp::Rcout << "Raw data...";

  // Raw Data
  y = Rcpp::as<arma::vec>(data["y"]);
  x = Rcpp::as<arma::vec>(data["x"]);

  // Dimensions
  n_i = y.size();
  dim_a = 2;
  dim_w = current_aw.size() - 2;
  dim_z = dim_w - 2;
  dim_alpha = common_pars->alpha.size(); // number of basis coefficients for the base curve

  // Basis Evaluation Matrix - warping function
  h_basis_mat = Rcpp::as<arma::mat>(curve_obj["h_basis_mat"]);

  // Rcpp::Rcout << "MCMC...";

  // MCMC working variables
  current_a = current_aw.subvec(0, dim_a - 1);
  current_w = current_aw.subvec(dim_a, current_aw.size()-1);
  current_dw = arma::diff(current_w);
  // Transform current_w back to euclidean space
  tmp_log_current_dw = arma::log(current_dw);
  tmp_log_centered_dw = tmp_log_current_dw - mean(tmp_log_current_dw);
  current_z = common_pars->chol_centering_mat.t() * tmp_log_centered_dw;
  current_warped_x = x;  // n_i x 1
  if((current_warped_x.min() < common_pars->f_left_bound) ||
     (current_warped_x.max() > common_pars->f_right_bound)){
    // Rcpp::Rcout << "Tugging for curve " << curve_id << std::endl;
    squish(&current_warped_x, common_pars->f_left_bound, common_pars->f_right_bound);
  }

  current_warped_f_basis_mat = arma::zeros(n_i,dim_alpha);

  proposed_w = arma::zeros(dim_w); // k_h x 1
  proposed_dw = arma::zeros(dim_w); // k_h - 1 x 1
  proposed_z = arma::zeros(dim_w); //dim_w x 1
  proposed_warped_x = arma::zeros(n_i); // n_i x 1
  proposed_warped_f_basis_mat = arma::zeros(n_i,dim_alpha);

  // Rcpp::Rcout << "Suff. Stat...";

  // Stochastic approximated sufficient statistics
  sapprox_a = arma::zeros(dim_a); // dim_a x 1
  sapprox_w = arma::zeros(dim_w); // dim_w x 1
  sapprox_warped_f_basis_mat = arma::zeros(n_i, dim_alpha); // n_i x dim_alpha

  sapprox_aug_warped_f_basis_mat = arma::zeros(n_i, dim_alpha + 1); // n_i x dim_alpha + 1
  sapprox_hat_mat = arma::zeros(dim_alpha + 1, dim_alpha + 1);      // (dim_alpha + 1) x (dim_alpha + 1)
  sapprox_sigma_a = arma::zeros(dim_a, dim_a);                      // dim_a x dim_a
  sapprox_log_dw = arma::zeros(dim_w - 1);                          // (dim_w - 1) x 1

  // Sufficient statistics based on the current MCMC draw
  current_aug_warped_f_basis_mat = arma::zeros(n_i, dim_alpha + 1); // n_i x dim_alpha + 1
  current_hat_mat = arma::zeros(dim_alpha + 1, dim_alpha + 1);      // (dim_alpha + 1) x (dim_alpha + 1)
  current_sigma_a = arma::zeros(dim_a, dim_a);                      // dim_a x dim_a
  current_log_dw = arma::zeros(dim_w - 1);                          // (dim_w - 1) x 1
}


// Initialize the f_basis_mat under current warpings
// Depends on: ...
// Changes: ...
void Curve::initialize_current_f_basis_mat(){
  gsl_vector *tmp_b_vec;
  gsl_bspline_workspace *tmp_bw;
  // Rcpp::Rcout << "allocate a cubic bspline workspace (k = 4)" << std::endl;
  // allocate a cubic bspline workspace (k = 4)
  tmp_b_vec = gsl_vector_alloc(dim_alpha);
  tmp_bw = gsl_bspline_alloc(common_pars->f_order,
                             common_pars->f_break_points.size());
  // Rcpp::Rcout << "evaluate current_warped_f_basis_mat" << std::endl;
  // evaluate current_warped_f_basis_mat
  gsl_bspline_knots(common_pars->f_break_points, tmp_bw);      // computes the knots associated with the given breakpoints and
  // stores them internally in bw->knots.
  for(int i = 0; i < n_i; ++i){                                // construct the basis evaluation matrix, warped_f_basis_mat
    // Rcpp::Rcout << "evaluating current_warped_f_basis_mat" << std::endl;
    gsl_bspline_eval(proposed_warped_x[i], tmp_b_vec, tmp_bw); // compute B_j(x_i) for all j
    for(int j = 0; j < dim_alpha; ++j){                        // fill in row i of X
      // Rcpp::Rcout << "updating proposed_warped_f_basis_mat" << std::endl;
      current_warped_f_basis_mat(i,j) = gsl_vector_get(tmp_b_vec, j);  // gsl_vector_get(B, j)
    }
  }
  // free GSL workspace
  gsl_bspline_free(tmp_bw);
  gsl_vector_free(tmp_b_vec);
}



// Draw a new proposed_w
// Depends on: common_pars, current_z, dim_w
// Changes: proposed_z, proposed_dw, proposed_w
void Curve::propose_new_w(){
  arma::vec tmp_dw(dim_w-1);
  // Update the random walk in R^dim_z from current_z
  proposed_z = current_z +
    arma::chol(common_pars->identity_cor_mat * common_pars->proposal_sigma).t() * arma::randn(dim_z);
  // update proposed_w
  tmp_dw = common_pars->chol_centering_mat * proposed_z;
  tmp_dw -= mean(tmp_dw);
  tmp_dw.transform(exp);
  proposed_dw = tmp_dw / arma::sum(tmp_dw);
  proposed_w(0) = 0;
  for(int idx = 0; idx < proposed_dw.size(); ++idx){
    proposed_w(idx + 1) = proposed_w(idx) + proposed_dw(idx);
  }
  return;
}



// Compute proposed_warped_f and  proposed_warped_y
// Depends on: h_basis_mat, proposed_w, common_pars
// Changes: proposed_warped_x, proposed_warped_f_basis_mat
void Curve::compute_proposed_warping_and_f_basis_mat(){
  // Squish in place to avoid out-of-bound error in gsl_bspline_eval
  proposed_warped_x = h_basis_mat * proposed_w;
  if((proposed_warped_x.min() < common_pars->f_left_bound) ||
     (proposed_warped_x.max() > common_pars->f_right_bound)){
    // Rcpp::Rcout << "Tugging for curve " << curve_id << std::endl;
    squish(&proposed_warped_x, common_pars->f_left_bound, common_pars->f_right_bound);
  }

  gsl_vector *tmp_b_vec;
  gsl_bspline_workspace *tmp_bw;

  // Rcpp::Rcout << "allocate a cubic bspline workspace (k = 4)" << std::endl;

  // allocate a cubic bspline workspace (k = 4)
  tmp_b_vec = gsl_vector_alloc(dim_alpha);
  tmp_bw = gsl_bspline_alloc(common_pars->f_order,
                         common_pars->f_break_points.size());

  // Rcpp::Rcout << "evaluate proposed_warped_f_basis_mat" << std::endl;

  // evaluate proposed_warped_f_basis_mat
  gsl_bspline_knots(common_pars->f_break_points, tmp_bw);      // computes the knots associated with the given breakpoints and
                                                               // stores them internally in bw->knots.

  for(int i = 0; i < n_i; ++i){           // construct the basis evaluation matrix, warped_f_basis_mat
    gsl_bspline_eval(proposed_warped_x[i], tmp_b_vec, tmp_bw); // compute B_j(x_i) for all j
    for(int j = 0; j < dim_alpha; ++j){                        // fill in row i of X
      // Rcpp::Rcout << "update proposed_warped_f_basis_mat" << std::endl;
      proposed_warped_f_basis_mat(i,j) = gsl_vector_get(tmp_b_vec, j);
    }
  }

  // free GSL workspace
  gsl_bspline_free(tmp_bw);
  gsl_vector_free(tmp_b_vec);
}



// Compute the Metropolis-Hasting ratio
// Depends on: current_a, common_pars,
//             proposed_warped_f_basis_mat, current_warped_f_basis_mat
// Changes: Nil
// Return: log metropolis-hasting ratio
double Curve::compute_log_mh_ratio(){
  double current_llk_data;
  double proposed_llk_data;
  double current_llk_w;
  double proposed_llk_w;
  double log_jacobian_term;

  arma::vec proposed_warped_f = proposed_warped_f_basis_mat * common_pars->alpha;
  arma::vec current_warped_f = current_warped_f_basis_mat * common_pars->alpha;
  arma::vec proposed_warped_y = current_a(0) + current_a(1) * proposed_warped_f;
  arma::vec current_warped_y = current_a(0) + current_a(1) * current_warped_f;

  // Compute the data log-likelihood (up to the common constant term)
  current_llk_data = -sum(square(y - current_warped_y)) / 2 / common_pars->sigma2;
  proposed_llk_data = -sum(square(y - proposed_warped_y)) / 2 / common_pars->sigma2;

  // Compute the dirichlet log-likelihood
  current_llk_w = compute_llk_dw(current_dw, common_pars->tau * common_pars->kappa);
  proposed_llk_w = compute_llk_dw(proposed_dw, common_pars->tau * common_pars->kappa);

  // Log jacobian term for the MH ratio
  log_jacobian_term = arma::sum(arma::log(proposed_dw) - arma::log(current_dw));

  return proposed_llk_data - current_llk_data + proposed_llk_w - current_llk_w + log_jacobian_term;
}



// Calls compute_log_mh_ratio() and accept/reject the propsed warping
// Depends on: see compute_log_mh_ratio()
// Change: current_z, current_dw, current_w, current_warped_x, current_warped_f_basis_mat, common_pars
// Note: Acceptances are tallied in the table of common_pars->mh_accept_rate_table
void Curve::mh_accept_reject(){
  double u = sum(arma::randu(1));
  if (std::log(u) < compute_log_mh_ratio()) {
    current_z = proposed_z;
    current_dw = proposed_dw;
    current_w = proposed_w;
    current_warped_x = proposed_warped_x;
    current_warped_f_basis_mat = proposed_warped_f_basis_mat;
    ++(common_pars->mh_accept_rate_table(curve_id, common_pars->mh_accept_rate_table_counter));
    return;
  }
}



// Draw a new amplitude effect (a) from the Gibbs sampler
// Depends on: current_warped_f_basis_mat, common_pars
// Changes: current_a
// Notes: Store updated current_a in common_pars->current_a_mat.col(curve_id);
void Curve::draw_new_a(){
  arma::mat tmp_f_mat = arma::ones(n_i, 2);
  arma::vec tmp_mu_post;
  arma::mat tmp_sigma_post;

  tmp_f_mat.col(1) = current_warped_f_basis_mat * common_pars->alpha;
  tmp_sigma_post = inv(tmp_f_mat.t() * tmp_f_mat / common_pars->sigma2 +
    common_pars->big_sigma_inverse);
  tmp_mu_post = tmp_sigma_post * (tmp_f_mat.t() * y / common_pars->sigma2 +
    common_pars->big_sigma_inverse * common_pars->mu);
  current_a = tmp_mu_post + arma::chol(tmp_sigma_post).t() * arma::randn(dim_a);
  common_pars->current_a_mat.col(curve_id) = current_a;
  // if(curve_id < 50){
  //   Rcpp::Rcout << "curve_id: " << curve_id << std::endl;
  //   Rcpp::Rcout << "tmp_mu_post: " << tmp_mu_post.t() << std::endl;
  //   Rcpp::Rcout << "tmp_sigma_post: " << std::endl << tmp_sigma_post << std::endl;
  //   Rcpp::Rcout << "sampled_a: " << std::endl << current_a.t() << std::endl;
  // }
}



// Wraper function to run the simulation step
void Curve::do_simulation_step(){
  for(int i = 0; i < common_pars->n_burn_mcmc; ++i){
    propose_new_w();
    compute_log_mh_ratio();
    compute_proposed_warping_and_f_basis_mat();
    mh_accept_reject();
    draw_new_a();
  }
}



// Centering step for current_a
void Curve::center_current_a(){
  if(!common_pars->need_centering)
    return;
  arma::mat centering_mat(dim_a, dim_a);
  arma::vec mean_current_a(dim_a);
  if(dim_a == 2){
    mean_current_a = mean(common_pars->current_a_mat, 1);
    centering_mat.eye();
    centering_mat(0, 1) = -mean_current_a(0) / mean_current_a(1);
    centering_mat(1, 1) = 1 / mean_current_a(1);
    current_a = centering_mat * current_a;
  }
  else {
    Rcpp::Rcout << "Warning! Centering only supporting for dim_a = 2.";
  }
}



// Update stochastic approximation of the sufficients statistics
// Depends on: ...
// Changes: ...
void Curve::update_sufficient_statistics_approximates(){
  double current_step_size = common_pars->saem_step_sizes(common_pars->saem_counter);
  arma::mat tmp_half_hat_mat(n_i, dim_alpha+1);

  // Compute sufficient statistics based on current MC state
  current_aug_warped_f_basis_mat.col(0) = current_a(0) * arma::ones(n_i);
  current_aug_warped_f_basis_mat.cols(1, dim_alpha) = current_a(1) * current_warped_f_basis_mat;
  tmp_half_hat_mat = current_aug_warped_f_basis_mat;
  tmp_half_hat_mat.col(0) = y - current_a(0);
  current_hat_mat = tmp_half_hat_mat.t() * tmp_half_hat_mat;
  current_sigma_a = (current_a - common_pars->mu) * (current_a - common_pars->mu).t();
  current_log_dw = arma::log(current_dw);

  // Update stochastic approximates
  sapprox_a = (1 - current_step_size) * sapprox_a +
    current_step_size * current_a;
  sapprox_w = (1 - current_step_size) * sapprox_w +
    current_step_size * current_w;
  sapprox_warped_f_basis_mat = (1 - current_step_size) * sapprox_warped_f_basis_mat +
    current_step_size * current_warped_f_basis_mat;
  sapprox_aug_warped_f_basis_mat = (1 - current_step_size) * sapprox_aug_warped_f_basis_mat +
    current_step_size * current_aug_warped_f_basis_mat;
  sapprox_hat_mat = (1 - current_step_size) * sapprox_hat_mat +
    current_step_size * current_hat_mat;
  sapprox_sigma_a = (1 - current_step_size) * sapprox_sigma_a +
    current_step_size * current_sigma_a;
  sapprox_log_dw = (1 - current_step_size) * sapprox_log_dw +
    current_step_size * current_log_dw;
}




Rcpp::List Curve::return_list(){
  arma::vec tmp_alpha_aug(dim_alpha + 1, arma::fill::ones);
  tmp_alpha_aug(arma::span(1, dim_alpha)) = common_pars->alpha;
  return Rcpp::List::create(
    Rcpp::Named("curve_id", Rcpp::wrap(curve_id)),
    Rcpp::Named("x", Rcpp::wrap(x)),
    Rcpp::Named("y", Rcpp::wrap(y)),
    Rcpp::Named("warped_x", Rcpp::wrap(h_basis_mat * sapprox_w)),
    Rcpp::Named("fitted_y", Rcpp::wrap(sapprox_aug_warped_f_basis_mat * tmp_alpha_aug)),
    Rcpp::Named("sapprox_a", Rcpp::wrap(sapprox_a)),
    Rcpp::Named("sapprox_w", Rcpp::wrap(sapprox_w)),
    Rcpp::Named("sapprox_warped_f_basis_mat", Rcpp::wrap(sapprox_warped_f_basis_mat)),
    Rcpp::Named("sapprox_aug_warped_f_basis_mat", Rcpp::wrap(sapprox_aug_warped_f_basis_mat)),
    Rcpp::Named("sapprox_hat_mat", Rcpp::wrap(sapprox_hat_mat)),
    Rcpp::Named("sapprox_sigma_a", Rcpp::wrap(sapprox_sigma_a)),
    Rcpp::Named("sapprox_log_dw", Rcpp::wrap(sapprox_log_dw))
  );
};

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include "class_unimodal_model.h"
#include "class_unimodal_curve.h"
#include "util.h"

// Constructor
Unimodal_Curve::Unimodal_Curve(Rcpp::List data, Unimodal_Model* pars, int id, int seed) : curve_id(id){

  // Temporary variables
  arma::vec tmp_log_current_dw;
  arma::vec tmp_log_centered_dw;

  // Point to common pars
  common_pars = pars;

  // Raw Data
  y = Rcpp::as<arma::vec>(data["y"]);
  x = Rcpp::as<arma::vec>(data["x"]);
  init_clust = Rcpp::as<arma::ivec>(data["init_clust"]).at(0) - 1;

  // Dimensions
  n_i = y.size();
  dim_a = common_pars->dim_a;
  dim_w = common_pars->dim_w;
  dim_z = dim_w - 2;
  n_cluster = common_pars->n_cluster; // number of clusters

  // Basis Evaluation Matrix - warping function
  h_basis_mat = arma::zeros(n_i, dim_w);


  // MCMC working variables
  // current_m = rand() % n_cluster;           // randomize initial configuration
  current_m = 0;                            // put every curve in group 0 at the beginning
  current_a = common_pars->mu_a;            // user-specified initial value -- expecting OLS
  current_dw = common_pars->kappa_id;       // initial with no warping
  current_w = arma::zeros(dim_w);
  current_w(arma::span(1,dim_w-1)) = arma::cumsum(current_dw);
  // Transform current_w back to euclidean space
  tmp_log_current_dw = arma::log(current_dw);
  tmp_log_centered_dw = tmp_log_current_dw - mean(tmp_log_current_dw);
  current_z = common_pars->chol_centering_mat.t() * tmp_log_centered_dw;
  current_warped_x = x;
  if((current_warped_x.min() < common_pars->h_left_bound) ||
     (current_warped_x.max() > common_pars->h_right_bound)){
    current_warped_x = arma::clamp(current_warped_x, common_pars->h_left_bound, common_pars->h_right_bound);
  }

  proposed_w = arma::zeros(dim_w);      // dim_w x 1
  proposed_dw = arma::zeros(dim_w - 1); // (dim_w - 1) x 1
  proposed_z = arma::zeros(dim_w - 2);  // (dim_w - 2) x 1
  proposed_warped_x = arma::zeros(n_i); // n_i x 1

  // Stochastic approximated sufficient statistics
  sapprox_a = arma::zeros(dim_a);                   // dim_a x 1
  sapprox_sq_a = arma::zeros(dim_a);                // dim_a x 1
  sapprox_w = arma::zeros(dim_w);                   // dim_w x 1
  sapprox_warped_f = arma::zeros(n_i);              // n_i x 1
  sapprox_fitted_y = arma::zeros(n_i);              // n_i x 1

  sapprox_log_dw = arma::zeros(dim_w - 1, n_cluster);            // (dim_w - 1) x 1
  sapprox_cluster_membership = arma::zeros(n_cluster);           // n_cluster x 1
  sapprox_residual_sum_of_squares = 0.0;

  // Sufficient statistics based on the current MCMC draw
  current_log_dw = arma::zeros(dim_w - 1, n_cluster);               // (dim_w - 1) x n_cluster
  current_cluster_membership = arma::zeros(n_cluster);              // n_cluster x 1
  current_cluster_membership(current_m) = 1;

  // Random number generator
  rng_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_gen, seed + curve_id);

  // Temporary or working variables
  // ... for draw_new_m()
  tmp_u = 0.0;
  tmp_pred_prob_clusters = arma::zeros(n_cluster);
  // ... for draw_new_a()
  tmp_f_mat = arma::ones(n_i, 2);
  tmp_mu_post = arma::zeros(dim_a);
  tmp_sigma_post = arma::zeros(dim_a, dim_a);
  //... for propose_new_w()
  tmp_dw = arma::zeros(dim_w - 1);
  // ... for compute_log_mh_ratio()
  proposed_warped_f = arma::zeros(n_i);
  proposed_fitted_y = arma::zeros(n_i);
  current_warped_f = arma::zeros(n_i);
  current_fitted_y = arma::zeros(n_i);
  proposed_residual_sum_of_squares = 0.0;
  current_residual_sum_of_squares = 0.0;
  proposed_minus_current_llk_data = 0.0;
  current_llk_w = 0.0;
  proposed_llk_w = 0.0;
  log_jacobian_term = 0.0;
  // ... for mh_accept_reject()
  mh_randu = 0.0;
  // ... for update_sufficient_statistics_approximates()
  current_step_size = 0.0;
  // ...for mh_accept_reject()
  tmp_mh_log_accept_prob = 0.0;

  // Feed back to pars
  common_pars->current_m_vec(curve_id) = current_m;
}



// Initialize the warping function basis evaluation matrix
// Depends on: x, common_pars
// Changes: h_basis_mat
void Unimodal_Curve::initialize_h_basis_mat(){
  if((x.min() < common_pars->h_left_bound) ||
     (x.max() > common_pars->h_right_bound)){
    x = arma::clamp(x, common_pars->h_left_bound, common_pars->h_right_bound);
  }

  gsl_vector *tmp_a_vec;
  gsl_bspline_workspace *tmp_aw;

  // allocate a cubic bspline workspace (k = 4)
  tmp_a_vec = gsl_vector_alloc(dim_w);
  tmp_aw = gsl_bspline_alloc(common_pars->h_order, common_pars->h_break_points.size());

  // evaluate current_warped_f_basis_mat
  gsl_bspline_knots(common_pars->h_break_points, tmp_aw);      // computes the knots associated with the given breakpoints and
  // stores them internally in tmp_aw->knots.
  for(int i = 0; i < n_i; ++i){                                // construct the basis evaluation matrix, warped_f_basis_mat
    gsl_bspline_eval(x(i), tmp_a_vec, tmp_aw);                 // compute B_j(x_i) for all j
    for(int j = 0; j < dim_w; ++j){                            // fill in row i of X
      h_basis_mat(i,j) = gsl_vector_get(tmp_a_vec, j);         // gsl_vector_get(B, j)
    }
  }
  // free GSL workspace
  gsl_bspline_free(tmp_aw);
  gsl_vector_free(tmp_a_vec);
  return;
}



// Draw a new proposed_w
// Depends on: common_pars, current_z, dim_w
// Changes: proposed_z, proposed_dw, proposed_w
void Unimodal_Curve::propose_new_w(){
  // Update the random walk in R^dim_z from current_z
  proposed_z = current_z +
    arma::chol(common_pars->identity_cor_mat * common_pars->proposal_sigma).t() * arma::randn(dim_z);
  // update proposed_w
  tmp_dw = common_pars->chol_centering_mat * proposed_z;
  tmp_dw -= mean(tmp_dw);
  tmp_dw.transform(exp);
  proposed_dw = tmp_dw / arma::sum(tmp_dw);
  proposed_w(0) = 0;
  proposed_w(arma::span(1, dim_w - 1)) = arma::cumsum(proposed_dw);
  return;
}



// Compute compute_proposed_warping_and_f_basis_mat
// Depends on: h_basis_mat, proposed_w, common_pars
// Changes: proposed_warped_x, proposed_warped_f_basis_mat
void Unimodal_Curve::compute_proposed_warping(){
  proposed_warped_x = h_basis_mat * proposed_w;
  // Squish in place to avoid out-of-bound error in gsl_bspline_eval (LEGACY)
  if((proposed_warped_x.min() < common_pars->h_left_bound) ||
     (proposed_warped_x.max() > common_pars->h_right_bound)){
    proposed_warped_x = arma::clamp(proposed_warped_x,
                                    common_pars->h_left_bound,
                                    common_pars->h_right_bound);
  }
  return;
}



// Compute the Metropolis-Hasting ratio
// Depends on: current_a, common_pars,
//             proposed_warped_f_basis_mat, current_warped_f_basis_mat
// Changes: Nil
// Return: log metropolis-hasting ratio
void Unimodal_Curve::compute_log_mh_ratio(){
  proposed_minus_current_llk_data = 0.0;
  current_llk_w = 0.0;
  proposed_llk_w = 0.0;
  log_jacobian_term = 0.0;

  proposed_warped_f = 4 * (arma::square(proposed_warped_x) - proposed_warped_x);
  proposed_fitted_y = current_a(0) + current_a(1) * proposed_warped_f;
  current_warped_f = 4 * (arma::square(current_warped_x) - current_warped_x);
  current_fitted_y = current_a(0) + current_a(1) * current_warped_f;

  proposed_residual_sum_of_squares = arma::sum(arma::square(y - proposed_fitted_y));
  current_residual_sum_of_squares = arma::sum(arma::square(y - current_fitted_y));

  // Compute the data log-likelihood (up to the common constant term)
  proposed_minus_current_llk_data = (current_residual_sum_of_squares -
    proposed_residual_sum_of_squares) / 2.0 / common_pars->sigma2;

  // Compute the dirichlet log-likelihood
  proposed_llk_w = compute_llk_dw(proposed_dw, common_pars->kappa_clusters.col(current_m));
  current_llk_w = compute_llk_dw(current_dw, common_pars->kappa_clusters.col(current_m));

  // Compute the log jacobian term for the MH ratio
  log_jacobian_term = arma::sum(arma::log(proposed_dw) - arma::log(current_dw));

  tmp_mh_log_accept_prob = proposed_minus_current_llk_data + proposed_llk_w - current_llk_w + log_jacobian_term;
}



// Calls compute_log_mh_ratio() and accept/reject the propsed warping
// Depends on: see compute_log_mh_ratio()
// Change: current_z, current_dw, current_w, current_warped_x, current_warped_f_basis_mat, common_pars
// Note: Acceptances are tallied in the table of common_pars->mh_accept_rate_table
void Unimodal_Curve::mh_accept_reject(){
  mh_randu = gsl_rng_uniform(rng_gen);
  if (std::log(mh_randu) < tmp_mh_log_accept_prob) {
    current_z = proposed_z;
    current_dw = proposed_dw;
    current_w = proposed_w;
    current_warped_x = proposed_warped_x;
    ++(common_pars->mh_accept_rate_table(curve_id, common_pars->mh_accept_rate_table_counter));
    return;
  }
  return;
}



// Draw a new amplitude effect (a) from the Gibbs sampler
// Depends on: current_warped_f_basis_mat, common_pars
// Changes: current_a
// Notes: Store updated current_a in common_pars->current_a_mat.col(curve_id);
void Unimodal_Curve::draw_new_a(){
  tmp_f_mat.col(1) = current_warped_f;
  tmp_sigma_post = inv(tmp_f_mat.t() * tmp_f_mat / common_pars->sigma2 +
    common_pars->sigma2_a_inv);
  tmp_mu_post = tmp_sigma_post * (tmp_f_mat.t() * y / common_pars->sigma2 +
    common_pars->sigma2_a_inv * common_pars->mu_a);
  current_a = tmp_mu_post + arma::chol(tmp_sigma_post).t() * arma::randn(dim_a);
  return;
}



// Draw a new cluster membership (m) from the Gibbs sampler
// Depends on: current_dw
// Changes: current_m;
void Unimodal_Curve::draw_new_m(){
  for(int cluster_idx = 0; cluster_idx < n_cluster; ++cluster_idx){
    tmp_pred_prob_clusters(cluster_idx) = compute_llk_dw(current_dw,
                           common_pars->kappa_clusters.col(cluster_idx));
  }
  tmp_pred_prob_clusters = exp(tmp_pred_prob_clusters) % common_pars->p_clusters;
  tmp_pred_prob_clusters = tmp_pred_prob_clusters / sum(tmp_pred_prob_clusters);

  tmp_u = gsl_rng_uniform(rng_gen);
  for(int cluster_idx = 0; cluster_idx < n_cluster; ++cluster_idx){
    if(tmp_u < tmp_pred_prob_clusters(cluster_idx)){
      current_m = cluster_idx;
      break;
    }
    else{
      tmp_u -= tmp_pred_prob_clusters(cluster_idx);
    }
  }
  current_cluster_membership.zeros();
  current_cluster_membership(current_m) = 1;
  common_pars->current_m_vec(curve_id) = current_m;
  return;
}



// Wraper function to run the simulation step
void Unimodal_Curve::do_simulation_step(){
  for(int mcmc_iter = 0; mcmc_iter < common_pars->n_burn_mcmc; ++mcmc_iter){
    propose_new_w();
    compute_proposed_warping();
    compute_log_mh_ratio();
    mh_accept_reject();
    draw_new_a();
    if(common_pars->saem_counter > common_pars->n_burn_saem){
      draw_new_m();
    }
  }
  return;
}



// Update stochastic approximation of the sufficients statistics
// Depends on: n_i, dim_alpha, current_a, current_dw, current_w,
//             current_warped_f_basis_mat, y, common_pars,
// Changes: current_hat_mat, current_sigma_a, current_log_dw,
//          sapprox_a, sapprox_w, sapprox_warped_f_basis_mat,
//          sapprox_aug_warped_f_basis_mat, sapprox_hat_mat,
//          sapprox_sigma_a, sapprox_log_dw
void Unimodal_Curve::update_sufficient_statistics_approximates(){
  current_step_size = common_pars->saem_step_sizes(common_pars->saem_counter);

  // Compute sufficient statistics based on current MC state
  current_log_dw.zeros();
  current_log_dw.col(current_m) = arma::log(current_dw);

  // Update stochastic approximates
  sapprox_residual_sum_of_squares = (1 - current_step_size) * sapprox_residual_sum_of_squares +
    current_step_size * current_residual_sum_of_squares;
  sapprox_a = (1 - current_step_size) * sapprox_a +
    current_step_size * current_a;
  sapprox_sq_a = (1 - current_step_size) * sapprox_sq_a +
    current_step_size * arma::square(current_a);
  sapprox_cluster_membership = (1 - current_step_size) * sapprox_cluster_membership +
    current_step_size * current_cluster_membership;
  sapprox_w = (1 - current_step_size) * sapprox_w +
    current_step_size * current_w;
  sapprox_log_dw = (1 - current_step_size) * sapprox_log_dw +
    current_step_size * current_log_dw;
  sapprox_warped_f = (1 - current_step_size) * sapprox_warped_f +
    current_step_size * current_warped_f;
  sapprox_fitted_y = (1 - current_step_size) * sapprox_fitted_y +
    current_step_size * current_fitted_y;
  return;
}



// Return fitted curve, predicted warping functions and sufficient statistics
Rcpp::List Unimodal_Curve::return_list(){
  return Rcpp::List::create(
    Rcpp::Named("curve_id", Rcpp::wrap(curve_id)),
    Rcpp::Named("x", Rcpp::wrap(x)),
    Rcpp::Named("y", Rcpp::wrap(y)),
    Rcpp::Named("init_clust", Rcpp::wrap(init_clust)),
    Rcpp::Named("warped_x", Rcpp::wrap(h_basis_mat * sapprox_w)),
    Rcpp::Named("fitted_y", Rcpp::wrap(sapprox_fitted_y)),
    Rcpp::Named("sapprox_residual_sum_of_squares", Rcpp::wrap(sapprox_residual_sum_of_squares)),
    Rcpp::Named("sapprox_a", Rcpp::wrap(sapprox_a)),
    Rcpp::Named("sapprox_w", Rcpp::wrap(sapprox_w)),
    Rcpp::Named("sapprox_log_dw", Rcpp::wrap(sapprox_log_dw)),
    Rcpp::Named("sapprox_cluster_membership", Rcpp::wrap(sapprox_cluster_membership))
  );
};



// Return fitted curve, predicted warping functions and sufficient statistics
Rcpp::List Unimodal_Curve::return_list(double y_scaling_factor){
  // Scaling matrices
  return Rcpp::List::create(
    Rcpp::Named("curve_id", Rcpp::wrap(curve_id)),
    Rcpp::Named("x", Rcpp::wrap(x)),
    Rcpp::Named("y", Rcpp::wrap(y * y_scaling_factor)),
    Rcpp::Named("init_clust", Rcpp::wrap(init_clust)),
    Rcpp::Named("warped_x", Rcpp::wrap(h_basis_mat * sapprox_w)),
    Rcpp::Named("fitted_y", Rcpp::wrap(y_scaling_factor * sapprox_fitted_y)),
    Rcpp::Named("sapprox_residual_sum_of_squares", Rcpp::wrap(y_scaling_factor * y_scaling_factor * sapprox_residual_sum_of_squares)),
    Rcpp::Named("sapprox_a", Rcpp::wrap(y_scaling_factor * sapprox_a)),
    Rcpp::Named("sapprox_w", Rcpp::wrap(sapprox_w)),
    Rcpp::Named("sapprox_log_dw", Rcpp::wrap(sapprox_log_dw)),
    Rcpp::Named("sapprox_cluster_membership", Rcpp::wrap(sapprox_cluster_membership))
  );
};



void Unimodal_Curve::print_random_number(){
  double u = gsl_rng_uniform(rng_gen);
  Rcpp::Rcout << "Curve: " << curve_id << " random number: " << u << std::endl;
}

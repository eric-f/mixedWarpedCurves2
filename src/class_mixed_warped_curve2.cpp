// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include "class_mixed_warped_model2.h"
#include "class_mixed_warped_curve2.h"
#include "util.h"

// Constructor
Mixed_Warped_Curve::Mixed_Warped_Curve(Rcpp::List data, Mixed_Warped_Model* pars, int id, int seed) : curve_id(id){

  // Point to common pars
  common_pars = pars;

  // Raw Data
  y = Rcpp::as<arma::vec>(data["y"]);
  x = Rcpp::as<arma::vec>(data["x"]);
  init_clust = Rcpp::as<arma::ivec>(data["init_clust"]).at(0) - 1;

  // Dimensions
  n_i = y.size();
  dim_f = common_pars->dim_f;
  dim_a = common_pars->dim_a;
  dim_w = common_pars->dim_w;
  dim_z = dim_w - 2;
  n_cluster = common_pars->n_cluster; // number of clusters

  // Basis Evaluation Matrix - warping function
  h_basis_mat = arma::zeros(n_i, dim_w);

  // MCMC working variables
  current_m = 0;                           // put every curve in group 0 at the beginning
  current_a = common_pars->mu_a;           // Initialize in R. Depending on model, OLS or fixed at (0, 1)
  current_dw = common_pars->kappa_id;      // Initial with no warping
  current_w = arma::zeros(dim_w);
  current_w(arma::span(1, dim_w-1)) = arma::cumsum(current_dw);

  // Linking current_w to random walk in R^dim_z
  tmp_log_current_dw = arma::log(current_dw);
  tmp_log_centered_dw = tmp_log_current_dw - mean(tmp_log_current_dw);
  current_z = common_pars->chol_centering_mat.t() * tmp_log_centered_dw;

  // Fitted object based on current warping
  current_warped_x = x;
  if((current_warped_x.min() < common_pars->h_left_bound) ||
     (current_warped_x.max() > common_pars->h_right_bound)){
    current_warped_x = arma::clamp(current_warped_x, common_pars->h_left_bound, common_pars->h_right_bound);
  }
  current_warped_f_basis_mat = arma::zeros(n_i, dim_f);

  // For MH proposal
  proposed_dw = arma::zeros(dim_w - 1);                 // (dim_w - 1) x 1
  proposed_w = arma::zeros(dim_w);                      // dim_w x 1
  proposed_z = arma::zeros(dim_w - 2);                  // (dim_w - 2) x 1
  proposed_warped_x = arma::zeros(n_i);                 // n_i x 1
  proposed_warped_f_basis_mat = arma::zeros(n_i, dim_f);// n_i x dim_f

  // Sufficient statistics based on the current MCMC draw
  current_log_dw = arma::zeros(dim_w - 1, n_cluster);           // (dim_w - 1) x n_cluster
  current_cluster_membership = arma::zeros(n_cluster);          // n_cluster x 1
  current_cluster_membership(current_m) = 1;                    // Initial to primary group
  current_aug_warped_f_basis_mat = arma::zeros(n_i, dim_f + 1); // n_i x (dim_f + 1)
  current_hat_mat = arma::zeros(dim_f + 1, dim_f + 1);          // (dim_f + 1) x (dim_f + 1)
  current_cov_a = arma::zeros(dim_a, dim_a);                    // dim_a x dim_a

  // Stochastic approximated sufficient statistics
  sapprox_a = arma::zeros(dim_a);                               // dim_a x 1
  sapprox_w = arma::zeros(dim_w);                               // dim_w x 1
  sapprox_log_dw = arma::zeros(dim_w - 1, n_cluster);           // (dim_w - 1) x n_cluster
  sapprox_cluster_membership = arma::zeros(n_cluster);          // n_cluster x 1
  sapprox_aug_warped_f_basis_mat = arma::zeros(n_i, dim_f + 1); // n_i x (dim_f + 1)          [Psi]
  sapprox_hat_mat = arma::zeros(dim_f + 1, dim_f + 1);          // (dim_f + 1) x (dim_f + 1)  [C_mat]
  sapprox_cov_a = arma::zeros(dim_a, dim_a);                    // dim_a x dim_a              [Sigma_a]

  // Private members

  // Random number generator
  rng_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_gen, seed + curve_id);

  // Temporary or working variables
  // ... for initialize_h_basis_mat()
  tmp_b_vec = gsl_vector_alloc(dim_f);   // Allocate a cubic bspline workspace
  tmp_bw = gsl_bspline_alloc(common_pars->f_order, common_pars->f_break_points.size());
  gsl_bspline_knots(common_pars->f_break_points, tmp_bw);   // Computes the knots associated with
                                                            // thegiven breakpoints and stores
                                                            // them internally in tmp_bw->knots.
  // ... for draw_new_m()
  tmp_u = 0.0;
  pred_prob_clusters = arma::zeros(n_cluster);
  cluster_idx = 0;
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
  mh_log_accept_prob = 0.0;
  // ... for mh_accept_reject()
  mh_randu = 0.0;
  // ... for center_current_a()
  centering_mat = arma::eye(dim_a, dim_a);
  mean_current_a = arma::zeros(dim_a);
  // ... for update_sufficient_statistics_approximates()
  current_step_size = 0.0;
  tmp_half_hat_mat = arma::zeros(n_i, dim_f + 1);
  // ... for do_simulation_step()
  mcmc_iter = 0;

  // Feed back to pars
  common_pars->current_m_vec(curve_id) = current_m;
}



// Initialize the warping function basis evaluation matrix
// Depends on: x, common_pars
// Changes: h_basis_mat
void Mixed_Warped_Curve::initialize_h_basis_mat(){
  if((x.min() < common_pars->h_left_bound) ||
     (x.max() > common_pars->h_right_bound)){
    x = arma::clamp(x, common_pars->h_left_bound, common_pars->h_right_bound);
  }

  gsl_vector *tmp_a_vec;
  gsl_bspline_workspace *tmp_aw;

  // allocate a cubic bspline workspace (k = 4)
  tmp_a_vec = gsl_vector_alloc(dim_w);
  tmp_aw = gsl_bspline_alloc(common_pars->h_order, common_pars->h_break_points.size());

  // evaluate h_basis_mat
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


// Initialize the f_basis_mat under current warpings
//
// Depends on: dim_f, current_warped_x
//
// Changes: current_warped_f_basis_mat
void Mixed_Warped_Curve::initialize_current_f_basis_mat(){
  // evaluate current_warped_f_basis_mat
  for(int i = 0; i < n_i; ++i){                               // construct the basis evaluation matrix, warped_f_basis_mat
    gsl_bspline_eval(current_warped_x[i], tmp_b_vec, tmp_bw); // compute B_j(x_i) for all j
    for(int j = 0; j < dim_f; ++j){                           // fill in row i of X
      current_warped_f_basis_mat(i,j) = gsl_vector_get(tmp_b_vec, j);
    }
  }
  return;
}


// Draw a new proposed_w
//
// Depends on: common_pars, current_z, dim_w
//
// Changes: proposed_z, proposed_dw, proposed_w
void Mixed_Warped_Curve::propose_new_w(){
  // Update the random walk in R^dim_z from current_z
  proposed_z = current_z +
    arma::chol(common_pars->identity_cor_mat * common_pars->proposal_sigma).t() * arma::randn(dim_z);
  // Update proposed_w
  tmp_dw = common_pars->chol_centering_mat * proposed_z;
  tmp_dw -= mean(tmp_dw);
  tmp_dw.transform(exp);
  proposed_dw = tmp_dw / arma::sum(tmp_dw);
  proposed_w(arma::span(1, dim_w - 1)) = arma::cumsum(proposed_dw);
  return;
}


// Compute compute_proposed_warping_and_f_basis_mat
//
// Depends on: h_basis_mat, proposed_w, common_pars
//
// Changes: proposed_warped_x, proposed_warped_f_basis_mat
void Mixed_Warped_Curve::compute_proposed_warping_and_f_basis_mat(){
  proposed_warped_x = h_basis_mat * proposed_w;
  // Squish in place to avoid out-of-bound error in gsl_bspline_eval (LEGACY)
  if((proposed_warped_x.min() < common_pars->h_left_bound) ||
     (proposed_warped_x.max() > common_pars->h_right_bound)){
    proposed_warped_x = arma::clamp(proposed_warped_x,
                                    common_pars->h_left_bound,
                                    common_pars->h_right_bound);
  }
  // evaluate proposed_warped_f_basis_mat
  for(int i = 0; i < n_i; ++i){                                 // construct the basis evaluation matrix, warped_f_basis_mat
    gsl_bspline_eval(proposed_warped_x(i), tmp_b_vec, tmp_bw);  // compute B_j(x_i) for all j
    for(int j = 0; j < dim_f; ++j){                             // fill in row i of X
      proposed_warped_f_basis_mat(i,j) = gsl_vector_get(tmp_b_vec, j);
    }
  }
  return;
}


// Compute the Metropolis-Hasting ratio
//
// Depends on: current_a, common_pars,
//             proposed_warped_f_basis_mat, current_warped_f_basis_mat,
//
// Changes: proposed_fitted_y, current_fitted_y,
//          proposed_residual_sum_of_squares,
//          current_residual_sum_of_squares,
//          proposed_llk_w, current_llk_w,
//          mh_log_accept_prob
void Mixed_Warped_Curve::compute_log_mh_ratio(){
  proposed_warped_f = proposed_warped_f_basis_mat * common_pars->alpha;
  proposed_fitted_y = current_a(0) + current_a(1) * proposed_warped_f;
  current_warped_f = current_warped_f_basis_mat * common_pars->alpha;
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

  mh_log_accept_prob = proposed_minus_current_llk_data + proposed_llk_w - current_llk_w + log_jacobian_term;
}



// Calls compute_log_mh_ratio() and accept/reject the propsed warping
//
// Depends on: see compute_log_mh_ratio()
//
// Change: current_z, current_dw, current_w,
//         current_warped_x, current_warped_f_basis_mat,
//         common_pars
//
// Note: Acceptances are tallied in the table of common_pars->mh_accept_rate_table
void Mixed_Warped_Curve::mh_accept_reject(){
  mh_randu = gsl_rng_uniform(rng_gen);
  if (std::log(mh_randu) < mh_log_accept_prob) {
    current_z = proposed_z;
    current_dw = proposed_dw;
    current_w = proposed_w;
    current_warped_x = proposed_warped_x;
    current_warped_f_basis_mat = proposed_warped_f_basis_mat;
    ++(common_pars->mh_accept_rate_table(curve_id, common_pars->mh_accept_rate_table_counter));
    return;
  }
  return;
}



// Draw a new amplitude effect (a) from the Gibbs sampler
//
// Depends on: current_warped_f, common_pars
//
// Changes: current_a
//
// Notes: Store updated current_a in common_pars->current_a_mat.col(curve_id);
void Mixed_Warped_Curve::draw_new_a(){
  tmp_f_mat.col(1) = current_warped_f;
  tmp_sigma_post = inv(tmp_f_mat.t() * tmp_f_mat / common_pars->sigma2 +
    common_pars->sigma2_a_inv);
  tmp_mu_post = tmp_sigma_post * (tmp_f_mat.t() * y / common_pars->sigma2 +
    common_pars->sigma2_a_inv * common_pars->mu_a);
  current_a = tmp_mu_post + arma::chol(tmp_sigma_post).t() * arma::randn(dim_a);
  common_pars->current_a_mat.col(curve_id) = current_a;
  return;
}



// Draw a new cluster membership (m) from the Gibbs sampler
//
// Depends on: current_dw
//
// Changes: current_m;
void Mixed_Warped_Curve::draw_new_m(){
  // Update posterior probability
  for(cluster_idx = 0; cluster_idx < n_cluster; ++cluster_idx){
    pred_prob_clusters(cluster_idx) = compute_llk_dw(current_dw,
                           common_pars->kappa_clusters.col(cluster_idx));
  }
  pred_prob_clusters = exp(pred_prob_clusters) % common_pars->p_clusters;
  pred_prob_clusters = pred_prob_clusters / sum(pred_prob_clusters);

  // Sample new current_m
  tmp_u = gsl_rng_uniform(rng_gen);
  for(cluster_idx = 0; cluster_idx < n_cluster; ++cluster_idx){
    if(tmp_u < pred_prob_clusters(cluster_idx)){
      current_m = cluster_idx;
      break;
    }
    else{
      tmp_u -= pred_prob_clusters(cluster_idx);
    }
  }
  current_cluster_membership.zeros();
  current_cluster_membership(current_m) = 1;
  common_pars->current_m_vec(curve_id) = current_m;
  return;
}


// Wraper function to run the simulation step
void Mixed_Warped_Curve::do_simulation_step(){
  for(mcmc_iter = 0; mcmc_iter < common_pars->n_burn_mcmc; ++mcmc_iter){
    propose_new_w();
    compute_proposed_warping_and_f_basis_mat();
    compute_log_mh_ratio();
    mh_accept_reject();
    draw_new_a();
    if(common_pars->saem_counter > common_pars->n_burn_saem){
      draw_new_m();
    }
  }
  return;
}


// Centering step for current_a
// Depends on: dim_a, common_pars, current_a
// Changes: current_a
void Mixed_Warped_Curve::center_current_a(){
  if(!common_pars->need_centering)
    return;
  if(dim_a == 2){
    mean_current_a = mean(common_pars->current_a_mat, 1);
    centering_mat(0, 1) = -mean_current_a(0) / mean_current_a(1);
    centering_mat(1, 1) = 1 / mean_current_a(1);
    current_a = centering_mat * current_a;
  }
  else {
    Rcpp::Rcout << "Warning! Centering only supporting for dim_a = 2.";
  }
  return;
}


// Update stochastic approximation of the sufficients statistics
//
// Depends on: n_i, dim_f, current_a, current_dw, current_w,
//             current_warped_f_basis_mat, y, common_pars,
//
// Changes: current_aug_warped_f_basis_mat,
//          current_hat_mat, current_cov_a, current_log_dw,
//          sapprox_a, sapprox_w,
//          sapprox_aug_warped_f_basis_mat, sapprox_hat_mat,
//          sapprox_cov_a, sapprox_log_dw
void Mixed_Warped_Curve::update_sufficient_statistics_approximates(){
  current_step_size = common_pars->saem_step_sizes(common_pars->saem_counter);

  // Compute sufficient statistics based on current MCMC state
  current_aug_warped_f_basis_mat.col(0) = -current_a(0) * arma::ones(n_i);
  current_aug_warped_f_basis_mat.cols(1, dim_f) = current_a(1) * current_warped_f_basis_mat;

  tmp_half_hat_mat = current_aug_warped_f_basis_mat;
  tmp_half_hat_mat.col(0) = y - current_a(0);
  current_hat_mat = tmp_half_hat_mat.t() * tmp_half_hat_mat;
  current_cov_a = (current_a - common_pars->mu_a) * (current_a - common_pars->mu_a).t();

  current_log_dw.zeros();
  current_log_dw.col(current_m) = arma::log(current_dw);

  // Update stochastic approximates
  sapprox_aug_warped_f_basis_mat = (1 - current_step_size) * sapprox_aug_warped_f_basis_mat +
    current_step_size * current_aug_warped_f_basis_mat;
  sapprox_hat_mat = (1 - current_step_size) * sapprox_hat_mat +
    current_step_size * current_hat_mat;
  sapprox_a = (1 - current_step_size) * sapprox_a +
    current_step_size * current_a;
  sapprox_cov_a = (1 - current_step_size) * sapprox_cov_a +
    current_step_size * current_cov_a;
  sapprox_w = (1 - current_step_size) * sapprox_w +
    current_step_size * current_w;
  sapprox_log_dw = (1 - current_step_size) * sapprox_log_dw +
    current_step_size * current_log_dw;
  sapprox_cluster_membership = (1 - current_step_size) * sapprox_cluster_membership +
    current_step_size * current_cluster_membership;
  return;
}



// Return fitted curve, predicted warping functions and sufficient statistics
Rcpp::List Mixed_Warped_Curve::return_list(){
  return Rcpp::List::create(
    Rcpp::Named("curve_id", Rcpp::wrap(curve_id)),
    Rcpp::Named("x", Rcpp::wrap(x)),
    Rcpp::Named("y", Rcpp::wrap(y)),
    Rcpp::Named("init_clust", Rcpp::wrap(init_clust)),
    Rcpp::Named("warped_x", Rcpp::wrap(h_basis_mat * sapprox_w)),
    Rcpp::Named("sapprox_a", Rcpp::wrap(sapprox_a)),
    Rcpp::Named("sapprox_w", Rcpp::wrap(sapprox_w)),
    Rcpp::Named("sapprox_log_dw", Rcpp::wrap(sapprox_log_dw)),
    Rcpp::Named("sapprox_cluster_membership", Rcpp::wrap(sapprox_cluster_membership))
  );
};



// Return fitted curve, predicted warping functions and sufficient statistics
Rcpp::List Mixed_Warped_Curve::return_list(double y_scaling_factor){
  // Scaling matrices
  return Rcpp::List::create(
    Rcpp::Named("curve_id", Rcpp::wrap(curve_id)),
    Rcpp::Named("x", Rcpp::wrap(x)),
    Rcpp::Named("y", Rcpp::wrap(y * y_scaling_factor)),
    Rcpp::Named("init_clust", Rcpp::wrap(init_clust)),
    Rcpp::Named("warped_x", Rcpp::wrap(h_basis_mat * sapprox_w)),
    Rcpp::Named("fitted_y", Rcpp::wrap(-sapprox_aug_warped_f_basis_mat * common_pars->alpha_aug * y_scaling_factor)),
    Rcpp::Named("sapprox_a", Rcpp::wrap(sapprox_a)),
    Rcpp::Named("sapprox_w", Rcpp::wrap(sapprox_w)),
    Rcpp::Named("sapprox_log_dw", Rcpp::wrap(sapprox_log_dw)),
    Rcpp::Named("sapprox_cluster_membership", Rcpp::wrap(sapprox_cluster_membership))
  );
};


void Mixed_Warped_Curve::print_random_number(){
  double u = gsl_rng_uniform(rng_gen);
  Rcpp::Rcout << "Curve: " << curve_id << " random number: " << u << std::endl;
}

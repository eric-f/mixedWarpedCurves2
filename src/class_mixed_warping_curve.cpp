// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include "class_mixed_warping_model.h"
#include "class_mixed_warping_curve.h"
#include "util.h"

// Constructor
Mixed_Warping_Curve::Mixed_Warping_Curve(Rcpp::List data, Mixed_Warping_Model* pars, int id, int seed) : curve_id(id){

  // Point to common pars
  common_pars = pars;

  // id for internal use

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
  // Allocate a cubic bspline workspace
  tmp_b_vec = gsl_vector_alloc(dim_f);
  tmp_bw = gsl_bspline_alloc(common_pars->f_order, common_pars->f_break_points.size());

  // Metropolis-Hasting
  current_a = common_pars->mu_a * arma::ones(n_cluster).t();      // dim_a x n_cluster
  current_dw = common_pars->kappa_id * arma::ones(n_cluster).t(); // initial with no warping
  current_w = arma::zeros(dim_w);
  current_w(arma::span(1,dim_w-1), arma::span(0, n_cluster)) = arma::cumsum(current_dw, 0);
  current_log_dw = arma::log(current_dw);                         // (dim_w - 1) x n_cluster
  // current_z: Transform current_w back to euclidean space
  current_z = common_pars->chol_centering_mat.t() * (current_log_dw - mean(current_log_dw));
  current_warped_x = x;
  if((current_warped_x.min() < common_pars->h_left_bound) ||
     (current_warped_x.max() > common_pars->h_right_bound)){
    current_warped_x = arma::clamp(current_warped_x, common_pars->h_left_bound, common_pars->h_right_bound);
  }
  current_warped_f_basis_mat = arma::zeros(n_i, dim_f, n_cluster);

  proposed_dw = arma::zeros(dim_w - 1); // (dim_w - 1) x 1
  proposed_w = arma::zeros(dim_w);      // dim_w x 1
  proposed_z = arma::zeros(dim_w - 2);  // (dim_w - 2) x 1
  proposed_warped_x = arma::zeros(n_i); // n_i x 1
  proposed_warped_f_basis_mat = arma::zeros(n_i, dim_f);

  // Sufficient Statistics
  sapprox_residual_sum_of_squares = arma::vec(n_cluster);          // n_cluster x 1
  sapprox_a = arma::zeros(dim_a, n_cluster);                       // dim_a x n_cluster
  sapprox_sq_a = arma::zeros(dim_a, n_cluster);                    // dim_a x n_cluster
  sapprox_w = arma::zeros(dim_w, n_cluster);                       // dim_w x n_cluster
  sapprox_dw = arma::zeros(dim_w-1, n_cluster);                    // (dim_w - 1) x n_cluster
  sapprox_log_dw = arma::zeros(dim_w - 1, n_cluster);              // (dim_w - 1) x n_cluster
  sapprox_warped_f_basis_mat = arma::zeros(n_i, dim_f, n_cluster); // n_i x dim_f x n_cluster
  sapprox_warped_f = arma::zeros(n_i, n_cluster);                  // n_i x n_cluster
  sapprox_fitted_y = arma::zeros(n_i, n_cluster);                  // n_i x n_cluster

  current_yy = arma::zeros(n_cluster);              // n_cluster x 1
  current_Xy = arma::zeros(dim_f, n_cluster);       // dim_f x n_cluster
  current_XX = arma::zeros(n_i, dim_f, n_cluster);  // dim_f x dim_f x n_cluster

  sapprox_yy = arma::zeros(n_cluster);              // n_cluster x 1
  sapprox_Xy = arma::zeros(dim_f, n_cluster);       // dim_f x n_cluster
  sapprox_XX = arma::zeros(n_i, dim_f, n_cluster);  // dim_f x dim_f x n_cluster

  pred_prob_clusters = arma::zeros(n_cluster);      // n_cluster x 1
  pred_prob_clusters(init_clust) = 1;

  // === Private ===

  // Random number generator
  rng_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_gen, seed + curve_id);

  // Temporary or working variables
  // ... for update_sufficient_statistics_approximates()
  current_step_size = 0.0;

  // Internal functions for the MH-within-Gibbs sampler
  // --------------------------------------------------
  // propose_new_w()
  tmp_dw = arma::zeros(dim_w - 1);
  // --------------------------------------------------
  // compute_proposed_warping_and_f_basis_mat()
  // --------------------------------------------------
  // compute_log_mh_ratio()
  proposed_warped_f = arma::zeros(n_i);
  proposed_fitted_y = arma::zeros(n_i);
  current_warped_f = arma::zeros(n_i, n_cluster);
  current_fitted_y = arma::zeros(n_i, n_cluster);
  proposed_residual_sum_of_squares = 0.0;
  current_residual_sum_of_squares = arma::zeros(n_cluster);
  proposed_minus_current_llk_data = 0.0;
  current_llk_w = arma::zeros(n_cluster);
  proposed_llk_w = 0.0;
  log_jacobian_term = 0.0;
  // --------------------------------------------------
  // mh_accept_reject()
  mh_randu = 0.0;
  tmp_mh_log_accept_prob = 0.0;
  // --------------------------------------------------
  // draw_new_a()
  tmp_f_mat = arma::ones(n_i, 2);
  tmp_mu_post = arma::zeros(dim_a);
  tmp_sigma_post = arma::zeros(dim_a, dim_a);
  // --------------------------------------------------
  // center_current_a()
  centering_mat = arma::mat(dim_a, dim_a);
  mean_current_a = arma::vec(dim_a);
  // --------------------------------------------------
  // update_post_prob()
  sapprox_llk_y = arma::zeros(n_cluster);
  sapprox_llk_a = arma::zeros(n_cluster);
  sapprox_llk_w = arma::zeros(n_cluster);

  // Feed back to pars
  // (delete) common_pars->current_m_vec(curve_id) = current_m;
}



// Initialize the warping function basis evaluation matrix
// Depends on: x, common_pars
// Changes: h_basis_mat
void Mixed_Warping_Curve::initialize_h_basis_mat(){
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
void Mixed_Warping_Curve::propose_new_w(int cluster_idx){
  // Update the random walk in R^dim_z from current_z
  proposed_z = current_z.col(cluster_idx) +
    arma::chol(common_pars->identity_cor_mat * common_pars->proposal_sigma).t() * arma::randn(dim_z);
  // update proposed_w
  tmp_dw = common_pars->chol_centering_mat * proposed_z;
  tmp_dw -= mean(tmp_dw);
  tmp_dw = arma::exp(tmp_dw);
  // tmp_dw.transform(exp);
  proposed_dw = tmp_dw / arma::sum(tmp_dw);
  proposed_w(0) = 0;
  proposed_w(arma::span(1, dim_w - 1)) = arma::cumsum(proposed_dw);
  return;
}



// Compute compute_proposed_warping_and_f_basis_mat
// Depends on: h_basis_mat, proposed_w, common_pars
// Changes: proposed_warped_x, proposed_warped_f_basis_mat
void Mixed_Warping_Curve::compute_proposed_warping_and_f_basis_mat(){
  // Evaluate proposed warped time
  proposed_warped_x = h_basis_mat * proposed_w;
  // Squish in place to avoid out-of-bound error in gsl_bspline_eval
  if((proposed_warped_x.min() < common_pars->h_left_bound) ||
     (proposed_warped_x.max() > common_pars->h_right_bound)){
    proposed_warped_x = arma::clamp(proposed_warped_x,
                                    common_pars->h_left_bound,
                                    common_pars->h_right_bound);
  }
  // Evaluate proposed_warped_f_basis_mat
  for(int i = 0; i < n_i; ++i){                                // construct the basis evaluation matrix, warped_f_basis_mat
    gsl_bspline_eval(proposed_warped_x(i), tmp_b_vec, tmp_bw); // compute B_j(x_i) for all j
    for(int j = 0; j < dim_f; ++j){                            // fill in row i of X
      proposed_warped_f_basis_mat(i,j) = gsl_vector_get(tmp_b_vec, j);
    }
  }
  return;
}



// Compute the Metropolis-Hasting ratio
// Depends on: current_a, common_pars,
//             proposed_warped_f_basis_mat, current_warped_f_basis_mat
// Changes: Nil
// Return: log metropolis-hasting ratio
void Mixed_Warping_Curve::compute_log_mh_ratio(int cluster_idx){

  proposed_warped_f = proposed_warped_f_basis_mat * common_pars->alpha;
  proposed_fitted_y = current_a(0, cluster_idx) + current_a(1, cluster_idx) * proposed_warped_f;

  current_warped_f.col(cluster_idx) = current_warped_f_basis_mat.slice(cluster_idx) * common_pars->alpha;
  current_fitted_y.col(cluster_idx) = current_a(0, cluster_idx) + current_a(1, cluster_idx) * current_warped_f.col(cluster_idx);

  proposed_residual_sum_of_squares = arma::sum(arma::square(y - proposed_fitted_y));
  current_residual_sum_of_squares(cluster_idx) = arma::sum(arma::square(y - current_fitted_y.col(cluster_idx)));

  // Compute the data log-likelihood (up to the common constant term)
  proposed_minus_current_llk_data = (current_residual_sum_of_squares(cluster_idx) -
    proposed_residual_sum_of_squares) / 2.0 / common_pars->sigma2;

  // Compute the dirichlet log-likelihood
  proposed_llk_w = compute_llk_dw(proposed_dw, common_pars->kappa_clusters.col(cluster_idx));
  current_llk_w = compute_llk_dw(current_dw, common_pars->kappa_clusters.col(cluster_idx));

  // Compute the log jacobian term for the MH ratio
  log_jacobian_term = arma::sum(arma::log(proposed_dw) - arma::log(current_dw.col(cluster_idx)));

  tmp_mh_log_accept_prob = proposed_minus_current_llk_data +
    proposed_llk_w - current_llk_w(cluster_idx) + log_jacobian_term;
}



// Calls compute_log_mh_ratio() and accept/reject the propsed warping
// Depends on: see compute_log_mh_ratio()
// Change: current_z, current_dw, current_w, current_warped_x, current_warped_f_basis_mat, common_pars
// Note: Acceptances are tallied in the table of common_pars->mh_accept_rate_table
void Mixed_Warping_Curve::mh_accept_reject(int cluster_idx){
  mh_randu = gsl_rng_uniform(rng_gen);
  if (std::log(mh_randu) < tmp_mh_log_accept_prob) {
    current_z.col(cluster_idx) = proposed_z;
    current_dw.col(cluster_idx) = proposed_dw;
    current_log_dw.col(cluster_idx) = arma::log(proposed_dw);
    current_w.col(cluster_idx) = proposed_w;
    current_warped_x.col(cluster_idx) = proposed_warped_x;
    current_warped_f_basis_mat.slice(cluster_idx) = proposed_warped_f_basis_mat;
    current_warped_f.col(cluster_idx) = proposed_warped_f;
    ++(common_pars->mh_accept_rate_table(curve_id, common_pars->mh_accept_rate_table_counter));
    return;
  }
  return;
}



// Draw a new amplitude effect (a) from the Gibbs sampler
// Depends on: current_warped_f_basis_mat, common_pars
// Changes: current_a
// Notes: Store updated current_a in common_pars->current_a_mat.col(curve_id);
void Mixed_Warping_Curve::draw_new_a(int cluster_idx){
  tmp_f_mat.col(1) = current_warped_f.col(cluster_idx);
  tmp_sigma_post = inv(tmp_f_mat.t() * tmp_f_mat / common_pars->sigma2 +
    common_pars->sigma2_a_inv);
  tmp_mu_post = tmp_sigma_post * (tmp_f_mat.t() * y / common_pars->sigma2 +
    common_pars->sigma2_a_inv * common_pars->mu_a);
  current_a.col(cluster_idx) = tmp_mu_post + arma::chol(tmp_sigma_post).t() * arma::randn(dim_a);
  return;
}


// Wraper function to run the simulation step
void Mixed_Warping_Curve::do_simulation_step(){
  for(int mcmc_iter = 0; mcmc_iter < common_pars->n_burn_mcmc; ++mcmc_iter){
    for(int cluster_iter = 0; cluster_iter < n_cluster; ++cluster_iter){
      propose_new_w(cluster_iter);
      compute_proposed_warping_and_f_basis_mat();
      compute_log_mh_ratio(cluster_iter);
      mh_accept_reject(cluster_iter);
      draw_new_a(cluster_iter);
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
void Mixed_Warping_Curve::update_sufficient_statistics_approximates(){
  current_step_size = common_pars->saem_step_sizes(common_pars->saem_counter);

  // Update stochastic approximates
  sapprox_residual_sum_of_squares = (1 - current_step_size) * sapprox_residual_sum_of_squares +
    current_step_size * current_residual_sum_of_squares;
  sapprox_a = (1 - current_step_size) * sapprox_a +
    current_step_size * current_a;
  sapprox_sq_a = (1 - current_step_size) * sapprox_sq_a +
    current_step_size * arma::square(current_a);

  sapprox_w = (1 - current_step_size) * sapprox_w +
    current_step_size * current_w;
  sapprox_dw = arma::diff(sapprox_w);
  sapprox_log_dw = (1 - current_step_size) * sapprox_log_dw +
    current_step_size * current_log_dw;

  for(int cluster_iter = 0; cluster_iter < n_cluster; ++cluster_iter){
    current_yy(cluster_iter) = arma::sum(arma::square(y - current_a(0, cluster_iter) * arma::ones(n_i)));
    current_Xy.col(cluster_iter) = current_a(1, cluster_iter) *
      (y - current_a(0, cluster_iter) * arma::ones(n_i)).t() *
      current_warped_f_basis_mat.slice(cluster_iter);
    current_XX.slice(cluster_iter) = current_a(1, cluster_iter) * current_a(1, cluster_iter) *
      current_warped_f_basis_mat.slice(cluster_iter).t() *
      current_warped_f_basis_mat.slice(cluster_iter);
  }

  sapprox_yy = (1 - current_step_size) * sapprox_yy +
    current_step_size * current_yy;
  sapprox_Xy = (1 - current_step_size) * sapprox_Xy +
    current_step_size * current_Xy;
  sapprox_XX = (1 - current_step_size) * sapprox_XX +
    current_step_size * current_XX;

  sapprox_warped_f_basis_mat = (1 - current_step_size) * sapprox_warped_f_basis_mat +
    current_step_size * current_warped_f_basis_mat;
  sapprox_warped_f = (1 - current_step_size) * sapprox_warped_f +
    current_step_size * current_warped_f;
  sapprox_fitted_y = (1 - current_step_size) * sapprox_fitted_y +
    current_step_size * current_fitted_y;

  if(common_pars->saem_counter > common_pars->n_burn_saem){
    update_post_prob();
  }

  return;
}


void Mixed_Warping_Curve::update_post_prob(){
  for(int cluster_idx = 0; cluster_idx < n_cluster; ++cluster_idx){

    sapprox_llk_y(cluster_idx) =
      - arma::as_scalar(sapprox_yy(cluster_idx) -
      2 * sapprox_Xy.col((cluster_idx)) * common_pars->alpha +
      common_pars->alpha.t() * sapprox_XX.slice(cluster_idx) * common_pars->alpha) /
        2 / common_pars->sigma2;

    sapprox_llk_a(cluster_idx) =
      - sapprox_sq_a(0, cluster_idx) / 2 / common_pars->sigma2_a(0) -
      (sapprox_sq_a(1, cluster_idx) - 2 * sapprox_a(1, cluster_idx) + 1) /
        2 / common_pars->sigma2_a(1);

    sapprox_llk_w(cluster_idx) =
      compute_llk_dw(sapprox_dw.col(cluster_idx),
                     common_pars->kappa_clusters.col(cluster_idx));
  }
  pred_prob_clusters = sapprox_llk_y + sapprox_llk_a + sapprox_llk_w;
  pred_prob_clusters = exp(pred_prob_clusters) % common_pars->p_clusters;
  pred_prob_clusters = pred_prob_clusters / sum(pred_prob_clusters);
  return;
}


// Return fitted curve, predicted warping functions and sufficient statistics
Rcpp::List Mixed_Warping_Curve::return_list(double y_scaling_factor){
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
    Rcpp::Named("sapprox_log_dw", Rcpp::wrap(sapprox_log_dw))
  );
};

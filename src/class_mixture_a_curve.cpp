// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include "class_mixture_a_model.h"
#include "class_mixture_a_curve.h"
#include "util.h"

// Constructor
Mixture_A_Curve::Mixture_A_Curve(Rcpp::List data, Mixture_A_Pars* pars, int id, int seed) : curve_id(id){

  // Temporary variables
  arma::vec tmp_log_current_dw;
  arma::vec tmp_log_centered_dw;

  // Point to common pars
  common_pars = pars;

  // Raw Data
  y = Rcpp::as<arma::vec>(data["y"]);
  x = Rcpp::as<arma::vec>(data["x"]);
  initial_m = Rcpp::as<arma::vec>(data["init_clust"]).at(0);
  // initial_m = rand() % num_clusters;

  // Dimensions
  n_i = y.size();
  dim_a = common_pars->dim_a;
  dim_alpha = common_pars->alpha.n_rows;    // number of basis coefficients for the base curve
  num_clusters = common_pars->alpha.n_cols; // number of clusters

  // Basis Evaluation Matrix - warping function
  f_basis_mat = arma::zeros(n_i, dim_alpha);

  // Allocate a cubic bspline workspace
  tmp_b_vec = gsl_vector_alloc(dim_alpha);
  tmp_bw = gsl_bspline_alloc(common_pars->f_order, common_pars->f_break_points.size());
  // Computes the knots associated with the given breakpoints and stores them internally in tmp_bw->knots.
  gsl_bspline_knots(common_pars->f_break_points, tmp_bw);

  // MCMC working variables
  current_m = initial_m; // Want to be able to initialize with truth
  current_cluster_membership = arma::zeros(num_clusters); // num_clusters x 1

  // Stochastic approximated sufficient statistics
  sapprox_cluster_membership = arma::zeros(num_clusters); // num_clusters x 1

  // Random number generator
  rng_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_gen, seed + curve_id);

  // Temporary or working variables
  // ... for compute_a_posterior()
  post_f = arma::zeros(n_i, num_clusters);
  post_1f = arma::ones(n_i, 2, num_clusters);
  post_1f_star = arma::ones(n_i, 2, num_clusters);
  post_amp_mu = arma::zeros(dim_a, num_clusters);
  tmp_post_amp_sigma = arma::zeros(dim_a, dim_a);
  post_amp_shsh = arma::zeros(num_clusters);
  post_amp_shsc = arma::zeros(num_clusters);
  post_amp_scsc = arma::zeros(num_clusters);
  post_y_sigma_inv = arma::zeros(n_i, n_i, num_clusters);
  post_y_sigma_log_det_scaled = arma::zeros(num_clusters);
  post_y_sigma_log_det_sign = arma::zeros(num_clusters);
  post_y_llk = arma::zeros(num_clusters);
  // ... for update_sufficient_statistics_approximates()
  current_step_size = 0.0;
  // Temporary variable
  YtY = arma::as_scalar(y.t() * y);
  Yt1 = arma::accu(y);
  Ft1 = arma::vec(num_clusters);
  FtY = arma::vec(num_clusters);
  FtF = arma::vec(num_clusters);
  SS_XtX = arma::zeros(dim_alpha, dim_alpha, num_clusters);
  SS_XtY = arma::zeros(dim_alpha, num_clusters);
  SS_post_prob = arma::zeros(num_clusters); // num_clusters x 1

  // Feed back to pars
  common_pars->current_m_vec(curve_id) = current_m;
}



// Initialize the warping function basis evaluation matrix
// Depends on: x, common_pars
// Changes: h_basis_mat
void Mixture_A_Curve::initialize_f_basis_mat(){
  // Rcpp::Rcout << "initialize_f_basis_mat" << std::endl;
  if((x.min() < common_pars->f_left_bound) ||
     (x.max() > common_pars->f_right_bound)){
    x = arma::clamp(x, common_pars->f_left_bound, common_pars->f_right_bound);
  }

  gsl_vector *tmp_a_vec;
  gsl_bspline_workspace *tmp_aw;

  // allocate a cubic bspline workspace (k = 4)
  tmp_a_vec = gsl_vector_alloc(dim_alpha);
  tmp_aw = gsl_bspline_alloc(common_pars->f_order,
                             common_pars->f_break_points.size());

  // evaluate current_warped_f_basis_mat
  gsl_bspline_knots(common_pars->f_break_points, tmp_aw);      // computes the knots associated with the given breakpoints and
  // stores them internally in tmp_aw->knots.
  for(int i = 0; i < n_i; ++i){                                // construct the basis evaluation matrix, warped_f_basis_mat
    gsl_bspline_eval(x(i), tmp_a_vec, tmp_aw);                 // compute B_j(x_i) for all j
    for(int j = 0; j < dim_alpha; ++j){                        // fill in row i of X
      f_basis_mat(i,j) = gsl_vector_get(tmp_a_vec, j);         // gsl_vector_get(B, j)
    }
  }
  // free GSL workspace
  gsl_bspline_free(tmp_aw);
  gsl_vector_free(tmp_a_vec);

  // Make BtB
  BtB = f_basis_mat.t() * f_basis_mat;

  return;
}



// Draw a new amplitude effect (a) from the Gibbs sampler
// Depends on: current_m, common_pars
// Changes: post_...
// Notes: Store updated current_a in common_pars->current_a_mat.col(curve_id);
void Mixture_A_Curve::compute_a_posterior(){

  post_f = f_basis_mat * common_pars->alpha;
  Ft1 = post_f.t() * arma::ones(n_i);
  FtY = post_f.t() * y;
  FtF = square(post_f).t() * arma::ones(n_i);

  for(int clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    post_1f.slice(clust_idx).col(1) = post_f.col(clust_idx);
    post_1f_star.slice(clust_idx) = post_1f.slice(clust_idx) * sqrt(common_pars->amp_sigma);
    // posterior
    tmp_post_amp_sigma = arma::inv_sympd(
      post_1f.slice(clust_idx).t() * post_1f.slice(clust_idx) / common_pars->sigma2 +
        common_pars->amp_sigma_inverse);
    post_amp_mu.col(clust_idx) = tmp_post_amp_sigma * (
      common_pars->amp_sigma_inverse * common_pars->mu0 +
        post_1f.slice(clust_idx).t() * y / common_pars->sigma2);
    post_amp_shsh(clust_idx) = tmp_post_amp_sigma(0, 0) + post_amp_mu(0, clust_idx) * post_amp_mu(0, clust_idx);
    post_amp_shsc(clust_idx) = tmp_post_amp_sigma(0, 1) + post_amp_mu(0, clust_idx) * post_amp_mu(1, clust_idx);
    post_amp_scsc(clust_idx) = tmp_post_amp_sigma(1, 1) + post_amp_mu(1, clust_idx) * post_amp_mu(1, clust_idx);

    // log-likelihood
    post_y_sigma_inv.slice(clust_idx) =
      arma::eye(n_i, n_i) / common_pars->sigma2 +
      post_1f.slice(clust_idx) * (
          common_pars->sigma2 * common_pars->amp_sigma_inverse +
            post_1f.slice(clust_idx).t()*post_1f.slice(clust_idx)) *
            post_1f.slice(clust_idx).t();
    arma::log_det(post_y_sigma_log_det_scaled(clust_idx),
                  post_y_sigma_log_det_sign(clust_idx),
                  arma::eye(2, 2) +
                    post_1f_star.slice(clust_idx).t() * post_1f_star.slice(clust_idx) / common_pars->sigma2);
    post_y_llk(clust_idx) = - post_y_sigma_log_det_scaled(clust_idx) -
      arma::as_scalar((y - post_f.col(clust_idx)).t() * post_y_sigma_inv.slice(clust_idx) * (y - post_f.col(clust_idx)));
  }

  return;
}



// Draw a new cluster membership (m) from the Gibbs sampler
// Depends on: post_y_llk
// Changes: SS_post_prob;
void Mixture_A_Curve::compute_post_prob(){

  //Rescale post_y_llk
  post_y_llk -= max(post_y_llk) * arma::ones(num_clusters);

  SS_post_prob = exp(post_y_llk) % common_pars->p_clusters;
  SS_post_prob = SS_post_prob / sum(SS_post_prob);

  // Rcpp::Rcout << "SS_post_prob: " << SS_post_prob.t() << std::endl;

  return;
}



// Draw a new cluster membership (m) from the Gibbs sampler
// Depends on: SS_post_prob
// Changes: current_m;
void Mixture_A_Curve::draw_new_m(){
  double u;

  u = gsl_rng_uniform(rng_gen);
  // Random shuffle
  current_m = rand() % num_clusters;
  for(int cluster_idx = 0; cluster_idx < num_clusters; ++cluster_idx){
    if(u < SS_post_prob(cluster_idx)){
      current_m = cluster_idx;
      break;
    }
    else{
      u -= SS_post_prob(cluster_idx);
    }
  }

  current_cluster_membership.zeros();
  current_cluster_membership(current_m) = 1;
  common_pars->current_m_vec(curve_id) = current_m;

  return;
}



// Update stochastic approximation of conditional pmf
void Mixture_A_Curve::update_sapprox_cluster_membership(){

  // Step size
  current_step_size = common_pars->saem_step_sizes(common_pars->saem_counter);

  // Update SA for cluster_membership
  sapprox_cluster_membership =
    (1 - current_step_size) * sapprox_cluster_membership +
    current_step_size * current_cluster_membership;

  // Transfer to SS_post_prob
  SS_post_prob = sapprox_cluster_membership;

  return;
}



// Wraper function to run the simulation step
void Mixture_A_Curve::do_simulation_step(){
  compute_a_posterior();
  compute_post_prob();
  for(int i = 0; i < common_pars->n_burn_mcmc; ++i){
      draw_new_m();
  }
  update_sapprox_cluster_membership();
  update_pieces_for_m_step();
  return;
}



// Wrapper function to run the E step
void Mixture_A_Curve::do_expectation_step(){
  compute_a_posterior();
  compute_post_prob();
  update_pieces_for_m_step();
  return;
}



// Update pieces for M-step
// Depends on:
// Changes:
void Mixture_A_Curve::update_pieces_for_m_step(){

  // Pieces alpha by clusters
  for(int clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    SS_XtX.slice(clust_idx) =
      SS_post_prob(clust_idx) * post_amp_scsc(clust_idx) * BtB;
    SS_XtY.col(clust_idx) =
      SS_post_prob(clust_idx) *
      f_basis_mat.t() *
      (post_amp_mu(1,clust_idx) * y - post_amp_shsc(clust_idx) * arma::ones(n_i));
  }

  // Pieces for sigma2
  SS_sigma2 = arma::as_scalar((
    YtY * arma::ones(num_clusters) +
      post_amp_shsh * n_i +
      post_amp_scsc % FtF -
      2 * post_amp_mu.row(0).t() * Yt1 +
      2 * post_amp_shsc % Ft1 -
      2 * post_amp_mu.row(1).t() % FtY).t() * SS_post_prob) /
        common_pars->n_total;

  // Pieces for sigma2_sh
  SS_sigma2_sh = arma::as_scalar(
    post_amp_shsh.t() * SS_post_prob) /
      common_pars->n_curve;

  // Pieces for sigma2_sc
  SS_sigma2_sc = arma::as_scalar(
    (post_amp_scsc - post_amp_mu.row(1).t() + arma::ones(num_clusters)).t() *
    SS_post_prob) /
    common_pars->n_curve;

  return;
}




// Return fitted curve, predicted warping functions and sufficient statistics
Rcpp::List Mixture_A_Curve::return_list(double y_scaling_factor){
  // Scaling matrices
  return Rcpp::List::create(
    Rcpp::Named("curve_id", Rcpp::wrap(curve_id)),
    Rcpp::Named("x", Rcpp::wrap(x)),
    Rcpp::Named("y", Rcpp::wrap(y * y_scaling_factor)),
    Rcpp::Named("fitted_y", Rcpp::wrap(f_basis_mat * common_pars->alpha * y_scaling_factor)),
    Rcpp::Named("post_amp_mu", Rcpp::wrap(post_amp_mu)),
    Rcpp::Named("post_amp_shsh", Rcpp::wrap(post_amp_shsh)),
    Rcpp::Named("post_amp_shsc", Rcpp::wrap(post_amp_shsc)),
    Rcpp::Named("post_amp_scsc", Rcpp::wrap(post_amp_scsc)),
    Rcpp::Named("sapprox_cluster_membership", Rcpp::wrap(SS_post_prob))
  );
};

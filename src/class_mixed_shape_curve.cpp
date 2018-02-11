// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include "class_mixed_shape_model.h"
#include "class_mixed_shape_curve.h"
#include "util.h"

// Constructor
Mixed_Shape_Curve::Mixed_Shape_Curve(Rcpp::List data,
                                     Mixed_Shape_Model* pars,
                                     int id,
                                     double y_scl_r,
                                     int seed) : curve_id(id), y_scl(y_scl_r){

  // Point to common pars
  common_pars = pars;

  // Raw Data
  y = Rcpp::as<arma::vec>(data["y"]);
  x = Rcpp::as<arma::vec>(data["x"]);

  // Dimensions
  n_i = y.size();
  dim_a = 2;
  dim_alpha = common_pars->alpha.n_rows;    // number of basis coefficients for the base curve
  num_clusters = common_pars->alpha.n_cols; // number of clusters

  // Basis Evaluation Matrix
  f_basis_mat = arma::zeros(n_i, dim_alpha);

  // Initialize with User Input
  initial_m = Rcpp::as<arma::ivec>(data["init_clust"]).at(0);
  // Random initialization
  // initial_m = rand() % num_clusters;

  // MCMC and Stochastic approximation
  current_m = initial_m;                                  // Want to be able to initialize with truth
  current_cluster_membership = arma::zeros(num_clusters); // num_clusters x 1
  sapprox_cluster_membership = arma::zeros(num_clusters); // num_clusters x 1
  current_a = arma::zeros(dim_a, num_clusters);           // dim_a x num_clusters

  // Feed back to pars
  common_pars->current_m_vec(curve_id) = current_m;

  // Sufficient Statistics
  SS_sigma2_sh = 0.0;
  SS_sigma2_sc = 0.0;
  SS_XtX = arma::zeros(dim_alpha, dim_alpha, num_clusters);
  SS_XtY = arma::zeros(dim_alpha, num_clusters);
  SS_sigma2 = 0.0;
  SS_post_prob = arma::zeros(num_clusters);

  // Random number generator
  rng_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_gen, seed + curve_id);

  // Working variables
  // ... B-spline workspace
  tmp_b_vec = gsl_vector_alloc(dim_alpha);
  tmp_bw = gsl_bspline_alloc(common_pars->f_order, common_pars->f_break_points.size());
  gsl_bspline_knots(common_pars->f_break_points, tmp_bw);
  // ... Fixed terms in the likelihood function
  YtY = arma::as_scalar(y.t() * y);
  Yt1 = arma::accu(y);
  BtB = arma::mat(dim_alpha, dim_alpha);
  Ft1 = arma::vec(num_clusters);
  FtY = arma::vec(num_clusters);
  FtF = arma::vec(num_clusters);
  // ... Fitted base shapes
  fitted_f = arma::zeros(n_i, num_clusters);
  fitted_1f_slice = arma::ones(n_i, 2);
  fitted_1f_star_slice = arma::ones(n_i, 2);
  // ... Posterior moments and likelihood
  post_amp_mu = arma::zeros(dim_a, num_clusters);
  post_amp_sigma_slice = arma::zeros(dim_a, dim_a);
  post_amp_sigma = arma::zeros(dim_a, dim_a, num_clusters);
  post_amp_shsh = arma::zeros(num_clusters);
  post_amp_shsc = arma::zeros(num_clusters);
  post_amp_scsc = arma::zeros(num_clusters);
  post_y_llk = arma::zeros(num_clusters);
  // ... Stochastic approximation step size
  current_step_size = 0.0;
}



// Initialize the warping function basis evaluation matrix
// Depends on: x, common_pars
// Changes: h_basis_mat
void Mixed_Shape_Curve::initialize_f_basis_mat(){
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
void Mixed_Shape_Curve::compute_a_posterior(){

  fitted_f = f_basis_mat * common_pars->alpha;

  for(clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    // Design matrix for the random effect model
    fitted_1f_slice.col(1) = fitted_f.col(clust_idx);

    // Posterior covariance matrix
    post_amp_sigma.slice(clust_idx) = arma::inv_sympd(
      fitted_1f_slice.t() * fitted_1f_slice / common_pars->sigma2 +
        common_pars->amp_sigma_inverse);
    // Posterior mean
    post_amp_mu.col(clust_idx) = post_amp_sigma.slice(clust_idx) * (
      common_pars->amp_sigma_inverse * common_pars->mu0 +
        fitted_1f_slice.t() * y / common_pars->sigma2);
    // Store as matrices for ease of collapsing in later step
    post_amp_shsh(clust_idx) = post_amp_sigma(0, 0, clust_idx) + post_amp_mu(0, clust_idx) * post_amp_mu(0, clust_idx);
    post_amp_shsc(clust_idx) = post_amp_sigma(0, 1, clust_idx) + post_amp_mu(0, clust_idx) * post_amp_mu(1, clust_idx);
    post_amp_scsc(clust_idx) = post_amp_sigma(1, 1, clust_idx) + post_amp_mu(1, clust_idx) * post_amp_mu(1, clust_idx);
  }

  return;
}



// Draw amplitude effect (a) from each mixture component
// Depends on: post_amp_mu
// Changes: current_m;
void Mixed_Shape_Curve::draw_new_a(){
  for(clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    // Draw a
    current_a.col(clust_idx) =   post_amp_mu.col(clust_idx) +
      arma::chol(post_amp_sigma.slice(clust_idx)).t() * arma::randn(dim_a);
  }
  return;
}



// Compute posterior probability of m conditional on a
// Depends on:
// Changes: SS_post_prob;
void Mixed_Shape_Curve::compute_conditional_post_prob(){

  for(clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    // Log-likelihood
    post_y_llk(clust_idx) =
      - arma::accu(
          arma::square(y -
          current_a(0, clust_idx) * arma::ones(n_i) -
          current_a(1, clust_idx) * fitted_f.col(clust_idx))
      ) / 2 / common_pars->sigma2;
  }

  //Rescale post_y_llk
  post_y_llk -= max(post_y_llk) * arma::ones(num_clusters);

  SS_post_prob = exp(post_y_llk) % common_pars->p_clusters;
  SS_post_prob = SS_post_prob / sum(SS_post_prob);

  return;
}



// Draw new cluster membership (m) from the Gibbs sampler
// Depends on: SS_post_prob
// Changes: current_m;
void Mixed_Shape_Curve::draw_new_m(){
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
void Mixed_Shape_Curve::update_sapprox_cluster_membership(){

  // Step size
  current_step_size = common_pars->saem_step_sizes(common_pars->saem_counter);

  // Update SA for cluster_membership
  sapprox_cluster_membership =
    (1 - current_step_size) * sapprox_cluster_membership +
    current_step_size * current_cluster_membership;

  return;
}



// Update pieces for M-step
// Depends on:
// Changes:
void Mixed_Shape_Curve::update_pieces_for_m_step(){

  // Pieces alpha by clusters
  for(clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    SS_XtX.slice(clust_idx) =
      SS_post_prob(clust_idx) * post_amp_scsc(clust_idx) * BtB;
    SS_XtY.col(clust_idx) =
      SS_post_prob(clust_idx) *
      f_basis_mat.t() *
      (post_amp_mu(1,clust_idx) * y - post_amp_shsc(clust_idx) * arma::ones(n_i));
  }

  // Pieces for sigma2
  Ft1 = fitted_f.t() * arma::ones(n_i);
  FtY = fitted_f.t() * y;
  FtF = square(fitted_f).t() * arma::ones(n_i);
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



// Wraper function to run the MH-based SA-E step
void Mixed_Shape_Curve::do_MH_simulation_step(){

  compute_a_posterior();
  for(int i = 0; i < common_pars->n_burn_mcmc; ++i){
    draw_new_a();
    compute_conditional_post_prob();
    draw_new_m();
  }
  update_sapprox_cluster_membership();
  update_pieces_for_m_step();

  return;
}



// Return fitted curve, predicted warping functions and sufficient statistics
Rcpp::List Mixed_Shape_Curve::return_list(){
  // Scaling matrices
  return Rcpp::List::create(
    Rcpp::Named("curve_id", Rcpp::wrap(curve_id)),
    Rcpp::Named("x", Rcpp::wrap(x)),
    Rcpp::Named("y", Rcpp::wrap(y * y_scl)),
    Rcpp::Named("fitted_base_shapes", Rcpp::wrap(f_basis_mat * common_pars->alpha * y_scl)),
    Rcpp::Named("y_scaling_factor", Rcpp::wrap(y_scl)),
    Rcpp::Named("scaled_post_amp_mu", Rcpp::wrap(post_amp_mu)),
    Rcpp::Named("scaled_post_amp_shsh", Rcpp::wrap(post_amp_shsh)),
    Rcpp::Named("scaled_post_amp_shsc", Rcpp::wrap(post_amp_shsc)),
    Rcpp::Named("scaled_post_amp_scsc", Rcpp::wrap(post_amp_scsc)),
    Rcpp::Named("sapprox_cluster_membership", Rcpp::wrap(SS_post_prob))
  );
};

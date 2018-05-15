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
  initial_m = Rcpp::as<arma::ivec>(data["init_clust"]).at(0);
  // initial_m = rand() % num_clusters;

  // Dimensions
  n_i = y.size();
  dim_a = common_pars->dim_a;
  dim_alpha = common_pars->alpha.n_rows;    // number of basis coefficients for the base curve
  num_clusters = common_pars->alpha.n_cols; // number of clusters

  // Basis Evaluation Matrix
  f_basis_mat = arma::zeros(n_i, dim_alpha);

  // Allocate a cubic bspline workspace
  tmp_b_vec = gsl_vector_alloc(dim_alpha);
  tmp_bw = gsl_bspline_alloc(common_pars->f_order, common_pars->f_break_points.size());
  // Computes the knots associated with the given breakpoints and stores them internally in tmp_bw->knots.
  gsl_bspline_knots(common_pars->f_break_points, tmp_bw);

  // MCMC working variables
  current_m = initial_m;                                  // Want to be able to initialize with truth
  current_cluster_membership = arma::zeros(num_clusters); // num_clusters x 1
  current_a = arma::zeros(dim_a, num_clusters);        // dim_a x num_clusters

  // Stochastic approximated sufficient statistics
  sapprox_cluster_membership = arma::zeros(num_clusters); // num_clusters x 1

  // Random number generator
  rng_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_gen, seed + curve_id);

  // Temporary or working variables
  // ... for compute_a_posterior()
  fitted_f = arma::zeros(n_i, num_clusters);
  fitted_1f_slice = arma::ones(n_i, 2);
  fitted_1f_star_slice = arma::ones(n_i, 2);
  post_amp_mu = arma::zeros(dim_a, num_clusters);
  post_amp_sigma_slice = arma::zeros(dim_a, dim_a);
  post_amp_sigma = arma::zeros(dim_a, dim_a, num_clusters);
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

  fitted_f = f_basis_mat * common_pars->alpha;
  Ft1 = fitted_f.t() * arma::ones(n_i);
  FtY = fitted_f.t() * y;
  FtF = square(fitted_f).t() * arma::ones(n_i);

  for(int clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    fitted_1f_slice.col(1) = fitted_f.col(clust_idx);

    // Posterior
    post_amp_sigma_slice = arma::inv_sympd(
      fitted_1f_slice.t() * fitted_1f_slice / common_pars->sigma2 +
        common_pars->amp_sigma_inverse);

    // if(curve_id == 1){
    //   Rcpp::Rcout << "post_amp_sigma_slice: " << std::endl << post_amp_sigma_slice << std::endl;
    // }

    post_amp_sigma.slice(clust_idx) = post_amp_sigma_slice;
    post_amp_mu.col(clust_idx) = post_amp_sigma_slice * (
      common_pars->amp_sigma_inverse * common_pars->mu0 +
        fitted_1f_slice.t() * y / common_pars->sigma2);
    post_amp_shsh(clust_idx) = post_amp_sigma_slice(0, 0) + post_amp_mu(0, clust_idx) * post_amp_mu(0, clust_idx);
    post_amp_shsc(clust_idx) = post_amp_sigma_slice(0, 1) + post_amp_mu(0, clust_idx) * post_amp_mu(1, clust_idx);
    post_amp_scsc(clust_idx) = post_amp_sigma_slice(1, 1) + post_amp_mu(1, clust_idx) * post_amp_mu(1, clust_idx);
  }

  return;
}



void Mixture_A_Curve::draw_new_a(){
  for(int clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    // Draw a
    current_a.col(clust_idx) =   post_amp_mu.col(clust_idx) +
      arma::chol(post_amp_sigma.slice(clust_idx)).t() * arma::randn(dim_a);
  }
  return;
}



// Compute posterior probability of m conditional on a
// Depends on:
// Changes: SS_post_prob;
void Mixture_A_Curve::compute_conditional_post_prob(){

  for(int clust_idx = 0; clust_idx < num_clusters; ++clust_idx){

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



// Compute marginal posterior probability of m
// Depends on:
// Changes: SS_post_prob;
void Mixture_A_Curve::compute_marginal_post_prob(){

  for(int clust_idx = 0; clust_idx < num_clusters; ++clust_idx){
    fitted_1f_slice.col(1) = fitted_f.col(clust_idx);
    fitted_1f_star_slice = fitted_1f_slice * sqrt(common_pars->amp_sigma);

    // log-likelihood
    // Inverse Covariance Matrix
    post_y_sigma_inv.slice(clust_idx) =
      arma::eye(n_i, n_i) / common_pars->sigma2 +
      fitted_1f_slice * (
          common_pars->sigma2 * common_pars->amp_sigma_inverse +
            fitted_1f_slice.t()*fitted_1f_slice) *
            fitted_1f_slice.t();
    // Determinant of Covariance Matrix
    arma::log_det(post_y_sigma_log_det_scaled(clust_idx),
                  post_y_sigma_log_det_sign(clust_idx),
                  arma::eye(2, 2) +
                    fitted_1f_star_slice.t() * fitted_1f_star_slice / common_pars->sigma2);
    // Log-likelihood
    post_y_llk(clust_idx) =
      - post_y_sigma_log_det_scaled(clust_idx) -
      arma::as_scalar(
        (y - fitted_f.col(clust_idx)).t() *
          post_y_sigma_inv.slice(clust_idx) *
          (y - fitted_f.col(clust_idx))) / 2;
  }

  //Rescale post_y_llk
  post_y_llk -= max(post_y_llk) * arma::ones(num_clusters);

  SS_post_prob = exp(post_y_llk) % common_pars->p_clusters;
  SS_post_prob = SS_post_prob / sum(SS_post_prob);

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



// Update stochastic approximation of amplitude effects
void Mixture_A_Curve::update_sapprox_amp_effect(){

  // Step size
  current_step_size = common_pars->saem_step_sizes(common_pars->saem_counter);

  // Update SA related to amplitude effects
  throw;
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

    // if(curve_id==1) {
    //   Rcpp::Rcout << "curve_id: " << curve_id << std::endl;
    //   Rcpp::Rcout << "clust_idx: " << clust_idx << std::endl;
    //   Rcpp::Rcout << "SS_XtX.slice(clust_idx): " << std::endl << SS_XtX.slice(clust_idx) << std::endl;
    //   Rcpp::Rcout << "SS_XtY.col(clust_idx): " << std::endl <<  SS_XtY.col(clust_idx) << std::endl;
    // }
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

  // if(curve_id==1) {
  //   Rcpp::Rcout << "curve_id: " << curve_id << std::endl;
  //   Rcpp::Rcout << "YtY: " << YtY << std::endl;
  //   Rcpp::Rcout << "post_amp_shsh * n_i: " << post_amp_shsh * n_i << std::endl;
  //   Rcpp::Rcout << "post_amp_scsc % FtF: " << post_amp_scsc % FtF << std::endl;
  //   Rcpp::Rcout << "- 2 * post_amp_mu.row(0).t() * Yt1: " << -2 * post_amp_mu.row(0).t() * Yt1 << std::endl;
  //   Rcpp::Rcout << "2 * post_amp_shsc % Ft1: " << 2 * post_amp_shsc % Ft1 << std::endl;
  //   Rcpp::Rcout << "- 2 * post_amp_mu.row(1).t() % FtY: " << -2 * post_amp_mu.row(1).t() % FtY << std::endl;
  //   Rcpp::Rcout << "SS_sigma2: " << SS_sigma2 << std::endl;
  // }

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
void Mixture_A_Curve::do_MH_simulation_step(){
  // if(curve_id==1){
  //   Rcpp::Rcout << "do_MH_simulation_step..." << std::endl;
  // }
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



// Wraper function to run the Gibbs-based SA-E step
void Mixture_A_Curve::do_Gibbs_expectation_step(){
  // if(curve_id==1){
  //   Rcpp::Rcout << "do_Gibbs_expectation_step..." << std::endl;
  // }
  compute_a_posterior();
  compute_marginal_post_prob();
  for(int i = 0; i < common_pars->n_burn_mcmc; ++i){
    draw_new_m();
  }
  update_sapprox_cluster_membership();
  update_pieces_for_m_step();
  return;
}



// Wrapper function to run the E step
void Mixture_A_Curve::do_expectation_step(){
  // if(curve_id==1){
  //   Rcpp::Rcout << "do_expectation_step..." << std::endl;
  // }
  compute_a_posterior();
  compute_marginal_post_prob();
  update_pieces_for_m_step();
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

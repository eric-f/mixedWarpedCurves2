// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_randist.h>
#include <boost/math/special_functions.hpp>
#include "class_unimodal_lite.h"
#include "util.h"

// Constructor
Unimodal_Model_Lite::Unimodal_Model_Lite(Rcpp::List pars_list,
                                         Rcpp::List aux_list,
                                         RcppGSL::Vector h_break_points_r) : h_break_points(h_break_points_r){
  sigma2 = Rcpp::as<double>(pars_list["sigma2"]);
  mu_a = Rcpp::as<arma::vec>(pars_list["mu_a"]);
  sigma2_a = Rcpp::as<arma::vec>(pars_list["sigma2_a"]);
  sigma2_a_mat = arma::diagmat(sigma2_a);
  sigma2_a_inv = arma::inv_sympd(sigma2_a_mat);
  p_clusters = Rcpp::as<arma::vec>(pars_list["p_clusters"]);
  cluster_sizes = p_clusters * n_curve;
  kappa_clusters = Rcpp::as<arma::mat>(pars_list["kappa_clusters"]);
  kappa_id = Rcpp::as<arma::vec>(pars_list["kappa_id"]);
  prop_tau = Rcpp::as<double>(pars_list["prop_tau"]);
  // Dimension
  dim_a = mu_a.size();
  n_cluster = kappa_clusters.n_cols;
  dim_kappa = kappa_clusters.n_rows;
  dim_w = dim_kappa + 1;
  dim_z = dim_w - 2.0;
  // Auxiliary variables
  n_total = Rcpp::as<int>(aux_list["n_total"]);
  n_curve = Rcpp::as<int>(aux_list["n_curve"]);
  h_order = Rcpp::as<int>(aux_list["h_order"]);
  h_left_bound = gsl_vector_min(h_break_points);
  h_right_bound = gsl_vector_max(h_break_points);
  chol_centering_mat = arma::zeros(dim_kappa, dim_z);
  identity_cor_mat = arma::eye(dim_z, dim_z);
  prop_sigma = Rcpp::as<double>(aux_list["prop_sigma"]);
};

// Initialization
void Unimodal_Model_Lite::generate_chol_centering_mat(){
  if(dim_w < 2)
    Rcpp::stop("Dimension of the warping function must be bigger than 2");
  for(int d = 0; d < dim_z; ++d){
    for(int idx_r = d; idx_r < dim_z + 1; ++idx_r){
      if(d == idx_r){
        chol_centering_mat(idx_r, d) = std::sqrt((dim_z - d) / (dim_z - d + 1));
      }
      else{
        chol_centering_mat(idx_r, d) = -std::sqrt(1 / (dim_z - d) / (dim_z - d + 1));
      }
    }
  }
  return;
}

// Return model
Rcpp::List Unimodal_Model_Lite::return_pars(){
  return Rcpp::List::create(
    Rcpp::Named("sigma2", Rcpp::wrap(sigma2)),
    Rcpp::Named("mu_a", Rcpp::wrap(mu_a)),
    Rcpp::Named("sigma2_a", Rcpp::wrap(sigma2_a)),
    Rcpp::Named("p_clusters", Rcpp::wrap(p_clusters)),
    Rcpp::Named("kappa_clusters", Rcpp::wrap(kappa_clusters))
  );
}

// Return auxiliary information
Rcpp::List Unimodal_Model_Lite::return_aux(){
  return Rcpp::List::create(
    Rcpp::Named("n_cluster", Rcpp::wrap(n_cluster)),
    Rcpp::Named("dim_a", Rcpp::wrap(dim_a)),
    Rcpp::Named("dim_kappa", Rcpp::wrap(dim_kappa)),
    Rcpp::Named("dim_w", Rcpp::wrap(dim_w)),
    Rcpp::Named("dim_z", Rcpp::wrap(dim_z)),
    Rcpp::Named("n_total", Rcpp::wrap(n_total)),
    Rcpp::Named("n_curve", Rcpp::wrap(n_curve)),
    Rcpp::Named("h_order", Rcpp::wrap(h_order)),
    Rcpp::Named("h_break_points", Rcpp::wrap(h_break_points)),
    Rcpp::Named("chol_centering_mat", Rcpp::wrap(chol_centering_mat)),
    Rcpp::Named("identity_cor_mat", Rcpp::wrap(identity_cor_mat)),
    Rcpp::Named("prop_sigma", Rcpp::wrap(prop_sigma))
  );
}

// Constructor
Unimodal_Curve_Lite::Unimodal_Curve_Lite(Rcpp::List data,
                                         Unimodal_Model_Lite* pars,
                                         int mc_mode_r,
                                         int n_mc_r,
                                         int id,
                                         int seed) : n_mc(n_mc_r), mc_mode(mc_mode_r){
  // Point to common pars
  common_pars = pars;

  // Raw Data
  curve_id = Rcpp::as<int>(data["curve_id"]) - 1;
  y = Rcpp::as<arma::vec>(data["y"]);
  x = Rcpp::as<arma::vec>(data["x"]);
  init_clust = Rcpp::as<arma::ivec>(data["init_clust"]).at(0) - 1;

  // Dimensions
  n_i = y.size();
  dim_a = common_pars->dim_a;
  dim_kappa = common_pars->dim_kappa;
  dim_w = common_pars->dim_w;
  dim_z = dim_w - 2;
  n_cluster = common_pars->n_cluster; // number of clusters

  // Basis Evaluation Matrix - warping function
  h_basis_mat = arma::zeros(n_i, dim_w);

  // Sufficient Statistics
  sapprox_a = Rcpp::as<arma::vec>(data["sapprox_a"]);
  sapprox_w = Rcpp::as<arma::vec>(data["sapprox_w"]);
  sapprox_dw = arma::diff(sapprox_w);
  sapprox_cluster_membership = Rcpp::as<arma::vec>(data["sapprox_cluster_membership"]);
  sapprox_warped_x = Rcpp::as<arma::vec>(data["warped_x"]);
  sapprox_fitted_y = Rcpp::as<arma::vec>(data["fitted_y"]);

  // Monte Carlo
  a_star = arma::zeros(dim_a);
  m_star = 0;
  u_star = 0.0;
  dw_star = arma::zeros(dim_kappa);
  w_star = arma::zeros(dim_w);
  warped_x_star = arma::zeros(n_i);
  warped_f_star = arma::zeros(n_i);
  fitted_y_star = arma::zeros(n_i);
  std_resid_star = arma::zeros(n_i);

  mc_rss = arma::zeros(n_mc);
  mc_logLik = arma::zeros(n_mc);
  mc_log_weight = arma::ones(n_mc);
  mc_log_weight_a = arma::ones(n_mc);
  mc_log_weight_w = arma::ones(n_cluster, n_mc);
  w = arma::zeros(dim_w, n_mc);
  a = arma::zeros(dim_a, n_mc);
  m = arma::zeros(n_mc);

  // Random number generator
  rng_gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_gen, seed + curve_id);
}

// Initializataion
void Unimodal_Curve_Lite::initialize_h_basis_mat(){
  if((x.min() < common_pars->h_left_bound) ||
     (x.max() > common_pars->h_right_bound)){
    x = arma::clamp(x, common_pars->h_left_bound, common_pars->h_right_bound);
  }
  // declare temporary variables
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

// Monte Carlo approx to logLik
void Unimodal_Curve_Lite::approx_logLik_full_mc(){
  // Calculate conditional likelihood
  for(idx = 0; idx < n_mc; ++idx){
    a_star = common_pars->mu_a + arma::chol(common_pars->sigma2_a_mat).t() * arma::randn(dim_a);
    m_star = 0;
    u_star = gsl_rng_uniform(rng_gen);
    for(idx_clust = 0; idx_clust < n_cluster-1; ++idx_clust){
      if(u_star < common_pars->p_clusters(idx_clust)){
        break;
      }
      else{
        ++m_star;
        u_star -= common_pars->p_clusters(idx_clust);
      }
    }
    // Sample Dirichlet
    for(idx_kappa = 0; idx_kappa < dim_kappa; ++idx_kappa){
      dw_star(idx_kappa) = gsl_ran_gamma(rng_gen,
              common_pars->kappa_clusters(idx_kappa, m_star), 1);
    }
    dw_star /= arma::sum(dw_star);
    w_star(arma::span(1, dim_kappa)) = arma::cumsum(dw_star);
    warped_x_star = h_basis_mat * w_star;

    warped_f_star = 4 * (arma::square(warped_x_star) - warped_x_star);
    fitted_y_star = a_star(0) + a_star(1) * warped_f_star;
    std_resid_star = (y - fitted_y_star) / std::sqrt(common_pars->sigma2);

    // Compute weight
    // mc_log_weight(idx) = 1;

    // mc_logLik(idx) = -std::log(2*M_PI) / 2.0 - arma::mean(arma::square(std_resid_star)) / 2.0;
    mc_rss(idx) = arma::sum(arma::square(y - fitted_y_star));
    mc_logLik(idx) = - arma::as_scalar(std_resid_star.t() * std_resid_star) / 2.0;
    w.col(idx) = w_star;
    a.col(idx) = a_star;
    m(idx) = m_star;
  }
  return;
}

// Markov Chain Monte Carlo approx to logLik
void Unimodal_Curve_Lite::approx_logLik_full_mcmc(){
  // Working variables
  int current_m=sapprox_cluster_membership.index_max();
  arma::vec current_a = sapprox_a;
  arma::vec current_w = sapprox_w;
  arma::vec current_dw = sapprox_dw;
  arma::vec current_z = common_pars->chol_centering_mat.t() *
    (arma::log(current_dw) - mean(arma::log(current_dw)));
  arma::vec current_warped_x = sapprox_warped_x;
  arma::vec current_warped_f = 4 * (arma::square(current_warped_x) - current_warped_x);
  arma::vec current_fitted_y = sapprox_fitted_y;

  arma::vec proposed_z = current_z;
  arma::vec proposed_dw = current_dw;
  arma::vec proposed_w = current_w;
  arma::vec proposed_warped_x = current_warped_x;
  arma::vec proposed_warped_f = 4 * (arma::square(proposed_warped_x) - proposed_warped_x);
  arma::vec proposed_fitted_y = current_fitted_y;

  double proposed_minus_current_llk_data = 0.0;
  double current_llk_w = 0.0;
  double proposed_llk_w = 0.0;
  double log_jacobian_term = 0.0;
  double mh_log_accept_prob = 0.0;
  double mh_randu = 0.0;

  arma::vec pred_prob_clusters(n_cluster);
  double mh_u = 0.0;

  arma::mat f_mat(n_i, dim_a, arma::fill::ones);
  arma::vec mu_post(dim_a, arma::fill::ones);
  arma::mat sigma_post(dim_a, dim_a, arma::fill::ones);

  // Calculate conditional likelihood
  for(idx = 0; idx < n_mc; ++idx){
    // Gibbs-step for m
    for(int cluster_idx = 0; cluster_idx < n_cluster; ++cluster_idx){
      pred_prob_clusters(cluster_idx) = compute_llk_dw(current_dw,
                         common_pars->kappa_clusters.col(cluster_idx));
    }
    pred_prob_clusters = exp(pred_prob_clusters) % common_pars->p_clusters;
    pred_prob_clusters = pred_prob_clusters / sum(pred_prob_clusters);
    mh_u = gsl_rng_uniform(rng_gen);
    for(int cluster_idx = 0; cluster_idx < n_cluster; ++cluster_idx){
      if(mh_u < pred_prob_clusters(cluster_idx)){
        current_m = cluster_idx;
        break;
      }
      else{
        mh_u -= pred_prob_clusters(cluster_idx);
      }
    }
    // MH-step for w
    // ... Propose
    proposed_z = current_z +
      arma::chol(common_pars->identity_cor_mat * common_pars->prop_sigma).t() * arma::randn(dim_z);
    proposed_dw = common_pars->chol_centering_mat * proposed_z;
    proposed_dw -= mean(proposed_dw);
    proposed_dw = arma::exp(proposed_dw);
    // proposed_dw.transform(exp);
    proposed_dw = proposed_dw / arma::sum(proposed_dw);
    proposed_w(0) = 0;
    proposed_w(arma::span(1, dim_w - 1)) = arma::cumsum(proposed_dw);
    proposed_warped_x = h_basis_mat * proposed_w;
    if((proposed_warped_x.min() < common_pars->h_left_bound) ||
       (proposed_warped_x.max() > common_pars->h_right_bound)){
      proposed_warped_x = arma::clamp(proposed_warped_x,
                                      common_pars->h_left_bound,
                                      common_pars->h_right_bound);
    }
    proposed_warped_f = 4 * (arma::square(proposed_warped_x) - proposed_warped_x);
    proposed_fitted_y = current_a(0) + current_a(1) * 4 * proposed_warped_f;
    current_fitted_y = current_a(0) + current_a(1) * 4 * current_warped_f;
    // ... Reject/Accept
    proposed_minus_current_llk_data =
      (arma::sum(arma::square(y - current_fitted_y)) -
       arma::sum(arma::square(y - proposed_fitted_y))) / 2.0 / common_pars->sigma2;
    proposed_llk_w = compute_llk_dw(proposed_dw, common_pars->kappa_clusters.col(current_m));
    current_llk_w = compute_llk_dw(current_dw, common_pars->kappa_clusters.col(current_m));
    log_jacobian_term = arma::sum(arma::log(proposed_dw) - arma::log(current_dw));
    mh_log_accept_prob = proposed_minus_current_llk_data + proposed_llk_w - current_llk_w + log_jacobian_term;
    mh_randu = gsl_rng_uniform(rng_gen);
    if (std::log(mh_randu) < mh_log_accept_prob) {
      current_z = proposed_z;
      current_dw = proposed_dw;
      current_w = proposed_w;
      current_warped_x = proposed_warped_x;
      current_warped_f = proposed_warped_f;
    }
    // Gibbs-step for a
    f_mat.col(1) = current_warped_f;
    sigma_post = inv(f_mat.t() * f_mat / common_pars->sigma2 +
      common_pars->sigma2_a_inv);
    mu_post = sigma_post * (f_mat.t() * y / common_pars->sigma2 +
      common_pars->sigma2_a_inv * common_pars->mu_a);
    current_a = mu_post + arma::chol(sigma_post).t() * arma::randn(dim_a);
    current_fitted_y = current_a(0) + current_a(1) * current_warped_f;

    // Compute log-likelihood (summand)
    std_resid_star = (y - current_fitted_y) / std::sqrt(common_pars->sigma2);
    // Compute weight
    // mc_log_weight(idx) = 1;

    // mc_logLik(idx) = -std::log(2*M_PI) / 2.0 - arma::mean(arma::square(std_resid_star)) / 2.0;
    mc_rss(idx) = arma::sum(arma::square(y - current_fitted_y));
    mc_logLik(idx) = - arma::as_scalar(std_resid_star.t() * std_resid_star) / 2.0;
    w.col(idx) = current_w;
    a.col(idx) = current_a;
    m(idx) = current_m;
  }
  return;
}

// Importance sampling approx to logLik
void Unimodal_Curve_Lite::approx_logLik_full_is(){
  // 0. Compute reuseable pieces
  // ... for a_star
  arma::mat f_mat(n_i, 2, arma::fill::ones);
  f_mat.col(1) = sapprox_fitted_y;
  arma::mat Sigma_q_inv = f_mat.t() * f_mat / common_pars->sigma2 +
    common_pars->sigma2_a_inv;
  arma::mat Sigma_q = arma::inv(Sigma_q_inv);
  double log_det_q;
  double log_det_hat;
  double diff_log_det;
  double sign;
  arma::log_det(log_det_q, sign, Sigma_q);
  arma::log_det(log_det_hat, sign, common_pars->sigma2_a_mat);
  diff_log_det = (log_det_hat - log_det_q) / 2.0;
  // ... for w_star
  // double tau_dw=common_pars->prop_tau * std::sqrt(n_i);
  double tau_dw=common_pars->prop_tau + n_i;
  arma::vec log_p = arma::log(common_pars->p_clusters);
  arma::vec diff_log_beta(n_cluster);
  arma::mat diff_kappa(dim_kappa, n_cluster);
  for(idx_clust=0; idx_clust < n_cluster; ++idx_clust){
    diff_kappa.col(idx_clust) = common_pars->kappa_clusters.col(idx_clust) - tau_dw*sapprox_dw;
    diff_log_beta(idx_clust) = boost::math::lgamma(arma::sum(common_pars->kappa_clusters.col(idx_clust)));
    diff_log_beta(idx_clust) -= boost::math::lgamma(tau_dw);
    for(idx_kappa=0; idx_kappa < dim_kappa; ++idx_kappa){
      diff_log_beta(idx_clust) -= boost::math::lgamma(common_pars->kappa_clusters(idx_kappa, idx_clust));
      diff_log_beta(idx_clust) += boost::math::lgamma(tau_dw*sapprox_dw(idx_kappa));
    }
  }

  // Rcpp::Rcout << "log_det_q: " << std::endl << log_det_q << std::endl;
  // Rcpp::Rcout << "log_det_hat: " << std::endl << log_det_hat << std::endl;
  // Rcpp::Rcout << "Sigma_q: " << std::endl << Sigma_q << std::endl;
  // Rcpp::Rcout << "Sigma_q_inv: " << std::endl << Sigma_q_inv << std::endl;
  // Rcpp::Rcout << "common_pars->sigma2_a_mat: " << std::endl << common_pars->sigma2_a_mat << std::endl;
  // Rcpp::Rcout << "common_pars->sigma2_a_inv: " << std::endl << common_pars->sigma2_a_inv << std::endl;

  // Loop for Importance sampling
  for(idx = 0; idx < n_mc; ++idx){
    // Rcpp::Rcout << "1. Generate sampling from instrumental distribution ... " << std::endl;
    // 1. Generate sampling from instrumental distribution
    // ... Gaussian a_star
    a_star = sapprox_a + arma::chol(Sigma_q).t() * arma::randn(dim_a);
    // ... Dirichlet w_star
    for(idx_kappa = 0; idx_kappa < dim_kappa; ++idx_kappa){
      dw_star(idx_kappa) = gsl_ran_gamma(rng_gen, tau_dw*sapprox_dw(idx_kappa), 1);
    }
    dw_star /= arma::sum(dw_star);
    w_star(arma::span(1, dim_kappa)) = arma::cumsum(dw_star);
    warped_x_star = h_basis_mat * w_star;

    // Rcpp::Rcout << "2. Compute log-likelihood ... " << std::endl;
    // 2. Compute log-likelihood (summand)
    warped_f_star = 4 * (arma::square(warped_x_star) - warped_x_star);
    fitted_y_star = a_star(0) + a_star(1) * warped_f_star;
    std_resid_star = (y - fitted_y_star) / std::sqrt(common_pars->sigma2);

    // Rcpp::Rcout << "3. Compute weight... " << std::endl;
    // 3. Compute weight
    // ... a_star part
    mc_log_weight_a(idx) =
      -diff_log_det / 2 +
      arma::as_scalar(-(a_star - common_pars->mu_a).t() * common_pars->sigma2_a_inv * (a_star - common_pars->mu_a))/2 -
      arma::as_scalar(-(a_star - sapprox_a).t() * Sigma_q_inv * (a_star - sapprox_a))/2;
    // ... w_star part for each m
    for(idx_clust=0; idx_clust < n_cluster; ++idx_clust){
      mc_log_weight_w(idx_clust, idx) =
        log_p(idx_clust) +
        diff_log_beta(idx_clust) +
        arma::as_scalar(arma::log(dw_star).t() * diff_kappa.col(idx_clust));
    }
    mc_log_weight(idx) = mc_log_weight_a(idx) + std::log(arma::sum(arma::exp(mc_log_weight_w.col(idx))));

    // 4. Store outputs
    // mc_logLik(idx) = -std::log(2*M_PI) / 2.0 - arma::mean(arma::square(std_resid_star)) / 2.0;
    mc_rss(idx) = arma::sum(arma::square(y - fitted_y_star));
    mc_logLik(idx) = - arma::as_scalar(std_resid_star.t() * std_resid_star) / 2.0;
    w.col(idx) = w_star;
    a.col(idx) = a_star;
    m(idx) = m_star;
  }
  return;
}

// Pack and output object
Rcpp::List Unimodal_Curve_Lite::return_obj(){
  // Scaling matrices
  return Rcpp::List::create(
    Rcpp::Named("curve_id", Rcpp::wrap(curve_id)),
    Rcpp::Named("x", Rcpp::wrap(x)),
    Rcpp::Named("y", Rcpp::wrap(y)),
    Rcpp::Named("init_clust", Rcpp::wrap(init_clust)),
    Rcpp::Named("sapprox_a", Rcpp::wrap(sapprox_a)),
    Rcpp::Named("sapprox_w", Rcpp::wrap(sapprox_w)),
    Rcpp::Named("sapprox_cluster_membership", Rcpp::wrap(sapprox_cluster_membership)),
    Rcpp::Named("mc_rss", Rcpp::wrap(mc_rss)),
    Rcpp::Named("mc_logLik", Rcpp::wrap(mc_logLik)),
    Rcpp::Named("mc_log_weight", Rcpp::wrap(mc_log_weight)),
    Rcpp::Named("mc_log_weight_a", Rcpp::wrap(mc_log_weight_a)),
    Rcpp::Named("mc_log_weight_w", Rcpp::wrap(mc_log_weight_w)),
    Rcpp::Named("w", Rcpp::wrap(w)),
    Rcpp::Named("a", Rcpp::wrap(a)),
    Rcpp::Named("m", Rcpp::wrap(m))
  );
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions.hpp>

#include <progress.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

// Function to compute Dirichlet log-likelihood
double dirichlet_llk(arma::vec w, arma::vec alpha) {
  int dim_w;
  bool any_degenerate(false);
  double log_d_terms = 0.0;
  double log_w_terms = 0.0;
  double density = 0.0;

  dim_w = w.size();
  // log reciprocal of generalized beta function
  log_d_terms = + boost::math::lgamma(arma::sum(alpha));
  for(int i=0; i<dim_w; i++){
    log_d_terms -= boost::math::lgamma(alpha(i));
  };
  // check degenerate case
  for(int i=0; i<dim_w; i++){
    if(w(i)==0 && alpha(i)==1){
      any_degenerate = true;
    }
  }
  // density kernel
  if(any_degenerate){
    log_w_terms = R_NegInf;
  } else{
    for(int i=0; i<dim_w; i++){
      log_w_terms += (alpha(i) - 1) * log(w(i));
    }
  }
  density = log_d_terms + log_w_terms;
  return density;
}


//' Internal function for fitting mixture of Dirichlet by EM
//'
//' @param data matrix of data point on simplex
//' @param init_clust initial cluster label
//' @param dim_m number of mixture components
//' @param maxit number of EM iterations
//' @export
// [[Rcpp::export]]
Rcpp::List em_mixture_of_dirichlet(
    Rcpp::NumericMatrix data,
    Rcpp::NumericVector init_clust,
    int dim_m,
    int maxit) {

  // Convert data to arma::mat
  arma::mat w = Rcpp::as<arma::mat>(data);
  int dim_n = w.n_cols;
  int dim_w = w.n_rows;

  // E-step
  arma::mat post_p(dim_n, dim_m, arma::fill::zeros);
  // M-step
  arma::vec p_hat(dim_m, arma::fill::zeros);
  arma::mat alpha_hat(dim_w, dim_m, arma::fill::ones);

  // Counter
  int i, j, m;
  // E-step temporary variable
  double data_llk, tmp_llk;
  arma::mat sum_log_w(dim_w, dim_m);
  // M-step temporary variable
  int newton_idx, newton_inner_idx;
  int newton_max_iter = 100;
  int newton_max_inner_iter = 100;
  arma::vec newton_update_step;
  double newton_abs_tol = 1.0e-8;
  arma::vec newton_g(dim_w, arma::fill::zeros);
  double newton_b = 0.0;
  double newton_z = 0.0;
  arma::vec newton_q(dim_w, arma::fill::zeros);
  arma::vec new_alpha_hat(dim_m);

  // Setup multi-threading
#ifdef _OPENMP
  omp_set_num_threads(1);
  REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif

  // Initialization
  p_hat = 1.0 / dim_m * arma::ones(dim_m);
  for(i = 0; i < dim_n; ++i){
    post_p(i, init_clust(i)-1) = 1;
  }
  alpha_hat *= 5;

  // EM loop
  for(i = 0; i < maxit; ++i){
    // Sufficient Statistics
    sum_log_w = arma::log(w) * post_p;
    // M-step
    p_hat = arma::mean(post_p, 0).t();
    for(m = 0; m < dim_m; ++m){
      if(p_hat(m) > 0){
        for(newton_idx = 0; newton_idx < newton_max_iter; ++newton_idx) {
          // compute hessian and gradient
          newton_z = dim_n * p_hat(m) * boost::math::trigamma(arma::sum(alpha_hat.col(m)));
          newton_g = sum_log_w.col(m) + dim_n * p_hat(m) * boost::math::digamma(arma::sum(alpha_hat.col(m)));
          for(j = 0; j < dim_w; ++j) {
            newton_q(j) = -dim_n * p_hat(m) * boost::math::trigamma(alpha_hat(j, m));
            newton_g(j) -= dim_n * p_hat(m) * boost::math::digamma(alpha_hat(j, m));
          }
          newton_b = arma::sum(newton_g / newton_q) / (1/newton_z + arma::sum(1 / newton_q));

          // step size
          newton_update_step = (newton_g - arma::ones(dim_w) * newton_b) / newton_q;

          // step halving if out-of-bound
          for(newton_inner_idx = 0; newton_inner_idx < newton_max_inner_iter; ++newton_inner_idx){
            new_alpha_hat = alpha_hat.col(m) - newton_update_step;
            if(arma::all(new_alpha_hat > 0)){
              break;
            }
            else{
              newton_update_step /= 2;
            }
          }

          // Check validity and update
          if(arma::all(new_alpha_hat > 0)){
            alpha_hat.col(m) = new_alpha_hat;
          }
          else{
            Rcpp::Rcout << "Maximum step-halving reached. Some estimated alpha still negative." << std::endl;
            throw;
          }
          // Check convergence
          if(arma::max(arma::abs(newton_update_step)) < newton_abs_tol) {
            break;
          }
        }
      }
      else{
        Rcpp::Rcout << "Empty cluster " << m << std::endl;
      }
    }
    // E-step
    for(i = 0; i < dim_n; ++i){
      for(m = 0; m < dim_m; ++m){
        tmp_llk = dirichlet_llk(w.col(i), alpha_hat.col(m));
        post_p(i,m) = p_hat(m) * exp(tmp_llk);
      }
      // normalizing step
      post_p.row(i) /= arma::sum(post_p.row(i));
    }
  }

  // Compute marginal log-likelihood
  data_llk = 0.0;
  for(i=0; i<dim_n; ++i){
    tmp_llk = 0.0;
    for(m=0; m<dim_m; ++m){
      tmp_llk += p_hat(m)*exp(dirichlet_llk(w.col(i), alpha_hat.col(m)));
    }
    data_llk += log(tmp_llk);
  }

  // Return to R
  return(Rcpp::List::create(
      Rcpp::Named("p_hat", p_hat),
      Rcpp::Named("alpha_hat", alpha_hat),
      Rcpp::Named("post_p", post_p),
      Rcpp::Named("llk", data_llk)));
}





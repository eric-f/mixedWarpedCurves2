// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>

// Lite weight model class for logLik
class Unimodal_Model_Lite{
public:
  // Model Parameters
  double sigma2;            // 1 x 1
  arma::vec mu_a;           // 2 x 1
  arma::vec sigma2_a;       // 2 x 1
  arma::mat sigma2_a_mat;   // 2 x 2
  arma::mat sigma2_a_inv;   // 2 x 2
  arma::vec p_clusters;     // n_cluster x 1
  arma::vec cluster_sizes;  // n_cluster x 1
  arma::mat kappa_clusters; // dim_kappa x n_cluster
  arma::mat kappa_id;       // dim_kappa x 1
  double prop_tau;

  // Dimension
  int dim_a;     // number of amplitude effects (default = 2)
  int n_cluster; // number of clusters
  int dim_kappa; // dim_kappa = dim_w - 1 = dimension of kappa;
  int dim_w;     // number of basis coefficients for the warping function
  double dim_z;  // dim_z = dim_w - 2;

  // Auxiliary variables
  int n_total;                    // total number of data points
  int n_curve;                    // number of curves
  int h_order;                    // order of the warping function splines
  RcppGSL::Vector h_break_points; // boundary and internal knot locations of the warping function splines
  double h_left_bound;
  double h_right_bound;
  arma::mat chol_centering_mat;   // (dim_w - 1) x (dim_w - 2), Cholesky decomposition of the I - 1/(dim_w - 1) * J matrix.
  arma::mat identity_cor_mat;     // dim_z x dim_z identity matrix (dim_z = dim_w - 2)
  double prop_sigma;

  // Constructor
  Unimodal_Model_Lite(Rcpp::List pars_list,
                      Rcpp::List aux_list,
                      RcppGSL::Vector h_break_points_r);

  // Initialization
  void generate_chol_centering_mat();

  // Return model
  Rcpp::List return_pars();

  // Return auxiliary information
  Rcpp::List return_aux();
};





// Lite weight data object class for logLik
class Unimodal_Curve_Lite{
public:
  Unimodal_Model_Lite* common_pars;

  int curve_id;
  int n_mc;
  int mc_mode;

  // Counter
  int idx, idx_clust, idx_kappa;

  // Raw Data
  arma::vec y;
  arma::vec x;
  int init_clust;

  // Dimensions
  int n_i;   // number of points per curve
  int dim_a; // number of amplitude effects (default = 2)
  int dim_kappa; // dim_w - 1
  int dim_w; // number of basis coefficients for the warping function
  int dim_z; // dim_w - 2
  int n_cluster;

  // Basis Evaluation Matrix - warping function
  arma::mat h_basis_mat;

  // Sufficient Statistics
  arma::vec sapprox_a;                       // dim_a x 1
  arma::vec sapprox_w;                       // dim_w x 1
  arma::vec sapprox_dw;
  arma::vec sapprox_cluster_membership;      // n_clsuter x 1
  arma::vec sapprox_warped_x;
  arma::vec sapprox_fitted_y;

  // Monte Carlo
  arma::vec a_star;
  int m_star;
  double u_star;
  arma::vec dw_star;
  arma::vec w_star;
  arma::vec warped_x_star;
  arma::vec warped_f_star;
  arma::vec fitted_y_star;
  arma::vec std_resid_star;

  arma::mat w;
  arma::mat a;
  arma::vec m;
  arma::vec mc_rss; // n_mc x 1
  arma::vec mc_logLik; // n_mc x 1
  arma::vec mc_log_weight; // n_mc x 1
  arma::vec mc_log_weight_a; // n_mc x 1
  arma::mat mc_log_weight_w; // n_cluster x n_mc

  // Constructor
  Unimodal_Curve_Lite(Rcpp::List data,
                      Unimodal_Model_Lite* pars,
                      int mc_mode_r,
                      int n_mc_r,
                      int id,
                      int seed);

  // Initializataion
  void initialize_h_basis_mat();

  // MC approximation to logLik
  void approx_logLik_full_mc();
  void approx_logLik_full_mcmc();
  void approx_logLik_full_is();

  // Pack and output object
  Rcpp::List return_obj();

private:
  // random number generator
  gsl_rng * rng_gen;
};

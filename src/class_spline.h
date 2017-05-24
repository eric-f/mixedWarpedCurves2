// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// class_pars.h
#ifndef CLASS_SPLINE_H
#define CLASS_SPLINE_H

class Spline_segment {
public:
  double left_bound;
  double right_bound;
  bool is_first_interval;
  bool is_last_interval;
  arma::vec polynomial_coefficients;

  //method
  double eval_spline(double x);
};

class Spline {
  public:

  int order; // order of the spline, 4 for a cubic spline
  arma::vec break_points;    // boundary and internal knot locations of the base curve
  int n_intervals; // break_points.size() - 1
  std::vector<Spline_segment> spline_segments; // length = n_intervals;

  // Constructor
  Spline(int spline_order, arma::vec spline_break_points);
};

#endif

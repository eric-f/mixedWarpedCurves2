library(gtools)
library(ggplot2)

simplex.y  <- function(x) {
  return( sqrt(0.75) *  x[,3] / rowSums(x) )
}
simplex.x  <- function(x) {
  return( (x[,2] + 0.5 * x[,3]) / rowSums(x) )
}

my_ddirichlet <- cppFunction(depends = c("RcppArmadillo", "BH"),
"
#include <boost/math/special_functions.hpp>
//
double my_ddirichlet(Rcpp::NumericVector x, Rcpp::NumericVector alpha_){
  arma::vec alpha = Rcpp::as<arma::vec>(alpha_);
  double log_gamma_sum_alpha = boost::math::lgamma(arma::sum(alpha));
  arma::vec log_gamma_alpha(alpha.size(), arma::fill::ones);
  for(int i = 0; i < alpha.size(); ++i){
    log_gamma_alpha(i) = boost::math::lgamma(alpha(i));
  }
  arma::vec log_x = log(x);
  double llk = log_gamma_sum_alpha - sum(log_gamma_alpha) +
  arma::as_scalar((alpha - arma::ones(alpha.size())).t() * log_x);
  return(llk);
  }
")


k <- 4
a <- runif(k) * 20
# a <- rep(20, k)

(1 + trigamma(sum(a))*sum(-1/trigamma(a)))*(n^k)*prod(-trigamma(a))

hess <- n * trigamma(sum(a)) - n * diag(trigamma(a))
grad <- digamma(sum(a)) - digamma(a)
det(hess)
det(solve(-hess))

eigen(hess)$value

range(diag(solve(-hess)))

### Plot density on simplex #################################################################
x0 <- expand.grid(seq(0, 1, length=200), seq(0, 1, length=200))
x0 <- x0[rowSums(x0)<=1,]
x0$Var3 = pmax(1 - x0$Var1 - x0$Var2, 0)
den <- ddirichlet(x0, a)

sy <- simplex.y(as.matrix(x0))
sx <- simplex.x(as.matrix(x0))

ggplot() +
  geom_point(aes(x=sx, y=sy, col=den))

m <- 100
a_grid <- seq(1,301, length=m)
a_grid_expand <- expand.grid(a_grid, a_grid)

### llk contour for k=3 #####################################################################
k <- 3
x0 <- rep(1/k, k)
# x0 <- rdirichlet(1, rep(10, k))
den3 <- array(NA, dim=rep(m, k))
for(i in 1:m){
  print(i)
  for(j in 1:m)
    for(l in 1:m)
      den3[i,j,l] <- my_ddirichlet(x0, cbind(a_grid[i], a_grid[j], a_grid[l]))
}

ggplot() +
  geom_point(aes(x=a_grid_expand$Var1, y=a_grid_expand$Var2, col=c(den3[50,,])))

contour(a_grid, a_grid, den3[20,,])

### llk contour for k=2 #####################################################################
k <- 2
x0 <- rep(1/k, k)
# x0 <- rdirichlet(1, rep(10, k))
den2 <- array(NA, dim=rep(m, k))
for(i in 1:m){
  print(i)
  for(j in 1:m)
    den2[i,j] <- my_ddirichlet(x0, c(a_grid[i], a_grid[j]))
}

ggplot() +
  geom_point(aes(x=a_grid_expand$Var1, y=a_grid_expand$Var2, col=c(den2)))

contour(a_grid, a_grid, den2)

curve(trigamma, 10, 20)

-trigamma(a)
max(-trigamma(a)) + length(a) * trigamma(sum(a))



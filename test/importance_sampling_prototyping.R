library(gtools)
library(mvtnorm)
library(splines)

out <- readRDS("~/org/project/Simulations/module/sim-20171113_unimodal_mixture_brillinger_start_3-200x1000/output/unimodal_mixture_brillinger_start_3-200x1000_set_000.rds")
crv <- out$saemOut$curves[[1]]
pars <- out$saemOut$pars
aux <- out$saemOut$aux
aux$h_inner_knots <- rev(aux$h_break_points[-1])[-1]
crv$h_basis <- bs(crv$x, knots = aux$h_inner_knots, degree = aux$h_order-1, intercept = TRUE)

## For simulating a
f_mat <- cbind(1, crv$fitted_y)
Sigma_a <- diag(c(pars$sigma2_a))
Sigma_q_inv <- t(f_mat) %*% f_mat / pars$sigma2 + diag(1/c(pars$sigma2_a))
Sigma_q <- solve(Sigma_q_inv)
log_det_hat = c(determinant(Sigma_a)$modulus)
log_det_q = c(determinant(Sigma_q)$modulus)

## For simulation w
tau_dw <- mean(colSums(pars$kappa_clusters)) + length(crv$y)*2
log_p <- log(pars$p_clusters)
kappa_q <- tau_dw*diff(c(crv$sapprox_w))

## For storage
n_mc <- 1000
std_resid <- numeric(n_mc)
llk <- numeric(n_mc)
w_a_q <- w_a_hat <- w_w_q <- w_w_hat <- weight <- numeric(n_mc)

layout(matrix(1:6, ncol=2))
plot(crv$y, type="l", col="red")
for(i in 1:n_mc){
  a_star <- rmvnorm(1, crv$sapprox_a, Sigma_q)
  dw_star <- rdirichlet(1, kappa_q)
  warped_x_star <- crv$h_basis %*% cumsum(c(0, dw_star))
  fitted_y_star <- a_star[1] + a_star[2] * 4 * (warped_x_star^2 - warped_x_star)
  lines(fitted_y_star, lty=3, col="grey80")
  llk[i] <- sum(dnorm(crv$y - fitted_y_star, sd=sqrt(pars$sigma2), log = TRUE))
  std_resid[i] <- sum(-(crv$y - fitted_y_star)^2) / pars$sigma2 / 2
  w_a_hat[i] <- dmvnorm(a_star, pars$mu_a, Sigma_a)
  w_a_q[i] <- dmvnorm(a_star, crv$sapprox_a, Sigma_q)
  w_w_hat[i] <- sum((exp(log_p) * apply(pars$kappa_clusters, 2, function(x){ddirichlet(dw_star, x)})))
  w_w_q[i] <- ddirichlet(dw_star, kappa_q)
}
lines(crv$y, type="l", col="red")
## Effective sample size
weight = w_a_hat / w_a_q * w_w_hat / w_w_q
range(weight)
range(std_resid)
sprintf("Effective sample size: %f",
        n_mc * mean(weight)^2 / mean(weight^2))
## IS approx to logLik
sprintf("Approximated obs. data. logLik: %f",
        log(mean(exp(std_resid) * weight)))
## Number of duds
sprintf("Number of zero L's: %d", sum(exp(std_resid)==0))
# boxplot(std_resid)
plot(w_a_hat ~ w_a_q)
plot(w_w_hat ~ w_w_q)
plot(exp(std_resid) ~ (w_w_q))
plot(exp(std_resid) ~ (w_a_q))
plot(exp(std_resid) ~ (w_w_q * w_a_q))

layout(1)

## At sapprox
sum(-(crv$y - crv$fitted_y)^2)/pars$sigma2/2
log(
  dmvnorm(c(crv$sapprox_a), pars$mu_a, Sigma_a) /
  dmvnorm(c(crv$sapprox_a), crv$sapprox_a, Sigma_q) *
  sum((exp(log_p) * apply(pars$kappa_clusters, 2, function(x){ddirichlet(diff(c(crv$sapprox_w)), x)}))) / ddirichlet(diff(c(crv$sapprox_w)), kappa_q)
  )

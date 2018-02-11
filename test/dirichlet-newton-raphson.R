rm(list=ls())
library(gtools)
k <- c(1,2,3,3,2,1)
n <- 10
S <- colSums(log(rdirichlet(n, k)))
tau0 <- sum(k)
kappa0 <- k / tau0
llk <- function(x){
  sum((x*kappa0 - 1) * S) -
    n * (sum(lgamma(x*kappa0)) - lgamma(x))
}

S <- c(-2.5068, -1.8054, -1.3960, -1.3965, -1.8016, -2.4910)
n <- 1
kappa0 <- c(0.0833, 0.1667, 0.2500, 0.2500, 0.1667, 0.0833)

s <- seq(1,200, length=100)
l <- sapply(s, llk)
plot(l~s, type="l")
abline(v=tau0, col="red")



tau <- 1000
kappa <- kappa0

s <- seq(1,200, length=100)
l <- sapply(s, llk)
plot(l~s, type="l")
abline(v=tau0, col="red")

iter=0
repeat{
  ## Gradient
  (g <- sum(kappa * S) - n * (sum(digamma(tau * kappa)*kappa) - digamma(tau)))
  if(abs(g) < 1e-8) break
    ## Hessian
    (h <- - n * (sum(trigamma(tau*kappa) * kappa^2) - trigamma(tau)))
  ## Newton's method
  step_size <- g / h
  repeat{
    new_tau <- tau - step_size
    if(new_tau > 0) break
    step_size <- step_size / 2
  }
  tau <- new_tau
  abline(v=tau, col="blue", lty=3)
  iter=iter+1
}
(iter)

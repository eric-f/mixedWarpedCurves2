library(gtools)
library(mixedWarpedCurves2)
library(ggplot2)

kappa1 <- c(0.5,1,1,0.5)
kappa1 <- kappa1 / sum(kappa1)

# for(i in 121:130){
  set.seed(123)
  dat0 <- sim_warping_mixture(100, rep(1, 1),
                              rbind(kappa1),
                              ni = 1001,
                              tau = 10,
                              mu_sh = -25, mu_sc = 500,
                              sd_sh = 10, sd_sc=50, sd_err = 10)
  w <- attr(dat0, "w")
  w_lst <- apply(w, 2, function(x){list(w=x)})
  out <- try(dirichlet_mle(w_lst))
  (out$alpha)
# }

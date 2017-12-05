library(gtools)
library(mixedWarpedCurves2)
library(ggplot2)

kappa0 <- c(1,2,2,1)
kappa0 <- kappa0 / sum(kappa0)
kappa1 <- c(0.5,1,1,0.5)
kappa1 <- kappa1 / sum(kappa1)
kappa2 <- c(1.2,1,1,1.2)
kappa2 <- kappa2 / sum(kappa2)
kappa3 <- c(0.5,1,2,1)
kappa3 <- kappa3 / sum(kappa3)
kappa4 <- c(1,1,0.5,0.5)
kappa4 <- kappa4 / sum(kappa4)

nsim <- 20
ns <- c(100, 1000, 5000)
alpha <- array(NA, c(length(ns), 4, nsim))
sd_alpha <- array(NA, c(length(ns), 4, nsim))

for(j in seq(along=ns))
  for(i in 1:nsim){
    dat0 <- sim_warping_mixture(ns[j], rep(1, 1),
                                rbind(kappa1),
                                ni = 11,
                                tau = 2,
                                mu_sh = -25, mu_sc = 500,
                                sd_sh = 10, sd_sc=50, sd_err = 10)
    w <- attr(dat0, "w")
    w_lst <- apply(w, 2, function(x){list(w=x)})
    out <- try(dirichlet_mle(w_lst))
    (out$alpha)
    alpha[j,,i] <- out$alpha
    sd_alpha[j,,i] <- sqrt(diag(solve(-out$hessian)))
  }


(mean_alpha <- apply(alpha, c(1, 2), mean))
(emp_sd_alpha <- apply(alpha, c(1, 2), sd))
(mod_sd_alpha <- apply(sd_alpha, c(1, 2), mean))


rm(list=ls())
gc()

# source("R/sim_warping_mixture.R")
library(mixedWarpedCurves2)
library(gtools)
library(plyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)


kappa1 <- c(0.3,1,1,0.3)
kappa1 <- kappa1 / sum(kappa1)
kappa2 <- c(1.2,1,1,1.2)
kappa2 <- kappa2 / sum(kappa2)
kappa3 <- c(0.1,1,2,1)
kappa3 <- kappa3 / sum(kappa3)

dat0 <- sim_warping_mixture(200, rep(1/3, 3),
                            rbind(kappa1,
                                  kappa2,
                                  kappa3),
                            ni = 1,
                            tau = 40,
                            mu_sh = -25, mu_sc = 500,
                            sd_sh = 10, sd_sc=50, sd_err = 10)


w <- attr(dat0, "w")
dw <- apply(w, 2, diff)
range(dw)
true_clust <- as.integer(dat0$clust[dat0$x==0])
set.seed(0)


## random start
out <- mixture_of_dirichlet(dw, 3, maxit = 200, nstart=1000)
out$llk
table(true_clust, apply(out$post_p, 1, which.max))

## True start
out <- mixture_of_dirichlet(dw, 3, true_clust)
out$llk
table(true_clust, apply(out$post_p, 1, which.max))

## Kmean start
kmean_clust <- kmeans(t(dw), 3, nstart = 20)$cluster
out <- mixture_of_dirichlet(dw, 3, kmean_clust)
out$llk
pred_clust <- apply(out$post_p, 1, which.max)
table(true_clust, pred_clust)
table(true_clust, kmean_clust)

purity(as.factor(true_clust), as.factor(pred_clust))
entropy(as.factor(true_clust), as.factor(pred_clust))

# dim_n <- ncol(dw)
# dim_w <- nrow(dw)
# dim_m <- ncol(out$alpha_hat)
# tmp_llk <- matrix(NA, dim_m, dim_n)
# dat_llk = 0.0
# for(i in 1:dim_n){
#   for(m in 1:dim_m){
#     tmp_llk[m,i] = out$p_hat[m,1] * ddirichlet(dw[,i], out$alpha_hat[,m])
#   }
# }
#
# sum(log(colSums(tmp_llk)))
#
# table(apply(tmp_llk, 2, which.max), apply(out$post_p, 1, which.max))

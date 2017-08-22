rm(list=ls())
gc()

library(mixedWarpedCurves2)
library(plyr)
library(ggplot2)

# dat0 <- readRDS("../data/test_data.rds")
dat0 <- readRDS("data/test_data.rds")
# dat0 <- readRDS("~/org/project/Simulations/data/w1-20x100/data_w1-20x100-set_000.rds")
my_fit <-
  mixedWarpedCurves2::fsim(y = dat0$y,
                           obs_time = dat0$x,
                           curve_id = dat0$id,
                           n_clust = 1,
                           saem_control = control_saem(n_saem_iter = 1000,
                                                       n_saem_burn = 500,
                                                       saem_step_seq_pow = 0.75,
                                                       n_mcmc_burn = 5,
                                                       n_core = 1,
                                                       accept_rate_window = 5,
                                                       # prop_sigma = 7.8125e-05,
                                                       prop_sigma = 1e-2,
                                                       need_centering = TRUE,
                                                       accept_rate_lb = 0.17,
                                                       accept_rate_ub = 0.33,
                                                       h_n_knots = 3+2,
                                                       f_n_knots = 3+2,
                                                       ind_amp = TRUE))



cbind(diag(solve(-my_fit$se_info$SA_H + my_fit$se_info$SA_C - my_fit$se_info$SA_G %*% t(my_fit$se_info$SA_G))),
      diag(my_fit$se_info$varcov))

my_fit$pars$tau_clusters
sigma_a_tmp <- sapply(my_fit$curves, function(x){x$sapprox_a})
plot(t(sigma_a_tmp))
plot(t(sigma_a_tmp), xlim=c(-1, 1), ylim=c(0.7, 1.3))
points(attributes(dat0)$RE[,1:2], col="red")
rowMeans(sigma_a_tmp)
sigma_a <- sigma_a_tmp %*% t(sigma_a_tmp)
big_sigma <- my_fit$pars$big_sigma
diag(sqrt(1/diag(big_sigma))) %*% big_sigma %*% diag(sqrt(1/diag(big_sigma)))

fisher <- my_fit$se_info$SA_H - my_fit$se_info$SA_C + my_fit$se_info$SA_G %*% t(my_fit$se_info$SA_G)
diag(1/diag(fisher))
diag(fisher)
smat <- my_fit$se_info$scaling_mat
round(sqrt(diag(smat %*% solve(-fisher) %*% smat)), 5)
se <- sqrt(diag(my_fit$se_info$varcov))
round(se, 4)
round(diag(1/se) %*% my_fit$se_info$varcov %*% diag(1/se), 3)

fisher <- my_fit$se_info$SA_H - my_fit$se_info$SA_C + my_fit$se_info$SA_G %*% t(my_fit$se_info$SA_G)
smat2 = diag(-sqrt(1/abs(diag(fisher))))
fisher2 = smat2 %*% fisher %*% smat2
round(sqrt(diag(smat2 %*% solve(-fisher2) %*% smat2)), 5)





out <- ldply(my_fit$curves, function(x){
  data.frame(id=x$curve_id,
             x=x$x,
             y=x$y,
             warped_x=x$warped_x,
             fitted_y=x$fitted_y)
})
out$id <- as.factor(out$id+1)
ggplot(out) +
  geom_line(aes(x=x, y=fitted_y, col=id)) +
  geom_line(aes(x=x, y=y, group=id)) +
  facet_wrap(~id)

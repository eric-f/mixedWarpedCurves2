rm(list=ls())
gc()

library(mixedWarpedCurves2)
# library(mixedWarpedCurves)
library(ggplot2)
library(tidyverse)
library(plyr)

dp <- readRDS("~/org/project/SESfda-simulation-and-analysis/data/00-03_seal-2010-18_interp_y.rds")
time <- as.numeric(rownames(dp))
dat0 <- as.data.frame(dp) %>%
  rownames_to_column("Time") %>%
  gather(key=id, value=y, -Time) %>%
  mutate(
    x = as.numeric(Time),
    id = as.factor(id))

# dat0 <- readRDS("data/test_data.rds")

system.time(my_fit <-
              mixedWarpedCurves2::fsim(y = dat0$y,
                                       obs_time = dat0$x,
                                       curve_id = dat0$id,
                                       n_clust = 1,
                                       saem_control = control_saem(seed = 103123,
                                                                   n_saem_iter = 100,
                                                                   n_saem_burn = 50,
                                                                   saem_step_seq_pow = 1,
                                                                   n_mcmc_burn = 5,
                                                                   n_core = 1,
                                                                   accept_rate_window = 10,
                                                                   prop_sigma = 1e-3,
                                                                   need_centering = TRUE,
                                                                   accept_rate_lb = 0.17,
                                                                   accept_rate_ub = 0.33,
                                                                   h_n_knots = 1+2,
                                                                   f_n_knots = 1+2)))

out <- ldply(my_fit$curves, function(x){
  data.frame(id=x$curve_id,
             x=x$x,
             y=x$y,
             warped_x=x$warped_x,
             fitted_y=x$fitted_y,
             cluster=which.max(x$sapprox_cluster_membership))
})
# out$sys_x <- dat0$system.time
out$id <- as.factor(out$id+1)
out %>%
  ggplot(data=.) +
  geom_line(aes(x=x, y=y, col=id), show.legend = FALSE) +
  geom_line(aes(x=x, y=fitted_y, group=id), show.legend = FALSE) +
  facet_wrap(~id)

my_fit$pars$sigma2
sa_a <- laply(my_fit$curves, function(x){
  x$sapprox_a
})

plot(sa_a[,1]);abline(h=0)
plot(sa_a[,2]);abline(h=1)
colMeans(sa_a)

out %>%
  ggplot(data=.) +
  geom_line(aes(x=x, y=sys_x, col=id))

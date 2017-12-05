rm(list=ls())
gc()

# source("R/sim_warping_mixture.R")
library(mixedWarpedCurves2)
library(plyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(gtools)

# dp <- readRDS("~/org/project/SESfda-simulation-and-analysis/data/00-03_seal-2010-18_interp_y.rds")
# time <- as.numeric(rownames(dp))
# dat0 <- as.data.frame(dp) %>%
#   rownames_to_column("Time") %>%
#   gather(key=id, value=y, -Time) %>%
#   mutate(
#     x = as.numeric(Time),
#     id = as.factor(id))

# dat0 <- readRDS("data/test_data.rds")
# dat0 <- readRDS("~/org/project/Simulations/data/w1-20x100/data_w1-20x100-set_000.rds")


kappa0 <- c(1,2,2,1)
kappa0 <- kappa0 / sum(kappa0)
kappa1 <- c(0,1,1,0)
kappa1 <- kappa1 / sum(kappa1)
kappa2 <- c(1.2,1,1,1.2)
kappa2 <- kappa2 / sum(kappa2)
kappa3 <- c(0,1,2,1)
kappa3 <- kappa3 / sum(kappa3)
kappa4 <- c(1,1,0,0)
kappa4 <- kappa4 / sum(kappa4)

dat0 <- sim_warping_mixture(200, rep(1/3, 3),
                            rbind(kappa1,
                                  kappa2,
                                  kappa3),
                            ni = 1000,
                            tau = 10,
                            mu_sh = -25, mu_sc = 500,
                            sd_sh = 10, sd_sc=50, sd_err = 10)




init_clust <- as.integer(dat0$clust[dat0$x==0]) - 1


x <- unique(dat0$x)
y <- dat0 %>%
  select(x, y, id) %>%
  spread(x, y) %>%
  remove_rownames %>%
  column_to_rownames("id") %>%
  t
clust_lbl <- dat0 %>%
  group_by(id) %>%
  summarise(clust=first(clust))

system.time(my_fit <- mixedWarpedCurves2::mixed_shape(y = y, x = x, nCat=3))
plot(my_fit$mlik_track)

clust_lbl$pred_clust_lbl <- apply(my_fit$p_jk, 1, which.max)
table(clust_lbl$clust, clust_lbl$pred_clust_lbl)
matplot(y, type="l", col=pred_clust_lbl)
matlines(my_fit$a_k, type="l", col="gold", lwd=2)


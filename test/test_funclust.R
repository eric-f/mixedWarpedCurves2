# rm(list=ls())
# gc()

library(plyr)
library(dplyr)
library(gtools)
library(ggplot2)
library(mixedWarpedCurves2)

kappa0 <- c(1,2,2,1)
kappa0 <- kappa0 / sum(kappa0)
kappa1 <- c(0,1,1,0)
kappa1 <- kappa1 / sum(kappa1)
kappa2 <- c(1.5,1,1,1.5)
kappa2 <- kappa2 / sum(kappa2)
kappa3 <- c(0,0,1,1)
kappa3 <- kappa3 / sum(kappa3)
kappa4 <- c(1,1,0,0)
kappa4 <- kappa4 / sum(kappa4)

ncurve <- 200
npt <- 201
dat0 <- sim_warping_mixture(ncurve,
                            rep(1/2, 2),
                            rbind(kappa3,
                                  # kappa2,
                                  kappa4),
                            ni = npt,
                            tau = 1000,
                            mu_sh = 0, mu_sc = 500,
                            sd_sh = 0, sd_sc=100, sd_err = 1)
dat0 %>%
  ggplot() +
  geom_line(aes(x=x, y=y, group=id, col=as.factor(clust)),
            show.legend = FALSE)

f0 <- create.bspline.basis(rangeval=c(0, 1),nbasis=20,norder=4)
fundata <- Data2fd(matrix(dat0$y, ncol=ncurve),
                   argvals=matrix(dat0$x, ncol=ncurve),
                   basisobj=f0);
true_clust <- dat0$clust[dat0$x==0]
plot(fundata, col=true_clust)

my_bs97 <- mixed_shape(y = matrix(dat0$y, ncol=ncurve),
                       x = dat0$x[dat0$id==1], nCat=2,
                       nTry = 30)
bs97_clust <- apply(my_bs97$p_jk,1,which.max)

# with varying dimensions (according to the results of the scree test)
res=funclust(fundata,K=2)
# res=funclust(fundata,K=2,hard=F,nbInit=30,nbIteration=100)

fun_clust <- res$cls

dat0$fun_clust <- rep(fun_clust, each=npt)
dat0$bs97_clust <- rep(bs97_clust, each=npt)

dat0 %>%
  ggplot() +
  geom_line(aes(x=x, y=y, group=id, col=as.factor(clust)),
            show.legend = FALSE) +
  facet_grid(~fun_clust)

res$V

library(mclust)
table(true_clust, bs97_clust)
table(true_clust, fun_clust)
plot(res$meanList[[1]][[1]])
lines(res$meanList[[1]][[2]], col="red")

rm(list=ls())
library(mixedWarpedCurves2)
library(mclust)
library(fda)

(DataID = sprintf("%03d", sample(0:199, 1)))
dat0 <- readRDS(paste0("~/org/project/Simulations/data/unimodal_mixture_3-200x1000/data_unimodal_mixture_3-200x1000-set_", DataID, ".rds"))


### Try another warping
warped_x <- matrix(dat0$warped_x, nrow=1000)
h2 <- function(x){2*x*(1-x)*(x-0.5)+x}
x0 = seq(0, 1, length=1000)
h2_x0 <- h2(x0)
h2_inv <- interpSpline(x0~h2_x0)
new_warped_x <- apply(warped_x, 2, function(wx){predict(h2_inv, wx)$y})
new_warped_x_diff <- new_warped_x - x0



## Change basis coefficients
nbasis = 9L
hbasis_diff <- create.bspline.basis(c(0, 1), nbasis = nbasis, dropind=c(1L, nbasis))
fd_dat0 <- Data2fd(x0, new_warped_x_diff, hbasis_diff)
hbasis <- create.bspline.basis(c(0, 1), nbasis = nbasis)
id <- Data2fd(x0, x0, hbasis)
id$coefs <- round(id$coefs / id$coefs[2])
id$coefs <- id$coefs / rev(id$coefs)[1]
fd_dat0$coefs <- rbind(0, fd_dat0$coefs, 0) + c(id$coefs)
fd_dat0$basis <- hbasis


## Check all points on simplex
new_w <- apply(fd_dat0$coefs, 2, diff)
all(new_w>0)
sum(new_w<0)

w <- apply(attr(dat0, "w"), 2, diff)
all(w>0)
(kmeans0_clust = kmeans(t(w), 3)$cluster)


# saveRDS(new_w, "test/plan_change_of_template3.rds")

# new_w <- readRDS("~/org/lib/mixedWarpedCurves2/test/plan_change_of_template3.rds")

init_clust = dat0$clust[dat0$x==0]
kmeans_clust = kmeans(t(new_w), 3)$cluster
if(all(new_w>0)){
  out <- mixture_of_dirichlet(new_w, nclust=3,
                              init_clust = kmeans_clust,
                              maxit=200, nstart = 50)
  out$llk
  pred_clust <- apply(out$post_p, 1, which.max)

  print(table(init_clust, pred_clust))
  print(adjustedRandIndex(init_clust, pred_clust))
}

table(init_clust, kmeans_clust)
adjustedRandIndex(init_clust, kmeans_clust)

table(init_clust, kmeans0_clust)
adjustedRandIndex(init_clust, kmeans0_clust)

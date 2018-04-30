library(mixedWarpedCurves2)
library(fda)

dat0 <- readRDS("~/org/project/Simulations/data/unimodal_mixture_3-200x1000/data_unimodal_mixture_3-200x1000-set_000.rds")


f <- function(x){
  4*(x-0.5)^2-1
}
f_star <- function(y, a=2){
  abs(2*y-1)^(2*a) - 1
}
h <- function(x, a=2){
  (sign(x-0.5) * abs(2*x-1)^a + 1)/2
}
h_inv <- function(y, a=2){
  (sign(y-0.5)*abs(2*y-1)^(1/a)+1)/2
}

### Un-smoothness at x=0.5 is problematic
a <- 0.8
new_warped_x <- h(dat0$warped_x, a=a)
new_warped_f <- f_star(new_warped_x, a=a)
new_f <- f_star(seq(0, 1, length=1000))

plot(dat0$y0[1:1000]~dat0$warped_x[1:1000], type="l")
curve(f, add=TRUE, col="blue", lty=3, lwd=2)
lines(new_warped_f[1:1000]~new_warped_x[1:1000], col="red")
curve(f_star(x, a), add=TRUE, col="blue", lty=3, lwd=2)


### Try another warping
h2 <- function(x){2*x*(1-x)*(x-0.5)+x}
x0 = seq(0, 1, length=1000)
h2_x0 <- h2(x0)
h2_inv <- interpSpline(x0~h2_x0)
new_warped_x <- apply(warped_x, 2, function(wx){predict(h2_inv, wx)$y})
new_warped_x_diff <- new_warped_x - x

matplot(warped_x, type="l")
matplot(new_warped_x, type="l")
matplot(new_warped_x_diff, type="l")


nbasis = 10L
hbasis_diff <- create.bspline.basis(c(0, 1), nbasis = nbasis, dropind=c(1L, nbasis))
fd_dat0 <- Data2fd(x0, new_warped_x_diff, hbasis_diff)
hbasis <- create.bspline.basis(c(0, 1), nbasis = nbasis)
id <- Data2fd(x0, x0, hbasis)
id$coefs <- round(id$coefs / id$coefs[2])
id$coefs <- id$coefs / rev(id$coefs)[1]
plot(id)
fd_dat0$coefs <- rbind(0, fd_dat0$coefs, 0) + c(id$coefs)
fd_dat0$basis <- hbasis

idx = sample(200, 1)
range(new_warped_x[,idx])
plot(new_warped_x[,idx]~x[,idx], type="l")
lines(warped_x[,idx]~x[,idx], type="l", col="blue")
lines(fd_dat0[idx], col="red")


new_w <- apply(fd_dat0$coefs, 2, diff)

all(new_w>0)


init_clust <- sample(3, size = 200, replace = T)
# em_mixture_of_dirichlet(new_warped_x, init_clust, 3, 1000)




curve(h(x, 1.2))
curve(h_inv(x, 1/1.2), add=T, col="red", lty=3, lwd=3)
curve(f_star(x, 1/1.3))

curve(-2*x^3+3*x^2, from=-.5, to=1.5)

library(mixedWarpedCurves2)

dat0 <- readRDS("../data/test_data.rds")
data <- data.frame(y = dat0$y,
                   x = dat0$x,
                   id = dat0$id)
saem_control <- control_saem(seed=NULL, n_core = 3)
saem_control$seed
pars <- NULL
pars$mu <- c(0, 1)
tmp_y <- tmp_x <- seq(0,1,length=1000)
h_knots <- sort(c(rep(range(saem_control$h_knots), saem_control$h_order-1), saem_control$h_knots))
bhx <- splineDesign(h_knots, tmp_x, saem_control$h_order, rep(0, 1000))
warping_ols <- lm(tmp_y~bhx-1, data=data)
pars$kappa <- diff(unname(warping_ols$coefficients))
f_knots <- sort(c(rep(range(saem_control$f_knots), saem_control$f_order-1), saem_control$f_knots))
bfx <- splineDesign(f_knots, data$x, saem_control$f_order, rep(0, length(data$x)))
shape_ols <- lm(y~bfx-1, data=data)
pars$alpha <- unname(shape_ols$coefficients)
pars$sigma2 <- var(shape_ols$residuals)
pars$num_clusters <- 3
saem_control$n_total = length(dat0$id)
saem_control$n_curve = length(unique(dat0$id))
data_lst <- split(data, data$id)


test_gsl_rng_multithread(data_lst, pars, saem_control, 1)

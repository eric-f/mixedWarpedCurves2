# rm(list=ls())
# gc()

library(mixedWarpedCurves2)
library(mixedWarpedCurves)
# library(profvis)
library(lineprof)
# library(ggplot2)
#
dat0 <- readRDS("~/org/lib/mixedWarpedCurves/data/test_data.rds")
system.time({
# lineprof({
  set.seed(1)
  saemOut <- fsim(
    y = dat0$y,
    obs_time = dat0$x,
    curve_id = dat0$id,
    init_pars = NULL,
    basis.control = attr(dat0, "basis.control"),
    pars.control = control_model_param(fixed.Sigma = FALSE),
    sa.control = control_sa(nIter = 500,
                            alphaSAEM = 0.75,
                            nCore = 1,
                            nBurnSAEM = 100,
                            nBurnMCMC = 5,
                            prop_sigma = 1e-2,
                            centering = TRUE,
                            accept_rate_lb = 0.15,
                            accept_rate_ub = 0.35),
    track_pars = FALSE)
})
# saveRDS(saemOut, "~/org/lib/mixedWarpedCurves/data/saemOut_for_test_memory_use.rds")

library(ggplot2)
saemOutDf <- predict_fsim(saemOut)
saemOutDf$id <- as.factor(saemOutDf$id)
ggplot(saemOutDf) +
  geom_line(aes(x=x, y=y, col=id)) +
  geom_line(aes(x=x, y=fitted_y, group=id), linetype=2) +
  facet_wrap(~id)






saemOut <- readRDS("~/org/lib/mixedWarpedCurves/data/saemOut_for_test_memory_use.rds")
saemOut$input$sa.control$nBurnSAEM <- 200
saemOut$input$sa.control$nIter <- 500
saemOut$input$sa.control$n_accept_rates <- 5
saemOut$input$sa.control$nBurnMCMC <- 5
saemOut$input$sa.control$accept_rate_lb
saemOut$input$sa.control$accept_rate_ub
# sink(file = "~/org/lib/mixedWarpedCurves/test/dump.txt")
system.time(my_fit <- test_oopc(saemOut$data_obj_lst, saemOut$pars, saemOut$input$sa.control))
# sink()



# saemOut$pars[c(2,5,7,4)]

library(plyr)
library(ggplot2)
out <- ldply(my_fit$curves, function(x){
  data.frame(id=x$curve_id, x=x$x, y=x$y, warped_x=x$warped_x, fitted_y=x$fitted_y)
})
ggplot(data=out) +
  geom_line(aes(x=x, y=y, col=as.factor(id))) +
  geom_line(aes(x=x, y=fitted_y, group=as.factor(id)), linetype=2) +
  facet_wrap(~id)

# rm(list=ls())
# gc()

library(mixedWarpedCurves)
# library(profvis)
# library(lineprof)
# library(ggplot2)
#
dat0 <- readRDS("~/org/lib/mixedWarpedCurves/data/test_data.rds")
# # system.time({
# lineprof({
#   set.seed(1)
#   saemOut <- fsim(
#     y = dat0$y,
#     obs_time = dat0$x,
#     curve_id = dat0$id,
#     init_pars = NULL,
#     basis.control = attr(dat0, "basis.control"),
#     pars.control = control_model_param(fixed.Sigma = FALSE),
#     sa.control = control_sa(nIter = 10,
#                             alphaSAEM = 0.75,
#                             nCore = 1,
#                             nBurnSAEM = 0,
#                             nBurnMCMC = 1,
#                             prop_sigma = 1e-2,
#                             centering = TRUE,
#                             accept_rate_lb = 0.15,
#                             accept_rate_ub = 0.35),
#     track_pars = FALSE)
# }) -> tmp
# saveRDS(saemOut, "~/org/lib/mixedWarpedCurves/data/saemOut_for_test_memory_use.rds")

saemOut <- readRDS("~/org/lib/mixedWarpedCurves/data/saemOut_for_test_memory_use.rds")
saemOut$input$sa.control$nBurnSAEM <- 10
saemOut$input$sa.control$nIter <- 20
saemOut$input$sa.control$n_accept_rates <- 5
saemOut$input$sa.control$nBurnMCMC <- 5
saemOut$input$sa.control$accept_rate_lb
saemOut$input$sa.control$accept_rate_ub
saemOut$input$sa.control$alphaSAEM <- 0.7
# sink(file = "~/org/lib/mixedWarpedCurves/test/dump.txt")
test_oopc(saemOut$data_obj_lst, saemOut$pars, saemOut$input$sa.control)
# sink()



# saemOut$pars[c(2,5,7,4)]

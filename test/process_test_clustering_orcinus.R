
out <- readRDS("test/out_saem_dp200_4clust-1496884671.rds")
my_fit <- readRDS("test/fit_saem_dp200_4clust-1496884671.rds")

out <- readRDS("test/out_saem_dp200_4clust-1496886564.rds")
my_fit <- readRDS("test/fit_saem_dp200_4clust-1496886564.rds")

out <- readRDS("test/out_saem_dp200_4clust-1496887085.rds")
my_fit <- readRDS("test/fit_saem_dp200_4clust-1496887085.rds")




my_fit$pars
str(my_fit$pars_track)
str(my_fit$aux)

ggplot(data=out) +
  geom_line(aes(x=x, y=y, col=cluster, group=id), alpha=0.4) +
  scale_y_reverse() +
  facet_grid(~cluster)

ggplot(data=out) +
  geom_line(aes(x=x, y=warped_x, group=id, col=cluster)) +
  geom_abline(intercept = 0, slope=1, col="gold")

out %>%
  dplyr::filter(cluster==4) %>%
  ggplot() +
  geom_line(aes(x=x, y=y, col=id, group=id), alpha=0.4) +
  scale_y_reverse() +
  facet_wrap(~id)

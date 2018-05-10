#' Simulate data - warped sinusodial curve
#'
#' A simulated dataset of 50 curves observed
#' on a common grid of 101 time points over
#' [0, 1]
#'
#' @format A data frame with 5050 rows and 8 variables:
#' \describe{
#'   \item{id}{curve id}
#'   \item{clust}{group label, all equal to 1 for this dataset}
#'   \item{sh}{amplitude shifting effect}
#'   \item{sc}{amplitude scaling effect}
#'   \item{x}{observation time}
#'   \item{warped_x}{warped time, h_i(x)}
#'   \item{y0}{warped curve without noise, a_i,sh + a_i,sc f(h_i(x))}
#'   \item{y}{warped curve with noise}
#' }
"data_ex1"

#' Simulate data - warped unimodal curve
#'
#' A simulated dataset of 50 curves observed
#' on a common grid of 101 time points over
#' [0, 1]
#'
#' @format A data frame with 5050 rows and 8 variables:
#' \describe{
#'   \item{id}{curve id}
#'   \item{clust}{true cluster label}
#'   \item{sh}{amplitude shifting effect}
#'   \item{sc}{amplitude scaling effect}
#'   \item{x}{observation time}
#'   \item{warped_x}{warped time, h_i(x)}
#'   \item{y0}{warped curve without noise, a_i,sh + a_i,sc f(h_i(x))}
#'   \item{y}{warped curve with noise}
#' }
"data_ex2"

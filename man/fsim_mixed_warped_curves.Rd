% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsim_mixed_warped_curves.R
\name{fsim_mixed_warped_curves}
\alias{fsim_mixed_warped_curves}
\title{Model-based curve registration and clustering with unknown base shape}
\usage{
fsim_mixed_warped_curves(y, obs_time, curve_id, init_clust = NULL,
  n_clust = 1, saem_control = control_saem(), trace = FALSE)
}
\arguments{
\item{y}{vector of observed curves}

\item{obs_time}{vector of the observation times}

\item{curve_id}{vector of curve IDs}

\item{init_clust}{vector of inital clutering label with length equals to the number of curves}

\item{n_clust}{integer, number of clusters (set to 1 if not clustering)}

\item{saem_control}{a list of values to control the MCMC and stochastic approximation. See control_saem().}

\item{trace}{if TRUE, tracing information of the estimated parameters are printed}
}
\value{
\describe{
  \item{pars}{List of estimated or fixed model parameters
    \describe{
      \item{alpha}{Estimated B-spline basis coefficient for the common base shape}
      \item{mu_a}{Fixed mean vector for the Gaussian amplitude effects}
      \item{sigma2}{Estimated error variance}
      \item{sigma2_a}{Estimated variance-covariance matrix of the Gaussian amplitude effect}
      \item{p_clusters}{Estimated mixing proportion}
      \item{kappa_id}{Fixed Dirichlet mean vector of the first component where the mean warping function is the identity function}
      \item{tau1}{Overall concentration of the Dirichlet distribution for the first component}
      \item{kappa_clusters}{Estimated Dirichlet concentration parameters}
    }
  }
  \item{curves}{List of curves with stochastic approximation to sufficient statistics, each curves has the following components
    \describe{
      \item{curve_id}{Curve ID. (Caution: this might be different from the inputted curve id, if the original id's is not a sequence from 1 to n.}
      \item{x}{Inputted observation time}
      \item{y}{Inputted observed curves}
      \item{y}{Inputted or random initial cluster label}
      \item{warped_x}{Estimated warped time}
      \item{fitted_y}{Fitted curve}
      \item{sapprox_residual_sum_of_squares}{Stochastic approximation to residual sum of squares}
      \item{sapprox_a}{Stochastic approximation to the conditional expectation of amplitude effects given data}
      \item{sapprox_w}{Stochastic approximation to the conditional expectation of warping coefficients given data}
      \item{sapprox_log_dw}{Stochastic approximation to sufficient statistics for SAEM}
      \item{sapprox_cluster_membership}{Stochastic approximation to predictive probabilities of cluster membership}
    }
  }
  \item{aux}{List of auxiliary information and intermediate variables for MCMC-SAEM}
  \item{pars_track}{Sequence of estimated parameters for convergence diagnostics}
  \item{se_info}{Not currently implemented}
  \item{y_scaling_factor}{Maximum absolute value of the observed curve}
}
}
\description{
This function fits the model \deqn{Y_i(t) = a_{i,sh} + a_{i,sc} f \circ h_i(y) + error}
by maximum likelihood via a stochastic approximation EM algorithm. In the model,
$f$ is a B-spline representing a common shape whereas \eqn{h_i:[0, 1] \to [0, 1]}
is a monotone B-spline representing a random time transformation, referred to as a
warping function or registration function. The vector of (warping) coefficients of $h_i$
follows a Dirichlet distributions. The function also supports a finite mixture of
Dirichlet distributions for the warping coefficient for simultaneous curve clustering and
registration.
}
\references{
Fu, E. and Heckman, N. (2017). Model-based curve registration via stochastic approximation EM algorithm. https://arxiv.org/abs/1712.07265
}
\seealso{
See https://github.com/eric-f/mixedWarpedCurves2 for examples
}

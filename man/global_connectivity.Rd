\name{global_connectivity}
\alias{global_connectivity}
\title{Global spatial-connectivity analysis}
\usage{
global_connectivity(
  x,
  graph,
  reference = "mean",
  scale = "sd",
  side = c("both", "superlevel", "sublevel"),
  n_perm = 999,
  seed = NULL,
  conf = 0.95,
  lambda_grid = NULL,
  n_grid = NULL,
  keep_permutations = FALSE
)
}
\arguments{
\item{x}{Numeric vector with one value per areal unit.}
\item{graph}{An \code{spdep} \code{nb}/\code{listw} object, a two-column edge list, or a symmetric adjacency matrix.}
\item{reference}{A numeric reference level or \code{"mean"}, \code{"median"}, or \code{"zero"}.}
\item{scale}{A positive numeric scale or \code{"sd"}, \code{"mad"}, or \code{"none"}.}
\item{side}{Analyse the superlevel direction, the sublevel direction, or both.}
\item{n_perm}{Number of random relabelings.}
\item{seed}{Optional random seed.}
\item{conf}{Level of the pointwise relabeling envelope.}
\item{lambda_grid}{Optional non-negative threshold grid for trapezoidal discrepancy integration.}
\item{n_grid}{Optional number of equally spaced thresholds from zero to the largest positive selected value.}
\item{keep_permutations}{Whether to retain all relabeled Betti-0 curves.}
}
\description{
Computes the Betti-0 connectivity curve for reference-centered areal data and assesses departures from spatial exchangeability using random relabeling.
}
\details{
The superlevel analysis uses the reference-centered vector, whereas the sublevel analysis applies the same superlevel construction to its sign-reversed version. The Monte Carlo calibration uses the common mean of the observed and all relabeled curves. The plotted reference mean and envelope are computed from the relabeled curves only.

If neither \code{lambda_grid} nor \code{n_grid} is supplied, the integrated absolute discrepancy is evaluated exactly over the distinct positive breakpoints. Supplying a grid uses trapezoidal integration; the Italian case-study reproduction uses \code{n_grid = 80}.
}
\value{
An object of class \code{spatconnect_global}, or an object containing \code{superlevel} and \code{sublevel} analyses when \code{side = "both"}.
}
\examples{
x <- c(3, 2, 1, -1)
edges <- data.frame(i = c(1, 2, 3), j = c(2, 3, 4))
global_connectivity(x, edges, reference = 0, scale = "none",
                    side = "superlevel", n_perm = 99, seed = 1)
}

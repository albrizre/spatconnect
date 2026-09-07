\name{local_connectivity}
\alias{local_connectivity}
\title{Local indicators of spatial connectivity}
\usage{
local_connectivity(
  x,
  graph,
  reference = "mean",
  scale = "sd",
  side = c("both", "superlevel", "sublevel"),
  n_perm = 999,
  seed = NULL,
  adjust = "BH",
  alpha = 0.05,
  activation_tail = 0.15,
  significance = c("unadjusted", "adjusted"),
  classify = TRUE,
  progress = interactive(),
  keep_permutations = FALSE
)
}
\arguments{
\item{x}{Numeric vector with one value per areal unit.}
\item{graph}{An \code{spdep} \code{nb}/\code{listw} object, a two-column edge list, or a symmetric adjacency matrix.}
\item{reference}{A numeric reference level or \code{"mean"}, \code{"median"}, or \code{"zero"}.}
\item{scale}{A positive numeric scale or \code{"sd"}, \code{"mad"}, or \code{"none"}.}
\item{side}{Analyse the superlevel direction, the sublevel direction, or both.}
\item{n_perm}{Number of conditional random relabelings for each focal area.}
\item{seed}{Optional random seed for the conditional relabelings.}
\item{adjust}{Method passed to \code{stats::p.adjust}; use \code{"none"} for no adjustment.}
\item{alpha}{Significance level used in profile classification.}
\item{activation_tail}{Upper-tail proportion of the selected vector used to define high activation.}
\item{significance}{Whether profiles use unadjusted or adjusted p-values.}
\item{classify}{Whether to construct connector profiles.}
\item{progress}{Whether to display a progress bar.}
\item{keep_permutations}{Whether to retain all conditional relabeling values of \eqn{M_i}.}
}
\description{
Computes integrated local activation \eqn{A_i}, merging \eqn{M_i}, and net contribution \eqn{L_i=A_i-M_i}, with conditional random-relabeling inference for the local merging contribution.
}
\details{
When an entering area touches a previously active connected component, one merging unit is generated. Half is assigned to the entering area and half is divided equally among its neighbours in that component. The normalization is applied separately to every previously active component touched by the entering area.

Large \eqn{M_i} values indicate strong cumulative participation in component fusions. Conditional inference fixes the focal selected value and tests whether \eqn{M_i} is unusually large. This is exactly equivalent to a lower-tail test for \eqn{L_i}. The current exact implementation assumes distinct positive values of the selected analysis vector; exact positive ties produce an informative error.
}
\value{
An object of class \code{spatconnect_local}. The \code{results} data frame contains \code{A}, \code{M}, \code{L}, entry merging counts, conditional null summaries, p-values, multiplicity-adjusted p-values, and optional activation--connectivity profiles. For \code{side = "both"}, \code{superlevel} and \code{sublevel} analyses are returned.
}
\examples{
x <- c(3, 2, 1, -1)
edges <- data.frame(i = c(1, 2, 3), j = c(2, 3, 4))
local_connectivity(x, edges, reference = 0, scale = "none",
                   side = "superlevel", n_perm = 99, seed = 1,
                   progress = FALSE)
}

# spatconnect

`spatconnect` is an R package for measuring **spatial connectivity** in data observed on a fixed spatial graph. It implements the global and local methodology developed in the manuscript *Global and local indicators of spatial connectivity for areal data*.

The current version is designed for scalar areal data, although the package name reflects the broader spatial-connectivity perspective of the framework.

## What does `spatconnect` measure?

Most classical spatial statistics focus on spatial association: whether neighboring areas tend to have similar values. `spatconnect` addresses a complementary question:

> **Do high- or low-valued areas form connected spatial regions, and which individual areas contribute most strongly to that connectivity?**

Starting from observations \(x_1,\ldots,x_n\) on an undirected adjacency graph, the package centers and optionally scales the data,

\[
z_i = \frac{x_i-b}{s},
\]

and studies the connected components of thresholded graphs across non-negative thresholds. The **superlevel analysis** uses \(\mathbf z\) to study areas above the reference level \(b\), while the **sublevel analysis** applies the same construction to \(-\mathbf z\) to study areas below it.

The number of connected components at each threshold is summarized by the Betti-0 curve.

## Main functions

The package has two primary user-facing functions:

- `global_connectivity()` computes the Betti-0 connectivity curve and performs a global random-relabeling test for departure from spatial exchangeability.
- `local_connectivity()` computes area-level activation, merging and net contributions, together with conditional random-relabeling inference and optional connectivity profiles.

The spatial graph may be supplied as:

- an `spdep::nb` object;
- an `spdep::listw` object;
- a binary adjacency matrix; or
- a two-column edge list.

Disconnected graphs and isolated areas are allowed.

## Installation

The package is currently available from GitHub:

```r
# install.packages("remotes")
remotes::install_github("albrizre/spatconnect")
```

Alternatively, after cloning the repository locally:

```r
remotes::install_local(".")
```

## Quick example

```r
library(spatconnect)

x <- c(2.4, 1.7, 0.6, -0.2, -1.1, -1.8)
edges <- data.frame(
  i = 1:5,
  j = 2:6
)
```

### Global spatial connectivity

```r
global <- global_connectivity(
  x,
  edges,
  reference = "mean",
  scale = "sd",
  side = "both",
  n_perm = 999,
  seed = 123
)

global$superlevel$p_value
global$sublevel$p_value

head(global$superlevel$curve)
```

For each analysis direction, the observed Betti-0 curve is compared with curves obtained after randomly relabeling the observed values over the fixed graph. The global discrepancy is the integrated absolute difference from the common mean of the observed and relabeled curves.

A small p-value indicates that the observed connectivity profile departs from spatial exchangeability. Because the discrepancy is two-sided with respect to the curves, its direction should be interpreted graphically: an observed Betti-0 curve below the relabeling mean indicates **fewer connected components and greater connectivity**, whereas a curve above the mean indicates **greater fragmentation**.

By default, the discrepancy integral is evaluated exactly over the positive data breakpoints. A fixed threshold grid can also be requested:

```r
global_grid <- global_connectivity(
  x,
  edges,
  side = "superlevel",
  n_perm = 999,
  n_grid = 80,
  seed = 123
)
```

The `n_grid = 80` option is used in the Italian COVID-19 case study accompanying the manuscript.

## Local indicators of spatial connectivity

```r
local <- local_connectivity(
  x,
  edges,
  reference = "mean",
  scale = "sd",
  side = "both",
  n_perm = 999,
  seed = 123,
  progress = FALSE
)

head(local$superlevel$results)
head(local$sublevel$results)
```

For each area \(i\), the returned table contains:

- `A`: integrated activation contribution \(A_i\);
- `M`: integrated merging contribution \(M_i\), the **local indicator of spatial connectivity (LISC)**;
- `L`: net contribution \(L_i=A_i-M_i\);
- `q_entry`: number of previously active components joined when the area enters the filtration;
- `M_null_mean` and `M_null_sd`: conditional random-relabeling summaries for \(M_i\);
- `M_excess`: \(M_i-E_0(M_i)\);
- `M_z_score`: standardized conditional excess;
- `p_value`: unadjusted conditional p-value;
- `p_adjusted`: multiplicity-adjusted p-value; and
- `profile`: optional activation-connectivity classification.

Large values of \(M_i(\mathbf z)\) indicate strong cumulative participation in the fusions through which **high-valued superlevel regions** become connected across thresholds. Large values of \(M_i(-\mathbf z)\) have the analogous interpretation for **low-valued sublevel regions**.

The raw magnitude of \(M_i\) is descriptive: areas that remain active over a wider range of thresholds have more opportunity to participate in merging. Local inference therefore conditions on the focal area's observed value. The resulting upper-tail test for \(M_i\) is exactly equivalent to a lower-tail test for \(L_i=A_i-M_i\).

## Activation--merging decomposition

When an area enters the filtration, it may touch one or more connected components that are already active. Each such contact with a previously disconnected component generates one non-redundant fusion.

For every fusion:

- one half of the merging unit is assigned to the entering area; and
- the remaining half is divided equally among the entering area's neighbors belonging to that previously active component.

The normalization is performed **separately for each connected component** touched by the entering area. This produces an exact local-to-global decomposition of the Betti-0 curve and avoids selecting an arbitrary maximum-spanning forest.

The current exact local implementation assumes that the positive values of the selected analysis vector are distinct, so areas enter the filtration one at a time. Exact positive ties trigger an informative error. Structural ties among edge activation levels do not require a tie-breaking rule under the component-specific allocation.

## Connectivity profiles

When `classify = TRUE`, `local_connectivity()` combines activation, conditional evidence of unusually large merging, and the sign of the net contribution into four descriptive-inferential categories:

- **Core connector**: high activation, significant conditional merging and \(L_i\leq 0\);
- **Bridge connector**: significant conditional merging without belonging to the selected activation tail;
- **Peak connector**: high activation and significant merging with \(L_i>0\) in the superlevel analysis;
- **Trough connector**: the corresponding sublevel category; and
- **Not significant** otherwise.

By default, high activation is defined using the upper 15% tail of the selected analysis vector (`activation_tail = 0.15`). Profiles use unadjusted p-values by default; set `significance = "adjusted"` to classify using multiplicity-adjusted values.

These labels summarize activation-connectivity configurations and should not be interpreted as causal mechanisms.

## Using an areal map

For polygon data, an adjacency graph can be generated with `spdep` and passed directly to `spatconnect`:

```r
library(sf)
library(spdep)
library(spatconnect)

# sf_object must contain one row for each entry of x.
nb <- poly2nb(sf_object, queen = TRUE)

res <- local_connectivity(
  x,
  nb,
  side = "superlevel",
  n_perm = 1999,
  seed = 2026,
  progress = TRUE
)

sf_object$M <- res$results$M
sf_object$p_connectivity <- res$results$p_value
sf_object$q_connectivity <- res$results$p_adjusted
sf_object$connectivity_profile <- res$results$profile
```

## Reference level and scale

The interpretation of connectivity is relative to a chosen reference value \(b\). The package defaults are:

```r
reference = "mean"
scale = "sd"
```

Other supported choices include `reference = "median"`, `reference = "zero"`, a numeric reference value, `scale = "mad"`, `scale = "none"`, or a positive numeric scale.

These choices affect the substantive meaning and units of the thresholds and integrated quantities and should therefore be reported in applications.

## Reproducing the Italian COVID-19 case study

The repository contains:

```text
analysis/reproduce_italy_case_study.R
```

which reproduces the two-wave Italian NUTS-3 COVID-19 analysis accompanying the methodological paper. The script downloads the required public data, constructs queen contiguity, calls `spatconnect` for the proposed global and local connectivity analyses, computes Local Moran's I as a benchmark, and produces the manuscript figures and tables.

The main settings are:

- 80 equally spaced thresholds for the global Betti-0 curves;
- 1999 global random relabelings;
- 1999 conditional local relabelings per focal area;
- 1999 conditional permutations for Local Moran's I;
- within-wave mean as the reference level;
- within-wave standard deviation as the scale; and
- activation-tail proportion \(\eta=0.15\).

The conditional local procedure is the computationally intensive part. With \(n\) areas and \(R\) relabelings, approximately \(nR\) local decompositions are recomputed for each analysis direction.

## Citation

If you use `spatconnect`, please cite the accompanying methodological paper:

> Briz-Redón, Á. and Ruiz-Valderrama, M. (2026). *Global and local indicators of spatial connectivity for areal data*. Manuscript.

Citation information is also available from R:

```r
citation("spatconnect")
```

The citation will be updated when the final bibliographic details of the article are available.

## Development status

`spatconnect` 0.1.0 is the initial implementation accompanying the methodological paper. The current package focuses on scalar areal data represented on a fixed undirected graph. The name `spatconnect` is intentionally broader than the present implementation to leave room for future extensions of the spatial-connectivity framework.

## License

MIT License. See `LICENSE` for details.

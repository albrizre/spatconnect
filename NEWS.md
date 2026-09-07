# spatconnect 0.1.0

- Initial GitHub-oriented release accompanying the methodological paper.
- `global_connectivity()` computes Betti-0 curves and a common-mean
  random-relabeling test using an exchangeable set of the observed curve and
  all relabeled curves.
- `global_connectivity()` supports exact breakpoint integration by default and
  optional fixed-grid trapezoidal integration for exact reproduction of the
  case study.
- `local_connectivity()` implements the component-specific symmetric
  activation--merging allocation used in the current manuscript.
- The local indicator `M_i` measures cumulative participation in component
  fusions; conditional random relabeling tests whether this participation is
  unusually large given the focal activation value and graph position.
- Public terminology is consistently `superlevel` / `sublevel`.
- Removed the earlier maximum-spanning-forest, edge tie-breaking, and
  `tie_sensitivity()` implementation.
- Supports `spdep` neighbour lists, `listw` objects, two-column edge lists, and
  symmetric adjacency matrices.

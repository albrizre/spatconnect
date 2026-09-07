# Italian COVID-19 case study

`reproduce_italy_case_study.R` reproduces the two-wave Italian NUTS-3 case
study accompanying the `spatconnect` methodological paper.

The script downloads the COVID-19 European Regional Tracker data and the NUTS
geometry, constructs queen contiguity, and calls the installed `spatconnect`
package for all proposed global and local connectivity calculations. Local
Moran's I and publication figures are produced in the script because they are
case-study components rather than package methodology.

Main analysis settings:

- 80 equally spaced thresholds for the global Betti-0 curves;
- 1999 global random relabelings;
- 1999 conditional local relabelings per focal area;
- 1999 conditional permutations for Local Moran's I;
- within-wave mean reference and standard-deviation scale;
- empirical activation tail proportion `eta = 0.15`.

The script writes wave-specific figures, tables, spatial outputs, and a combined
summary to its output directory.

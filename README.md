# Reproduction code for the main simulation study

This directory contains the R code used to produce the main simulation results
reported in Section 4 of the manuscript. The study compares eight anomaly
detection procedures for fixed-effects panel models and evaluates common-slope
estimation after the detected units are removed.

## Files

- `simulation_comparison.R` contains the data-generating mechanisms and the
  implementations of all eight procedures.
- `run_expanded_main_comparison.R` runs the detection experiment underlying
  Tables 1--4.
- `run_estimation_comparison.R` runs the common-slope estimation experiment
  underlying Tables 5--6.

The two runner scripts load functions from `simulation_comparison.R`. The latter
should therefore remain in the same directory as the runner scripts. It is not
necessary to run `simulation_comparison.R` separately.

## Simulation design

The formal experiments use `n = 200` panel units, `T = 100` observations per
unit, and `K = m = 4`. Predictor processes are generated under i.i.d., AR, and MA
designs, and the errors follow either a Gaussian or a mixture-normal
distribution. All reported settings use 100 Monte Carlo replications.

The detection experiment uses contamination proportions
`rho_sim = 0.10, 0.20` and anomaly magnitudes
`zeta = 0.10, 0.15, 0.20, 0.25`. The estimation experiment fixes
`rho_sim = 0.20` and uses `zeta = 0.10, 0.30, 0.50`.

The procedures are reported in the following order:

1. Proposed
2. Oracle LM--BH
3. Pooled LM--BH
4. Batch Conformal--BH (Oracle)
5. Batch Conformal--BH (Screened)
6. VSOM
7. Cook's distance
8. DFFITS

The estimation experiment also reports the raw estimator computed from all
units before anomaly removal.

## Software requirements

The final results were produced with R 4.6.0 on Windows 11. The formal runner
scripts use only packages distributed with R, primarily `stats` and `parallel`;
no additional CRAN packages are required. The scripts can also be run on Linux
or macOS with a current R installation.


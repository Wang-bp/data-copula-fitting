# Reproducibility code: Copula-Based Joint Probability and Return-Period Analysis of Wind–Wave Conditions

This repository contains the main calculation code used for the revised manuscript. It reproduces the site-level statistical workflow from the four hourly hindcast time series: marginal models, temporally consistent monthly extreme pairs, single/mixed copulas, model-selection diagnostics, OR/AND return periods, marginal-model sensitivity, seasonal sensitivity, and key diagnostic figures.

## 1. Input data

Put these files in `data/`:

- `time_station_20yrs_2.npy`
- `u10_station_20yrs_2.npy`
- `hs_station_20yrs_2.npy`
- `tp_station_20yrs_2.npy`

The time array is interpreted as **hours since 1970-01-01 00:00:00**, proleptic Gregorian. After unavailable months are excluded, the manuscript analysis contains 161,400 common-valid hourly samples and 221 complete calendar-month extreme pairs.

## 2. Environment

Python 3.10+ is recommended.

```bash
pip install -r requirements.txt
```

The manuscript uses the classical **Sheather-Jones direct plug-in bandwidth** for the log-KDE. To reproduce that selector exactly, install R and make `Rscript` available on the PATH. `R/sj_dpi_bandwidth.R` calls `stats::bw.SJ(..., method="dpi")`. The monthly SJ-DPI bandwidths reported in the Supplementary Material are also stored in `paper_parameters.json` so the return-period calculations do not depend on re-estimating them.

## 3. Run

```bash
python run_all.py
```

Outputs are written to `outputs/`.

## 4. Script-to-manuscript map

- `01_prepare_monthly_extremes.py`: Section 2.7; creates the 221 temporally consistent monthly pairs. Companion variables are retained at the timestamp of the governing monthly maximum.
- `02_fit_marginals.py`: Sections 2.4 and 3.3; log-SJ-DPI KDE, 95th-percentile Kernel-Pareto, GEV, KS/AD and upper-tail diagnostics, plus GPD threshold diagnostics.
- `03_fit_copulas.py`: Sections 2.5–2.6 and Tables 5–6; common maximum pseudo-likelihood fitting for Gaussian, Student-t, Gumbel, Clayton, Frank and the Gumbel-Clayton-Frank mixture; AIC/BIC, temporal holdout, Cramer-von Mises-type `S_n`, and parametric bootstrap.
- `04_return_periods.py`: Section 2.8 and Tables 7–9 / Supplementary Table S5; monthly sampling interval `mu_I = 1/12 year`, OR and AND definitions, fixed equal-univariate thresholds, and balanced-probability physical-coordinate sensitivity.
- `05_seasonality.py`: Supplementary Tables S7–S8; cool-season (Nov–Mar) versus warm-season (Apr–Oct) diagnostics and the two-regime annual exceedance-rate sensitivity.
- `06_make_key_figures.py`: reproducible versions of the marginal-PDF and Kernel-Pareto threshold-diagnostic figures.

## 5. Important implementation notes

### Monthly return period
The monthly extreme-event sequence uses `mu_I = 1/12 year`. For a threshold pair mapped to copula coordinates `(u,v)`, the code uses

- OR exceedance probability: `1 - C(u,v)`
- AND exceedance probability: `1 - u - v + C(u,v)`
- return period: `T = mu_I / p_exceed`

For Tables 7 and 8, the two variables are first assigned their **separate equal univariate design levels** (e.g., both are 10-year univariate values). That same physical threshold pair is then evaluated under each competing copula.

### Mixed copula
All models are fitted by the same maximum pseudo-likelihood principle. The mixed density is

`a*c_Gumbel + b*c_Clayton + c*c_Frank`, with non-negative weights summing to one.

Finite-mixture likelihoods are non-convex and may have local/boundary optima. The code therefore uses multiple starts. `paper_parameters.json` contains the parameter snapshots reported in the revised manuscript and is used by the fixed-threshold return-period script for exact table reconstruction. Re-running the optimizer from different software versions or starts can give slightly different mixture decompositions while producing almost indistinguishable fitted densities; this is why the manuscript treats the mixture as a sensitivity model rather than an intrinsically superior model.

### Kernel-Pareto
The revised threshold is the 95th percentile for all variables. The body is represented by a Gaussian KDE and the exceedances by a GPD. Threshold diagnostics are evaluated from the 90th to 99.5th percentiles.

### Copula-model evaluation
The primary code reports log pseudo-likelihood, AIC/BIC and temporal holdout performance. For monthly extreme pairs, it additionally computes `S_n` and a 100-replicate parametric bootstrap p-value. Random seeds are fixed in `analysis_config.json`.

## 6. Reproducibility scope

The package contains the main statistical calculation code requested during peer review. COAWST/WRF/ROMS/WW3 model integration source code and model namelists are not included here; those belong to the numerical hindcast production workflow rather than the post-processing/statistical analysis reported in this repository.

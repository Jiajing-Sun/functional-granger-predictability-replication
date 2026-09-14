# Relevant Functional Granger Predictability

R code for **Self-Normalized Inference for Relevant Functional Granger Predictability**, by Jiajing Sun, Abderrahim Taamouti, and Zhuo Lin.

Updated 14 September 2026. The empirical analysis uses calendar 2025, constant-30-day Bitcoin option curves, and next-calendar-day realized volatility. Reported relevance thresholds are 0.01, 0.025, and 0.05. The two smaller thresholds are exploratory sensitivity benchmarks.

## Requirements

Tested with R 4.5.1 and `data.table` 1.17.8. Simulation scripts use `parallel`; the diagnostic scripts use forked workers and should run on macOS or Linux.

```r
install.packages("data.table", repos = "https://cloud.r-project.org")
```

## Files

| Script | Purpose |
|---|---|
| `R/FGC_simulation_functions_lambda20_plugin.R` | FARX simulation, population moments, relevant tests, and the separate approximate point-null benchmark. |
| `R/FGC_cv_lambda20_fixed.R` | Critical values for the fixed 19-point grid, `0.25 + 0.75 * (1:19)/20`. |
| `R/run_simulation_analysis.R` | Main Monte Carlo study and simulation figures. |
| `R/build_smiles.R` | Constant-30-day IV and trading-intensity curves and construction variants. |
| `R/run_calendar_analysis.R` | Full-sample and rolling tests, forecast comparisons, and sampling sensitivity. |
| `R/run_finite_sample_diagnostics.R` | Dimension, standardization, and diagonal-correction diagnostics. |
| `R/run_matched_diagnostics.R` | Boundary diagnostics with HAR controls, reverse designs, and calendar masking. |

## Simulations

Run from the repository root:

```bash
SN_MAIN_R=1000 SN_NCORES=4 Rscript R/run_simulation_analysis.R
```

The main study has 72 configurations and 1,000 replications per configuration. Seeds are fixed in the script. Results and figures are written to `processed/simulation/`. To check execution with fewer replications:

```bash
SN_MAIN_R=2 SN_NCORES=2 SN_OUTPUT_DIR=/tmp/fgc_smoke \
  Rscript R/run_simulation_analysis.R
```

The 5% critical values are 8.694475491210449 for the quadratic normalizer and 3.1993661030719607 for the adjusted-range normalizer. The point-null benchmark has a separate approximate calibration.

## Empirical inputs

This repository contains R source and this README. It does not include data, manuscript files, or generated results. Raw-data acquisition and figure/table formatting are outside this R release.

Place these prepared inputs in `processed/empirical/`:

| File | Required columns |
|---|---|
| `option_day_expiry_bins.csv` | `date`, `tau`, `u_bin`, `iv_dec`, `usd_notional`, `n_trades`, `min_u`, `max_u` |
| `spot_rv_alternative_sampling.csv` | `date`; `rv_1min`, `rv_5min`, `rv_10min`, `rv_15min`, `rv_5min_subsampled`; the corresponding `spot_vol_*` columns |

Dates are UTC, from 2025-01-01 through 2025-12-31. Option bins use Deribit BTC option trades: `tau` is time to expiry in days measured from the observation day's end; `u_bin` is log strike/index-price rounded to 0.001; `iv_dec` is BTC-amount-weighted implied volatility in decimal units; and `usd_notional` sums BTC amount times the contemporaneous USD index price. `min_u` and `max_u` record observed support within each bin. Keep maturities of 7–90 days and log-moneyness within [-0.35, 0.35].

Realized variance is constructed from Coinbase Exchange BTC-USD minute candles; `spot_vol_*` is the square root of the corresponding `rv_*` measure. The five-minute measure uses 288 intraday returns, taking the close of the preceding completed minute at each boundary. A previous close may be carried for at most five missing minutes. Longer gaps invalidate the daily measure. Keep missing dates on the calendar and encode missing values as `NA`.

Run from the repository root:

```bash
Rscript R/build_smiles.R
Rscript R/run_calendar_analysis.R
```

The first script creates the main curves and all four construction variants. The second writes full-sample results, 180-calendar-day rolling results with a five-day step, threshold-gain calculations, forecast comparisons, and sensitivity results to `processed/empirical/`.

The main specification requires supported expiries on both sides of 30 days, uses 21 log-moneyness coordinates from -0.25 to 0.25, and conditions forward tests on 1/5/22-calendar-day HAR history. It applies the diagonal correction and retains every fixed self-normalization prefix. It uses no nearest-maturity fallback. Input preparation must preserve the daily coverage and units used in the paper to reproduce its numerical results.

## Additional diagnostics

The diagnostic scripts require `processed/simulation/critical_values.csv` from the study's calibration outputs, with columns `trim`, `method`, `alpha`, and `quantile`. It must include `trim` values 0, 0.25, and 0.5; methods `quad` and `range`; and `alpha=0.05`. The main simulation writes critical values for trim 0.25 only and does not generate this multi-grid input.

From the repository root:

```bash
cd R
Rscript run_finite_sample_diagnostics.R
SN_DIAG_REVERSE_ONLY=1 Rscript run_finite_sample_diagnostics.R
Rscript run_matched_diagnostics.R
```

These scripts use 1,000 replications per configuration and write to `processed/diagnostics/`.

# ACTL90007 Assignment 2 — Question 3

## Run in RStudio

1. Extract the complete ZIP and keep the folder structure unchanged.
2. Open `Q3_analysis.R` in RStudio.
3. Use Session > Set Working Directory > To Source File Location.
4. Install packages once: `install.packages(c("StMoMo", "reshape2", "forecast"))`.
5. Click **Source** to execute the script. Figures and result tables are written to `output/`.

Alternatively, from a terminal in this folder: `Rscript Q3_analysis.R`.

## Data and method

- The four original HMD text files are in `data/`, renamed with AUS/UK prefixes to avoid duplicate filenames.
- The files have two header lines before the column names; the script uses `skip=2`.
- Both sexes (`Total`), ages **40–90 inclusive**, years **1961–2005 inclusive**.
- Each population has a complete 51-by-45 grid (2,295 cells), with positive deaths/exposures and no missing selected observations. No ages are silently dropped.
- The original HMD files match the previously uploaded Excel data exactly over the selected window.
- LC: `StMoMo::lc(link="log")`, fitted to deaths and central exposures (Poisson log link).
- The UK's stronger cohort effect is identified from the diagonal residual patterns, not from comparing raw deviance totals between countries.
- RH: lecturer's `rh(link="log", cohortAgeFun="NP", approxConst=TRUE)`, fitted to UK data.
- All residual figures use **unscaled Poisson deviance residuals**, `residuals(fit, scale=FALSE)`, with matching colour limits within each before/after pair. Residual signs are plotted separately to avoid confusing larger UK exposures with stronger cohort structure.
- Blue: observed deaths below fitted deaths; red: above fitted. With year horizontal and age vertical, constant-birth-cohort paths slope upwards to the right.
- No mortality forecasting, smoothing-parameter search or replacement SVD estimator is used; these are not requested in Q3.

## Classroom example and changes

`LC_APC_teacher_original.R` is preserved unchanged for reference; run `Q3_analysis.R`, not the original demonstration.
The classroom example uses Male, ages 20–89 and years 1981–2021; the assignment version changes these to Total, 40–90 and 1961–2005. It also analyses both countries and adds the required residual comparisons. The original file references an undefined `RH.simple.fit` near the end; the assignment script consistently uses `RH.fit`.

## Outputs

- `Q3_Report.docx` / `Q3_Report.pdf`: written answer with figures and conclusions.
- `output/Q3a_LC_parameters.*`: both countries' ax, bx and kt.
- `output/Q3a_LC_residual_heatmaps.*` / `...signs.*`: Australia versus UK.
- `output/Q3b_UK_LC_RH_residual_heatmaps.*` / `...signs.*`: UK before/after cohort term.
- `output/Q3b_UK_RH_cohort_parameters.*`: cohort sensitivity and factor.
- `output/*parameters.csv`, `*fitted_rates.csv`, `*deviance_residuals.csv`: numeric results.
- `output/model_diagnostics.csv`, `model_summaries.txt`, `sessionInfo.txt`: convergence and reproducibility information.
- `.rds` files preserve the fitted models and can be read using `readRDS()`.

The deviance reduction is descriptive in-sample evidence, not a formal chi-square likelihood-ratio test or proof of improved future forecasts. Residual structure remains and should not be described as entirely eliminated.

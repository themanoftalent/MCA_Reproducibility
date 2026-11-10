# Technical Documentation: Grouped Vehicle Emission Regression

This technical document provides detailed methodological and implementation documentation for the computational analysis accompanying:

> **Article:** Weighting and Mileage Coding in Grouped Vehicle Emission Regression  
> **Journal:** *Mathematical and Computational Applications* (MDPI)  
> **Manuscript ID:** `mca-4666564`

---

## 1. Regression Model Specifications

The analytical models fit grouped vehicle emission observations across 18 age–mileage cells for three regulated pollutants: carbon monoxide (CO), non-methane hydrocarbons (NMHC), and nitrogen oxides ($\text{NO}_x$). Predictions are evaluated at the designated reference condition, cell $g_{17}$ (vehicle age = 5 years, mileage category = 2).

The following specifications are formulated:

- **Model A (Age-Only):**
  $$y = \beta_0 + \beta_1 \cdot \text{age}$$
  A two-parameter baseline linear model assuming emission rates vary solely as a function of vehicle age.
  
- **Model AM (Age and Continuous/Ordinal Mileage with Interaction):**
  $$y = \beta_0 + \beta_1 \cdot \text{age} + \beta_2 \cdot m + \beta_3 \cdot (\text{age} \cdot m)$$
  A four-parameter specification where mileage category $m \in \{1, 2, 3\}$ is treated as an ordinal continuous covariate with an interaction term between age and mileage.

- **Model CAT (Categorical Mileage with Full Interaction):**
  $$y = \beta_0 + \beta_1 \cdot \text{age} + \beta_2 \cdot \mathbb{I}(m=2) + \beta_3 \cdot \mathbb{I}(m=3) + \beta_4 \cdot (\text{age} \cdot \mathbb{I}(m=2)) + \beta_5 \cdot (\text{age} \cdot \mathbb{I}(m=3))$$
  A six-parameter specification allowing distinct intercepts and distinct age slopes for each of the three mileage categories.

- **Model ACAT (Additive Categorical Mileage):**
  In comparative model discussions, an additive categorical formulation ($y = \beta_0 + \beta_1 \cdot \text{age} + \beta_2 \cdot \mathbb{I}(m=2) + \beta_3 \cdot \mathbb{I}(m=3)$) assumes common age slope across mileage classes without interaction. In the present computational codebase, the primary models implemented and evaluated are `A`, `AM`, and the fully interacted `CAT` model.

Design matrices are generated in `src/analyze.py` via `design(a, m, kind)`.

---

## 2. Weighting Methods

Four weighting schemes are evaluated across all regression models:

1. **Equal Weighting (`equal`):**
   $$w_i = 1$$
   Ordinary least squares on cell means, treating all cells with equal weight regardless of sample size or dispersion.

2. **Inverse Standard Deviation (`inverse_sd`):**
   $$w_i = \frac{1}{s_i}$$
   Weights inversely proportional to the sample standard deviation $s_i$ of the cell.

3. **Inverse Variance (`inverse_variance`):**
   $$w_i = \frac{1}{s_i^2}$$
   Weights inversely proportional to cell sample variance.

4. **Count-Informed Inverse Variance (`count_informed`):**
   $$w_i = \frac{n_i}{s_i^2}$$
   Weights proportional to cell count $n_i$ and inversely proportional to cell variance, corresponding to minimum-variance weighting under classical independent sampling assumptions.

---

## 3. Covariance Estimation and Standard Errors

Uncertainty around coefficient estimates $\hat{\beta}$ and reference cell predictions $\hat{y}_{\text{ref}} = c^T \hat{\beta}$ is quantified through three covariance formulations:

1. **Working Covariance:**
   $$\operatorname{Cov}_{\text{working}}(\hat{\beta}) = s_w^2 (X^T W X)^{-1}, \quad s_w^2 = \frac{\sum_{i=1}^N w_i e_i^2}{N - p}$$
   where $e_i = y_i - x_i^T \hat{\beta}$, $N = 18$ cells, and $p = \operatorname{rank}(X)$.

2. **Heteroskedasticity-Consistent HC3 Covariance:**
   $$\operatorname{Cov}_{\text{HC3}}(\hat{\beta}) = (X^T W X)^{-1} \left( \sum_{i=1}^N \left(\frac{w_i e_i}{1 - h_{ii}}\right)^2 x_i x_i^T \right) (X^T W X)^{-1}$$
   where $h_{ii}$ represents the $i$-th diagonal element of the weighted hat matrix $H = W^{1/2} X (X^T W X)^{-1} X^T W^{1/2}$.

3. **Count-Informed Plug-in Covariance:**
   $$\operatorname{Cov}_{\text{plugin}}(\hat{\beta}) = (X^T W X)^{-1}$$
   Evaluated specifically under count-informed weights $w_i = n_i / s_i^2$, reflecting the theoretical sampling covariance when within-cell variances are known.

Critical values for working and HC3 intervals employ Student's $t$-distribution with $18 - p$ degrees of freedom. Plug-in intervals use the standard normal critical value ($z_{0.975} = 1.959964$).

---

## 4. Cross-Validation and Influence Diagnostics

- **Leave-One-Cell-Out (LOCO) Validation:**
  Evaluates out-of-sample prediction error by omitting each of the 18 cells in turn. The verification suite confirms exact numerical equivalence with the closed-form PRESS relationship:
  $$y_g - \hat{y}_{(g)} = \frac{e_g}{1 - h_{gg}}$$
  LOCO mean squared error is recorded across all model and weighting combinations.

- **Age-Level Validation:**
  Partitions the 18 cells into 6 folds corresponding to the 6 distinct age levels, omitting all three mileage categories at a given age simultaneously. For the `CAT` model, age-level validation errors coincide with cell-level errors because separate regression lines are fitted independently for each mileage class. This identity is an algebraic property of block-diagonal design structures.

- **Influence Diagnostics:**
  Maximum single-cell deletion changes on reference predictions, leverage values ($h_{ii}$), and Cook's distance statistics are computed and tabulated.

---

## 5. Robust Regression and Weight Sensitivity

- **Fixed-Scale Huber M-Estimation:**
  Robust regression estimates are computed using a fixed scale with Huber tuning parameter $c = 1.345$. The optimization is independently validated in `src/verify.py` by verifying first-order gradient conditions ($\|\nabla \mathcal{L}\|_\infty < 10^{-6}$) and confirming objective agreement against BFGS optimization.

- **Bounded Weight Sensitivity Analysis:**
  To evaluate the stability of reference cell predictions against arbitrary reweighting, bounded sensitivity intervals are computed via Sequential Least Squares Programming (SLSQP). For each pollutant and direction (minimization and maximization), eight random starts are executed. Attained extrema and solver convergence diagnostics are archived in `results/table06_sensitivity.csv` and `results/sensitivity_solver_runs.csv`.

---

## 6. Monte Carlo Simulation Studies

- **Baseline Simulation Study:**
  - 6 data-generating scenarios: equal counts & equal SD, unequal counts, unequal SD, counts and SD combined, contaminated cell, and categorical misspecification.
  - 5 estimators evaluated per scenario over 1,000 replications (30,000 total fits).
  - Target: uncontaminated latent mean at reference cell $g_{17}$.
  - Random seed: `20261017` (NumPy generator sequence preserving the original simulation draw order).

- **Extension Simulation Study:**
  - 11 scenarios evaluating sensitivity to non-normality (standardized $t_5$ distributions), structural departures, and variance uncertainty under both known and estimated dispersion.
  - 2 regression models (`AM` and `CAT`) over 1,000 replications (44,000 total fits).
  - Random seed: `20261030`.
  - In each replication, 18 cell-mean errors are drawn first, followed by 18 independent chi-square dispersion draws with $n_i - 1$ degrees of freedom. Pairing heavy-tailed cell means with independent chi-square variance estimates serves as a deliberate stress test.

---

## 7. Required Input Files and Directory Structure

The computational code does not redistribute original vehicle data or supplementary spreadsheets. To execute `run_pipeline.py`, the following input files must be provided in a `data/` subdirectory:

```
data/
├── grouped_inputs.csv          # 18-cell summary data (age, mileage, count, mean, SD)
├── resampled_assignments.csv   # Cell assignment frequencies
├── reference_table_01.csv      # Manuscript comparison reference for Table 1
├── reference_table_02.csv      # Manuscript comparison reference for Table 2
├── ...
└── reference_table_11.csv      # Manuscript comparison reference for Table 11
```

The study inputs are available from Table 1 of the article and the journal supplementary materials (see [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md)).

---

## 8. Verification and Computational Integrity

The repository includes a comprehensive verification suite:

1. **Algebraic and Numerical Checks (`validation/verification.json`):**
   - 83 formal checks covering normal equations, scale invariance, leverage traces ($\sum h_{ii} = p$), covariance matrix symmetry, closed-form PRESS identities, Huber stationarity, and simulation consistency.
   - All 83 checks pass.

2. **Manuscript Value Comparison (`results/manuscript_comparison.csv`):**
   - 622 distinct values compared against the manuscript tables.
   - All 622 values match within rounding tolerance.

3. **Manifest Verification (`verify_integrity.py`):**
   - Confirms SHA-256 checksums across all tracked repository files against `SHA256_MANIFEST.json`.

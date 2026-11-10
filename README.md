# Weighting and Mileage Coding in Grouped Vehicle Emission Regression

This repository contains the computational code, results, figures, and verification records accompanying the research article:

> **Article:** Weighting and Mileage Coding in Grouped Vehicle Emission Regression  
> **Journal:** *Mathematical and Computational Applications* (MDPI)  
> **Manuscript ID:** `mca-4666564`

---

## Study Overview and Computational Methods

The study evaluates how alternative weighting specifications, mileage coding conventions, and standard error estimators influence estimated emission levels at a reference operating point and affect out-of-sample prediction performance. The empirical analysis utilizes grouped vehicle emission data across 18 age–mileage cells for carbon monoxide (CO), non-methane hydrocarbons (NMHC), and nitrogen oxides ($\text{NO}_x$), and is supported by Monte Carlo simulation experiments.

The computational procedures implement:

1. **Regression Model Specifications:**
   - **Model A:** Linear regression with vehicle age ($y = \beta_0 + \beta_1 \cdot \text{age}$).
   - **Model AM:** Regression with vehicle age, ordinal mileage class ($m \in \{1, 2, 3\}$), and an interaction term ($y = \beta_0 + \beta_1 \cdot \text{age} + \beta_2 \cdot m + \beta_3 \cdot (\text{age} \cdot m)$).
   - **Model CAT:** Categorical mileage specification with class-specific intercepts and slopes ($y = \beta_0 + \beta_1 \cdot \text{age} + \beta_2 \cdot \mathbb{I}(m=2) + \beta_3 \cdot \mathbb{I}(m=3) + \beta_4 \cdot (\text{age} \cdot \mathbb{I}(m=2)) + \beta_5 \cdot (\text{age} \cdot \mathbb{I}(m=3))$).
   - Reference predictions evaluate cell $g_{17}$ (age 5 years, mileage category 2).

2. **Weighting Schemes:**
   - Equal weights ($w_i = 1$)
   - Inverse standard deviation ($w_i = 1/s_i$)
   - Inverse variance ($w_i = 1/s_i^2$)
   - Count-informed inverse variance ($w_i = n_i/s_i^2$)

3. **Uncertainty and Covariance Estimation:**
   - Standard working covariance
   - Heteroskedasticity-consistent HC3 covariance
   - Count-informed plug-in covariance

4. **Validation and Diagnostics:**
   - Leave-one-cell-out (LOCO) error evaluation verified against closed-form PRESS residuals
   - Leave-one-age-out validation folds (omitting all mileage classes at each age)
   - Cook's distance, leverage, and single-cell omission shifts

5. **Robust Estimation and Sensitivity Analysis:**
   - Fixed-scale Huber M-estimation ($c = 1.345$)
   - Multistart SLSQP constrained optimization for bounded weight perturbation sensitivity

6. **Monte Carlo Simulations:**
   - Baseline study: 6 data-generating scenarios, 5 estimation methods, 1,000 replications (30,000 fits)
   - Extended study: 11 scenarios under known and estimated dispersion, 2 model specifications (44,000 fits)

---

## Repository Structure

```
├── DATA_AVAILABILITY.md        # Data availability and source attribution statement
├── README.md                   # Repository overview and execution guidance
├── README_TECHNICAL.md         # Detailed technical notes and methodology
├── requirements.txt            # Python dependencies
├── run_pipeline.py             # Main pipeline orchestration script
├── verify_integrity.py         # SHA-256 hash verification script
├── figures/                    # High-resolution figures (PNG and vector PDF)
│   ├── figure1.{pdf,png}
│   ├── figure2.{pdf,png}
│   └── figure3.{pdf,png}
├── results/                    # Full numerical results, tables, and diagnostics
│   ├── coefficients.csv
│   ├── table01_*.csv ... table11_*.csv
│   ├── manuscript_comparison.csv
│   └── ...
├── src/                        # Core Python analytical and simulation modules
│   ├── analyze.py
│   ├── extension.py
│   ├── figures.py
│   ├── report.py
│   ├── tables.py
│   └── verify.py
└── validation/                 # Computational checks and execution logs
    ├── execution.json
    ├── run.log
    └── verification.json
```

---

## Data Availability

The numerical inputs used in this study are presented in Table 1 of the article. The supplementary Excel workbooks, `samplecounts.xlsx` and `18-cell_data.xlsx`, are provided separately through the journal as supplementary materials. 

This repository contains the computational code, results, figures, and verification materials. The European Commission Joint Research Centre Euro 7 dataset is available from its original source at [https://doi.org/10.2905/JRC.8897YXR](https://doi.org/10.2905/JRC.8897YXR).

For further details, refer to [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

---

## Environment and Installation

The code is developed for Python 3.12 (and compatible versions). To configure a dedicated virtual environment:

```bash
# Create a virtual environment
python3 -m venv .venv

# Activate the virtual environment
# On macOS / Linux:
source .venv/bin/activate
# On Windows:
# .venv\Scripts\Activate.ps1

# Install dependencies
python -m pip install -r requirements.txt
```

---

## Verification and Execution

### 1. Verifying Archived File Integrity

To verify the integrity of all repository files against the SHA-256 manifest:

```bash
python verify_integrity.py
```

### 2. Inspecting Computational Results

All numerical outputs, tables, and verification records are tracked directly within `results/` and `validation/`. Specifically:
- `validation/verification.json` documents 83 formal algebraic and computational checks (all passing).
- `results/manuscript_comparison.csv` documents 622 comparisons against the manuscript numerical values (all matching within rounding tolerance).

### 3. Re-executing the Computational Pipeline

The pipeline is designed to be re-executable when supplementary input files are placed in the expected `data/` directory:
- `data/grouped_inputs.csv` (18 cell summaries)
- `data/resampled_assignments.csv` (resampled assignment records)
- `data/reference_table_01.csv` through `data/reference_table_11.csv` (manuscript reference tables)

When the inputs are in place, the complete analysis and simulation suite can be rerun with:

```bash
python run_pipeline.py
```

A quick diagnostic run with reduced simulation replications is also supported:

```bash
python run_pipeline.py --reps 50
```

If run without the required supplementary files, `run_pipeline.py` prints an informative message detailing the necessary inputs and their expected locations.

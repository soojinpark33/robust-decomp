# When Do Traditional and Causal Decomposition Methods Diverge? A Practical Guide for Studying Health Disparities
R codes for "When Do Traditional and Causal Decomposition Methods Diverge? A Practical Guide for Studying Health Disparities"
Soojin Park<sup>1</sup>, Su Yeon Kim<sup>1</sup>, and Chioun Lee<sup>2</sup>

<sup>1</sup> School of Education, University of California, Riverside  
<sup>2</sup> Department of Sociology, University of California, Riverside


## Overview

A central objective among researchers across disciplines is to identify malleable factors that can reduce social disparities. Traditionally, researchers have relied on the difference-in-coefficients and Kitagawa-Oaxaca-Blinder frameworks. More recently, methods grounded in the potential outcomes framework, such as causal decomposition analysis or causal mediation analysis, have emerged. While these methods share the same goal of identifying drivers of disparity, they frequently yield divergent results depending on the underlying confounding structures and settings. Despite these significant differences, applied researchers lack clear guidance on selecting appropriate methods for their specific research contexts. To address this gap, this study provides a systemic review and offer an intuitive guidance through Directed Acyclic Graphs and comparative simulation studies. We begin by reviewing each method assuming no unmeasured confounding, which is often violated in observational settings. Consequently, we extend our analysis to two realistic scenarios: 1) unmeasured confounding exists in the relationship between intermediate confounders and the mediator, and 2) unmeasured confounding exists in the relationship between the mediator and the outcome. Finally, we illustrate these recommendations through a case study examining the role of educational attainment in explaining racial disparities in later-life cognition.

For more details of our proposed methods, see [our paper](https://www.degruyter.com/document/doi/10.1515/jci-2022-0031/html). 
Here, we provide `R` code to (1) apply the decomposition methods to your own data and (2) reproduce our simulation study.

## Repository Contents

| File | Purpose |
|---|---|
| `causal_decomposition_analysis.R` | **Start here to analyze your own data.** Compares decomposition methods and runs a sensitivity analysis. |
| `synthetic_data.dta` | Synthetic example data for trying out `causal_decomposition_analysis.R` |
| `Simulation study.R` | Reproduces the simulation results in the paper (for replication only) |
| `Comp_source.R` | Helper functions used only by `Simulation study.R` |

> **Applying the methods to your own data?** Use `causal_decomposition_analysis.R`. The simulation files are written for the simulated data in the paper and are not intended for applied analyses.

---

# Using the Methods on Your Own Data

`causal_decomposition_analysis.R` estimates how much of a group disparity in an outcome would be reduced by equalizing a mediator, using:

1. **Difference in Coefficients (DC)**
2. **Kitagawa-Oaxaca-Blinder (KOB)**, with each group in turn as the reference
3. **Causal Decomposition Analysis (CDA)** via sequential mediation imputation (`smi()` from the `causal.decomp` package)

It then runs a **sensitivity analysis** for unmeasured mediator–outcome confounding, using Cinelli and Hazlett's (2020) covariate benchmarking (Equations 54 and 56 of the supplementary material).

## Step 1. Install packages

```r
install.packages(c("haven", "causal.decomp", "ggplot2", "ggrepel"))
```

## Step 2. Try it on the synthetic data

Set your working directory to the folder where you cloned this repository, then run:

```r
source("causal_decomposition_analysis.R")
```

With the default settings, this analyzes `synthetic_data.dta`. Running it once first is a good way to see what the output looks like.

## Step 3. Prepare your data

Your data need one row per person and the following variables:

| Role | Description | Format |
|---|---|---|
| Group | The two groups being compared (e.g., Black vs. non-Black) | 0/1. The script converts it to a factor; **1 is the comparison group**, 0 the reference group. |
| Outcome | The outcome with the disparity | numeric |
| Mediator | The factor you would intervene on to reduce the disparity | numeric |
| Baseline covariate | **One** pre-group covariate such as age. The script centers it and adds a squared term. | numeric |
| Intermediate confounders | Variables affected by group membership that also affect the mediator and outcome (e.g., childhood SES) | numeric (dummy-code categorical variables) |

Remove or impute missing values before running. The script warns you if it finds any.

## Step 4. Edit the configuration section

Open `causal_decomposition_analysis.R` and change the **CONFIGURATION SECTION** at the top. You don't need to edit anything below it.

```r
data_file         <- "my_data.dta"                         # your data file
group             <- "black"                               # group variable (0/1)
mediator_var      <- "educy"                               # mediator
outcome_var       <- "cog27"                               # outcome
base_cov          <- "age"                                 # one baseline covariate
intermediate_covs <- c("chdSES", "RTHLTHCH", "Pdivorce16") # intermediate confounders

benchmark_covariates <- c("age_centered", "chdSES", "RTHLTHCH", "Pdivorce16")

n_boot   <- 1000   # bootstrap iterations for DC and KOB intervals
smi_sims <- 1000   # simulations for CDA
max_rsq  <- 0.3    # largest partial R-squared shown in the sensitivity plot
```

Notes:

* **`benchmark_covariates`**: the observed covariates used as reference points in the sensitivity plot. The centered baseline covariate is always called **`age_centered`**, whatever `base_cov` is named, so list it as `"age_centered"`.
* **Not using Stata?** Replace `read_dta(data_file)` in the data-loading section with, for example, `read.csv(data_file)`.
* For a quick test run, lower `n_boot` and `smi_sims` (e.g., to 100).

## Step 5. Run and read the results

```r
source("causal_decomposition_analysis.R")
```

The script prints a **comparison table** (also stored as `comparison_table`) with these rows:

| Row | Meaning |
|---|---|
| Total Disparity | Initial difference in the outcome between groups |
| Explained Component (Disparity Reduction) | Portion of the disparity that would be removed by equalizing the mediator |
| Unexplained Component (Disparity Remaining) | Portion that would remain |
| % Explained | Disparity reduction ÷ total disparity × 100 |

Each method has an estimate and a 95% confidence interval:

| Columns | Method |
|---|---|
| `DC_*` | Difference in Coefficients (bootstrap CI) |
| `KOB_NB_*` | KOB using the reference group's (group = 0) coefficients |
| `KOB_B_*` | KOB using the comparison group's (group = 1) coefficients |
| `SMI_*` | Causal decomposition analysis |

It also saves a **sensitivity contour plot** (`sensitivity_reduction.png`). The plot shows how the disparity reduction would change under unmeasured confounding of a given strength, with your benchmark covariates marked for comparison. Other useful objects left in your R session:

* `res.1a1`: full `smi()` result
* `sens_results`: full sensitivity analysis result

## Models used

For transparency, these are the models the script fits, using the names in your configuration:

* **DC:** `outcome ~ group + base + base² + intermediates`, compared with the same model plus `mediator`
* **KOB:** `outcome ~ base + base² + intermediates + mediator`, fitted separately within each group
* **CDA:** mediator model `mediator ~ group + base + base²`, and outcome model `outcome ~ group * mediator + base + base² + intermediates`

See [our paper](https://www.degruyter.com/document/doi/10.1515/jci-2022-0031/html) for when each method gives a valid estimate.

---

# Reproducing the Simulation Study

`Simulation study.R` (with `Comp_source.R`) reproduces the simulation tables in the paper. It generates six data-generating scenarios and compares DIC, KOB, CDA, and the modified DIC/KOB estimators against the true values.

Before running:

1. Change the `source(...)` line near the top to `source("Comp_source.R")`.
2. Change `out_dir` to a folder on your computer (e.g., `out_dir <- "sim_results"`).
3. The full run uses a population of 1,000,000 and `n_iter <- 500`, which takes a long time. For a quick test, set `n_iter` to a small number such as 10.

Additional packages needed: `install.packages(c("dplyr", "oaxaca", "truncnorm", "causal.decomp"))`.

---

These supplementary materials are provided solely for the purpose of reproducibility and must be used in compliance with academic ethical guidelines. If you reference these materials in your own work, please ensure proper citation of the original sources.

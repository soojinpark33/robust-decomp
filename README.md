# A Systematic Comparison of Traditional and Causal Decomposition Methods: Using Directed Acyclic Graphs and Simulation Studies
R codes for "A Systematic Comparison of Traditional and Causal Decomposition Methods: Using Directed Acyclic Graphs and Simulation Studies"
Soojin Park<sup>1</sup>, Su Yeon Kang<sup>1</sup>, and Chioun Lee<sup>2</sup>

<sup>1</sup> School of Education, University of California, Riverside  
<sup>2</sup> Department of Sociology, University of California, Riverside


## Overview

A central objective among researchers across disciplines is to identify malleable factors that can reduce social disparities. Traditionally, researchers have relied on the difference-in-coefficients and Kitagawa-Oaxaca-Blinder frameworks. More recently, methods grounded in the potential outcomes framework, such as causal decomposition analysis or causal mediation analysis, have emerged. While these methods share the same goal of identifying drivers of disparity, they frequently yield divergent results depending on the underlying confounding structures and settings. Despite these significant differences, applied researchers lack clear guidance on selecting appropriate methods for their specific research contexts. To address this gap, this study provides a systemic review and offer an intuitive guidance through Directed Acyclic Graphs and comparative simulation studies. We begin by reviewing each method assuming no unmeasured confounding, which is often violated in observational settings. Consequently, we extend our analysis to two realistic scenarios: 1) unmeasured confounding exists in the relationship between intermediate confounders and the mediator, and 2) unmeasured confounding exists in the relationship between the mediator and the outcome. Finally, we illustrate these recommendations through a case study examining the role of educational attainment in explaining racial disparities in later-life cognition.

For more details of our proposed methods, see [our paper](https://www.degruyter.com/document/doi/10.1515/jci-2022-0031/html). 
Here, we provide `R` codes to reproduce our simulation study. 


## Simulation Study

* `Simulation_study.R`  

   This `R` file contains the simulation codes for comparing three methods. This code replicates our results in Tables 1 and 2 of our paper.

* `Comp_source.R` 
 
   This `R` file includes source functions required to run our simulation codes. 

These supplementary materials are provided solely for the purpose of reproducibility and must be used in compliance with academic ethical guidelines. If you reference these materials in your own work, please ensure proper citation of the original sources.

## Case Study

* `synthetic_data.dta`
  
  Synthetic data created to run the code.

* `causal_decomposition_analysis.R`

  This code implements:
1. **Multiple Decomposition Methods**:
   - Difference in Coefficients (DC)
   - Oaxaca-Blinder Decomposition (KOB) with two reference groups
   - Causal Decomposition using Sequential Mediation Imputation (SMI)

2. **Sensitivity Analysis**:
   - Implements Cinelli and Hazlett's (2020) benchmarking procedure
   - Uses equations 54 and 56 from the supplementary material
   - Creates contour plots showing sensitivity to unmeasured confounding
     2. **Prepare your data**: If using your own data, ensure it's in Stata format (.dta) or modify the data loading section.

## Requirements

### R Packages
```r
install.packages(c("haven", "causal.decomp", "ggplot2", "ggrepel"))
```

### Required Packages:
- `haven`: For reading Stata files (.dta)
- `causal.decomp`: For causal decomposition analysis
- `ggplot2`: For plotting
- `ggrepel`: For label positioning in plots
  
## Usage
1. **Configure variables**: Edit the configuration section at the top of `causal_decomposition_analysis.R` to match your data:
   ```r
   # Configuration
   data_file <- "synthetic_data.dta"  # or your data file
   treatment_var <- "black"      # Treatment/exposure variable
   mediator_var <- "educy"        # Mediator variable
   outcome_var <- "cog27"         # Outcome variable
   covariates <- c("chdSES", "RTHLTHCH", "Pdivorce16")
   ```

2. **Run the analysis**:
   ```r
   source("Case Study/causal_decomposition_analysis.R")
   ```

### Output

The script produces:
- **Comparison Table**: Estimates and 95% confidence intervals for all decomposition methods
- **Sensitivity Plots**: Contour plots showing sensitivity to unmeasured confounding with covariate benchmarks
- **Benchmark Summary**: Partial R² values for each benchmark covariate


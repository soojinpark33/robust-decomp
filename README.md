# Identifying Causally-Robust Mediators of Health Disparities: A Review and Simulation Studies With Directed Acyclic Graphs
R codes for "Identifying Causally-Robust Mediators of Health Disparities: \\A Review and Simulation Studies With Directed Acyclic Graphs"
Soojin Park<sup>1</sup>, Su Yeon Kang<sup>1</sup>, and Chioun Lee<sup>2</sup>

<sup>1</sup> School of Education, University of California, Riverside  
<sup>2</sup> Department of Sociology, University of California, Riverside


## Overview

Background. Traditionally, researchers have used linear approaches, such as difference-in-coefficients (DIC) and Kitagawa-Oaxaca-Blinder (KOB) decomposition, to identify risk factors or resources (referred to as `mediators') underlying health disparities. More recently, causal decomposition analysis (CDA) has gained popularity by defining clear causal effects of interest and estimating them without any modeling restrictions. 
Methods. We start with a brief review of each method, assuming no unmeasured confounders. We then move to two more realistic scenarios: 1) unmeasured confounders affect the relationship between intermediate confounders and the mediator, and 2) unmeasured confounders affect the relationship between the mediator and the outcome. For each scenario, we generate simulated data, apply the three methods, compare their estimates, and interpret the results using Directed Acyclic Graphs. 
Results. While the DIC approach is appropriate when no intermediate confounders are present, it is unrealistic to assume the absence of intermediate confounders as health disparities arise from multiple factors over the life-course. The KOB decomposition is appropriate when controlling for baseline covariates (such as age) is not required. When unmeasured confounding exists, 1) the DIC method yields biased estimates in both scenarios, 2) both KOB and CDA produce biased results in the second scenario; however, CDA accompanied with sensitivity analysis can help assess the robustness of those estimates. 
Conclusions. We advise against using the DIC method when investigating drivers of health disparities. We recommend CDA combined with sensitivity analysis as a robust strategy for identifying mediators of health disparities.

For more details of our proposed methods, see [our paper](https://www.degruyter.com/document/doi/10.1515/jci-2022-0031/html). 
Here, we provide `R` codes to reproduce our simulation study. 


## Simulation Study

* `Simulation_study.R`  

   This `R` file contains the simulation codes for comparing three methods. This code replicates our results in Tables 1 and 2 of our paper.

* `Comp_source.R` 
 
   This `R` file includes source functions required to run our simulation codes. 

These supplementary materials are provided solely for the purpose of reproducibility and must be used in compliance with academic ethical guidelines. If you reference these materials in your own work, please ensure proper citation of the original sources.

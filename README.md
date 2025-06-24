# 'Inference procedures in sequential trial emulation with survival outcomes: comparing confidence intervals based on the sandwich variance estimator, bootstrap and jackknife'
Juliette M. Limozin, Shaun R. Seaman, Li Su
https://arxiv.org/abs/2407.08317

This repository contains the code scripts for running the simulation study (Section 5) and HERS data analysis (Section 6).

## Dependency
**TrialEmulation** package https://cran.r-project.org/web/packages/TrialEmulation/index.html
## Code/
- **simulate_MSM_simplified.R**: data-generating algorithm function script
- **simulation_low_jackknife.R**, **simulation_med_jackknife.R**, **simulation_high_jackknife.R**: simulation study scripts for low, medium and high outcome event rate scenarios respectively
- **simu_*.slurm**: slurm scripts for HPC to run the simulation study
- **true_value_gen_newsimus.R**: script to generate the true values of the MRDs for each simulation scenario by generating data for a very large randomized controlled trial, as proposed by Keogh et al. (2023). The true marginal risks in trial 0 when all patients were always treated or all patients were never treated were approximated by Kaplan-Meier estimates from two extremely large datasets (n = 1, 000, 000). The results are saved in **true_value_surv1.rda** and **true_value_surv0.rda**.
- **weight_func_efficient.R**: script for weight estimation function that can recalculate new inverse probability weights based on a new set of treatment and censoring models' coefficients, and/or a bootstrap sample. Used for generating non-parametric and LEF bootstrap CIs.
- **CI_SEs_generator.R**: script to generate the standard error estimates from each CI method after simulation scripts are run. 
- **results_with_jackknife_iter_seed.R**: script to generate simulation results as seen in Section 5.2
- **HERS_modelling_fixed_boot.R**: script to generate HERS data analysis results as seen in Section 6
- **Archive/**: additional scrips and .rda objects used for the Supplementary Materials

## Results and images/
**Accepted plots/** contains all the images included in the main text and supplementary materials of the accepted manuscript.

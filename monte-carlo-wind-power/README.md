# Monte Carlo Methods for Wind Power Estimation

This project uses **Monte Carlo simulation and variance reduction techniques** to estimate electricity production from a Vestas V164 (9.5 MW) wind turbine. Wind speeds are modelled using Weibull distributions, with monthly variation in the distribution parameters. The analysis is then extended to two turbines exposed to **dependent wind speeds**.

The simulations were implemented in **MATLAB** as part of the course *Monte Carlo and Empirical Methods for Stochastic Inference* (MASM11/FMSN50) at Lund University.

## Objectives

- Estimate monthly expected turbine output and construct approximate **95% confidence intervals**.
- Compare standard Monte Carlo simulation with variance reduction methods.
- Estimate the probability of nonzero power production, the average power coefficient, and the capacity factor.
- Investigate combined power production, covariance, variance, and threshold probabilities for two turbines with dependent wind speeds.

## Methods

| Method | Application |
| --- | --- |
| Standard Monte Carlo | Baseline estimation of expected power output |
| Conditional (truncated) sampling | Restrict sampling to wind speeds at which the turbine generates power, with a probability adjustment |
| Control variates | Use wind speed, whose expectation is known, to reduce estimator variance |
| Importance sampling | Reweight samples from normal proposal distributions |
| Antithetic sampling | Pair draws from a conditional Weibull distribution to reduce variance |
| Bivariate Weibull model | Analyse dependent wind speeds and joint turbine production |

## Selected findings

The variance reduction methods produced narrower confidence intervals than the standard Monte Carlo baseline in the reported comparisons. Under the assumptions of the model, the report estimates an average **capacity factor of approximately 39.4%** and a **fraction of wind conditions with positive production of approximately 85.0%** (each calculated by averaging the monthly estimates equally).

These figures are model-based estimates, **not a complete technical or economic assessment** of a proposed wind farm. Detailed calculations, tables, and discussion are provided in the report.

## Files

```text
monte-carlo-wind-power/
├── README.md
├── proj1.m
├── powercurve_V164.mat
└── report/
    └── Monte_Carlo_Assignment_1.pdf
```

- **`proj1.m`** — MATLAB script containing the simulations, confidence interval calculations, and plots.
- **`powercurve_V164.mat`** — supplied turbine power-curve data/function used by the script.
- **`report/Monte_Carlo_Assignment_1.pdf`** — mathematical derivations, methods, numerical results, and discussion.

## Requirements and usage

You need **MATLAB** with the **Statistics and Machine Learning Toolbox** (for functions such as `wblrnd`, `wblpdf`, `wblcdf`, `wblinv`, `normpdf`, `mvnrnd`, and `mvnpdf`).

1. Download or clone the repository.
2. In MATLAB, set the current folder to `monte-carlo-wind-power`, where both `proj1.m` and `powercurve_V164.mat` are located.
3. Run:

   ```matlab
   proj1
   ```

The script prints estimates and confidence intervals to the Command Window and generates comparison figures. A fixed random seed is used for reproducibility. The PDF contains numerical results from a particular version of the analysis; keep its tables synchronized with the script when making changes.

## Academic context and attribution

The project originated as a university group assignment. The original report was co-authored by **Ruolin Wang and Leonard**, and the material has subsequently been revised for portfolio presentation. The wind-speed parameters and turbine power-curve data were supplied in the course assignment.

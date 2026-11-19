# Sequential Monte Carlo: Self-Avoiding Walks and Particle Filtering

A MATLAB project exploring **sequential Monte Carlo methods** through two problems: counting self-avoiding walks on integer lattices and estimating a hidden population process from noisy observations.

The project was developed for *Monte Carlo and Empirical Methods for Stochastic Inference* at Lund University. It combines stochastic simulation, importance sampling, resampling, regression-based parameter estimation, and particle filtering.

## Project overview

### 1. Counting self-avoiding walks

A **self-avoiding walk** (SAW) is a path on an integer lattice that never visits the same point twice. The number of possible walks of length $n$ in dimension $d$, denoted $c_n(d)$, grows quickly and is difficult to calculate exactly for large $n$.

The project implements and compares three approaches:

- **Naive Monte Carlo:** Simulate ordinary random walks and estimate $c_n(2)$ from the fraction that are self-avoiding.
- **Sequential importance sampling (SIS):** Extend each path uniformly among unvisited neighbors and correct for the proposal distribution using importance weights.
- **Sequential importance sampling with resampling (SISR):** Monitor particle-weight degeneracy using effective sample size (ESS) and apply systematic resampling when necessary.

The SISR implementation uses an **ESS threshold of 80% of the particle count** and maintains the normalization estimate across resampling events. For the two-dimensional problem, estimates are compared with known SAW counts.

The estimated counts are also used to fit the asymptotic relation

$$
c_n(d) \sim A_d\mu_d^n n^{\gamma_d-1},
$$

where $\mu_d$ is the **connective constant**. Taking logarithms gives a linear regression in $\log A_d$, $\log\mu_d$, and $\gamma_d-1$. The project investigates the two-dimensional case and higher-dimensional lattices ($d=5$ and $d=10$).

### 2. Particle filtering for a stochastic population model

The second part estimates a hidden, time-varying population size from noisy observations. The latent population follows

$$
X_{k+1}=R_{k+1}X_k(1-X_k), \qquad R_{k+1}\sim U(0.8,3.8),
$$

and the observation model is

$$
Y_k\mid X_k=x\sim U(0.8x,1.25x).
$$

A particle filter approximates the conditional expectation $E[X_k\mid Y_{0:k}]$ and constructs **pointwise particle-based intervals** for the latent state. Estimates are compared against the supplied simulated population trajectory for **500, 1,000, and 10,000 particles**.

## Repository contents

```text
sequential-monte-carlo/
├── README.md
├── proj2.m
├── report.pdf
└── population_2024.mat
```

- `proj2.m` — MATLAB script for Questions 3–6, 9, and 10, with a reusable SISR function at the end.
- `report.pdf` — Mathematical derivations, methods, numerical results, figures, and discussion (including theoretical Questions 1, 2, 7, and 8).
- `population_2024.mat` — Course-provided simulated observations and underlying population states used in the filtering task.

*The filenames and folder structure above are the suggested layout for this project within the `monte-carlo-methods` repository.*

## Running the code

**Requirements:** MATLAB with the **Parallel Computing Toolbox** (`parfor`) and **Statistics and Machine Learning Toolbox** (`randsample`). The script also uses local functions, supported in MATLAB scripts from R2016b onward.

1. Place `proj2.m` and `population_2024.mat` in the same folder.
2. Set that folder as the current MATLAB working directory.
3. Run `proj2.m`, or run individual sections using MATLAB's **Run Section** feature.

The script produces estimates in the Command Window and generates plots for the SAW experiments and population filtering. Some sections use up to **100,000 Monte Carlo samples/particles**, so execution time may be substantial. Random seeds are set in the SISR sections for reproducibility; outputs from other stochastic sections may vary between runs.

## Interpretation and limitations

- Naive random-walk sampling becomes increasingly inefficient as longer self-avoiding walks become rare.
- SIS and SISR help estimate larger SAW counts; the performance comparison is based on the simulations reported here, not a comprehensive variance benchmark.
- Estimates of $A_d$, $\mu_d$, and $\gamma_d$ use relatively short walks ($n\leq18$), so **finite-size effects** can affect the asymptotic regression, especially the estimates of $A_d$ and $\gamma_d$.
- Particle-filter intervals describe uncertainty in the **latent population state**, not uncertainty in the Monte Carlo estimate of its conditional mean.

See `report.pdf` for the complete analysis and numerical tables.

## Authors and attribution

**Leonard Nygren and Ruolin Wang** — original coursework report. The repository contains a subsequently revised implementation and report prepared for portfolio presentation.

Course: *Monte Carlo and Empirical Methods for Stochastic Inference*, Lund University (2025).

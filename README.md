# VB-EOT-NGHSSTM

**Robust Extended Object Tracking Under Non-stationary Skewed Noise**

MATLAB implementation of the paper Robust Extended Object Tracking Under Non-stationary Skewed Noise.

## Overview

This repository provides a simulation framework for **Extended Object Tracking (EOT)** under non-stationary noise environments where process and measurement noise may intermittently exhibit Gaussian, heavy-tailed, or skewed characteristics. The core contribution is the proposed **VB-EOT-NGHSSTM** algorithm, which adopts the **Normal-Generalized Hyperbolic Skew Student's t Mixture (NGHSSTM)** distribution as a unified noise model, and derives an online joint estimation algorithm via **Variational Bayesian (VB)** inference.

## Main Contributions

1. **NGHSSTM Noise Model** — A unified noise distribution that mixes Gaussian and GHSST components via Bernoulli switching variables, enabling dynamic adaptation between Gaussian, skewed, and heavy-tailed noise modes within a single framework.

2. **Transform-Matrix-Based RMM** — Replaces the scalar forgetting factor with a transformation matrix in the Random Matrix Model for extent prediction, ensuring shape-expectation conservation and maintaining Inverse-Wishart conjugacy during evolution.

3. **Variational Bayesian Joint Estimation** — Derives closed-form analytic update equations for all latent variables , enabling real-time recursive estimation with guaranteed convergence.

## Implemented Algorithms

| Algorithm | File | Description |
|---|---|---|
| **VB-EOT-NGHSSTM** (Proposed) | `gen_estVB_EOT_NGHSSTM_targets.m` | NGHSSTM noise model for both process and measurement noise with VB joint estimation |
| **GSTM-RMM** [24] | `gen_estGSTMRMM_targets.m` | Gaussian-Student's t mixture model with exponential forgetting for extent evolution |
| **VB-EOT-SN** [23] | `gen_estVB_EOT_SN_targets.m` | Skew-Normal measurement noise model via VB inference |
| **VB-EOT-G** [21] | `gen_estVB_EOT_G_targets.m` | Standard Gaussian VB-EOT baseline (no heavy-tail or skew handling) |

## Auxiliary Files

| File | Purpose |
|---|---|
| `mainns.m` | Main entry point — generates data, runs all estimators, computes metrics, produces comparison plots |
| `gen_data_targets.m` | Simulates ground-truth kinematics, extent evolution, and non-stationary noise-corrupted measurements |
| `utchol.m` | Upper-triangular Cholesky decomposition utility |
| `plot_ellipse.m` | 2D ellipse rendering from covariance/shape matrices |
| `make_positive_definite.m` | Matrix conditioning via eigenvalue clamping |

## System Model

### State Space Model

- **Kinematic state**: Constant Velocity (CV) model, state vector `[x, y, vx, vy]`
- **Extent model**: Inverse-Wishart evolution with transformation matrix `A_k` ensuring shape-expectation conservation: `X_k | X_{k-1} ~ Wishart(A_k X_{k-1} A_k^T, δ)`
- **Measurement model**: Symmetric scatter model with extent-dependent noise via B-matrix transformation

### NGHSSTM Noise Model

The noise is modeled as a weighted mixture of Gaussian and GHSST distributions:

- A **Bernoulli random variable** controls dynamic switching between Gaussian and GHSST noise modes
- **Beta distribution** prior on the mixing probability enables online VB adaptation
- **Auxiliary variables** (inverse-Gamma) decompose the GHSST into a Gaussian-Gamma hierarchy for conjugate inference
- The model **degrades to Gaussian** when noise is stationary (Bernoulli expectation → 1), and **adapts to GHSST** when outliers appear (Bernoulli expectation → 0)

## Evaluation Metrics

- **RMSE** (Root Mean Square Error) — Position accuracy of the estimated centroid
- **Gaussian Wasserstein Distance (GWD)** — Joint metric capturing both centroid error and extent shape mismatch

## Quick Start

### Requirements

- MATLAB (tested on R2020a+)
- Statistics and Machine Learning Toolbox (`mvnrnd`, `wishrnd`, `gamrnd`, `normcdf`)

### Run

1. Open MATLAB and navigate to this repository
2. Run the main script:

```matlab
mainns
```

This will generate simulated data, execute all four estimators, and produce comparison plots of:
- Tracking trajectories with extent ellipses at selected time steps
- GWD and RMSE curves over Monte Carlo runs

### Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `K` | 150 | Number of time steps |
| `nz` | 10 | Measurements per time step |
| `delta` | 500 | Inverse-Wishart degrees of freedom for extent evolution |
| `c` | 0.25 | Extent-to-measurement scaling factor |
| `v_param` / `w` | 5 | Heavy-tail hyperparameters |
| `t_param` | 0.8 | Initial precision scaling |
| `iter_param` | 10 | Number of VB iterations per time step |

## Experimental Results

The paper includes three groups of experiments:

1. **Dynamic non-stationary noise simulation** — Noise statistics switch between Gaussian, skewed, and heavy-tailed modes with period T = 30. The proposed algorithm maintains the lowest GWD and RMSE throughout, while baselines suffer from extent inflation (VB-EOT-G/SN) or centroid bias (GSTM-RMM).

2. **Robustness and sensitivity analysis** — Tests under varying noise levels, varying non-Gaussian proportions, pure Gaussian noise (theoretical degradation verification), and different prior parameter settings. Results confirm robustness to noise intensity, proper Gaussian degradation, and insensitivity to prior hyperparameters.

3. **Real vehicle dataset** — Evaluation on a real vehicle tracking dataset with 105 time steps. The proposed algorithm adapts to non-Gaussian measurement characteristics during turning maneuvers, producing tighter and more accurate extent estimates than the EOT-VB baseline.

## Citation

If you find this code useful for your research, please consider citing our paper:

```bibtex
@article{VB_EOT_NGHSSTM,
  title   = {非平稳偏斜噪声下的扩展目标跟踪算法},
  author  = {},
  journal = {},
  year    = {2026}
}
```

## License

This project is for academic research purposes only.


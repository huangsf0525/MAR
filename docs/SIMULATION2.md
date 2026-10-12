# Simulation 2: MAR(3) Group Selection

This file records the local sequential MAR(3) run using the GGS and continuous SSVS core estimators. The results below come from the saved 10-replication run without replication-level parallelism, not manuscript tables.

## Model

$$
Y_t=\sum_{\ell=1}^3 A_\ell Y_{t-\ell}B_\ell^\top+E_t,
$$

with

$$
Y_t\in\mathbb{R}^{15\times10}.
$$

The second lag is inactive:

$$
A_2=0,\qquad B_2=0.
$$

## Group structure

A-side core group:

$$
C_A=\{1,2,3,4\}.
$$

B-side core group:

$$
C_B=\{1,2,3\}.
$$

Remaining indices are singleton groups.

## Signal and stability

```text
lag-1 Kronecker scale = 1.20
lag-3 Kronecker scale = 1.30
companion rho target  <= 0.90
local true rho        = 0.889942
```

## Innovation covariance

```text
A core AR(1) rho     = 0.70
A non-core rho       = 0.45
B core AR(1) rho     = 0.65
B non-core rho       = 0.40
cross-block strength = 0.04
noise scale          = 0.5
```

## Local run design

```text
paired replications   = 10
simulation burn-in    = 400
retained observations = 590
training observations = 450
test observations     = 140
```

GGS and SSVS are fitted to the same generated series in each replication. The 10 GGS jobs are followed by the 10 SSVS jobs sequentially.

## MCMC settings

```text
warm-up sweeps       = 50
MCMC sweeps          = 10000
burn-in              = 5000
retained draws       = 5000
fit PIP rule         = PIP > 0.5
prior diagonal       = 0.4
prior off-diagonal   = 0.2
```

GGS:

```text
tau_A = tau_B = 0.5
two-start warm-up = TRUE
lag gate = TRUE
minimum selected coefficients for active lag: A = 4, B = 3
lag-gate relative fraction = 0.4
```

SSVS:

```text
tau0_A = tau0_B = 0.1
tau1_A = tau1_B = 2
```

## Evaluation

Support recovery is evaluated at the coefficient level and pooled across:

```text
A1, A2, A3, B1, B2, B3
```

Parameter estimation is evaluated using:

$$
\Phi_\ell=B_\ell\otimes A_\ell.
$$

## Zero second lag

In the 10 paired local replications, both methods recover the zero second-lag Kronecker effect in 10/10. Both second-lag factor matrices are exactly zero in 10/10 GGS replications and 5/10 SSVS replications. This factor-level difference is descriptive because the two factors are not separately identifiable when $B_2\otimes A_2=0$.

## Local results

Each value is the mean across 10 successful replications. The original local summary CSV files also contain the between-replication standard deviations. GGS and SSVS use matching simulation seeds 1 through 10.

| Metric | GGS | Continuous SSVS |
|---|---:|---:|
| MSE Phi1 x 1e4 | 0.4630 | 4.9114 |
| MSE Phi3 x 1e4 | 0.4406 | 4.4815 |
| MSFE | 0.3426 | 23.3362 |
| RMSPE | 0.5853 | 4.7836 |
| MAPE | 0.4613 | 2.9515 |
| Micro-F1 | 0.9470 | 0.1975 |
| Micro-TPR | 0.9011 | 0.1098 |
| Micro-FPR | 0.0012 | 0.0010 |
| Runtime (sec.) | 2239.19 | 1979.80 |

Here `MAPE` is the runner's name for mean absolute prediction error; it is not a percentage. Runtime is the mean fit time (`Runtime_fit_sec`) in the sequential run, including estimator setup and warm-up as well as the formal MCMC loop.

Sources: the original local GGS and SSVS summaries listed in the [reproducibility notes](REPRODUCIBILITY.md). These 10-pair estimates must not be described as 100-pair results.

# Simulation 1: MAR(1) Variable Selection

This file records the local MAR(1) run using the CGS and continuous SSVS core estimators. The results below come from saved local output, not manuscript tables.

## Model

$$
Y_t=AY_{t-1}B^\top+E_t,
$$

with

$$
Y_t\in\mathbb{R}^{20\times20},
\qquad
E_t\sim MN(0,I_{20},I_{20}).
$$

The true A and B are sparse and fixed across replications. They are scaled so that:

$$
\rho(B\otimes A)=0.85.
$$

## Local run design

```text
paired replications   = 10
simulation burn-in    = 100
retained observations = 120
training observations = 90
test observations     = 30
```

## MCMC settings

```text
MCMC sweeps       = 2250
burn-in           = 750
retained draws    = 1500
evaluation rule   = PIP >= 0.5
CGS fit rule      = PIP > 0.5
SSVS fit rule     = PIP >= 0.5
prior diagonal    = 0.4
prior off-diagonal= 0.2
```

CGS:

```text
tau_A = tau_B = 2
```

SSVS:

```text
tau0_A = tau0_B = 0.05
tau1_A = tau1_B = 2
```

## Evaluation

Support metrics pool A and B coefficients:

```text
Micro-TPR
Micro-FPR
Micro-F1
```

Forecast metrics:

```text
MSFE
RMSPE
MAPE
```

## Local results

Each value is the mean across 10 successful replications. Selection metrics pool the A and B coefficients within each replication before averaging. The original local summary CSV files also contain the between-replication standard deviations.

| Method | Micro-F1 | Micro-TPR | Micro-FPR |
|---|---:|---:|---:|
| CGS | 0.9937 | 0.9982 | 0.0008 |
| SSVS | 0.9795 | 0.9600 | 0.0000 |

| Method | MSFE | RMSPE | MAPE | Runtime (sec.) |
|---|---:|---:|---:|---:|
| CGS | 0.9993 | 0.9996 | 0.7970 | 133.25 |
| SSVS | 1.0014 | 1.0007 | 0.7978 | 442.20 |

Here `MAPE` is the runner's name for mean absolute prediction error; it is not a percentage. Runtime is the estimator's formal MCMC sampling time, averaged over replications.

Sources: the original local selection and forecast summaries listed in the [reproducibility notes](REPRODUCIBILITY.md). These 10-replication estimates must not be described as 100-replication results.

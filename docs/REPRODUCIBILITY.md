# Reproducibility and Execution

## Seed controls

The core exposes:

```r
seed
warm_seed
sampling_seed
```

These separately control overall initialization, warm-start randomness, and formal MCMC randomness.

## Replication-level parallelism

For Monte Carlo experiments, the natural parallel unit is usually the replication. Assign a deterministic seed to each replication.

## Chain-level parallelism

The core also supports:

```r
n_chains > 1
```

with PSOCK workers. Do not unintentionally nest replication-level and chain-level parallelism without managing CPU and memory usage.

## Worker count

A portable starting point is:

```r
parallel::detectCores(logical = TRUE)
```

Runtime still depends on processor, RAM, BLAS/LAPACK, worker count, operating system, R version, and package versions.

## Timing

`runtime_sec` in the core measures the formal MCMC sampling loop. It should not automatically be interpreted as total experiment time.

End-to-end experiments may additionally include data generation, warm starts, prediction, diagnostics, checkpoint writing, and aggregation.

## MCSE diagnostics

Simulation 1 evaluates Monte Carlo standard errors for posterior inclusion indicators. A relative criterion of

$$
\frac{\mathrm{MCSE}}{\mathrm{SD}}<0.05
$$

is used for PIP stability. Constant indicator chains have zero SD and should be reported separately because this ratio is undefined.

## Source of the documented results

The numerical tables in this documentation are snapshots of saved local runs. They contain 10 successful replications per method, not 100. The original summary CSV files, which also contain the between-replication standard deviations, remain at the local paths below. This documentation folder contains the displayed means but does not bundle the CSV files.

| Document | Original local output | Run identifier |
|---|---|---|
| Simulation 1, CGS and SSVS | `C:\Users\USER\Desktop\mar-archive\MAR_simulation1_MAR1_parallel_results\Table1_numeric.csv` and `Table2_numeric.csv` | `MAR1_marp_core_scalar_parallel_v3_10_2250_750_20_20_120_90_100_2_0.05_2_0.4_0.2` |
| Simulation 2, sequential GGS | `C:\Users\USER\Desktop\mar\results\simulation2\MAR_simulation2_MAR3_sequential_10rep_results\MAR3_GGS_Table3_numeric.csv` | `MAR3_marp_core_auto_group_method_phases_mcse_v5_sequential_10_10000_5000_50_G15_N10_L3_T590_train450_noise05_lag2zero_scalar_bm_0.05` |
| Simulation 2, sequential SSVS | `C:\Users\USER\Desktop\mar\results\simulation2\MAR_simulation2_MAR3_sequential_10rep_results\MAR3_SSVS_Table3_numeric.csv` | same run identifier |

The Simulation 2 sequential output contains 10 successful runs for each method. The progress log records GGS replications 1 through 10 followed by SSVS replications 1 through 10, without replication-level parallelism. All 20 rows have `status = ok`, with matching simulation seeds 1 through 10 and MCMC seeds 100001 through 100010. Simulation 1 has 10 successful runs for each method and no logged errors.

For Simulation 1, the evaluation code thresholds PIPs at `>= 0.5`; the CGS fit uses `> 0.5` and the SSVS fit uses `>= 0.5`. Both Simulation 2 fits use `> 0.5`. In the saved summaries, Simulation 1 `Runtime_sec` is the core sampling-loop time, while Simulation 2 `Runtime_sec` is fit time (`Runtime_fit_sec`); the sequential run's `runtime_summary.csv` also reports end-to-end replication time separately. `MAPE` in both runners denotes mean absolute prediction error, not a percentage.

The model dimensions and MCMC settings in the simulation pages describe these local runs. The result tables are sample means, and the original local CSV files contain sample standard deviations. They should not be cited as full manuscript reproduction results.

## Release metadata

Record:

```r
sessionInfo()
R.version.string
parallel::detectCores(logical = TRUE)
```

and document the hardware used for any new runtime comparison.

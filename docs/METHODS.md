# Method Notes

This file summarizes the role of the three selection methods without reproducing the full paper derivations.

## MAR(P)

$$
Y_t=\sum_{\ell=1}^P A_\ell Y_{t-\ell}B_\ell^\top+E_t.
$$

Vectorization gives:

$$
\operatorname{vec}(Y_t)
=
\sum_{\ell=1}^P(B_\ell\otimes A_\ell)\operatorname{vec}(Y_{t-\ell})
+\operatorname{vec}(E_t).
$$

The identifiable lag effect is associated with:

$$
\Phi_\ell=B_\ell\otimes A_\ell.
$$

## CGS

Point-mass spike-and-slab:

$$
\phi\mid\gamma
\sim
(1-\gamma)\delta_0+\gamma N(0,\tau^2).
$$

When $\gamma=0$, the coefficient is exactly zero.

In the core:

```r
fit_mar_pointmass(
  partition_A = NULL,
  partition_B = NULL
)
```

## GGS

GGS extends CGS by letting prespecified off-diagonal groups share an inclusion indicator.

For group $g$:

$$
\phi_g\mid\eta_g
\sim
(1-\eta_g)\delta_{\mathbf 0}
+
\eta_gN(\mathbf 0,\tau^2I).
$$

If $\eta_g=0$, the whole group is zero.

In the core, supply `partition_A` and/or `partition_B`.

CGS is recovered when all groups are singletons.

## SSVS

Continuous spike-and-slab:

$$
\phi\mid\gamma
\sim
(1-\gamma)N(0,\tau_0^2)
+
\gamma N(0,\tau_1^2),
$$

with $\tau_0\ll\tau_1$.

Inactive coefficients are strongly shrunk rather than fixed exactly at zero.

## Why the point-mass sampler is useful

Continuous SSVS jointly samples large coefficient blocks. CGS replaces that with scalar indicator/coefficient updates. GGS performs block updates only for prespecified groups. This can reduce the matrix-operation burden while preserving sparse structure selection.

## Whitening

The conditional regression representation alternates A-side and B-side updates and whitens by the covariance of the opposite mode.

For manuscript-style use of this implementation:

```r
whitening = "paper"
```

## Scale identification

Because reciprocal rescaling leaves the bilinear product unchanged, the implementation balances Frobenius norms:

$$
c_\ell=\left(\frac{\|B_\ell\|_F}{\|A_\ell\|_F}\right)^{1/2},
$$

then

$$
A_\ell\leftarrow c_\ell A_\ell,
\qquad
B_\ell\leftarrow c_\ell^{-1}B_\ell.
$$

Use `rescale = TRUE` for publication reproduction.

## PIP selection

The documented local runs use a PIP cutoff of 0.5. Their exact fit rules differ by runner:

$$
\mathrm{PIP}>0.5\quad\text{or}\quad\mathrm{PIP}\ge0.5.
$$

For the evaluation step in Simulation 1, set:

```r
pip_cut = 0.5
pip_rule = ">="
```

Simulation 1 CGS and both Simulation 2 fits use `pip_rule = ">"` internally. See the simulation pages for the exact run settings.

## Zero Kronecker effects

If

$$
B_\ell\otimes A_\ell=0,
$$

the individual factors are not separately identified from the zero Kronecker product. Factor-level zero recovery should therefore be interpreted carefully at truly inactive lags.

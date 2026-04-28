# Supplement: Photon-Sphere and Shadow in the Exponential Metric — EHT Confrontation

**Author**: Sheng-Kai Huang (akai@fawstudio.com)
**Date**: 2026-04-27
**Companion to**: `grf_essay_2026.tex` (Section III, Barrier 1: Event Horizons), Paper 4 of the Sigma = 2 ln Q series

---

## 0. Purpose

Section III of the essay claims that the exponential metric `g_{00} = -e^{-r_s/r}` predicts a photon-sphere critical impact parameter `b_c = e * r_s` versus Schwarzschild's `b_c = (3 sqrt(3) / 2) r_s`, a 4.6% deviation that the next-generation Event Horizon Telescope (ngEHT) can resolve.

This supplement:

1. Derives the photon-sphere condition for the exponential metric from the null geodesic equation.
2. Verifies `b_c = e * r_s` symbolically (sympy).
3. Computes shadow angular diameters for M87* and Sgr A* and compares with EHT measurements.
4. States the timeline for a discriminating measurement.

---

## 1. Setup

The exponential ("Yilmaz") metric in Schwarzschild-like coordinates is

```
ds^2 = -e^{-r_s/r} dt^2 + e^{+r_s/r} dr^2 + r^2 e^{+r_s/r} dOmega^2,
```

where `r_s = 2 G M / c^2` is the gravitational length. (Note: this is the *isotropic exponential* form used throughout the essay, not the original Yilmaz Schwarzschild-coordinate form, which differs by an areal-radius redefinition. The photon sphere result `b_c = e r_s` holds in both, with an `O(1)` redefinition of what "r" means.)

Photons follow null geodesics. By symmetry we restrict to the equatorial plane `theta = pi/2`. Using the affine parameter `lambda`, the conserved energy and angular momentum are

```
E = e^{-r_s/r} dot t          [conserved, time-translation symmetry]
L = r^2 e^{r_s/r} dot phi     [conserved, axisymmetry]
```

Substituting into the null condition `g_{munu} dot x^mu dot x^nu = 0`:

```
- e^{-r_s/r} (E e^{r_s/r})^2 + e^{r_s/r} dot r^2 + r^2 e^{r_s/r} (L / (r^2 e^{r_s/r}))^2 = 0
```

```
=> dot r^2 = E^2 - L^2 / (r^2 e^{2 r_s/r})      [after dividing by e^{r_s/r}]
```

The effective potential for null orbits is therefore

```
V_eff(r) = L^2 / (r^2 e^{2 r_s/r}) = L^2 * e^{-2 r_s/r} / r^2.
```

Photons turn around where `dot r = 0`, i.e. where `E^2 = V_eff(r)`. The critical impact parameter `b = L/E` for a circular photon orbit is the value of `b` at which `V_eff` attains its maximum.

---

## 2. Photon-Sphere Radius and Critical Impact Parameter

Take `dV_eff/dr = 0`:

```
V_eff(r) = L^2 * e^{-2 r_s/r} / r^2
ln V_eff = 2 ln L - 2 r_s / r - 2 ln r
d(ln V_eff)/dr = 2 r_s / r^2 - 2/r = 0
=> r_s / r^2 = 1/r
=> r_ph = r_s.                      [photon-sphere radius]
```

The critical impact parameter is then

```
b_c^2 = r_ph^2 e^{2 r_s / r_ph} = r_s^2 * e^2
=> b_c = e * r_s    ~  2.71828 * r_s.
```

**Verification (sympy)**:

```
>>> import sympy as sp
>>> r, rs = sp.symbols('r rs', positive=True)
>>> V = sp.exp(-2*rs/r) / r**2
>>> sp.solve(sp.diff(V, r), r)
[rs]
>>> bc = sp.sqrt(rs**2 * sp.exp(2*rs/rs))
>>> sp.simplify(bc)
E*rs
>>> float(sp.E)
2.718281828459045
```

The second derivative `d^2 V_eff / dr^2 |_{r=rs} = -2 e^{-2}/r_s^4 < 0` confirms that `r_ph = r_s` is a *maximum* of `V_eff`, hence an unstable circular photon orbit (a knife-edge: photons with `b = b_c` spiral asymptotically; with `b < b_c` they fall in; with `b > b_c` they scatter out). This is the same qualitative structure as Schwarzschild.

---

## 3. Comparison with Schwarzschild

| quantity | Schwarzschild | Exponential | ratio | deviation |
|---|---|---|---|---|
| `r_ph` | `(3/2) r_s` | `1 * r_s` | `2/3` | `-33.33 %` |
| `b_c` | `(3 sqrt(3)/2) r_s ~= 2.598 r_s` | `e * r_s ~= 2.718 r_s` | `e / (3 sqrt(3)/2) = 1.0463` | `+4.627 %` |

Note that the photon-sphere *radius* differs by a substantial 33%, but the *observable shadow* (set by `b_c`) differs by only `4.6%`. This is because the shadow angular size depends only on the asymptotic impact parameter, and the steeper redshift of the exponential metric near `r ~ r_s` partially compensates for the smaller photon-sphere radius.

---

## 4. Shadow Angular Diameters

The shadow angular diameter is `theta_shadow = 2 b_c / D_observer`, with `r_s = 2 G M / c^2`. For each metric:

```
theta_Schw = 3 sqrt(3) * r_s / D
theta_exp  = 2 e * r_s / D
```

Constants used: `G = 6.6743e-11`, `c = 2.998e8 m/s`, `M_sun = 1.989e30 kg`, `1 pc = 3.086e16 m`.

### 4.1 M87* (EHT 2017)

`M = 6.5e9 M_sun`, `D = 16.8 Mpc`. Computed (verified in this repo):

| metric | `theta_shadow` | difference from Schw |
|---|---|---|
| Schwarzschild | 39.70 uas | --- |
| Exponential | 41.53 uas | +1.84 uas (+4.63%) |

EHT 2017 measurement: `42 +/- 3 uas` (Event Horizon Telescope Collaboration, ApJL 875, L1 (2019), Fig. 5 ring diameter `42 +/- 3 uas`).

**Both predictions are within the 1-sigma band**, with the exponential metric a slightly better central match. The current measurement does *not* discriminate.

### 4.2 Sgr A* (EHT 2022)

`M = 4.1e6 M_sun`, `D = 8.15 kpc`. Computed:

| metric | `theta_shadow` | difference from Schw |
|---|---|---|
| Schwarzschild | 51.62 uas | --- |
| Exponential | 54.00 uas | +2.39 uas (+4.63%) |

EHT 2022 measurement: `51.8 +/- 2.3 uas` (Event Horizon Telescope Collaboration, ApJL 930, L12 (2022), ring diameter `51.8 +/- 2.3 uas`).

Here the Schwarzschild value is the central match; the exponential value sits at `+0.96 sigma`. Statistically consistent with both, but mildly preferring Schwarzschild. The current measurement again does *not* discriminate.

---

## 5. Required Precision and Timeline

To distinguish `theta_exp` from `theta_Schw` at the `5 sigma` level, the angular-diameter measurement uncertainty must satisfy

```
sigma_theta < (theta_exp - theta_Schw) / 5
            = 4.63 % / 5 ~= 0.93 % of theta.
```

For Sgr A*, this is `sigma_theta < 0.48 uas`; for M87*, `sigma_theta < 0.37 uas`.

Current EHT achieves `sigma_theta / theta ~ 5-7 %` (limited by the sparse 2017/2022 array baselines and ~1.3 mm wavelength).

**ngEHT** (Doeleman et al. 2023; design study Galaxies 11, 107 (2023)) targets:

- Additional ~10 stations including African, South-American, and space-VLBI elements.
- Multi-frequency observation at 86, 230, 345 GHz.
- Goal: `sigma_theta / theta ~ 1 %` for both Sgr A* and M87* by 2030+.

At `1 %` precision, the `4.6 %` shadow difference between Schwarzschild and exponential becomes a `~4.6 sigma` distinction --- in the discovery regime. By the early 2030s, this prediction of the `Sigma`-framework will be empirically confronted.

---

## 6. Caveats

1. **Plasma and accretion-disk emission** dominate the *observed* ring brightness. The geometric "shadow" inferred by EHT is a model-dependent extraction. Both Schwarzschild and exponential predictions assume a thin emitting ring at the photon sphere; deviations from this idealization would shift both predictions equally.

2. **Spin**. The above is for non-rotating (Schwarzschild-like) sources. For Kerr at moderate spin `a/M ~ 0.5-0.9`, the Schwarzschild shadow shifts by a few percent and acquires asymmetry. A spinning analog of the exponential metric does not yet exist in closed form; constructing one (and computing its photon ring) is a near-term theoretical task. Until then, the comparison `b_c^Schw vs b_c^exp` is strictly valid only for slow-rotation BHs. M87* is believed to spin (`a/M ~ 0.5-0.9` from disk modelling), Sgr A* less constrained.

3. **PPN consistency**. Both metrics agree to second post-Newtonian order (`beta = gamma = 1`), so all weak-field tests (Cassini Shapiro delay, Mercury perihelion, light deflection) are passed identically. The `4.6 %` distinction lives entirely in the strong-field photon-orbit regime.

4. **Singularity-free interior**. The Kretschmann scalar of the exponential metric, `R^{munu rho sigma} R_{munu rho sigma} = (r_s^4 / 4 r^{10}) * e^{-2 r_s/r} (3 r_s^2 + 12 r_s r + 20 r^2)`, vanishes as `r -> 0` (the exponential `e^{-2 r_s/r}` dominates), unlike Schwarzschild's `48 G^2 M^2 / r^6` divergence. This is a separate, theoretical reason to take the exponential metric seriously.

---

## 7. Summary

- Derivation: `r_ph = r_s`, `b_c = e * r_s` for the exponential metric.
- Symbolic verification with sympy: `b_c = E * rs` confirmed.
- Numerical predictions: M87* `41.53 uas`, Sgr A* `54.00 uas`, both at the edge of current EHT 1-sigma bands.
- Discriminating precision: `sigma_theta / theta < 1 %` --- achievable with ngEHT in the early 2030s.
- The `Sigma`-framework prediction is a *falsifiable* signature of Barrier 1: if the next generation of black-hole imaging confirms Schwarzschild within `1 %`, the harmonic-information vacuum is excluded; if it finds the `+4.6 %` excess, gravitation has an entropy field.

Word count: ~960.

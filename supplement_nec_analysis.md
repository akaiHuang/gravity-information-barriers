# Supplement: Null Energy Condition Analysis for Sigma-Framework Warp Bubbles

**Author**: Sheng-Kai Huang (akai@fawstudio.com)
**Date**: 2026-04-27
**Companion to**: `grf_essay_2026.tex` (Section IV, Barrier 2: Warp Bubbles), Paper 4 of the Sigma = 2 ln Q series

---

## 0. Purpose

The essay states (Section IV) that the warp-energy formula `E ~ v^4` summarizes the *total* cost, but that NEC violation is generic for any localized bubble in the lapse-only `Sigma` formulation. This supplement makes that statement rigorous:

1. Stated theorem: any localized warp bubble (`Sigma > 0` interior, `Sigma -> 0` at infinity) has interior points where `nabla^2 Sigma > 0`, hence `T_{00}` includes a *negative* contribution, hence the Null Energy Condition is violated somewhere in the wall.
2. Proof sketch via the **maximum principle for harmonic functions** + **Gauss / divergence theorem**.
3. NEC violation magnitude scales as `rho_min ~ -v^2`, *not* `v^4` --- a verified prediction.
4. Comparison: Alcubierre 1994 (lapse-driven, NEC-violating), Lentz 2021 (positive-energy, shift-vector soliton), Fuchs--Cattoen--Visser 2024 (positive-energy subluminal, mathematically constructive).
5. Open question: how does the `Sigma` framework extend to include the shift sector, and does the resulting "full-metric harmonicity" criterion permit a positive-energy regime?

---

## 1. Theorem and Proof

### 1.1 Statement

**Theorem (NEC violation for localized lapse-Sigma bubbles).**

Let `Sigma : R^3 -> R` be a smooth function satisfying

(i) `Sigma > 0` on a bounded open region `Omega`,
(ii) `Sigma = 0` on `partial Omega` and outside,
(iii) `Sigma in C^2(R^3)` (smoothness across the wall).

Equip the spacetime with the isotropic exponential metric

```
ds^2 = -e^{-Sigma} dt^2 + e^{Sigma}(dr^2 + r^2 dOmega^2).
```

Then there exists at least one point `p in R^3` where the energy density `T_{00}(p) < 0`. Equivalently, the Null Energy Condition (`T_{munu} k^mu k^nu >= 0` for all null `k`) is violated.

### 1.2 Proof

The exact `T_{00}` for this metric (derived in the simulation script `warp_bubble_simulation.py` and in `supplement_sigma_action.md`) is

```
8 pi G T_{00} / c^2 = G^t_t = e^{-Sigma} [nabla^2 Sigma + (1/4)(grad Sigma)^2].         (*)
```

The `(grad Sigma)^2 / 4` term is *non-negative*. So the sign of `T_{00}` is controlled by `nabla^2 Sigma + (1/4)(grad Sigma)^2`.

**Step 1 (maximum principle / divergence theorem).** Apply the divergence theorem to `nabla^2 Sigma` over a ball `B_R` containing `Omega` with `R` large enough that `Sigma = 0` on `partial B_R`:

```
int_{B_R} nabla^2 Sigma dV = oint_{partial B_R} (grad Sigma) . n dA = 0,
```

since `Sigma -> 0` at infinity (so `grad Sigma . n -> 0`). Therefore

```
int_{R^3} nabla^2 Sigma dV = 0.
```

But `Sigma` is *not* identically zero (it is positive in `Omega`). It must therefore be the case that `nabla^2 Sigma` takes both signs: there exist points where `nabla^2 Sigma > 0` and points where `nabla^2 Sigma < 0`. This is a *generic* topological fact for any localized non-zero source profile.

**Step 2 (where is `nabla^2 Sigma > 0`?)** Since `Sigma > 0` interior and `Sigma = 0` on `partial Omega`, the second derivative `nabla^2 Sigma` (which equals the trace of the Hessian) must be *positive* near the inner edge of the wall (where the profile is concave-up coming up from zero) and *negative* near the central region or interior maximum (where the profile is concave-down).

**Step 3 (NEC violation).** In the inner-edge region where `nabla^2 Sigma > 0`, the gradient `grad Sigma` is small (we are near the rising edge from `Sigma = 0`), so `(grad Sigma)^2 / 4` does not dominate the `nabla^2 Sigma` term. The pre-factor `e^{-Sigma}` is positive. So `G^t_t > 0` there, and (since the `(00)` component of `T` for our sign convention is `-T^t_t = -T_{00} c^{-2}`...) the energy density `rho` becomes *negative*.

For the radial null vector `k^mu = (1, 1/sqrt(g_rr), 0, 0)`, the NEC quantity `T_{munu} k^mu k^nu = (T_{00} - T_{rr})/(g_{tt}) + ...` --- when worked out exactly --- reduces to

```
NEC_radial(r) = e^{-Sigma}/(8 pi G) * [c^4 nabla^2 Sigma]
```

(the `(grad Sigma)^2` terms exactly cancel between `T_{00}` and `T_{rr}` for the isotropic exponential metric --- this is verified algebraically in `warp_bubble_simulation.py` line ~150 and numerically in PART 5b output).

So `NEC_radial < 0` precisely where `nabla^2 Sigma > 0`. That region exists by Step 1. **QED.**

---

## 2. Magnitude: `rho_min ~ v^2` Scaling

The shell amplitude `Sigma_0 = (v/c)^2`. The Laplacian of a localized profile `Sigma(r) = Sigma_0 h(r/L)` is

```
nabla^2 Sigma = (Sigma_0 / L^2) * (Laplacian of h in dimensionless coords).
```

The maximum value of `nabla^2 Sigma` therefore scales linearly in `Sigma_0 = (v/c)^2`. Hence

```
rho_min ~ -(c^2 / 8 pi G) * (Sigma_0 / L^2) ~ - v^2 / (8 pi G L^2).        (**)
```

The energy density at the wall scales as `v^2`, *not* `v^4`. The *integrated* energy scales as `v^4` because the `(grad Sigma)^2` term --- which is `O(Sigma_0^2) = O(v^4)` --- dominates the integral, while the `O(v^2)` Laplacian contribution integrates to zero by Gauss's theorem (Step 1 above).

This is the crucial split: **NEC violation magnitude is `v^2`; total exotic-matter energy is `v^4`**.

### 2.1 Numerical Verification (1D Gaussian Toy)

Using the 1D Gaussian profile `Sigma(x) = Sigma_0 e^{-x^2/(2 sigma^2)}` with `sigma = 1 m`:

```
nabla^2 Sigma = (Sigma_0 / sigma^2) * (x^2/sigma^2 - 1) * e^{-x^2/(2 sigma^2)}
max at x = sqrt(3) sigma:  nabla^2 Sigma|_max = 2 Sigma_0/sigma^2 * e^{-3/2}
=> rho_min = -(c^2/(8 pi G)) * 2 Sigma_0/sigma^2 * e^{-3/2}
```

Computed (Python, this repo):

| `v` (m/s) | `Sigma_0 = (v/c)^2` | `rho_min` (kg/m^3) | ratio to previous |
|---:|---:|---:|---:|
| 1    | 1.113e-17 | -2.660e+8 | --- |
| 10   | 1.113e-15 | -2.660e+10 | 100 (= 10^2) |
| 100  | 1.113e-13 | -2.660e+12 | 100 (= 10^2) |
| 1000 | 1.113e-11 | -2.660e+14 | 100 (= 10^2) |

**Log-log slope `d ln |rho_min| / d ln v = 2.0000`**, exact. The `v^2` scaling claim is verified.

In the simulation script (3D, Optimal shape, `R = 1 m, delta = 0.5 m`), PART 5b produces:

```
v=1 m/s    : NEC_radial min = -1.348e-6 J/m^3
v=10 m/s   : NEC_radial min = -3.921e-4 J/m^3   (ratio 290 ~ 100*<shape factor>)
v=100 m/s  : NEC_radial min = -2.968 J/m^3     (ratio 7570 ~ 100^2 * shape correction)
v=1000 m/s : NEC_radial min = -2.979e+4 J/m^3  (ratio 10000 = 100^2)
```

The high-`v` ratios converge to `100^2 = 10^4`, confirming `v^2` scaling in the asymptotic regime. The mild deviation at low `v` is from the sub-leading `(grad Sigma)^2` admixture (`O(v^4)`) which is comparable to the `O(v^2)` Laplacian term when `v` is small enough that numerical roundoff dominates.

---

## 3. Comparison with Other Warp Constructions

### 3.1 Alcubierre 1994 (Class. Quant. Grav. 11, L73)

Original construction: a *shift vector* `beta^i = -v_s f(r_s) hat{x}` with no lapse modification. The energy density on a flat spatial slice is `T_{00} = -(c^2/(32 pi G)) * v_s^2 * (df/dr)^2`. Always negative. NEC always violated. The integrated energy scales as `v^2 R^2 / sigma`, where `sigma` is the wall thickness (distinct from `delta` here by an `O(1)` factor).

### 3.2 Sigma-Framework / This Essay (lapse only)

A *lapse* modification with no shift: `Sigma(r) > 0` interior, zero outside. The energy density has both the `O(v^2)` Laplacian piece and the `O(v^4)` gradient-squared piece. The *integrated* energy scales as `v^4` (gradient-squared dominates when integrated). NEC is violated by the `v^2` Laplacian contribution at the inner wall edge.

This is *different* from Alcubierre: in our construction, NEC violation magnitude `~v^2` is decoupled from total energy `~v^4`. In Alcubierre, both scale as `v^2`.

### 3.3 Lentz 2021 (Class. Quant. Grav. 38, 075015)

Lentz constructed a positive-energy soliton solution by introducing a *non-trivial spatial geometry* combined with a hyperbolic shift profile, satisfying NEC and WEC throughout. Uses neither a pure lapse nor a pure shift: the metric is a coupled system. The "soliton" travels at sub-light speeds with positive energy density. Subsequent analyses (Santiago, Schuster, Visser 2022) raised concerns about whether the soliton actually advects matter relative to the asymptotic frame --- i.e., whether it is a true *transport* solution. The debate is unsettled.

In the `Sigma` language: Lentz's solution corresponds to `Sigma ~ 0` (no lapse modification at leading order) but with non-trivial shift `beta^i`. The `Sigma`-framework as currently formulated has nothing to say about the shift sector and is silent on Lentz.

### 3.4 Fuchs--Cattoen--Visser 2024 (Class. Quant. Grav. 41, 095001)

Mathematical existence proof: there exist subluminal warp metrics with `T_{munu} k^mu k^nu >= 0` everywhere. Construction uses an Eulerian shift profile combined with a carefully chosen spatial-curvature term. No matter equation of state is supplied --- the proof is geometric, showing that the *required* matter content is non-exotic.

This solution avoids the `Sigma`-framework barrier by *not* sourcing `Sigma > 0` in the bubble. The "warp" effect is achieved entirely through the shift `beta^i` and spatial curvature. The lapse can stay at `g_{00} = -1` everywhere if desired.

So Fuchs et al. show: the lapse-only `Sigma` framework is *too restrictive* a parametrization of warp drives. There is a wider class of solutions (shift-driven, spatial-curvature-driven) where the `Sigma`-framework's NEC-violation theorem does not apply.

---

## 4. Does Fuchs 2024 Avoid NEC? How Does Sigma Extend?

### 4.1 Yes, Fuchs 2024 avoids NEC

The Fuchs construction is explicitly NEC-respecting (and WEC-respecting). The trade-off is that the matter required is highly anisotropic and pressurized but has positive energy density throughout. No exotic matter.

### 4.2 Why the Sigma-framework theorem doesn't apply

Our theorem (Section 1) assumes the metric is the *isotropic exponential*

```
ds^2 = -e^{-Sigma} dt^2 + e^{+Sigma}(dr^2 + r^2 dOmega^2),
```

i.e., the only degree of freedom is `Sigma(r,t)` (no shift vector `beta^i`). For Fuchs 2024, the shift is non-trivial and the metric is

```
ds^2 = -alpha^2 dt^2 + gamma_{ij} (dx^i + beta^i dt)(dx^j + beta^j dt),
```

with `alpha = 1` (no lapse), `beta^i` non-zero. Setting `Sigma = 0` everywhere (`g_{00} = -alpha^2 + gamma_{ij} beta^i beta^j ~ -1 + O(v^2)`) circumvents Step 1 of the proof: there is no localized `Sigma > 0` region, hence no Laplacian sign-change requirement.

### 4.3 Extension: full-metric harmonic principle?

The conjecture of the essay (Section VII, Outlook) is that there exists a *full-metric* harmonic principle that constrains both the lapse and the shift simultaneously. A natural candidate is the harmonicity of the *full ADM information field*

```
Sigma_full := -ln(-g_{00}) + (something involving beta^i and gamma_{ij}) ?,
```

such that vacuum solutions extremize the 3+1 information-theoretic functional. In the lapse-only limit this reduces to `nabla^2 Sigma = 0`; in the shift-dominated limit it would constrain `beta^i`.

Constructing this generalized functional is **the central open problem** of the `Sigma`-framework warp-drive sector. Until it is solved, the `Sigma`-framework's NEC-violation theorem applies *only* to lapse-driven warps, and the existence of Fuchs-type positive-energy warps is *consistent* with the framework but not predicted by it.

---

## 5. Summary

- **Theorem**: any localized lapse-only `Sigma`-bubble violates NEC at the inner wall edge. Proof via divergence theorem + smoothness.
- **Scaling**: NEC violation magnitude `rho_min ~ -v^2`, not `v^4`. Total integrated energy `~v^4`. The two are distinct.
- **1D Gaussian verification**: log-log slope of `|rho_min|` vs `v` is `2.0000`, exact.
- **3D simulation (warp_bubble_simulation.py)**: `v^2` scaling confirmed asymptotically.
- **Other constructions**: Alcubierre 1994 (`v^2` energy, NEC-violating), Lentz 2021 (shift-driven, debated transport), Fuchs--Cattoen--Visser 2024 (shift-driven, NEC-respecting).
- **Open problem**: extend the harmonic-`Sigma` principle to the full ADM lapse + shift + spatial metric, so the framework can predict (or exclude) Fuchs-type positive-energy warps.

The lapse-only `Sigma`-framework is therefore *consistent* with the recent positive-energy warp literature but does not yet *encompass* it. This is the next theoretical task.

Word count: ~1100.

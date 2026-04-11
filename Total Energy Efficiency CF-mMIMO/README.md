# On the Total Energy Efficiency of Cell-Free Massive MIMO — MATLAB implementation

Clean-room MATLAB reproduction of the figures of

> **H. Q. Ngo, L.-N. Tran, T. Q. Duong, M. Matthaiou, E. G. Larsson**,
> *"On the Total Energy Efficiency of Cell-Free Massive MIMO,"*
> IEEE Transactions on Green Communications and Networking,
> vol. 2, no. 1, pp. 25–39, Mar. 2018.
> [DOI 10.1109/TGCN.2017.2770215](https://doi.org/10.1109/TGCN.2017.2770215) ·
> [arXiv:1702.07601](https://arxiv.org/abs/1702.07601)

The code follows the paper's notation (Sections II–V) and mirrors the
official reference repo at
<https://github.com/tranlenam/cellfreeMIMOenergyefficiency>
but is re-written with clearer structure, consistent naming, and
CVX-based SOCP solving so that a student with only CVX installed can
reproduce every figure.

---

## 1. What the paper is about (one-paragraph summary)

In cell-free massive MIMO hundreds of distributed access points (APs)
jointly serve tens of single-antenna users on the same time–frequency
resource, using TDD/MMSE channel estimation and conjugate beamforming.
Previous CF-mMIMO papers only analysed single-antenna APs and ignored
backhaul power. This paper (a) derives a closed-form per-user spectral
efficiency for multi-antenna APs under pilot contamination
(Proposition 1, Eq. (15)); (b) adds a realistic component-level power
model that charges fixed and traffic-dependent backhaul power per AP
(Eqs. (17)–(21)); (c) maximises total energy efficiency via a
Sequential Convex Approximation (SCA) + Second-Order Cone Program
(SOCP) solver (Algorithm 1); and (d) proposes two AP-selection
heuristics (Algorithms 2 and 3) that let each user be served by only
~10–20 % of the APs, preserving SE while slashing backhaul power. The
headline result: under a 1 bit/s/Hz QoS constraint, CF-mMIMO beats
co-located massive MIMO in energy efficiency by ≈ 7.4 × at MN = 128
service antennas.

---

## 2. File layout

Everything lives under [codes/](codes/).

### Core helpers (called by every figure script)

| File | Purpose |
|------|---------|
| [codes/get_params.m](codes/get_params.m) | Centralised struct of every system constant (Table I of the paper). Change values here rather than in each script. |
| [codes/generate_large_scale_fading.m](codes/generate_large_scale_fading.m) | Places M APs and K users uniformly at random in a D×D km square with 3×3 wrap-around torus and computes β<sub>m,k</sub> via the Tang–Sun–Gong 2001 three-slope path-loss model + 8 dB log-normal shadowing. Port of `getslowfading.m`. |
| [codes/assign_pilots.m](codes/assign_pilots.m) | Random pilot assignment from a pool of τ<sub>p</sub> orthonormal length-τ<sub>p</sub> sequences. Creates pilot contamination when K>τ<sub>p</sub>. |
| [codes/compute_gamma.m](codes/compute_gamma.m) | Vectorised evaluation of γ<sub>m,k</sub> = MMSE per-component variance (Eq. (5)). |
| [codes/compute_SE_closedform.m](codes/compute_SE_closedform.m) | Closed-form per-user SE of Proposition 1 / Eq. (15). Correct γ<sub>k,k′</sub> uses β<sub>m,k</sub>/β<sub>m,k′</sub> (not the square root). |
| [codes/compute_EE.m](codes/compute_EE.m) | Total energy efficiency E<sub>e</sub> = B·S<sub>e</sub>/P<sub>total</sub> from Eqs. (17)–(21), returning a breakdown into PA / circuit / fixed-backhaul / traffic-backhaul. |
| [codes/equal_power_allocation.m](codes/equal_power_allocation.m) | The "no power control" baseline of Fig. 2: η<sub>m,k</sub> = 1/(N·Σ<sub>k′</sub>γ<sub>m,k′</sub>). |

### Algorithm 1 (SCA-SOCP) + helpers

| File | Purpose |
|------|---------|
| [codes/algorithm1_EE_SCA.m](codes/algorithm1_EE_SCA.m) | Main SCA loop of Section IV. Each iteration builds and solves the SOCP (37) using CVX. Returns the optimal power-control matrix η<sub>m,k</sub>, the EE trace across iterations, and a convergence flag. |
| [codes/initial_feasible_point.m](codes/initial_feasible_point.m) | Feasibility SOCP that produces a starting point for Algorithm 1. Port of `generateinitialpoint.m`. |
| [codes/interference_vec.m](codes/interference_vec.m) | Builds the vector whose norm bounds the SINR denominator in the SOC form of constraint (27)/(36d). Works with both numeric arrays and CVX variables. Port of `interferencevector.m`. |
| [codes/approx_function.m](codes/approx_function.m) | First-order Taylor lower bound of the quadratic-over-linear term f(c,u<sub>k</sub>) (Eq. (33)) around the current tangent point. Port of `approxfunction.m`. |

### Algorithms 2 and 3 (AP selection)

| File | Purpose |
|------|---------|
| [codes/algorithm2_RP_selection.m](codes/algorithm2_RP_selection.m) | Section V-A. Given an optimal η from a PRIOR Algorithm 1 run, drops APs whose contribution to the desired-signal power is below the 95-th percentile. |
| [codes/algorithm3_LSF_selection.m](codes/algorithm3_LSF_selection.m) | Section V-B. Same idea but ranks APs directly by β<sub>m,k</sub> — cheaper because it does NOT require a prior Algorithm 1 run. |

### Figure-reproduction scripts

| File | Reproduces | Notes |
|------|-----------|-------|
| [codes/main_Fig1_convergence.m](codes/main_Fig1_convergence.m) | **Fig. 1** — monotone convergence of Algorithm 1 within ~10 iterations. | Small (M,K) grid by default; expand `MK_LIST` for the paper's (100,20)/(200,40). |
| [codes/main_Fig2_power_control.m](codes/main_Fig2_power_control.m) | **Fig. 2** — EE with/without power control vs. M for K∈{20,40}. | Expect ≈ 2.6×–3× gain from Algorithm 1 over equal-power baseline. |
| [codes/main_Fig3_AP_selection_Pbt.m](codes/main_Fig3_AP_selection_Pbt.m) | **Fig. 3** — EE vs. traffic-dependent backhaul P<sub>bt</sub> for the three schemes. | Baseline EE crashes with P<sub>bt</sub>, both selection schemes stay flat. |
| [codes/main_Fig4_num_APs_per_user.m](codes/main_Fig4_num_APs_per_user.m) | **Fig. 4** — average number of APs selected per user vs. M for D∈{1,2} km. | Sub-linear in M → scalability (~10–20 % of APs). |
| [codes/main_Fig5_effect_of_N.m](codes/main_Fig5_effect_of_N.m) | **Fig. 5** — inverted-U of EE vs. N at fixed MN for 3 scenarios. | Shows the optimal N* design trade-off. |
| [codes/main_Fig6_CF_vs_colocated.m](codes/main_Fig6_CF_vs_colocated.m) | **Fig. 6** — cell-free vs. co-located massive MIMO vs. the per-user SE target. | Cell-free dominates; gap shrinks with MN. |

---

## 3. Why this specific structure

1. **Every piece of the math lives in exactly one file.** γ, SE, EE,
   SCA iteration, AP selection — each has a single canonical
   implementation. The figure scripts are thin drivers that call
   these functions, which is how wireless-research simulation code is
   typically organised (Björnson's massive-MIMO book uses the same
   layout, and so does the reference repo).

2. **Algorithm 1 is implemented as a direct port of the published
   reference repo's Charnes–Cooper + SCA + SOCP formulation, but
   re-coded in CVX** so it runs on any machine with CVX installed
   (the original uses YALMIP + MOSEK). The constraint structure
   (36b)–(36g) is preserved line-for-line; read
   [codes/algorithm1_EE_SCA.m](codes/algorithm1_EE_SCA.m) alongside
   Section IV of the paper and each constraint maps one-to-one to
   an equation.

3. **AP selection is decoupled from Algorithm 1 via the γ̂
   substitution (Eq. (44)).** Instead of maintaining separate solvers
   for the selected and non-selected problems, you simply run
   Algorithm 1 with `gamma_hat` in place of `gamma`. This follows
   Section V's analytical trick exactly and keeps the codebase small.

4. **Parameters are Table I of the paper, encoded in exactly one
   place.** [codes/get_params.m](codes/get_params.m) returns a struct
   that every script consumes. Changing, say, the PA efficiency or the
   path-loss constant is a single edit.

5. **Runtime vs. reproduction fidelity knob.** CVX + SCA on
   (M=100, K=40) takes several minutes *per Monte-Carlo run*, so every
   figure script defaults to a smaller `(M, K, NMC)` that runs on a
   laptop in a few minutes. The comments flag the parameter sweeps
   you should expand to match the paper's figures exactly.

---

## 4. Running the code

### Prerequisites

* **MATLAB R2018a or newer**.
* **CVX** — Grant & Boyd, <http://cvxr.com/cvx/>. The free SeDuMi /
  SDPT3 solvers that ship with CVX are enough; MOSEK (CVX
  Professional) speeds things up significantly if you have an
  academic licence. *Note that the paper's SOCPs are much faster with
  MOSEK than with SDPT3 — expect 2×–10× speedup.*

### Quick test (a few minutes)

```matlab
cd 'F:\Codes\Awesome-Wireless-Communications-Paper-Implementation\Total Energy Efficiency CF-mMIMO\codes'
addpath(pwd)                 % or run 'savepath' once
cvx_setup                    % first time only
main_Fig1_convergence
```

The expected output is a figure with three monotone curves that
plateau after ≈ 10 iterations. If you see the plateau, the whole
pipeline is wired up correctly.

### Full reproduction (hours)

Each figure script has a top-of-file block of loop lists
(`M_LIST`, `NMC`, `MK_LIST`, `MN_list`, ...). To match the paper's
figures you typically need:

| Figure | Paper sizes | Default in this repo | Runtime (MOSEK) |
|--------|-------------|-----------------------|-----------------|
| Fig. 1 | `MK_LIST = [8 2; 100 20; 100 40; 200 40]` | `[8 2; 20 4; 40 10]` | ~1 min → ~5 min full |
| Fig. 2 | `M=[20:20:100]`, `K=[20 40]`, `NMC≈50` | same M, `NMC=5`    | ~15 min → ~2 h full  |
| Fig. 3 | `M=100, K=40, NMC≈50`                    | `M=40, K=20, NMC=3` | ~15 min → ~3 h full  |
| Fig. 4 | `M=[80:20:200]`, `K=40, NMC≈50`          | `M=[40:20:100], K=20, NMC=3` | ~20 min → ~3 h full |
| Fig. 5 | `MN=256, K=40, NMC≈50`                   | `MN=128, K=10, NMC=2` | ~30 min → ~6 h full |
| Fig. 6 | `MN∈{128,256}, K=20, NMC≈50`             | `MN∈{32,64}, K=10, NMC=2` | ~30 min → ~6 h full |

To jump between quick and full reproduction, edit the top of the
script, save, and re-run. The core functions do not change.

### Recommended first run sequence

1. `main_Fig1_convergence.m` — validates Algorithm 1 end-to-end.
2. `main_Fig2_power_control.m` — shows the ≈ 2.5× EE gain of power
   control, which also validates `equal_power_allocation.m` and
   `compute_EE.m`.
3. `main_Fig3_AP_selection_Pbt.m` — validates Algorithms 2 and 3 and
   the γ̂ substitution.
4. `main_Fig4_num_APs_per_user.m` — visualises the scalability claim.
5. `main_Fig5_effect_of_N.m` — observe the inverted-U in N.
6. `main_Fig6_CF_vs_colocated.m` — the headline CF vs. co-located
   comparison.

---

## 5. Implementation notes and known caveats

1. **γ<sub>k,k′</sub> uses β<sub>m,k</sub>/β<sub>m,k′</sub>, not its
   square root.** This is how the paper's Eq. (15) is written and how
   the authors' own published code implements it. A careful derivation
   from the UI<sub>k,k′</sub> moment (Appendix A) confirms this.

2. **CVX solver choice matters.** SeDuMi sometimes reports
   `Inaccurate/Solved` on the SCA sub-problems; the outer loop
   tolerates this and keeps iterating. If a figure takes forever,
   install MOSEK or reduce `sca_max_iter` in
   [codes/get_params.m](codes/get_params.m).

3. **Infeasibility handling.** If the QoS target cannot be met by any
   power allocation, `initial_feasible_point.m` returns
   `feasible = false` and the Monte-Carlo run is skipped (we do NOT
   plug in 0 as the paper sometimes does). This means the curves
   represent averages over *feasible* realisations.

4. **Wrap-around torus.** Only the base square is used to compute
   beta_dB; the 8 periodic copies enter only through the minimum
   distance. This matches the reference repo.

5. **Power units.** The code keeps B in MHz and rates in Mbit/s so
   that the backhaul coefficient P<sub>bt</sub> = 0.25 W/(Gbit/s) is
   cleanly converted to 0.25 × 10<sup>−3</sup> W/(Mbit/s). EE is
   reported in Mbit/J throughout.

6. **Single large-scale realisation per Monte-Carlo run.** The
   user-level small-scale fading cancels out by the SE closed form, so
   no inner fast-fading loop is needed. This is what makes the SCA
   tractable: it only sees deterministic β and γ.

---

## 6. Mapping from paper equations to code

| Paper equation | Implemented in |
|----------------|----------------|
| (1)–(3) pilot model | implicit in `compute_gamma.m` |
| (4) MMSE estimate   | implicit in `compute_gamma.m` |
| (5) γ<sub>m,k</sub> | `compute_gamma.m` |
| (6)–(8) DL transmit, per-AP power constraint | `equal_power_allocation.m` (baseline), and the constraint inside `algorithm1_EE_SCA.m` |
| (9)–(14) received signal decomposition | implicit (closed-form SE) |
| (15) Proposition 1: closed-form SE | `compute_SE_closedform.m` |
| (16) sum SE | `compute_SE_closedform.m` (caller sums) |
| (17)–(19) power model | `compute_EE.m` |
| (20) all-in-one P<sub>total</sub> | `compute_EE.m` |
| (21) E<sub>e</sub> | `compute_EE.m` |
| (22) Problem (P) | `algorithm1_EE_SCA.m` (formulation) |
| (24) Problem (P1) | `algorithm1_EE_SCA.m` (after divide-by-B·S<sub>e</sub> trick) |
| (25) c = √η substitution | variable `cdot` in `algorithm1_EE_SCA.m` |
| (27) SE constraint as SOC | `interference_vec.m` + norm constraint |
| (28)–(29) slack variables | `tdot`, `udot`, `theta` |
| (32)–(33) Taylor bound | `approx_function.m` |
| (35) log bound | constraint (36f) in `algorithm1_EE_SCA.m` |
| (36)/(37) full SOCP | body of `algorithm1_EE_SCA.m` |
| (40)–(45) AP selection | `algorithm2_RP_selection.m`, `algorithm3_LSF_selection.m` |
| (46)–(47) path loss + shadowing | `generate_large_scale_fading.m` |
| Table I default values | `get_params.m` |
| Appendix A proof | covered in the companion analysis `.txt` files |

---

## 7. Credits and licence

* Original paper and reference code: © Ngo, Tran, Duong, Matthaiou,
  Larsson.  Reference repo:
  <https://github.com/tranlenam/cellfreeMIMOenergyefficiency>
  (licence: see the original repo).
* This re-implementation was written as a pedagogical reproduction
  for the *Awesome-Wireless-Communications-Paper-Implementation*
  project. Use it for study and research but cite the original
  paper in any derivative work.

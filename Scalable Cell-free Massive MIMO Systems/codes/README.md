# Scalable Cell-Free Massive MIMO Systems - MATLAB Implementation

## Paper Reference

**E. Bjornson and L. Sanguinetti**, "Scalable Cell-Free Massive MIMO Systems,"
*IEEE Transactions on Communications*, vol. 68, no. 7, pp. 4247-4261, Jul. 2020.

- **DOI:** [10.1109/TCOMM.2020.2987091](https://doi.org/10.1109/TCOMM.2020.2987091)
- **arXiv:** [1908.03119](https://arxiv.org/abs/1908.03119)
- **Original code:** [github.com/emilbjornson/scalable-cell-free](https://github.com/emilbjornson/scalable-cell-free)

---

## Overview

This implementation reproduces the three simulation figures (Figures 4, 5, and 6) from the paper.
The paper proposes a scalable framework for Cell-Free Massive MIMO using Dynamic Cooperation
Clustering (DCC), where each AP serves only a bounded subset of UEs.

### Figures Reproduced

| Figure | Description | Script |
|--------|-------------|--------|
| **Fig. 4** | PDF of normalized channel gain variations (channel hardening) | `main_figure4.m` |
| **Fig. 5** | CDF of UL SE per UE (scalable vs. non-scalable combining) | `main_figure5.m` |
| **Fig. 6** | CDF of DL SE per UE (scalable precoding, hardening vs. perfect CSI) | `main_figure6.m` |

Figures 1-3 in the paper are conceptual illustrations (not simulation results).

---

## File Descriptions

### Main Scripts (run these to generate figures)

| File | Purpose |
|------|---------|
| `main_figure4.m` | **Figure 4**: Computes the PDF of `\|h_k^H w_k\|^2 / E{\|h_k^H w_k\|^2}` for MR and LP-MMSE precoding, showing channel hardening. LP-MMSE concentrates near 1 (strong hardening); MR has a long tail (weak hardening). |
| `main_figure5.m` | **Figure 5**: Compares UL SE CDFs for 5 schemes: scalable P-MMSE, scalable LP-MMSE, non-scalable MMSE (All), non-scalable L-MMSE (All), and non-scalable MR (All). Two subplots: (a) L=400, N=1 and (b) L=100, N=4. |
| `main_figure6.m` | **Figure 6**: Compares DL SE CDFs for 3 scalable schemes (P-MMSE, LP-MMSE, MR), each with hardening bound and perfect-CSI bound. Shows that P-MMSE achieves 98% of the genie-aided SE. |

### Helper Functions

| File | Purpose |
|------|---------|
| `generateSetup.m` | Generates random network layouts with L APs and K UEs in a 2x2 km area with wrap-around. Implements the **three-step algorithm** (Section V-A) for: (1) Master AP selection, (2) pilot assignment minimizing pilot contamination, (3) DCC cluster formation with threshold-based AP inclusion. |
| `functionChannelEstimates.m` | Generates correlated Rayleigh fading channels `h_{kl} ~ CN(0, R_{kl})` and computes **MMSE channel estimates** (Eq. 3-4). Returns channel estimates `Hhat`, true channels `H`, estimate correlation `B`, and error correlation `C`. |
| `functionComputeSE_uplink.m` | Computes UL SE for four combining schemes: **MR** (closed-form UatF bound), **LP-MMSE** (UatF bound via Monte Carlo), **P-MMSE** (instantaneous SINR, Prop. 1), and **MMSE** (non-scalable benchmark). |
| `functionComputeSE_downlink.m` | Computes DL SE for three precoding schemes (**MR**, **LP-MMSE**, **P-MMSE**), each with both the **hardening bound** (Prop. 3) and a **perfect-CSI** genie-aided bound. Uses two-pass normalization for correct precoding power. |
| `functionRlocalscattering.m` | Generates the `N x N` spatial correlation matrix for a ULA using the **local scattering model** with Gaussian angular distribution and a specified angular standard deviation (ASD). |

---

## Simulation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Coverage area | 2 x 2 km | Square area with wrap-around |
| Setup (a) | L=400, N=1 | 400 single-antenna APs |
| Setup (b) | L=100, N=4 | 100 four-antenna APs |
| K | 100 | Number of UEs |
| tau_c | 200 | Coherence block length (samples) |
| tau_p | 10 | Pilot sequence length |
| p_k | 100 mW | UL transmit power per UE |
| rho | 1 W (1000 mW) | DL transmit power per AP |
| Bandwidth | 20 MHz | Communication bandwidth |
| Noise figure | 7 dB | At the APs |
| Pathloss exponent | 3.76 | Three-slope model from [2, Sec. 4.1.3] |
| Shadow fading | 10 dB std | Log-normal |
| AP height | 10 m | Above UE level |
| ASD | 20 degrees | Angular standard deviation |
| Cluster threshold | -40 dB | AP serves UE if gain within 40 dB of master |
| nbrOfSetups | 25 | Random network realizations |
| nbrOfRealizations | 1000 | Channel realizations per setup |

---

## How to Run

### Prerequisites
- **MATLAB** R2019a or later (uses `histcounts`, `xline`, `linspace`)
- No additional toolboxes required

### Quick Test (reduced parameters)
For a quick test run (~5-10 minutes per figure), edit the main scripts and reduce:
```matlab
nbrOfSetups = 5;         % Instead of 25
nbrOfRealizations = 100; % Instead of 1000
```

### Full Reproduction (matches paper)
1. Open MATLAB and navigate to the `codes` folder:
   ```matlab
   cd('F:\Codes\Awesome-Wireless-Communications-Paper-Implementation\Scalable Cell-free Massive MIMO Systems\codes')
   ```

2. Run each figure script:
   ```matlab
   % Figure 4: Channel hardening PDF
   main_figure4

   % Figure 5: UL SE CDF comparison
   main_figure5

   % Figure 6: DL SE CDF comparison
   main_figure6
   ```

3. Each script generates two subplot figures (one per setup) and prints summary statistics.

### Expected Runtime
- **Figure 4**: ~30-60 minutes (10 setups x 1000 realizations)
- **Figure 5**: ~2-4 hours (25 setups x 1000 realizations, includes MMSE matrix inversions)
- **Figure 6**: ~3-6 hours (25 setups x 1000 realizations, two-pass precoding)

Runtime depends on hardware. The P-MMSE and MMSE schemes require matrix inversions at each realization, which is the computational bottleneck.

---

## Key Results Expected

### Figure 5 (Uplink)
- Scalable P-MMSE achieves ~89% of the average SE of non-scalable MMSE (All)
- Distributed LP-MMSE average SE is ~2.7x higher than MR (All)
- LP-MMSE (scalable) matches L-MMSE (All) with negligible gap

### Figure 6 (Downlink)
- P-MMSE hardening bound achieves ~98% of the perfect-CSI genie-aided SE
- LP-MMSE achieves ~90% of genie-aided SE
- MR only achieves ~60% of genie-aided SE (weak channel hardening)
- Centralized P-MMSE outperforms distributed LP-MMSE despite using 40x less total power

---

## Mathematical Key Equations

- **Channel estimation (Eq. 3):** `hhat_{kl} = sqrt(p_k * tau_p) * R_{kl} * Psi_{tl}^{-1} * y_pilot`
- **LP-MMSE combining (Eq. 29):** `v_{kl} = p_k * (sum_{i in D_l} p_i*(hhat_{il}*hhat_{il}' + C_{il}) + I)^{-1} * hhat_{kl}`
- **P-MMSE combining (Eq. 23):** Like MMSE but restricted to UEs with overlapping clusters
- **DL precoding (Eq. 33):** `w_i = v_i / sqrt(E{||v_i||^2})` (UL-DL duality)
- **DL power allocation (Eq. 43):** `rho_{kl} = rho * sqrt(beta_{kl}) / sum_j sqrt(beta_{jl})`

---

## References

[1] E. Bjornson, L. Sanguinetti, "Scalable Cell-Free Massive MIMO Systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247-4261, 2020.

[2] E. Bjornson, J. Hoydis, L. Sanguinetti, "Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency," Foundations and Trends in Signal Processing, vol. 11, no. 3-4, pp. 154-655, 2017.

[5] H. Q. Ngo et al., "Cell-Free Massive MIMO Versus Small Cells," IEEE Trans. Wireless Commun., vol. 16, no. 3, pp. 1834-1850, 2017.

[7] E. Bjornson, L. Sanguinetti, "Making Cell-Free Massive MIMO Competitive With MMSE Processing," IEEE Trans. Wireless Commun., vol. 19, no. 1, pp. 77-90, 2020.

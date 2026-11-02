# Noise-Limited Phase Estimation with Entangled Probes: A Qiskit Aer Study

This project compares GHZ and product probes for phase estimation under local phase-damping noise. All results come from Qiskit Aer simulations. No hardware is used.

## Question

GHZ probes can in principle estimate a phase with precision scaling as 1/N, compared with 1/√N for product probes. We ask how the GHZ advantage changes with noise strength and probe size.

## Method

A phase φ is encoded on N qubits in one of two ways:

- **Product probe:** each qubit is prepared in |+⟩, phased, and measured in the X basis independently. The signal is the qubit-averaged probability of outcome 0.
- **GHZ probe:** the qubits are entangled into a GHZ state, phased collectively, then disentangled and measured. The signal is the probability of the all-zeros outcome.

Noise is local phase damping with strength λ, applied to each qubit's phase gate. For each probe and each N, the phase is scanned over 50 points evenly spaced in [0, 2π), with 20000 shots per point. The signal is fitted by linear least squares to a + b cos(kφ) + c sin(kφ), with k = 1 for the product probe and k = N for the GHZ probe.

The classical Fisher information comes from the fitted visibility V = 2√(b² + c²), clipped to at most 1:

F_single = k²(1 − √(1 − V²))

The product probe's Fisher information is N times the single-qubit value. The phase uncertainty per shot is Δφ = 1/√F. Each measured Fisher value carries sampling scatter of a few percent.

The reported Fisher values are classical Fisher information from a fixed X-basis readout. They are not the quantum Fisher information. The QFI is computed separately in closed form for this channel: F_Q = N(1−λ) for the product probe and F_Q = N²(1−λ)^N for the GHZ probe. These closed forms are checked against a full density-matrix calculation for N = 2, 3, 4 at each λ tested.

## Validation

- Single-qubit Ramsey fringes match cos²(φ/2). The noisy fringe contrast matches 0.5·√(1−λ) at λ = 0.1, 0.3, and 0.5.
- Noiseless product Fisher information matches N within about 6 percent for all N from 2 to 8.
- Noiseless GHZ Fisher information matches N² exactly at N = 2, 3, 4, 6, and 7, and is within about 5 percent at N = 5 and 8.
- Repeated noiseless GHZ runs at N = 5 (8 runs, 20000 shots each) give a mean of about 24.6 against the target of 25, with a standard deviation of about 0.5. Five of the eight runs are capped at 25 by the visibility clip.
- Numerical QFI matches the closed forms for N = 2, 3, 4 at every λ tested.

## Repository contents

| File | Purpose |
|------|---------|
| `probes.py` | Circuits, noise model, estimator, and the full sweep. Writes `scaling_results.csv`. |
| `check_single.py` | Single-qubit and two-qubit Fisher check, plus noiseless GHZ fringe amplitude at N = 5. |
| `repeat_ghz.py` | Repeated noiseless GHZ runs at N = 5, to test for estimator bias. |
| `check_fisher.py` | Compares measured and ideal GHZ fringes at N = 5 using a finite-difference Fisher estimate. |
| `plots.py` | Produces the classical Fisher curves, GHZ-to-product ratio, and GHZ slope figures. |
| `qfi.py` | Computes the QFI bounds, checks them numerically, writes `qfi_bounds.csv`, and produces `fisher_vs_qfi.png`. |
| `scaling_results.csv` | Measured classical Fisher information, with columns lambda, probe, N, fisher, dphi_per_shot. |
| `qfi_bounds.csv` | Closed-form QFI, with columns lambda, probe, N, qfi. |

## Requirements

- Python 3.13.5 (tested)
- qiskit 2.6.1 and qiskit-aer 0.18.2 (tested)
- numpy, scipy, matplotlib

Install with: pip install qiskit qiskit-aer matplotlib numpy scipy

## How to run

The main pipeline is:

1. `python probes.py` (full sweep; writes `scaling_results.csv`)
2. `python plots.py` (classical Fisher figures)
3. `python qfi.py` (QFI bounds and comparison figure)

Optional diagnostics:

- `python check_single.py` (verifies the single-qubit and GHZ pipeline; not needed to reproduce the results)
- `python repeat_ghz.py`
- `python check_fisher.py`

`probes.py` is a long run. It executes 3500 circuits at 20000 shots each (5 noise levels × 2 probes × 7 values of N × 50 phases). The CSV is rewritten after each noise level, so partial results survive an interruption.

## Results

**Scaling exponents.** The slope of Δφ versus N, fitted over N = 2 to 8:

| λ | Product slope | GHZ slope |
|---|---|---|
| 0.00 | −0.49 | −0.995 |
| 0.05 | −0.50 | −0.82 |
| 0.10 | −0.50 | −0.70 |
| 0.20 | −0.50 | −0.46 |
| 0.30 | −0.50 | −0.18 (not a power law; see Limitations) |

**Product probe.** The measured classical Fisher information is about 1/(1 + √λ) times the QFI at every λ tested, to within about 1 to 2 percent.

**GHZ probe.** At N = 8 the measured GHZ Fisher information is about 0.53 of its QFI bound at λ = 0.2 and about 0.50 at λ = 0.3.

**Crossover.**
- At λ = 0.3, the measured GHZ curve is above the product curve at N = 4 (2.05 against 1.81) and below it at N = 5 (2.22 against 2.26). The measured GHZ curve peaks at N = 5 and then declines. The QFI ratio N(1−λ)^(N−1) falls below 1 between N = 6 and N = 7.
- At λ = 0.2, no crossover is observed up to N = 8. The measured GHZ-to-product ratio is about 1.29 at N = 8. The QFI ratio is 8(0.8)⁷ ≈ 1.68 at the same N.
- At λ ≤ 0.1, the measured GHZ curve stays above the product curve throughout the range studied.

## Limitations

- Probe sizes only go up to N = 8, so any crossover at λ = 0.2 lies outside the data.
- One noise model (local, uncorrelated phase damping) and one readout (X basis). The results describe this readout and do not characterize optimal measurements for either probe.
- Classical Fisher information is computed from the fitted visibility. It carries a few percent sampling scatter per point and a small low bias from the visibility cap, up to about 5 percent at some N.
- The GHZ curve at λ = 0.3 is non-monotonic, so its fitted slope is a summary, not a power law.

## References

- Huelga, S. F., et al. (1997). Improvement of frequency standards with quantum entanglement. *Physical Review Letters*, 79, 3865.
- Demkowicz-Dobrzański, R., Kołodyński, J., & Guţă, M. (2012). The elusive Heisenberg limit in quantum-enhanced metrology. *Nature Communications*, 3, 1063.
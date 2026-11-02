import numpy as np
from qiskit import transpile
from qiskit_aer import AerSimulator
from probes import (make_noise, run_probe, ghz_circuit, fit_fringe,
                    PHASES, BASIS, SHOTS)

# Single-qubit and two-qubit product Fisher information with and without noise
for lam in [0.0, 0.1]:
    backend = AerSimulator(noise_model=make_noise(lam))
    F1 = run_probe(backend, "product", 1)
    F2 = run_probe(backend, "product", 2)
    print(f"lambda={lam:.2f}  F(N=1)={F1:.3f}  F(N=2)={F2:.3f}  "
          f"expected N=1 ~ {1 - lam:.3f}  expected N=2 ~ {2 * (1 - lam):.3f}")

# Noiseless GHZ fringe amplitude at N=5
backend = AerSimulator()
n = 5
s = []
for phi in PHASES:
    tqc = transpile(ghz_circuit(n, phi), basis_gates=BASIS, optimization_level=0)
    counts = backend.run(tqc, shots=SHOTS).result().get_counts()
    s.append(counts.get("0" * n, 0) / SHOTS)
coef = fit_fringe(np.array(s), k=n)
print(f"GHZ N=5 fit a,b,c = {np.round(coef, 4)}  amplitude = {np.hypot(coef[1], coef[2]):.4f}")
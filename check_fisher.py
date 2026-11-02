import numpy as np
from qiskit import transpile
from qiskit_aer import AerSimulator
from probes import ghz_circuit, BASIS, SHOTS, PHASES

n = 5
backend = AerSimulator()

s = []
for phi in PHASES:
    tqc = transpile(ghz_circuit(n, phi), basis_gates=BASIS, optimization_level=0)
    counts = backend.run(tqc, shots=SHOTS).result().get_counts()
    s.append(counts.get("0" * n, 0) / SHOTS)
s = np.array(s)

ideal = 0.5 * (1 + np.cos(n * PHASES))
dphi = PHASES[1] - PHASES[0]
mask = (ideal > 0.05) & (ideal < 0.95)   # same phase points for both curves

def fisher_fd(p):
    dp = (np.roll(p, -1) - np.roll(p, 1)) / (2 * dphi)
    p = np.clip(p, 1e-3, 1 - 1e-3)
    return np.mean(dp[mask] ** 2 / (p[mask] * (1 - p[mask])))

print(f"max |data - ideal| = {np.max(np.abs(s - ideal)):.4f}")
print(f"Fisher from data  = {fisher_fd(s):.3f}")
print(f"Fisher from ideal = {fisher_fd(ideal):.3f}")
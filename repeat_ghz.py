import numpy as np
from qiskit_aer import AerSimulator
from probes import run_probe

backend = AerSimulator()
vals = [run_probe(backend, "ghz", 5) for _ in range(8)]
print("GHZ N=5 noiseless Fisher over 8 runs")
print("values =", np.round(vals, 3))
print("mean =", round(np.mean(vals), 3), " std =", round(np.std(vals), 3), " target = 25")
import csv
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, phase_damping_error

SHOTS = 20000
PHASES = np.linspace(0, 2 * np.pi, 50, endpoint=False)
NS = [2, 3, 4, 5, 6, 7, 8]
LAMBDAS = [0.0, 0.05, 0.1, 0.2, 0.3]
BASIS = ["h", "p", "cx", "measure"]


def make_noise(lam):
    noise = NoiseModel()
    if lam > 0:
        noise.add_all_qubit_quantum_error(phase_damping_error(lam), ["p"])
    return noise


def ghz_circuit(n, phi):
    qc = QuantumCircuit(n, n)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)          # create GHZ
    for i in range(n):
        qc.p(phi, i)             # encode phase on every qubit
    for i in reversed(range(n - 1)):
        qc.cx(i, i + 1)          # undo entangling
    qc.h(0)
    qc.measure(range(n), range(n))
    return qc


def product_circuit(n, phi):
    qc = QuantumCircuit(n, n)
    for i in range(n):
        qc.h(i)
        qc.p(phi, i)
        qc.h(i)
    qc.measure(range(n), range(n))
    return qc


def ghz_signal(counts, n):
    # probability of all-zeros, which oscillates as cos^2(N phi / 2)
    return counts.get("0" * n, 0) / SHOTS


def product_signal(counts, n):
    # average per-qubit probability of 0, which oscillates as cos^2(phi / 2)
    zeros = np.zeros(n)
    for key, c in counts.items():
        val = int(key, 2)
        for i in range(n):
            if not (val >> i) & 1:
                zeros[i] += c
    return zeros.mean() / SHOTS


def fit_fringe(s, k):
    X = np.column_stack([np.ones_like(PHASES), np.cos(k * PHASES), np.sin(k * PHASES)])
    coef, *_ = np.linalg.lstsq(X, s, rcond=None)
    return coef  # a, b, c for model a + b cos(k phi) + c sin(k phi)


def fisher_from_fit(coef, k):
    a, b, c = coef
    V = min(2 * np.hypot(b, c), 1.0)           # fitted fringe visibility, clipped to physical range
    return k**2 * (1 - np.sqrt(max(0.0, 1 - V**2)))


def run_probe(backend, name, n):
    s_vals = []
    for phi in PHASES:
        if name == "ghz":
            qc = ghz_circuit(n, phi)
        else:
            qc = product_circuit(n, phi)
        tqc = transpile(qc, basis_gates=BASIS, optimization_level=0)
        counts = backend.run(tqc, shots=SHOTS).result().get_counts()
        s_vals.append(ghz_signal(counts, n) if name == "ghz" else product_signal(counts, n))
    s_vals = np.array(s_vals)

    if name == "ghz":
        coef = fit_fringe(s_vals, k=n)
        F = fisher_from_fit(coef, k=n)
    else:
        coef = fit_fringe(s_vals, k=1)
        F = n * fisher_from_fit(coef, k=1)     # independent qubits add Fisher information
    return F


if __name__ == "__main__":
    rows = []
    for lam in LAMBDAS:
        backend = AerSimulator(noise_model=make_noise(lam))
        for name in ["product", "ghz"]:
            dphis = []
            for n in NS:
                F = run_probe(backend, name, n)
                dphi = 1 / np.sqrt(F)
                dphis.append(dphi)
                rows.append([lam, name, n, F, dphi])
                print(f"lambda={lam:.2f}  {name:8s}  N={n}  F={F:8.3f}  dphi/shot={dphi:.4f}")
            slope = np.polyfit(np.log(NS), np.log(dphis), 1)[0]
            print(f"  --> lambda={lam:.2f} {name} log-log slope of dphi vs N: {slope:.3f}\n")

        with open("scaling_results.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["lambda", "probe", "N", "fisher", "dphi_per_shot"])
            w.writerows(rows)
        print(f"saved after lambda={lam:.2f}\n")
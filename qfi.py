import csv
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce

LAMBDAS = [0.0, 0.05, 0.1, 0.2, 0.3]
NS = list(range(2, 9))
NUM_CHECK_NS = [2, 3, 4]


def kraus_pd(lam):
    return [np.diag([1.0, np.sqrt(1 - lam)]), np.diag([0.0, np.sqrt(lam)])]


def apply_channel(rho, lam, n):
    if lam == 0:
        return rho
    Ks = kraus_pd(lam)
    for q in range(n):
        out = np.zeros_like(rho)
        for K in Ks:
            ops = [np.eye(2)] * n
            ops[q] = K
            full = reduce(np.kron, ops)
            out += full @ rho @ full.conj().T
        rho = out
    return rho


def ghz_state(n):
    psi = np.zeros(2**n, dtype=complex)
    psi[0] = psi[-1] = 1 / np.sqrt(2)
    return psi


def plus_state(n):
    return np.ones(2**n, dtype=complex) / np.sqrt(2**n)


def phase_unitary(n, phi):
    d = np.diag([1.0, np.exp(1j * phi)])
    return reduce(np.kron, [d] * n)


def rho_of_phi(psi0, n, lam, phi):
    psi = phase_unitary(n, phi) @ psi0
    return apply_channel(np.outer(psi, psi.conj()), lam, n)


def qfi_numeric(psi0, n, lam, phi=0.3, d=1e-5):
    rho = rho_of_phi(psi0, n, lam, phi)
    drho = (rho_of_phi(psi0, n, lam, phi + d) - rho_of_phi(psi0, n, lam, phi - d)) / (2 * d)
    p, V = np.linalg.eigh(rho)
    dV = V.conj().T @ drho @ V
    F = 0.0
    for i in range(len(p)):
        for j in range(len(p)):
            s = p[i] + p[j]
            if s > 1e-10:
                F += 2 * abs(dV[i, j]) ** 2 / s
    return F.real


def qfi_analytic(probe, n, lam):
    if probe == "product":
        return n * (1 - lam)
    return n**2 * (1 - lam) ** n


if __name__ == "__main__":
    print("Numerical check of QFI against closed forms (N = 2, 3, 4)")
    for lam in LAMBDAS:
        for n in NUM_CHECK_NS:
            Fp = qfi_numeric(plus_state(n), n, lam)
            Fg = qfi_numeric(ghz_state(n), n, lam)
            print(f"lambda={lam:.2f} N={n}  product num={Fp:.4f} ana={qfi_analytic('product', n, lam):.4f}"
                  f"  |  ghz num={Fg:.4f} ana={qfi_analytic('ghz', n, lam):.4f}")

    rows = []
    for lam in LAMBDAS:
        for n in NS:
            rows.append([lam, "product", n, qfi_analytic("product", n, lam)])
            rows.append([lam, "ghz", n, qfi_analytic("ghz", n, lam)])
    with open("qfi_bounds.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lambda", "probe", "N", "qfi"])
        w.writerows(rows)

    measured = {}
    with open("scaling_results.csv") as f:
        for row in csv.DictReader(f):
            measured[(float(row["lambda"]), row["probe"], int(row["N"]))] = float(row["fisher"])

    fig, axes = plt.subplots(1, len(LAMBDAS), figsize=(3.2 * len(LAMBDAS), 3.6), sharey=False)
    for ax, lam in zip(axes, LAMBDAS):
        for probe, color in [("product", "C0"), ("ghz", "C1")]:
            ax.plot(NS, [measured[(lam, probe, n)] for n in NS], "o-", color=color,
                    label=f"{probe} measured")
            ax.plot(NS, [qfi_analytic(probe, n, lam) for n in NS], "--", color=color,
                    label=f"{probe} QFI bound")
        ax.set_title(f"λ = {lam:.2f}")
        ax.set_xlabel("N")
        ax.set_ylabel("Fisher information")
        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig("fisher_vs_qfi.png", dpi=200)
    print("saved qfi_bounds.csv and fisher_vs_qfi.png")
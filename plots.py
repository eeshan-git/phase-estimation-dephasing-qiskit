import csv
import numpy as np
import matplotlib.pyplot as plt

rows = []
with open("scaling_results.csv") as f:
    for row in csv.DictReader(f):
        rows.append((float(row["lambda"]), row["probe"], int(row["N"]), float(row["fisher"])))

lambdas = sorted({x[0] for x in rows})
Ns = sorted({x[2] for x in rows})


def get(lam, probe):
    d = {n: F for l, p, n, F in rows if l == lam and p == probe}
    return np.array([d[n] for n in Ns])


# Figure 1: classical Fisher information vs N, one panel per lambda
fig, axes = plt.subplots(1, len(lambdas), figsize=(3.2 * len(lambdas), 3.4))
for ax, lam in zip(axes, lambdas):
    ax.loglog(Ns, get(lam, "product"), "o-", label="product")
    ax.loglog(Ns, get(lam, "ghz"), "s-", label="GHZ")
    ax.set_title(f"λ = {lam:.2f}")
    ax.set_xlabel("N")
    ax.set_ylabel("classical Fisher info")
    ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig("fisher_vs_N.png", dpi=200)

# Figure 2: GHZ-to-product Fisher ratio vs N
plt.figure(figsize=(5, 3.8))
for lam in lambdas:
    ratio = get(lam, "ghz") / get(lam, "product")
    plt.plot(Ns, ratio, "o-", label=f"λ = {lam:.2f}")
plt.axhline(1, color="k", lw=0.8, ls="--")
plt.xlabel("N")
plt.ylabel("F_GHZ / F_product")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("ghz_over_product.png", dpi=200)

# Figure 3: fitted GHZ log-log slope vs lambda
slopes = []
for lam in lambdas:
    dphi = 1 / np.sqrt(get(lam, "ghz"))
    slopes.append(np.polyfit(np.log(Ns), np.log(dphi), 1)[0])

plt.figure(figsize=(5, 3.8))
plt.plot(lambdas, slopes, "o-")
plt.axhline(-1, color="k", lw=0.8, ls="--", label="Heisenberg (-1)")
plt.axhline(-0.5, color="gray", lw=0.8, ls=":", label="SQL (-0.5)")
plt.xlabel("λ")
plt.ylabel("GHZ log-log slope of Δφ vs N")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("ghz_slope_vs_lambda.png", dpi=200)

print("saved fisher_vs_N.png, ghz_over_product.png, ghz_slope_vs_lambda.png")
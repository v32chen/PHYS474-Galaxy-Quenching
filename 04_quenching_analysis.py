#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt

# Load the combined galaxy measurements
with np.load("galaxy_combined_preclean.npz") as d:
    log_mass = d["log_stellar_mass"].copy()
    log_ssfr = d["log_ssfr"].copy()
    z = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()

# Load the quality mask
clean = np.load("galaxy_quality_mask.npy")

# Keep the same exploratory redshift range
sample = clean & (z >= 0.02) & (z < 0.08)

# Quenched galaxy classification
quenched = log_ssfr < -11.0

# Five narrower mass bins
mass_bins = np.arange(10.0, 10.51, 0.1)
mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2

# Wilson 95% confidence interval
def wilson_interval(k, n):
    p = k / n
    t = 1.96
    denom = 1 + t**2 / n
    center = (p + t**2 / (2*n)) / denom
    margin = t * np.sqrt(
        p*(1-p)/n + t**2/(4*n**2)
    ) / denom
    return center - margin, center + margin

results = {}

for name, classification in [
    ("Centrals", p_sat < 0.5),
    ("Satellites", p_sat > 0.5)
]:
    fractions = []
    errors_low = []
    errors_high = []

    print(f"\n{name}")
    print("Mass bin | Total | Quenched | f_Q")

    for lo, hi in zip(mass_bins[:-1], mass_bins[1:]):
        selected = (
            sample &
            classification &
            (log_mass >= lo) &
            (log_mass < hi)
        )

        n = np.count_nonzero(selected)
        k = np.count_nonzero(selected & quenched)

        if n == 0:
            fractions.append(np.nan)
            errors_low.append(np.nan)
            errors_high.append(np.nan)
            continue

        f = k / n
        lower, upper = wilson_interval(k, n)

        fractions.append(f)
        errors_low.append(f - lower)
        errors_high.append(upper - f)

        print(f"{lo:.1f}-{hi:.1f} | {n} | {k} | {f:.4f}")

    results[name] = (
        np.array(fractions),
        np.array([errors_low, errors_high])
    )

# Create the scientific figure
plt.figure(figsize=(9, 6))

for name, marker in [
    ("Centrals", "o"),
    ("Satellites", "s")
]:
    fractions, errors = results[name]

    plt.errorbar(
        mass_centers,
        fractions,
        yerr=errors,
        marker=marker,
        capsize=4,
        linewidth=1.5,
        label=name
    )

plt.xlabel(r"$\log_{10}(M_*/M_\odot)$")
plt.ylabel("Quenched fraction")
plt.ylim(0, 1)
plt.legend()
plt.grid(alpha=0.2)

plt.tight_layout()
plt.savefig(
    "quenched_fraction_vs_stellar_mass_exploratory.pdf",
    bbox_inches="tight"
)
plt.show()


# In[ ]:





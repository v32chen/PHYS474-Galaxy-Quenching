#!/usr/bin/env python
# coding: utf-8

# In[1]:
import numpy as np
import matplotlib.pyplot as plt

# Load the combined dataset
with np.load("galaxy_combined_preclean.npz") as d:
    log_mass = d["log_stellar_mass"].copy()
    log_ssfr = d["log_ssfr"].copy()
    redshift = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()
    halo_mass = d["halo_mass"].copy()

# Load the quality mask
clean = np.load("galaxy_quality_mask.npy")

# Match the stellar-mass and redshift ranges
base = (
    clean &
    (redshift >= 0.02) &
    (redshift < 0.08) &
    (log_mass >= 10.2) &
    (log_mass < 10.3)
)

quenched = log_ssfr < -11.0

centrals = base & (p_sat < 0.5)
satellites = base & (p_sat > 0.5)

# Central baseline at this stellar mass
f_central = (
    np.count_nonzero(centrals & quenched) /
    np.count_nonzero(centrals)
)

print("Central quenched fraction:", round(f_central, 4))
print("Satellites in selected mass bin:",
      np.count_nonzero(satellites))

# Halo mass in log10(M_sun/h)
log_halo = np.full(len(halo_mass), np.nan)
valid_halo = np.isfinite(halo_mass) & (halo_mass > 0)
log_halo[valid_halo] = np.log10(halo_mass[valid_halo])

halo_bins = np.arange(11.5, 15.51, 0.5)
centers = (halo_bins[:-1] + halo_bins[1:]) / 2

fractions = []
errors_low = []
errors_high = []

print("\nSATELLITE QUENCHING BY HALO MASS")
print("Halo bin | Total | Quenched | f_Q")

for lo, hi in zip(halo_bins[:-1], halo_bins[1:]):
    selected = (
        satellites &
        (log_halo >= lo) &
        (log_halo < hi)
    )

    n = np.count_nonzero(selected)
    k = np.count_nonzero(selected & quenched)

    # Skip bins with too few galaxies
    if n < 30:
        fractions.append(np.nan)
        errors_low.append(np.nan)
        errors_high.append(np.nan)
        print(f"{lo:.1f}-{hi:.1f} | {n} | skipped")
        continue

    f = k / n

    # Wilson 95% confidence interval
    t = 1.96
    denom = 1 + t**2 / n
    mid = (f + t**2 / (2*n)) / denom
    margin = t * np.sqrt(
        f*(1-f)/n + t**2/(4*n**2)
    ) / denom

    lower = mid - margin
    upper = mid + margin

    fractions.append(f)
    errors_low.append(f - lower)
    errors_high.append(upper - f)

    print(f"{lo:.1f}-{hi:.1f} | {n} | {k} | {f:.4f}")

# Plot results
plt.figure(figsize=(9, 6))

plt.errorbar(
    centers,
    fractions,
    yerr=[errors_low, errors_high],
    fmt="o-",
    capsize=4,
    label="Satellite galaxies"
)

plt.axhline(
    f_central,
    linestyle="--",
    label="Central baseline (same stellar-mass bin)"
)

plt.xlabel(r"$\log_{10}[M_{\rm halo}/(M_\odot/h)]$")
plt.ylabel("Quenched fraction")
plt.ylim(0, 1)

plt.legend()
plt.grid(alpha=0.2)
plt.tight_layout()

plt.savefig(
    "satellite_quenching_vs_halo_mass_exploratory.pdf",
    bbox_inches="tight"
)

plt.show()


# In[2]:
import numpy as np

# Keep the same stellar-mass selection
mass_selected = (
    clean &
    (log_mass >= 10.2) &
    (log_mass < 10.3)
)

# Two redshift intervals
redshift_ranges = [
    (0.02, 0.05),
    (0.05, 0.08)
]

# Same host-halo mass bins
halo_edges = np.arange(11.5, 15.01, 0.5)

for z_low, z_high in redshift_ranges:

    print(f"\nREDSHIFT: {z_low:.2f} - {z_high:.2f}")
    print("--------------------------------")

    z_selected = (
        mass_selected &
        (redshift >= z_low) &
        (redshift < z_high)
    )

    cen = z_selected & (p_sat < 0.5)
    sat = z_selected & (p_sat > 0.5)

    n_cen = np.count_nonzero(cen)
    n_sat = np.count_nonzero(sat)

    print("Central galaxies:", n_cen)
    print("Satellite galaxies:", n_sat)

    if n_cen > 0:
        f_cen = np.mean(log_ssfr[cen] < -11)
        print(f"Central quenched fraction: {f_cen:.3f}")

    print("\nHalo mass bin | N satellites | f_Q")

    for lo, hi in zip(
        halo_edges[:-1], halo_edges[1:]
    ):
        selected = (
            sat &
            (log_halo >= lo) &
            (log_halo < hi)
        )

        n = np.count_nonzero(selected)

        if n < 30:
            print(f"{lo:.1f}-{hi:.1f} | {n} | Too few")
            continue

        f_q = np.mean(log_ssfr[selected] < -11)

        print(f"{lo:.1f}-{hi:.1f} | {n} | {f_q:.3f}")


# In[3]:
import numpy as np
import csv

# Load our combined measurements
with np.load("galaxy_combined_preclean.npz") as d:
    log_mass = d["log_stellar_mass"].copy()
    log_ssfr = d["log_ssfr"].copy()
    redshift = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()
    halo_mass = d["halo_mass"].copy()

clean = np.load("galaxy_quality_mask.npy")

# Same exploratory redshift selection
base = clean & (redshift >= 0.02) & (redshift < 0.08)
quenched = log_ssfr < -11.0

# Five stellar-mass bins, 0.1 dex wide
mass_edges = np.round(np.arange(10.0, 10.51, 0.1), 2)

# Seven host-halo mass bins, 0.5 dex wide
halo_edges = np.arange(11.5, 15.01, 0.5)

log_halo = np.full(len(halo_mass), np.nan)
valid = np.isfinite(halo_mass) & (halo_mass > 0)
log_halo[valid] = np.log10(halo_mass[valid])

shape = (len(mass_edges)-1, len(halo_edges)-1)

sat_fraction = np.full(shape, np.nan)
efficiency = np.full(shape, np.nan)
sat_counts = np.zeros(shape, dtype=int)
central_fraction = np.full(shape[0], np.nan)

rows = []

for i, (m_lo, m_hi) in enumerate(zip(
    mass_edges[:-1], mass_edges[1:]
)):
    # Hold stellar mass approximately fixed
    in_mass = (
        base &
        (log_mass >= m_lo) &
        (log_mass < m_hi)
    )

    central = in_mass & (p_sat < 0.5)
    satellite = in_mass & (p_sat > 0.5)

    n_cen = np.count_nonzero(central)

    if n_cen < 30:
        continue

    f_cen = np.mean(quenched[central])
    central_fraction[i] = f_cen

    print(f"\nSTELLAR MASS: {m_lo:.1f}-{m_hi:.1f}")
    print(f"Centrals: {n_cen}, f_Q = {f_cen:.3f}")
    print("Halo mass | N sat | f_Q sat | Efficiency")

    for j, (h_lo, h_hi) in enumerate(zip(
        halo_edges[:-1], halo_edges[1:]
    )):
        selected = (
            satellite &
            (log_halo >= h_lo) &
            (log_halo < h_hi)
        )

        n = np.count_nonzero(selected)
        sat_counts[i, j] = n

        if n < 30 or f_cen >= 1:
            print(f"{h_lo:.1f}-{h_hi:.1f} | {n} | skipped")
            continue

        f_sat = np.mean(quenched[selected])

        eps = (f_sat - f_cen) / (1 - f_cen)

        sat_fraction[i, j] = f_sat
        efficiency[i, j] = eps

        rows.append([
            m_lo, m_hi, h_lo, h_hi,
            n_cen, n, f_cen, f_sat, eps
        ])

        print(
            f"{h_lo:.1f}-{h_hi:.1f} | {n:5d} | "
            f"{f_sat:.3f} | {eps:.3f}"
        )

# Save numerical results for the report
with open("quenching_by_mass_and_halo.csv",
          "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "stellar_mass_low", "stellar_mass_high",
        "halo_mass_low", "halo_mass_high",
        "n_central", "n_satellite",
        "fQ_central", "fQ_satellite",
        "satellite_quenching_efficiency"
    ])
    writer.writerows(rows)

print("\nAnalysis completed and results saved!")


# In[5]:
import matplotlib.pyplot as plt

mass_centers = (mass_edges[:-1] + mass_edges[1:]) / 2
halo_centers = (halo_edges[:-1] + halo_edges[1:]) / 2

plots = [
    (
        sat_fraction,
        "Satellite quenched fraction",
        "satellite_fraction_mass_halo.pdf"
    ),
    (
        efficiency,
        "Satellite quenching efficiency",
        "satellite_efficiency_mass_halo.pdf"
    )
]

for values, colorbar_label, filename in plots:

    plt.figure(figsize=(10, 6))

    image = plt.imshow(
        np.ma.masked_invalid(values),
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        extent=[
            halo_edges[0], halo_edges[-1],
            mass_edges[0], mass_edges[-1]
        ],
        cmap="viridis"
    )

    plt.colorbar(image, label=colorbar_label)

    plt.xticks(halo_centers)
    plt.yticks(mass_centers)

    plt.xlabel(
        r"$\log_{10}[M_{\rm halo}/(M_\odot/h)]$"
    )
    plt.ylabel(
        r"$\log_{10}(M_*/M_\odot)$"
    )

    plt.tight_layout()
    plt.savefig(filename, bbox_inches="tight")
    plt.show()

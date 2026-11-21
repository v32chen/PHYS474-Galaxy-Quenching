#!/usr/bin/env python
# coding: utf-8

# In[2]:
import numpy as np
import matplotlib.pyplot as plt
import csv

# Load the measurements and final sample
with np.load("galaxy_combined_preclean.npz") as d:
    mass = d["log_stellar_mass"].copy()
    ssfr = d["log_ssfr"].copy()
    z = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()
    halo_mass = d["halo_mass"].copy()
    group_id = d["group_id"].copy()

primary = np.load("galaxy_final_sample_mask.npy")

assert len(primary) == len(mass)
print("Final sample size:", np.count_nonzero(primary))

quenched = ssfr < -11.0
mass_edges = np.round(np.arange(10.0, 10.51, 0.1), 2)
mass_centers = (mass_edges[:-1] + mass_edges[1:]) / 2

rng = np.random.default_rng(474)
n_boot = 500

f_cen, f_sat, efficiency = [], [], []
cen_ci, sat_ci, eff_ci = [], [], []
table_rows = []

for lo, hi in zip(mass_edges[:-1], mass_edges[1:]):

    selected = primary & (mass >= lo) & (mass < hi)

    # Work with the galaxies in this bin
    g = group_id[selected]
    q = quenched[selected]
    c = p_sat[selected] < 0.5
    s = p_sat[selected] > 0.5

    groups, inverse = np.unique(g, return_inverse=True)
    ng = len(groups)

    # Group-level counts for bootstrap resampling
    nc = np.bincount(inverse, weights=c, minlength=ng)
    ns = np.bincount(inverse, weights=s, minlength=ng)
    kc = np.bincount(inverse, weights=c & q, minlength=ng)
    ks = np.bincount(inverse, weights=s & q, minlength=ng)

    n_c = int(nc.sum())
    n_s = int(ns.sum())
    k_c = int(kc.sum())
    k_s = int(ks.sum())

    fc = k_c / n_c
    fs = k_s / n_s
    eps = (fs - fc) / (1 - fc)

    boot_fc, boot_fs, boot_eps = [], [], []

    for _ in range(n_boot):
        draw = rng.integers(0, ng, size=ng)

        nc_b = nc[draw].sum()
        ns_b = ns[draw].sum()

        if nc_b == 0 or ns_b == 0:
            continue

        fc_b = kc[draw].sum() / nc_b
        fs_b = ks[draw].sum() / ns_b

        if fc_b >= 1:
            continue

        boot_fc.append(fc_b)
        boot_fs.append(fs_b)
        boot_eps.append((fs_b-fc_b)/(1-fc_b))

    ci_c = np.percentile(boot_fc, [2.5, 97.5])
    ci_s = np.percentile(boot_fs, [2.5, 97.5])
    ci_e = np.percentile(boot_eps, [2.5, 97.5])

    f_cen.append(fc)
    f_sat.append(fs)
    efficiency.append(eps)
    cen_ci.append(ci_c)
    sat_ci.append(ci_s)
    eff_ci.append(ci_e)

    table_rows.append([
        lo, hi, n_c, n_s, k_c, k_s,
        fc, fs, fs-fc, eps, ci_e[0], ci_e[1]
    ])

    print(
        f"{lo:.1f}-{hi:.1f}: "
        f"Ncen={n_c}, Nsat={n_s}, "
        f"fQcen={fc:.3f}, fQsat={fs:.3f}, "
        f"eff={eps:.3f} "
        f"[{ci_e[0]:.3f}, {ci_e[1]:.3f}]"
    )

f_cen = np.array(f_cen)
f_sat = np.array(f_sat)
efficiency = np.array(efficiency)
cen_ci = np.array(cen_ci)
sat_ci = np.array(sat_ci)
eff_ci = np.array(eff_ci)

# Save the final numerical results
with open("final_quenching_results.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "mass_low", "mass_high", "n_central",
        "n_satellite", "quenched_central",
        "quenched_satellite", "fQ_central",
        "fQ_satellite", "excess_fraction",
        "satellite_efficiency",
        "efficiency_CI_low", "efficiency_CI_high"
    ])
    writer.writerows(table_rows)

print("\nFinal results table saved.")


# In[3]:
plt.figure(figsize=(8, 5.5))

for values, ci, label, marker in [
    (f_cen, cen_ci, "Centrals", "o"),
    (f_sat, sat_ci, "Satellites", "s")
]:
    errors = [
        values - ci[:, 0],
        ci[:, 1] - values
    ]

    plt.errorbar(
        mass_centers, values,
        yerr=errors,
        fmt=marker + "-",
        capsize=4,
        label=label
    )

plt.xlabel(r"$\log_{10}(M_*/M_\odot)$")
plt.ylabel("Quenched fraction")
plt.ylim(0, 1)
plt.legend()
plt.grid(alpha=0.2)
plt.tight_layout()

plt.savefig("FINAL_FIGURE_1.pdf", bbox_inches="tight")
plt.show()


# In[4]:
halo_edges = np.arange(11.5, 15.01, 0.5)
halo_centers = (halo_edges[:-1] + halo_edges[1:]) / 2

log_halo = np.full(len(halo_mass), np.nan)
valid = np.isfinite(halo_mass) & (halo_mass > 0)
log_halo[valid] = np.log10(halo_mass[valid])

eff_grid = np.full((5, 7), np.nan)
halo_rows = []

# Recalculate environmental results for the final sample
for i in range(5):
    mlo, mhi = mass_edges[i:i+2]

    in_mass = primary & (mass >= mlo) & (mass < mhi)

    for j in range(7):
        hlo, hhi = halo_edges[j:j+2]

        sat = (
            in_mass &
            (p_sat > 0.5) &
            (log_halo >= hlo) &
            (log_halo < hhi)
        )

        n = np.count_nonzero(sat)

        if n < 30:
            continue

        k = np.count_nonzero(sat & quenched)
        fs = k / n

        eps = (fs - f_cen[i]) / (1 - f_cen[i])

        eff_grid[i, j] = eps
        halo_rows.append([mlo, mhi, hlo, hhi, n, k, fs, eps])

with open("final_environmental_results.csv",
          "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "mass_low", "mass_high", "halo_low",
        "halo_high", "n_satellite", "n_quenched",
        "satellite_fraction", "efficiency"
    ])
    writer.writerows(halo_rows)

# FIGURE 2: Fixed stellar mass, varying halo mass
i = 2    # Stellar mass 10.2-10.3
sat_values = np.full(7, np.nan)
sat_low = np.full(7, np.nan)
sat_high = np.full(7, np.nan)

for j in range(7):
    selected = (
        primary &
        (mass >= mass_edges[i]) &
        (mass < mass_edges[i+1]) &
        (p_sat > 0.5) &
        (log_halo >= halo_edges[j]) &
        (log_halo < halo_edges[j+1])
    )

    n = np.count_nonzero(selected)

    if n < 30:
        continue

    sat_values[j] = np.mean(quenched[selected])

    # Resample entire groups for satellite uncertainty
    g = group_id[selected]
    q = quenched[selected]

    groups, inverse = np.unique(g, return_inverse=True)
    ng = len(groups)
    counts = np.bincount(inverse, minlength=ng)
    qcounts = np.bincount(
        inverse, weights=q, minlength=ng
    )

    estimates = []

    for _ in range(n_boot):
        draw = rng.integers(0, ng, size=ng)
        estimates.append(
            qcounts[draw].sum() / counts[draw].sum()
        )

    sat_low[j], sat_high[j] = np.percentile(
        estimates, [2.5, 97.5]
    )

plt.figure(figsize=(8, 5.5))
plt.errorbar(
    halo_centers, sat_values,
    yerr=[
        sat_values - sat_low,
        sat_high - sat_values
    ],
    fmt="o-", capsize=4,
    label="Satellites"
)

plt.axhline(
    f_cen[i], linestyle="--",
    label="Central baseline"
)

plt.axhspan(
    cen_ci[i, 0], cen_ci[i, 1],
    alpha=0.12
)

plt.xlabel(r"$\log_{10}[M_{\rm halo}/(M_\odot/h)]$")
plt.ylabel("Quenched fraction")
plt.ylim(0, 1)
plt.legend()
plt.grid(alpha=0.2)
plt.tight_layout()

plt.savefig("FINAL_FIGURE_2.pdf", bbox_inches="tight")
plt.show()

# FIGURE 3: Efficiency across all mass bins
plt.figure(figsize=(9, 5.5))

image = plt.pcolormesh(
    halo_edges, mass_edges,
    np.ma.masked_invalid(eff_grid),
    cmap="viridis",
    shading="auto"
)

plt.colorbar(image, label="Satellite quenching efficiency")
plt.xticks(halo_centers)
plt.yticks(mass_centers)
plt.xlabel(r"$\log_{10}[M_{\rm halo}/(M_\odot/h)]$")
plt.ylabel(r"$\log_{10}(M_*/M_\odot)$")

plt.tight_layout()
plt.savefig("FINAL_FIGURE_3.pdf", bbox_inches="tight")
plt.show()

print("All three final figures saved.")


# In[5]:
print("FINAL SAMPLE: THRESHOLD SENSITIVITY")
print("Mass bin | -11.2 | -11.0 | -10.8")

for lo, hi in zip(mass_edges[:-1], mass_edges[1:]):
    selected = primary & (mass >= lo) & (mass < hi)
    cen = selected & (p_sat < 0.5)
    sat = selected & (p_sat > 0.5)

    values = []

    for threshold in [-11.2, -11.0, -10.8]:
        fc = np.mean(ssfr[cen] < threshold)
        fs = np.mean(ssfr[sat] < threshold)
        values.append((fs - fc) / (1 - fc))

    print(
        f"{lo:.1f}-{hi:.1f} | "
        + " | ".join(f"{v:.3f}" for v in values)
    )

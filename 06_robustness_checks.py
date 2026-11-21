#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt

# Load galaxy measurements
with np.load("galaxy_combined_preclean.npz") as d:
    log_mass = d["log_stellar_mass"].copy()
    log_ssfr = d["log_ssfr"].copy()
    z = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()

clean = np.load("galaxy_quality_mask.npy")

assert len(clean) == len(z)

# Keep the same quenching threshold
quenched = log_ssfr < -11.0

# Stellar-mass bins
mass_edges = np.round(np.arange(10.0, 10.51, 0.1), 2)
mass_centers = (mass_edges[:-1] + mass_edges[1:]) / 2

# Two redshift intervals
redshift_bins = [
    (0.02, 0.05),
    (0.05, 0.08)
]

plt.figure(figsize=(9, 6))

for z_low, z_high in redshift_bins:

    efficiencies = []

    print(f"\nREDSHIFT {z_low:.2f} - {z_high:.2f}")
    print("Mass | N_cen | N_sat | fQ_cen | fQ_sat | Efficiency")

    for m_low, m_high in zip(
        mass_edges[:-1], mass_edges[1:]
    ):

        selected = (
            clean &
            (z >= z_low) &
            (z < z_high) &
            (log_mass >= m_low) &
            (log_mass < m_high)
        )

        central = selected & (p_sat < 0.5)
        satellite = selected & (p_sat > 0.5)

        n_cen = np.count_nonzero(central)
        n_sat = np.count_nonzero(satellite)

        if n_cen < 30 or n_sat < 30:
            efficiencies.append(np.nan)
            continue

        f_cen = np.mean(quenched[central])
        f_sat = np.mean(quenched[satellite])

        if f_cen >= 1:
            efficiencies.append(np.nan)
            continue

        epsilon = (f_sat - f_cen) / (1 - f_cen)
        efficiencies.append(epsilon)

        print(
            f"{m_low:.1f}-{m_high:.1f} | "
            f"{n_cen:5d} | {n_sat:5d} | "
            f"{f_cen:.3f} | {f_sat:.3f} | "
            f"{epsilon:.3f}"
        )

    plt.plot(
        mass_centers,
        efficiencies,
        marker="o",
        linewidth=2,
        label=f"{z_low:.2f} <= z < {z_high:.2f}"
    )

plt.xlabel(r"$\log_{10}(M_*/M_\odot)$")
plt.ylabel("Satellite quenching efficiency")

plt.legend()
plt.grid(alpha=0.2)
plt.tight_layout()

plt.savefig(
    "quenching_efficiency_redshift_comparison.pdf",
    bbox_inches="tight"
)

plt.show()


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import csv

# Load data
with np.load("galaxy_combined_preclean.npz") as d:
    log_mass = d["log_stellar_mass"].copy()
    log_ssfr = d["log_ssfr"].copy()
    z = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()
    halo_mass = d["halo_mass"].copy()

clean = np.load("galaxy_quality_mask.npy")

quenched = log_ssfr < -11.0

# Stellar-mass and halo-mass bins
mass_edges = np.round(np.arange(10.0, 10.51, 0.1), 2)
halo_edges = np.arange(11.5, 15.01, 0.5)

# Six narrow redshift bins
z_edges = np.linspace(0.02, 0.08, 7)

log_halo = np.full(len(halo_mass), np.nan)
valid_halo = np.isfinite(halo_mass) & (halo_mass > 0)
log_halo[valid_halo] = np.log10(halo_mass[valid_halo])

shape = (len(mass_edges)-1, len(halo_edges)-1)

efficiency_adjusted = np.full(shape, np.nan)
efficiency_original = np.full(shape, np.nan)

rows = []

for i, (m_lo, m_hi) in enumerate(
    zip(mass_edges[:-1], mass_edges[1:])
):
    mass_sample = (
        clean &
        (z >= 0.02) &
        (z < 0.08) &
        (log_mass >= m_lo) &
        (log_mass < m_hi)
    )

    central = mass_sample & (p_sat < 0.5)
    satellite = mass_sample & (p_sat > 0.5)

    if np.count_nonzero(central) < 30:
        continue

    f_cen_original = np.mean(quenched[central])

    print(f"\nSTELLAR MASS {m_lo:.1f}-{m_hi:.1f}")
    print("Halo bin | N_sat | Adjusted efficiency")

    for j, (h_lo, h_hi) in enumerate(
        zip(halo_edges[:-1], halo_edges[1:])
    ):
        sat = (
            satellite &
            (log_halo >= h_lo) &
            (log_halo < h_hi)
        )

        n_sat = np.count_nonzero(sat)

        if n_sat < 30:
            continue

        f_sat = np.mean(quenched[sat])

        # Original, unadjusted efficiency
        eps_original = (
            (f_sat - f_cen_original) /
            (1 - f_cen_original)
        )

        # Reweight central quenched fraction
        # to match the satellites' redshift distribution
        f_cen_adjusted = 0.0
        sufficient_data = True

        for z_lo, z_hi in zip(
            z_edges[:-1], z_edges[1:]
        ):
            in_z = (z >= z_lo) & (z < z_hi)

            sat_z = sat & in_z
            cen_z = central & in_z

            n_sat_z = np.count_nonzero(sat_z)
            n_cen_z = np.count_nonzero(cen_z)

            if n_sat_z == 0:
                continue

            if n_cen_z < 30:
                sufficient_data = False
                break

            weight = n_sat_z / n_sat

            f_cen_adjusted += (
                weight * np.mean(quenched[cen_z])
            )

        if not sufficient_data or f_cen_adjusted >= 1:
            continue

        eps_adjusted = (
            (f_sat - f_cen_adjusted) /
            (1 - f_cen_adjusted)
        )

        efficiency_adjusted[i, j] = eps_adjusted
        efficiency_original[i, j] = eps_original

        rows.append([
            m_lo, m_hi, h_lo, h_hi,
            n_sat, f_sat,
            f_cen_adjusted, eps_adjusted,
            eps_original
        ])

        print(
            f"{h_lo:.1f}-{h_hi:.1f} | "
            f"{n_sat:4d} | "
            f"{eps_adjusted:.3f}"
        )

# Save numerical results
with open(
    "quenching_redshift_adjusted.csv",
    "w",
    newline=""
) as file:
    writer = csv.writer(file)

    writer.writerow([
        "mass_low", "mass_high",
        "halo_low", "halo_high",
        "n_satellite",
        "fQ_satellite",
        "fQ_central_redshift_adjusted",
        "efficiency_adjusted",
        "efficiency_original"
    ])

    writer.writerows(rows)

# Compare the two methods
difference = efficiency_adjusted - efficiency_original

print("\nMAXIMUM ABSOLUTE CHANGE IN EFFICIENCY:")
if np.isfinite(difference).any():
    print(np.nanmax(np.abs(difference)))

print("\nResults saved.")


# In[3]:


plt.figure(figsize=(10, 6))

plt.imshow(
    np.ma.masked_invalid(efficiency_adjusted),
    origin="lower",
    aspect="auto",
    extent=[
        halo_edges[0],
        halo_edges[-1],
        mass_edges[0],
        mass_edges[-1]
    ],
    cmap="viridis"
)

plt.colorbar(label="Redshift-adjusted quenching efficiency")

plt.xticks((halo_edges[:-1] + halo_edges[1:]) / 2)
plt.yticks((mass_edges[:-1] + mass_edges[1:]) / 2)

plt.xlabel(
    r"$\log_{10}[M_{\rm halo}/(M_\odot/h)]$"
)
plt.ylabel(
    r"$\log_{10}(M_*/M_\odot)$"
)

plt.tight_layout()

plt.savefig(
    "quenching_efficiency_redshift_adjusted.pdf",
    bbox_inches="tight"
)

plt.show()


# In[4]:


import numpy as np
import matplotlib.pyplot as plt

# Load the combined catalogue
with np.load("galaxy_combined_preclean.npz") as d:
    mass = d["log_stellar_mass"].copy()
    ssfr = d["log_ssfr"].copy()
    z = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()
    tinker_idx = d["tinker_row"].copy()

clean = np.load("galaxy_quality_mask.npy")

# Read Vmax from the Tinker catalogue (column index 5)
tinker = np.loadtxt("sdss_kdgroups_v1.0.dat")
vmax = tinker[tinker_idx, 5]

valid_vmax = np.isfinite(vmax) & (vmax > 0)

sample = clean & (z >= 0.02) & (z < 0.08)

print("GALAXIES WITH VALID VMAX")
print("-----------------------")
print("Selected:", np.count_nonzero(sample))
print("Valid Vmax:", np.count_nonzero(sample & valid_vmax))
print("Invalid Vmax:", np.count_nonzero(sample & ~valid_vmax))

# Construct inverse-volume weights
weights = np.zeros(len(vmax), dtype=float)
weights[valid_vmax] = 1.0 / vmax[valid_vmax]

# Keep the same galaxies for weighted and unweighted results
sample = sample & valid_vmax
quenched = ssfr < -11.0

def fractions(mask):
    n = np.count_nonzero(mask)

    if n == 0:
        return np.nan, np.nan, np.nan

    f_raw = np.mean(quenched[mask])
    w = weights[mask]
    f_weighted = np.average(quenched[mask], weights=w)

    # Effective number of observations under the weights
    n_eff = w.sum()**2 / np.sum(w**2)

    return f_raw, f_weighted, n_eff

mass_edges = np.round(np.arange(10.0, 10.51, 0.1), 2)
mass_centers = (mass_edges[:-1] + mass_edges[1:]) / 2

eps_raw = []
eps_weighted = []

print("\nINVERSE-VMAX SENSITIVITY TEST")
print("Mass bin | Raw efficiency | Weighted efficiency | N_eff cen/sat")

for lo, hi in zip(mass_edges[:-1], mass_edges[1:]):

    in_mass = sample & (mass >= lo) & (mass < hi)

    cen = in_mass & (p_sat < 0.5)
    sat = in_mass & (p_sat > 0.5)

    fc, fc_w, nc_eff = fractions(cen)
    fs, fs_w, ns_eff = fractions(sat)

    e_raw = (fs - fc) / (1 - fc)
    e_weighted = (fs_w - fc_w) / (1 - fc_w)

    eps_raw.append(e_raw)
    eps_weighted.append(e_weighted)

    print(
        f"{lo:.1f}-{hi:.1f} | "
        f"{e_raw:.3f} | {e_weighted:.3f} | "
        f"{nc_eff:.0f}/{ns_eff:.0f}"
    )


# In[5]:


plt.figure(figsize=(9, 6))

plt.plot(
    mass_centers, eps_raw,
    marker="o", linewidth=2,
    label="Unweighted"
)

plt.plot(
    mass_centers, eps_weighted,
    marker="s", linewidth=2,
    label="Inverse-Vmax weighted"
)

plt.xlabel(r"$\log_{10}(M_*/M_\odot)$")
plt.ylabel("Satellite quenching efficiency")

plt.legend()
plt.grid(alpha=0.2)
plt.tight_layout()

plt.savefig(
    "vmax_weighting_diagnostic.pdf",
    bbox_inches="tight"
)

plt.show()


# In[6]:


import numpy as np

with np.load("galaxy_combined_preclean.npz") as d:
    mass = d["log_stellar_mass"].copy()
    ssfr = d["log_ssfr"].copy()
    z = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()

clean = np.load("galaxy_quality_mask.npy")

# Keep the same sample as our previous analyses
base = clean & (z >= 0.02) & (z < 0.08)

mass_edges = np.round(np.arange(10.0, 10.51, 0.1), 2)
thresholds = [-11.2, -11.0, -10.8]

print("QUENCHING THRESHOLD SENSITIVITY")
print("Mass bin | Threshold | Efficiency")

for m_lo, m_hi in zip(mass_edges[:-1], mass_edges[1:]):

    selected = (
        base &
        (mass >= m_lo) &
        (mass < m_hi)
    )

    central = selected & (p_sat < 0.5)
    satellite = selected & (p_sat > 0.5)

    for threshold in thresholds:

        quenched = ssfr < threshold

        f_cen = np.mean(quenched[central])
        f_sat = np.mean(quenched[satellite])

        epsilon = (f_sat - f_cen) / (1 - f_cen)

        print(
            f"{m_lo:.1f}-{m_hi:.1f} | "
            f"{threshold:.1f} | "
            f"{epsilon:.3f}"
        )


# In[7]:


import numpy as np

with np.load("galaxy_combined_preclean.npz") as d:
    mass = d["log_stellar_mass"].copy()
    ssfr = d["log_ssfr"].copy()
    z = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()
    group_id = d["group_id"].copy()

clean = np.load("galaxy_quality_mask.npy")

base = clean & (z >= 0.02) & (z < 0.08)
quenched = ssfr < -11.0

mass_edges = np.round(np.arange(10.0, 10.51, 0.1), 2)

rng = np.random.default_rng(474)
n_bootstrap = 400

print("GROUP-BOOTSTRAP UNCERTAINTIES")
print("Mass bin | Efficiency | 95% bootstrap CI")

for lo, hi in zip(mass_edges[:-1], mass_edges[1:]):

    selected = base & (mass >= lo) & (mass < hi)

    # Galaxies in this stellar-mass bin
    g = group_id[selected]
    q = quenched[selected]
    cen = p_sat[selected] < 0.5
    sat = p_sat[selected] > 0.5

    # Identify distinct host groups
    unique_groups, inverse = np.unique(
        g, return_inverse=True
    )

    ng = len(unique_groups)

    # Number of central/satellite galaxies per group
    nc = np.bincount(inverse, weights=cen, minlength=ng)
    ns = np.bincount(inverse, weights=sat, minlength=ng)

    # Number quenched per group
    kc = np.bincount(
        inverse, weights=(cen & q), minlength=ng
    )
    ks = np.bincount(
        inverse, weights=(sat & q), minlength=ng
    )

    def calc_efficiency(indices):
        ncen = nc[indices].sum()
        nsat = ns[indices].sum()

        if ncen == 0 or nsat == 0:
            return np.nan

        fcen = kc[indices].sum() / ncen
        fsat = ks[indices].sum() / nsat

        if fcen >= 1:
            return np.nan

        return (fsat - fcen) / (1 - fcen)

    # Original efficiency
    original = calc_efficiency(np.arange(ng))

    # Resample whole groups with replacement
    bootstrap = []

    for b in range(n_bootstrap):
        draw = rng.integers(0, ng, size=ng)
        bootstrap.append(calc_efficiency(draw))

    bootstrap = np.array(bootstrap)
    bootstrap = bootstrap[np.isfinite(bootstrap)]

    lower, upper = np.percentile(
        bootstrap, [2.5, 97.5]
    )

    print(
        f"{lo:.1f}-{hi:.1f} | "
        f"{original:.3f} | "
        f"{lower:.3f} to {upper:.3f}"
    )


# In[8]:


import numpy as np

with np.load("galaxy_combined_preclean.npz") as d:
    mass = d["log_stellar_mass"].copy()
    z = d["redshift"].copy()
    p_sat = d["satellite_probability"].copy()

clean = np.load("galaxy_quality_mask.npy")

# Published SDSS DR7 empirical mixture limit
# Zu & Mandelbaum (2015), using h = 0.7
h = 0.7
z_for_limit = np.maximum(z, 0.025)

log_mass_limit = (
    8.0
    + 5.4 * (z_for_limit - 0.025)**0.33
    - 2 * np.log10(h)
)

# Primary sample: nearby galaxies in five mass bins
primary = (
    clean &
    (z >= 0.02) &
    (z < 0.05) &
    (mass >= 10.0) &
    (mass < 10.5) &
    (mass >= log_mass_limit)
)

print("FINAL PRIMARY SAMPLE")
print("--------------------")
print("Total:", np.count_nonzero(primary))
print("Centrals:", np.count_nonzero(primary & (p_sat < 0.5)))
print("Satellites:", np.count_nonzero(primary & (p_sat > 0.5)))

print("\nPublished approximate mass boundary at z=0.05:")
limit_005 = (
    8 + 5.4 * (0.05 - 0.025)**0.33
    - 2 * np.log10(h)
)
print(f"log stellar mass limit = {limit_005:.3f}")

# Save this mask without altering previous results
np.save("galaxy_final_sample_mask.npy", primary)

print("\nFinal sample mask saved.")


# In[ ]:





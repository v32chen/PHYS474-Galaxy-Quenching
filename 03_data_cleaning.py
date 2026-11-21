#!/usr/bin/env python
# coding: utf-8

# In[1]:
import numpy as np

# Load the combined galaxy dataset
with np.load("galaxy_combined_preclean.npz") as d:
    log_mass = d["log_stellar_mass"]
    log_sfr = d["log_sfr"]
    log_ssfr = d["log_ssfr"]
    z = d["redshift"]
    p_sat = d["satellite_probability"]
    halo_mass = d["halo_mass"]

    sfr_flag = d["sfr_flag"]
    ssfr_flag = d["ssfr_flag"]
    z_warning = d["redshift_warning"]

print("DATA QUALITY INSPECTION")
print("------------------------")
print("Total records:", len(log_mass))

# Inspect measurement ranges
for name, values in [
    ("log stellar mass", log_mass),
    ("log SFR", log_sfr),
    ("log sSFR", log_ssfr),
    ("Redshift", z)
]:
    finite = values[np.isfinite(values)]

    print(f"\n{name}")
    print("Minimum:", np.min(finite))
    print("Median:", np.median(finite))
    print("Maximum:", np.max(finite))

# Count suspicious measurements
print("\nSUSPICIOUS MEASUREMENTS")
print("------------------------")

print("Stellar mass <= 0:",
      np.count_nonzero(log_mass <= 0))

print("Stellar mass exactly -1:",
      np.count_nonzero(log_mass == -1))

print("log SFR below -20:",
      np.count_nonzero(log_sfr < -20))

print("log sSFR below -20:",
      np.count_nonzero(log_ssfr < -20))

print("Redshift <= 0:",
      np.count_nonzero(z <= 0))

# Count galaxies passing initial quality checks
quality = (
    np.isfinite(log_mass) &
    np.isfinite(log_sfr) &
    np.isfinite(log_ssfr) &
    np.isfinite(z) &
    np.isfinite(halo_mass) &
    np.isfinite(p_sat) &
    (log_mass > 0) &
    (z > 0) &
    (halo_mass > 0) &
    (p_sat >= 0) &
    (p_sat <= 1) &
    (z_warning == 0) &
    (sfr_flag == 0) &
    (ssfr_flag == 0)
)

print("\nPRELIMINARY QUALITY FILTER")
print("---------------------------")
print("Records passing:", np.count_nonzero(quality))

print("Centrals passing:",
      np.count_nonzero(quality & (p_sat < 0.5)))

print("Satellites passing:",
      np.count_nonzero(quality & (p_sat > 0.5)))


# In[2]:
import numpy as np

# Check the extreme values more carefully
print("SFR exactly -99:",
      np.count_nonzero(log_sfr == -99))

print("sSFR exactly -99:",
      np.count_nonzero(log_ssfr == -99))

print("SFR < -20 but not -99:",
      np.count_nonzero((log_sfr < -20) &
                       (log_sfr != -99)))

print("sSFR < -20 but not -99:",
      np.count_nonzero((log_ssfr < -20) &
                       (log_ssfr != -99)))

# Refine the preliminary quality mask
# The -20 limits exclude clearly extreme estimates,
# not galaxies classified as quenched.
clean = (
    quality &
    (log_sfr > -20) &
    (log_ssfr > -20)
)

print("\nQUALITY-CLEANED SAMPLE")
print("----------------------")
print("Total:", np.count_nonzero(clean))
print("Centrals:", np.count_nonzero(clean & (p_sat < 0.5)))
print("Satellites:", np.count_nonzero(clean & (p_sat > 0.5)))

# Keep the original data and save which rows passed
np.save("galaxy_quality_mask.npy", clean)

print("\nQuality mask saved.")


# In[3]:
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plt.figure(figsize=(9, 6))

plt.hist2d(
    z[clean],
    log_mass[clean],
    bins=80,
    norm=LogNorm(),
    cmap="viridis"
)

plt.colorbar(label="Number of galaxies")

plt.xlabel("Redshift (z)", fontsize=12)
plt.ylabel(r"$\log_{10}(M_*/M_\odot)$", fontsize=12)

plt.tight_layout()
plt.savefig("stellar_mass_vs_redshift.pdf",
            bbox_inches="tight")
plt.show()


# In[4]:
import numpy as np
import matplotlib.pyplot as plt

# Preliminary redshift and stellar-mass selection
sample = (
    clean &
    (z >= 0.02) &
    (z < 0.08) &
    (log_mass >= 10.0) &
    (log_mass < 10.5)
)


# In[5]:
# Separate central and satellite galaxies
central = sample & (p_sat < 0.5)
satellite = sample & (p_sat > 0.5)

print("EXPLORATORY GALAXY SAMPLE")
print("-------------------------")
print("Total galaxies:", np.count_nonzero(sample))
print("Central galaxies:", np.count_nonzero(central))
print("Satellite galaxies:", np.count_nonzero(satellite))


# In[6]:
import numpy as np
import matplotlib.pyplot as plt

# Specific star-formation rate bins
bins = np.arange(-14, -8, 0.15)

plt.figure(figsize=(9, 6))

# Central galaxies
plt.hist(
    log_ssfr[central],
    bins=bins,
    density=True,
    histtype="step",
    linewidth=2,
    label="Central galaxies"
)

# Satellite galaxies
plt.hist(
    log_ssfr[satellite],
    bins=bins,
    density=True,
    histtype="step",
    linewidth=2,
    label="Satellite galaxies"
)

# Proposed quenching threshold
plt.axvline(
    -11,
    color="red",
    linestyle="--",
    linewidth=2,
    label="Quenching threshold"
)

plt.xlabel(
    r"$\log_{10}(\mathrm{sSFR}/\mathrm{yr}^{-1})$"
)
plt.ylabel("Probability density")

plt.legend()
plt.tight_layout()

plt.savefig(
    "ssfr_central_vs_satellite_exploratory.pdf",
    bbox_inches="tight"
)

plt.show()


# In[7]:
import numpy as np

# Proposed quenching threshold
threshold = -11.0
quenched = log_ssfr < threshold

# Count galaxies
n_cen = np.count_nonzero(central)
n_sat = np.count_nonzero(satellite)

q_cen = np.count_nonzero(central & quenched)
q_sat = np.count_nonzero(satellite & quenched)

# Calculate quenched fractions
f_cen = q_cen / n_cen
f_sat = q_sat / n_sat

# Wilson 95% confidence intervals
def wilson_interval(k, n):
    p = k / n
    z95 = 1.96
    denom = 1 + z95**2 / n
    center = (p + z95**2 / (2*n)) / denom
    margin = z95 * np.sqrt(
        p*(1-p)/n + z95**2/(4*n**2)
    ) / denom
    return center - margin, center + margin

ci_cen = wilson_interval(q_cen, n_cen)
ci_sat = wilson_interval(q_sat, n_sat)

# Excess satellite quenching
delta_f = f_sat - f_cen

print("QUENCHING ANALYSIS")
print("------------------")

print(f"Central quenched: {q_cen} / {n_cen}")
print(f"Satellite quenched: {q_sat} / {n_sat}")

print(f"\nCentral quenched fraction: {f_cen:.4f}")
print(f"95% CI: {ci_cen[0]:.4f} to {ci_cen[1]:.4f}")

print(f"\nSatellite quenched fraction: {f_sat:.4f}")
print(f"95% CI: {ci_sat[0]:.4f} to {ci_sat[1]:.4f}")

print(f"\nExcess satellite quenching: {delta_f:.4f}")
print(f"Difference: {100*delta_f:.2f} percentage points")

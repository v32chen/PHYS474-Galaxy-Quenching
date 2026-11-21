#!/usr/bin/env python
# coding: utf-8

# In[10]:
import numpy as np
from astropy.io import fits

# Load the matching indices from our previous analysis
with np.load("galaxy_crossmatch.npz") as matches:
    mpa_idx = matches["mpa_row"].copy()
    tinker_idx = matches["tinker_row"].copy()
    separation = matches["separation_arcsec"].copy()
    delta_z = matches["redshift_difference"].copy()

print("Matched galaxies:", len(mpa_idx))

# Helper: read only the required rows from an MPA-JHU file
def read_mpa(filename, column):
    with fits.open(filename, memmap=True) as hdul:
        return np.array(hdul[1].data[column][mpa_idx])

# MPA-JHU galaxy properties
ra = read_mpa("gal_info_dr7_v5_2.fit", "RA")
dec = read_mpa("gal_info_dr7_v5_2.fit", "DEC")
z = read_mpa("gal_info_dr7_v5_2.fit", "Z")
z_warning = read_mpa("gal_info_dr7_v5_2.fit", "Z_WARNING")

log_mass = read_mpa("totlgm_dr7_v5_2.fit", "MEDIAN")

log_sfr = read_mpa(
    "gal_totsfr_dr7_v5_2.fits", "MEDIAN"
)

log_ssfr = read_mpa(
    "gal_totspecsfr_dr7_v5_2.fits", "MEDIAN"
)

sfr_flag = read_mpa(
    "gal_totsfr_dr7_v5_2.fits", "FLAG"
)

ssfr_flag = read_mpa(
    "gal_totspecsfr_dr7_v5_2.fits", "FLAG"
)

# Tinker environmental properties
tinker = np.loadtxt("sdss_kdgroups_v1.0.dat")

p_sat = tinker[tinker_idx, 6]
halo_mass = tinker[tinker_idx, 7]
group_id = tinker[tinker_idx, 10].astype(int)

# Classify centrals and satellites
is_satellite = p_sat > 0.5
is_central = p_sat < 0.5

print("All measurements successfully extracted.")
print("Number of combined records:", len(log_mass))


# In[11]:
print("COMBINED GALAXY SAMPLE")
print("----------------------")
print("Total matched records:", len(z))

print("Central galaxies:", np.count_nonzero(is_central))
print("Satellite galaxies:", np.count_nonzero(is_satellite))

print("\nMEASUREMENT CHECKS")
print("------------------")
print("Finite stellar masses:",
      np.count_nonzero(np.isfinite(log_mass)))

print("Finite SFR values:",
      np.count_nonzero(np.isfinite(log_sfr)))

print("Finite sSFR values:",
      np.count_nonzero(np.isfinite(log_ssfr)))

print("SFR FLAG = 0:",
      np.count_nonzero(sfr_flag == 0))

print("sSFR FLAG = 0:",
      np.count_nonzero(ssfr_flag == 0))

print("Redshift warning = 0:",
      np.count_nonzero(z_warning == 0))

print("\nFIRST FIVE MATCHED GALAXIES")
print("--------------------------")

for i in range(5):
    print(
        f"Galaxy {i+1}: "
        f"z={z[i]:.4f}, "
        f"logM={log_mass[i]:.2f}, "
        f"logSSFR={log_ssfr[i]:.2f}, "
        f"P_sat={p_sat[i]:.3f}"
    )


# In[12]:
np.savez_compressed(
    "galaxy_combined_preclean.npz",
    mpa_row=mpa_idx,
    tinker_row=tinker_idx,
    ra=ra,
    dec=dec,
    redshift=z,
    redshift_warning=z_warning,
    log_stellar_mass=log_mass,
    log_sfr=log_sfr,
    log_ssfr=log_ssfr,
    sfr_flag=sfr_flag,
    ssfr_flag=ssfr_flag,
    satellite_probability=p_sat,
    halo_mass=halo_mass,
    group_id=group_id,
    separation_arcsec=separation,
    redshift_difference=delta_z
)

print("Combined dataset saved successfully!")

#!/usr/bin/env python
# coding: utf-8

# In[30]:


import sys
import site

# Find where your personal Python packages are installed
user_site = site.getusersitepackages()

# Tell Jupyter to look there
site.addsitedir(user_site)

# Try importing Astropy
import astropy

print("Astropy version:", astropy.__version__)
print("Astropy location:", astropy.__file__)


# In[31]:


from astropy.io import fits

# Open your SDSS galaxy catalogue
filename = "gal_info_dr7_v5_2.fit"

with fits.open(filename) as hdul:
    data = hdul[1].data

    print("Number of galaxies:", len(data))
    print("Columns:", data.columns.names)

    print("\nFirst galaxy:")
    print("RA:", data["RA"][0])
    print("Dec:", data["DEC"][0])
    print("Redshift:", data["Z"][0])


# In[32]:


from astropy.io import fits
from pathlib import Path

files = [
    "totlgm_dr7_v5_2",
    "gal_totsfr_dr7_v5_2",
    "gal_totspecsfr_dr7_v5_2"
]

for name in files:
    matches = list(Path(".").glob(name + ".*"))

    if not matches:
        print(f"\n{name}: FILE NOT FOUND")
        continue

    filename = matches[0]

    with fits.open(filename) as hdul:
        data = hdul[1].data

        print(f"\nFile: {filename}")
        print("Number of galaxies:", len(data))
        print("Columns:", data.columns.names)


# In[24]:


from astropy.io import fits

with fits.open("totlgm_dr7_v5_2.fit") as h:
    mass = h[1].data["MEDIAN"][:5].copy()

with fits.open("gal_totsfr_dr7_v5_2.fits") as h:
    sfr = h[1].data["MEDIAN"][:5].copy()

with fits.open("gal_totspecsfr_dr7_v5_2.fits") as h:
    ssfr = h[1].data["MEDIAN"][:5].copy()

print("Galaxy | log(M*) | log(SFR) | log(sSFR)")
print("---------------------------------------")

for i in range(5):
    print(f"{i:6d} | {mass[i]:7.2f} | {sfr[i]:8.2f} | {ssfr[i]:9.2f}")


# In[25]:


import numpy as np

filename = "sdss_kdgroups_v1.0.dat"

# Load the galaxy group catalogue
group = np.loadtxt(filename)

# Examine its structure
print("Number of galaxies:", group.shape[0])
print("Number of columns:", group.shape[1])

# Check that galaxy IDs follow the expected order
ids = group[:, 0].astype(int)

print("IDs are sequential:",
      np.array_equal(ids, np.arange(len(group))))

# Show the first five galaxies
print("\nFirst five galaxy records:")
print(group[:5])


# In[26]:


import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u

# Load MPA-JHU galaxy information
with fits.open("gal_info_dr7_v5_2.fit") as hdul:
    mpa_ra = np.array(hdul[1].data["RA"], dtype=float)
    mpa_dec = np.array(hdul[1].data["DEC"], dtype=float)
    mpa_z = np.array(hdul[1].data["Z"], dtype=float)

# Load the Tinker group catalogue
tinker = np.loadtxt("sdss_kdgroups_v1.0.dat")

tinker_ra = tinker[:, 1]
tinker_dec = tinker[:, 2]
tinker_z = tinker[:, 3]

print("MPA-JHU observations:", len(mpa_ra))
print("Tinker galaxies:", len(tinker_ra))


# In[27]:


# Keep only valid positions and redshifts
valid_mpa = (
    np.isfinite(mpa_ra) &
    np.isfinite(mpa_dec) &
    np.isfinite(mpa_z) &
    (mpa_z > 0) &
    (mpa_dec >= -90) &
    (mpa_dec <= 90)
)

valid_tinker = (
    np.isfinite(tinker_ra) &
    np.isfinite(tinker_dec) &
    np.isfinite(tinker_z) &
    (tinker_z > 0) &
    (tinker_dec >= -90) &
    (tinker_dec <= 90)
)

# Preserve original row numbers
mpa_rows = np.flatnonzero(valid_mpa)
tinker_rows = np.flatnonzero(valid_tinker)

# Convert coordinates to astronomical sky positions
mpa_coords = SkyCoord(
    ra=mpa_ra[valid_mpa] * u.deg,
    dec=mpa_dec[valid_mpa] * u.deg
)

tinker_coords = SkyCoord(
    ra=tinker_ra[valid_tinker] * u.deg,
    dec=tinker_dec[valid_tinker] * u.deg
)

# Find every position pair within 1 arcsecond
i_t, i_m, separations, _ = search_around_sky(
    tinker_coords,
    mpa_coords,
    1 * u.arcsec
)

# Convert back to original file row numbers
i_t = tinker_rows[i_t]
i_m = mpa_rows[i_m]

print("Position-matching candidate pairs:", len(i_t))


# In[28]:


# Calculate redshift differences
dz = np.abs(tinker_z[i_t] - mpa_z[i_m])

# Keep candidate pairs with consistent redshifts
good_z = dz < 0.001

i_t = i_t[good_z]
i_m = i_m[good_z]
sep_arcsec = separations.arcsec[good_z]
dz = dz[good_z]

print("Candidate pairs passing both cuts:", len(i_t))

# Count how many candidates point to each record
count_t = np.bincount(i_t, minlength=len(tinker))
count_m = np.bincount(i_m, minlength=len(mpa_ra))

# Require a unique match in BOTH catalogues
unique = (
    (count_t[i_t] == 1) &
    (count_m[i_m] == 1)
)

matched_tinker = i_t[unique]
matched_mpa = i_m[unique]

matched_sep = sep_arcsec[unique]
matched_dz = dz[unique]

print("\nUnambiguous one-to-one matches:",
      len(matched_tinker))

print("Fraction of Tinker catalogue matched:",
      f"{100*len(matched_tinker)/len(tinker):.2f}%")

if len(matched_tinker) > 0:
    print("Median sky separation:",
          f"{np.median(matched_sep):.3f} arcsec")

    print("Median redshift difference:",
          f"{np.median(matched_dz):.6f}")

# Save row correspondences for later analysis
np.savez_compressed(
    "galaxy_crossmatch.npz",
    tinker_row=matched_tinker,
    mpa_row=matched_mpa,
    separation_arcsec=matched_sep,
    redshift_difference=matched_dz
)

print("\nCross-match indices saved.")


# In[29]:


print("First five unambiguous matches:\n")

print("Tinker ID | MPA-JHU row | Separation | Delta z")
print("------------------------------------------------")

for k in range(min(5, len(matched_tinker))):
    print(
        f"{int(tinker[matched_tinker[k], 0]):9d} | "
        f"{matched_mpa[k]:11d} | "
        f"{matched_sep[k]:7.3f}\" | "
        f"{matched_dz[k]:.6f}"
    )


# In[ ]:





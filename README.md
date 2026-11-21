# PHYS 474 — Environmental Dependence of Galaxy Quenching

## Project Overview

This project investigates how galaxy environment influences the suppression of star formation, with a focus on comparing central and satellite galaxies.

The main research question is:

**At fixed stellar mass, are satellite galaxies more likely to be quenched than central galaxies, and how does this excess quenching depend on host dark-matter halo mass?**

The analysis uses observational galaxy catalogues from Sloan Digital Sky Survey Data Release 7 (SDSS DR7), combining galaxy physical properties from MPA-JHU with environmental classifications from the Tinker galaxy group catalogue.

This project was developed for PHYS 474 (Fall 2026).

## Data Sources

The analysis uses two primary data resources.

### 1. MPA-JHU SDSS DR7 Galaxy Catalogues

Website: https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/

Four files are used:

- `gal_info_dr7_v5_2.fit` — galaxy coordinates, redshifts, and observational information
- `totlgm_dr7_v5_2.fit` — stellar-mass estimates
- `gal_totsfr_dr7_v5_2.fits` — total star-formation-rate estimates
- `gal_totspecsfr_dr7_v5_2.fits` — total specific star-formation-rate estimates

### 2. Tinker Galaxy Group Catalogue

Website: https://www.galaxygroupfinder.net/catalogs

Catalogue: `sdss_kdgroups_v1.0.dat`

This catalogue provides central/satellite classifications, satellite probabilities, estimated host-halo masses, and group membership information.

The version used in this analysis contains 559,028 records. The catalogue webpage describes 559,038 records, so this discrepancy is noted as a data-version limitation.

The original large astronomical catalogues are not included in this repository and must be obtained separately.

## Analysis Method

### Catalogue Cross-Matching

Galaxies from the MPA-JHU and Tinker catalogues are matched using sky coordinates and redshifts.

The matching criteria are:

- Angular separation less than 1 arcsecond
- Absolute redshift difference less than 0.001
- A unique one-to-one correspondence between catalogue records

This procedure produces 498,431 unambiguous matched records.

### Data Cleaning and Sample Selection

Measurements are screened for invalid values and quality flags.

The primary sample is restricted to:

- Redshift: \(0.02 \leq z < 0.05\)
- Stellar mass: \(10.0 \leq \log_{10}(M_*/M_\odot) < 10.5\)

An empirical stellar-mass completeness guideline is also applied. This is an approximate selection rather than a guarantee of complete sampling.

The primary sample contains **16,602 galaxies**:

- 10,397 central galaxies
- 6,205 satellite galaxies

### Galaxy Quenching Classification

Galaxies are classified as quenched when:

\[
\log_{10}(\mathrm{sSFR}/\mathrm{yr}^{-1}) < -11
\]

The quenched fraction is calculated as:

\[
f_Q = \frac{N_{\mathrm{quenched}}}{N_{\mathrm{total}}}
\]

The excess satellite quenched fraction is:

\[
\Delta f_Q = f_{Q,\mathrm{sat}} - f_{Q,\mathrm{cen}}
\]

Satellite-quenching efficiency is defined as:

\[
\epsilon_{\mathrm{sat}} =
\frac{f_{Q,\mathrm{sat}} - f_{Q,\mathrm{cen}}}
{1 - f_{Q,\mathrm{cen}}}
\]

Central and satellite galaxies are compared within narrow stellar-mass bins. Host-halo-mass bins are used to examine the environmental dependence of satellite quenching.

## Python Analysis Scripts

The scripts should be executed in numerical order.

| File | Purpose |
|---|---|
| `01_catalogue_crossmatch.py` | Matches MPA-JHU and Tinker galaxies |
| `02_build_galaxy_sample.py` | Combines galaxy physical and environmental measurements |
| `03_data_cleaning.py` | Checks data quality and examines sample distributions |
| `04_quenching_analysis.py` | Calculates quenched fractions as functions of stellar mass |
| `05_environmental_analysis.py` | Investigates host-halo-mass dependence |
| `06_robustness_checks.py` | Tests observational selections, thresholds, and uncertainties |
| `07_final_results.py` | Generates the primary-sample results, tables, and final figures |

## Requirements

Python 3 with the following scientific packages:

- NumPy
- Astropy
- Matplotlib

Install them using:

```bash
python -m pip install numpy astropy matplotlib
```

The Python standard-library `csv`, `sys`, `site`, and `pathlib` modules are also used.

## Running the Analysis

1. Download the four MPA-JHU catalogues and the Tinker group catalogue from the sources above.
2. Extract compressed FITS files where necessary and use the exact filenames expected by the Python scripts.
3. Place the data files in the same working directory as the scripts, or update the file paths.
4. Run scripts `01` through `07` in numerical order.

The scripts generate intermediate files including:

- `galaxy_crossmatch.npz`
- `galaxy_combined_preclean.npz`
- `galaxy_quality_mask.npy`
- `galaxy_final_sample_mask.npy`

These intermediate files are used by subsequent scripts and are not necessary to store on GitHub if the original analysis can be reproduced.

## Main Results

The primary sample shows that satellite galaxies have higher quenched fractions than central galaxies in all five investigated stellar-mass bins.

Satellite-quenching efficiency ranges from approximately 0.425 to 0.493 across these bins.

The results also show a general increase in satellite quenched fraction with host-halo mass at fixed stellar mass, consistent with an environmental contribution to galaxy quenching.

### Final Scientific Figures

- `FINAL_FIGURE_1.pdf` — Quenched fraction versus stellar mass for central and satellite galaxies
- `FINAL_FIGURE_2.pdf` — Satellite quenched fraction versus host-halo mass at fixed stellar mass
- `FINAL_FIGURE_3.pdf` — Satellite-quenching efficiency as a function of stellar mass and host-halo mass

### Numerical Results

- `final_quenching_results.csv` — Quenched fractions, efficiencies, and confidence intervals
- `final_environmental_results.csv` — Environmental quenching measurements in stellar-mass and halo-mass bins

## Uncertainty Analysis

Statistical uncertainty is estimated using bootstrap resampling of galaxy groups, with 95% percentile confidence intervals.

Robustness tests include:

- Alternative quenching thresholds of −11.2, −11.0, and −10.8
- Different redshift selections
- Redshift adjustment of the central-galaxy comparison
- Exploratory inverse-\(V_{\max}\) weighting

The qualitative finding of excess satellite quenching persists across the tested quenching thresholds.

## Limitations

The analysis uses archival SDSS DR7 observations and inferred galaxy-group assignments. The principal limitations include approximate stellar-mass completeness, uncertainty in halo-mass and central/satellite classifications, and observational selection effects.

Quenched fractions alone cannot establish whether ram-pressure stripping, starvation, or another physical mechanism caused an individual galaxy to quench.

The Tinker catalogue record-count discrepancy is also documented and should be considered when reproducing the analysis.

## Scientific Background

Relevant literature includes:

- Gunn & Gott (1972) — Ram-pressure stripping
- Peng et al. (2012) — Environmental quenching and satellite-quenching efficiency
- Wetzel et al. (2013) — Satellite-quenching timescales
- Wright et al. (2022) — Environmental quenching in simulations
- Goubert et al. (2024) — Central and satellite quenching in observations and simulations

## Acknowledgements

This project uses publicly available SDSS DR7 data products, the MPA-JHU galaxy-property catalogues, and the Tinker galaxy group catalogue.

The original observational and catalogue data remain the work of their respective survey and research teams.

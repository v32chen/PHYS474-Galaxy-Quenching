# PHYS 474 — Environmental Dependence of Galaxy Quenching

## Project Overview

This project investigates how galaxy environment influences the suppression of star formation, with a focus on comparing central and satellite galaxies.

The main research question is:

**At fixed stellar mass, are satellite galaxies more likely to be quenched than central galaxies, and how does this excess quenching depend on the mass of their host dark-matter haloes?**

The analysis uses observational data from the Sloan Digital Sky Survey Data Release 7 (SDSS DR7). Galaxy physical properties from the MPA-JHU catalogues are combined with environmental information from the Tinker galaxy group catalogue.

The project uses Python for catalogue cross-matching, data cleaning, statistical analysis, uncertainty estimation, and scientific visualization.

This research was conducted as part of PHYS 474 (Fall 2026).

---

## Data Sources

Two publicly available astronomical data resources were used.

### 1. MPA-JHU SDSS DR7 Galaxy Catalogues

Website: https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/

The following four files provide galaxy physical properties:

| Catalogue | Description |
|---|---|
| `gal_info_dr7_v5_2.fit` | Galaxy positions, redshifts, and spectroscopic information |
| `totlgm_dr7_v5_2.fit` | Stellar-mass estimates |
| `gal_totsfr_dr7_v5_2.fits` | Total star-formation-rate estimates |
| `gal_totspecsfr_dr7_v5_2.fits` | Total specific star-formation-rate estimates |

The MPA-JHU measurements provide the quantities required to determine the star-formation activity of individual galaxies.

### 2. Tinker Galaxy Group Catalogue

Website: https://www.galaxygroupfinder.net/catalogs

Catalogue: `sdss_kdgroups_v1.0.dat`

This catalogue provides environmental information, including:

- Satellite probabilities
- Central and satellite classifications
- Estimated host-halo masses
- Galaxy group membership identifiers
- Maximum observable survey volumes ($V_{\max}$)

Galaxies are classified as satellites when their satellite probability satisfies $P_{\mathrm{sat}}>0.5$, and as centrals when $P_{\mathrm{sat}}<0.5$.

The version of the catalogue used in this analysis contains 559,028 records. The catalogue website describes 559,038 records; this discrepancy is documented as a potential data-version limitation.

**Note:** The original astronomical catalogues are not included in this repository because of their large file sizes. They can be obtained from the websites above.

---

## Research Methodology

### 1. Galaxy Catalogue Cross-Matching

The MPA-JHU and Tinker catalogues use different galaxy indexing systems. Therefore, galaxies are matched using their sky coordinates and redshifts rather than assuming that corresponding entries share the same row numbers.

The matching criteria are:

- Angular separation less than 1 arcsecond
- Absolute redshift difference less than 0.001
- A unique one-to-one correspondence between the catalogues

The cross-matching procedure identifies **498,431 unambiguous galaxy matches**.

### 2. Data Cleaning and Sample Selection

The combined data are screened for invalid measurements, unreliable redshifts, missing values, and relevant catalogue quality flags.

The final primary sample uses the redshift selection:

$$
0.02 \leq z < 0.05
$$

and the stellar-mass range:

$$
10.0 \leq \log_{10}(M_*/M_\odot) < 10.5
$$

An empirical stellar-mass completeness guideline is also applied to reduce observational selection effects. The resulting sample is treated as approximately mass-complete within the selected range, rather than assuming perfect completeness.

The final primary sample contains **16,602 galaxies**:

- **10,397 central galaxies**
- **6,205 satellite galaxies**

### 3. Galaxy Quenching Classification

The specific star-formation rate is defined as:

$$
\mathrm{sSFR} = \frac{\mathrm{SFR}}{M_*}
$$

where $\mathrm{SFR}$ is the star-formation rate and $M_*$ is the galaxy's stellar mass.

The specific star-formation-rate estimates are obtained from the MPA-JHU catalogue.

Galaxies are classified as quenched when:

$$
\log_{10}\left(\frac{\mathrm{sSFR}}{\mathrm{yr}^{-1}}\right) < -11
$$

Galaxies above this threshold are classified as star-forming for the purposes of this analysis.

The sensitivity of the results to this definition is tested using alternative thresholds of −11.2 and −10.8.

### 4. Quenched Fraction

The quenched fraction of a galaxy population is defined as:

$$
f_Q = \frac{N_{\mathrm{quenched}}}{N_{\mathrm{total}}}
$$

where $N_{\mathrm{quenched}}$ is the number of quenched galaxies and $N_{\mathrm{total}}$ is the total number of galaxies in the selected population.

Quenched fractions are calculated separately for central and satellite galaxies within narrow stellar-mass bins of width 0.1 dex.

### 5. Excess Satellite Quenching

The excess quenched fraction among satellite galaxies is:

$$
\Delta f_Q = f_{Q,\mathrm{sat}} - f_{Q,\mathrm{cen}}
$$

where $f_{Q,\mathrm{sat}}$ and $f_{Q,\mathrm{cen}}$ represent the satellite and central quenched fractions, respectively.

A positive value indicates a higher quenched fraction among satellites compared with centrals of similar stellar mass.

### 6. Satellite-Quenching Efficiency

The satellite-quenching efficiency is defined as:

$$
\epsilon_{\mathrm{sat}} =
\frac{f_{Q,\mathrm{sat}} - f_{Q,\mathrm{cen}}}
{1 - f_{Q,\mathrm{cen}}}
$$

This quantity measures the excess satellite quenched fraction relative to the fraction of central galaxies that remain star-forming.

The analysis follows the general environmental-quenching framework discussed by Peng et al. (2012).

### 7. Environmental Dependence

Satellite galaxies are divided into bins of host dark-matter halo mass.

The analysis examines how the satellite quenched fraction and satellite-quenching efficiency vary with host-halo mass while controlling approximately for stellar mass.

Host-halo masses are taken from the Tinker catalogue and expressed in units of $M_\odot/h$.

---

## Python Analysis Scripts

The repository contains seven Python analysis scripts.

| File | Description |
|---|---|
| `01_catalogue_crossmatch.py` | Matches galaxy records between the MPA-JHU and Tinker catalogues |
| `02_build_galaxy_sample.py` | Combines galaxy physical properties with environmental information |
| `03_data_cleaning.py` | Inspects data quality, removes unreliable measurements, and explores sample distributions |
| `04_quenching_analysis.py` | Calculates central and satellite quenched fractions as functions of stellar mass |
| `05_environmental_analysis.py` | Investigates satellite quenching as a function of host-halo mass |
| `06_robustness_checks.py` | Performs redshift, weighting, threshold, and uncertainty tests; defines the final primary sample |
| `07_final_results.py` | Calculates final results and generates the publication figures and numerical tables |

The scripts are designed to be run sequentially, since later stages use intermediate files generated by earlier stages.

---

## Software Requirements

The analysis uses Python 3 with the following scientific libraries:

- NumPy
- Astropy
- Matplotlib

Install these packages using:

`python -m pip install numpy astropy matplotlib`

The Python standard-library modules `csv`, `sys`, `site`, and `pathlib` are also used.

---

## Reproducing the Analysis

To reproduce the results:

1. Download the required SDSS DR7 MPA-JHU FITS catalogues.
2. Download the Tinker galaxy group catalogue.
3. Place the original data files in the working directory, using the filenames expected by the scripts.
4. Install the required Python packages.
5. Execute the analysis scripts in numerical order, from `01_catalogue_crossmatch.py` to `07_final_results.py`.

The analysis generates several intermediate files:

| File | Purpose |
|---|---|
| `galaxy_crossmatch.npz` | Matched galaxy row indices and matching diagnostics |
| `galaxy_combined_preclean.npz` | Combined physical and environmental galaxy measurements |
| `galaxy_quality_mask.npy` | Initial data-quality selection |
| `galaxy_final_sample_mask.npy` | Primary scientific sample selection |

These intermediate files are generated during the analysis and are not required to be stored in the repository when the original pipeline can be rerun.

---

## Main Scientific Results

The final sample contains 16,602 galaxies selected from SDSS DR7.

The results show that satellite galaxies have higher quenched fractions than central galaxies across all five investigated stellar-mass bins.

### Quenched Fraction and Satellite-Quenching Efficiency

| Stellar-mass bin | Central quenched fraction | Satellite quenched fraction | Satellite-quenching efficiency |
|---|---|---|---|
| 10.0–10.1 | 0.200 | 0.594 | 0.493 |
| 10.1–10.2 | 0.239 | 0.604 | 0.480 |
| 10.2–10.3 | 0.294 | 0.642 | 0.493 |
| 10.3–10.4 | 0.358 | 0.654 | 0.460 |
| 10.4–10.5 | 0.475 | 0.698 | 0.425 |

Both central and satellite quenched fractions generally increase with stellar mass.

However, the satellite population remains more strongly quenched at every stellar mass examined.

### Host-Halo-Mass Dependence

The results also show a general increase in satellite quenched fraction with host-halo mass, particularly toward group and cluster environments.

The trend remains visible when galaxies are compared within narrow stellar-mass intervals.

These findings are consistent with environmental processes contributing to the suppression of star formation in satellite galaxies.

However, the analysis does not directly establish which physical quenching mechanism is responsible.

---

## Scientific Figures

The repository contains three principal figures generated from the final sample.

### Figure 1 — Quenched Fraction Versus Stellar Mass

**File:** `FINAL_FIGURE_1.pdf`

Compares central and satellite quenched fractions in five stellar-mass bins. Error bars represent 95% group-bootstrap confidence intervals.

### Figure 2 — Satellite Quenching Versus Host-Halo Mass

**File:** `FINAL_FIGURE_2.pdf`

Shows the satellite quenched fraction as a function of host-halo mass for galaxies with:

$$
10.2 \leq \log_{10}(M_*/M_\odot) < 10.3
$$

A central-galaxy reference fraction is included for comparison.

### Figure 3 — Satellite-Quenching Efficiency Heatmap

**File:** `FINAL_FIGURE_3.pdf`

Illustrates satellite-quenching efficiency as a function of both stellar mass and host-halo mass.

Blank cells correspond to bins containing fewer than 30 satellite galaxies.

---

## Numerical Results

The repository includes two CSV files containing the final quantitative results.

### `final_quenching_results.csv`

Contains:

- Stellar-mass bins
- Central and satellite galaxy counts
- Quenched galaxy counts
- Quenched fractions
- Excess satellite quenched fractions
- Satellite-quenching efficiencies
- 95% bootstrap confidence intervals for quenching efficiency

### `final_environmental_results.csv`

Contains:

- Stellar-mass and host-halo-mass bins
- Satellite galaxy counts
- Quenched satellite counts
- Satellite quenched fractions
- Satellite-quenching efficiencies

These tables provide the numerical values used to construct the final scientific figures.

---

## Statistical Uncertainties and Robustness Tests

### Group-Bootstrap Uncertainties

Statistical uncertainties are estimated through bootstrap resampling of galaxy groups rather than individual galaxies.

This method accounts approximately for correlations between galaxies assigned to the same dark-matter halo.

The analysis uses 500 bootstrap realizations for the final primary-sample results and reports percentile-based 95% confidence intervals.

### Quenching-Threshold Sensitivity

The analysis tests three quenching thresholds:

- −11.2
- −11.0
- −10.8

Satellite-quenching efficiency remains positive across all five stellar-mass bins under each tested threshold.

### Redshift Dependence

Quenching measurements are compared across different redshift intervals to investigate sensitivity to observational selection.

Additional redshift standardization is performed to assess the effects of differences between central and satellite redshift distributions.

### Survey-Volume Weighting

An exploratory inverse-$V_{\max}$ analysis examines how the quenching measurements respond to survey-volume weighting.

The weighted and unweighted estimates are compared as a diagnostic of observational selection effects.

This weighting test is not considered a complete correction for survey incompleteness.

---

## Limitations

Several limitations should be considered when interpreting the findings.

**Survey completeness:** The stellar-mass completeness selection is approximate. Magnitude limits and differing galaxy mass-to-light ratios may affect the observed quenched fractions.

**Central/satellite classifications:** The galaxy group catalogue provides inferred environmental classifications, which may contain uncertainties.

**Host-halo masses:** Halo masses are catalogue estimates rather than direct measurements and may have systematic uncertainties.

**Quenching definition:** The exact satellite-quenching efficiency depends on the adopted specific star-formation-rate threshold.

**Catalogue consistency:** The analyzed Tinker group catalogue contains 559,028 records, whereas the documentation describes 559,038 records.

**Physical interpretation:** The observed quenching trends are consistent with environmental effects such as ram-pressure stripping and starvation, but the measurements alone cannot uniquely identify the physical mechanism.

These limitations should be considered when comparing the results with previous observational and simulation studies.

---

## Scientific References

The scientific background and methodology draw on the following literature:

- Gunn & Gott (1972). Foundational work on ram-pressure stripping.
- Peng et al. (2012). Environmental quenching and satellite-quenching efficiency.
- Wetzel et al. (2013). Satellite galaxy quenching histories and timescales.
- Zu & Mandelbaum (2015). Galaxy stellar-mass selection and completeness.
- Wright et al. (2022). Environmental quenching processes in galaxy simulations.
- Goubert et al. (2024). Comparison of central and satellite quenching in observations and simulations.

Full bibliographic references and comparisons with previous studies are provided in the accompanying PHYS 474 conference paper.

---

## Acknowledgements

This project uses publicly available observational data from the Sloan Digital Sky Survey (SDSS), galaxy physical-property measurements from the MPA-JHU collaboration, and environmental information from the Tinker galaxy group catalogue.

The observational surveys and original catalogue measurements are credited to their respective research teams.

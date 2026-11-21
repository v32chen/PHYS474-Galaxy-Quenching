# PHYS 474 — Environmental Dependence of Galaxy Quenching

## About This Project

This project explores how a galaxy's environment affects its ability to form new stars. In particular, I wanted to investigate whether satellite galaxies are more likely to stop forming stars than central galaxies, even when they have similar stellar masses.

The main question I am trying to answer is:

**At fixed stellar mass, are satellite galaxies more likely to be quenched than central galaxies, and how does this depend on the mass of their host dark-matter halo?**

To investigate this, I used observational data from the Sloan Digital Sky Survey (SDSS DR7). I combined measurements of galaxy stellar masses and star-formation rates with information about whether each galaxy is a central or satellite and the mass of its surrounding dark-matter halo.

The analysis was carried out in Python, from matching the original galaxy catalogues to calculating quenched fractions, estimating uncertainties, and producing the final graphs.

This research was completed as part of **PHYS 474 (Fall 2026)**.

---

## Where the Data Came From

I used two main astronomical data sources.

### 1. MPA-JHU SDSS DR7 Catalogues

Website: https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/

The MPA-JHU catalogues contain measurements of the physical properties of galaxies observed by SDSS. For this project, I used four files:

| File | Information provided |
|---|---|
| `gal_info_dr7_v5_2.fit` | Galaxy coordinates, redshifts, and spectroscopic information |
| `totlgm_dr7_v5_2.fit` | Stellar masses |
| `gal_totsfr_dr7_v5_2.fits` | Total star-formation rates |
| `gal_totspecsfr_dr7_v5_2.fits` | Specific star-formation rates |

These measurements let me determine how massive each galaxy is and how actively it is forming stars.

### 2. Tinker Galaxy Group Catalogue

Website: https://www.galaxygroupfinder.net/catalogs

File: `sdss_kdgroups_v1.0.dat`

The Tinker catalogue provides information about the environments in which galaxies live.

I mainly used it to identify:

- Whether a galaxy is classified as a central or satellite
- The probability that a galaxy is a satellite
- The estimated mass of its host dark-matter halo
- The group to which it belongs

For this analysis, galaxies with satellite probabilities above 0.5 are classified as satellites, while those below 0.5 are classified as centrals.

The catalogue version used here contains 559,028 records, although the website describes 559,038 records. I have noted this difference as a limitation of the dataset.

The original astronomical catalogues are not uploaded to this repository because of their large sizes, but they can be accessed through the websites linked above.

---

## How I Analyzed the Data

### 1. Matching the Galaxy Catalogues

The first challenge was identifying which galaxies in the MPA-JHU catalogues corresponded to galaxies in the Tinker catalogue.

Since the catalogues do not use the same row ordering, I matched the galaxies using their sky coordinates and redshifts.

I required:

- An angular separation of less than 1 arcsecond
- A redshift difference smaller than 0.001
- A unique one-to-one match between the two catalogues

This produced **498,431 matched galaxy records**.

### 2. Cleaning and Selecting the Galaxies

After matching the catalogues, I removed records with invalid or unreliable measurements using the available quality flags.

I also needed to limit the sample to reduce observational selection effects. More distant galaxies generally need to be brighter to be included in a magnitude-limited survey, which can make lower-mass galaxies underrepresented.

For the main analysis, I selected galaxies with:

$$
0.02 \leq z < 0.05
$$

and stellar masses satisfying:

$$
10.0 \leq \log_{10}(M_*/M_\odot) < 10.5
$$

I also applied an empirical stellar-mass completeness guideline. The resulting sample is treated as approximately mass-complete, although some selection effects may remain.

The final sample contains **16,602 galaxies**, consisting of:

- 10,397 central galaxies
- 6,205 satellite galaxies

### 3. Identifying Quenched Galaxies

To determine whether a galaxy is quenched, I used its specific star-formation rate (sSFR).

This is the star-formation rate divided by stellar mass:

$$
\mathrm{sSFR} = \frac{\mathrm{SFR}}{M_*}
$$

A lower sSFR means that a galaxy is forming fewer new stars relative to the mass of stars it already contains.

For this project, galaxies were classified as quenched when:

$$
\log_{10}\left(\frac{\mathrm{sSFR}}{\mathrm{yr}^{-1}}\right) < -11
$$

The sSFR measurements were obtained from the MPA-JHU catalogue.

I also tested slightly different thresholds to check whether this choice affected the conclusions.

### 4. Calculating Quenched Fractions

The quenched fraction tells us what proportion of galaxies in a given population are no longer forming stars at a significant rate.

It is calculated using:

$$
f_Q = \frac{N_{\mathrm{quenched}}}{N_{\mathrm{total}}}
$$

I calculated this separately for central and satellite galaxies in five stellar-mass bins, each 0.1 dex wide.

Comparing galaxies within narrow mass ranges is important because stellar mass itself affects how likely a galaxy is to be quenched.

To measure the difference between the two populations, I calculated:

$$
\Delta f_Q = f_{Q,\mathrm{sat}} - f_{Q,\mathrm{cen}}
$$

I also calculated the satellite-quenching efficiency:

$$
\epsilon_{\mathrm{sat}} =
\frac{f_{Q,\mathrm{sat}} - f_{Q,\mathrm{cen}}}
{1 - f_{Q,\mathrm{cen}}}
$$

This measures the excess satellite quenched fraction relative to the fraction of central galaxies that are still star-forming.

### 5. Studying the Effect of Environment

After comparing central and satellite galaxies, I investigated whether the quenched fraction of satellites changes with their host-halo mass.

To do this, I divided satellite galaxies into halo-mass bins while keeping their stellar masses within narrow ranges.

This allowed me to examine how the environmental quenching signal changes between galaxies living in relatively small haloes and those in more massive galaxy groups and clusters.

---

## Python Files

I organized the analysis into seven Python scripts.

| Script | What it does |
|---|---|
| `01_catalogue_crossmatch.py` | Matches galaxies between the MPA-JHU and Tinker catalogues |
| `02_build_galaxy_sample.py` | Combines the galaxy measurements into one dataset |
| `03_data_cleaning.py` | Checks the data quality and examines the galaxy distributions |
| `04_quenching_analysis.py` | Calculates quenched fractions for central and satellite galaxies |
| `05_environmental_analysis.py` | Investigates how quenching depends on host-halo mass |
| `06_robustness_checks.py` | Tests the reliability of the results and selects the primary sample |
| `07_final_results.py` | Produces the final results, uncertainty estimates, and figures |

The scripts should be run in numerical order because later scripts use the outputs of earlier ones.

### Software Used

The project was completed using Python 3 and the following libraries:

- NumPy
- Astropy
- Matplotlib

They can be installed using:

```bash
python -m pip install numpy astropy matplotlib
```

### Reproducing the Analysis

To reproduce the analysis, download the original MPA-JHU and Tinker catalogues from the websites above and place them in the working directory using the filenames expected by the scripts.

Run the seven Python scripts in numerical order.

The analysis generates intermediate files such as `galaxy_crossmatch.npz`, `galaxy_combined_preclean.npz`, `galaxy_quality_mask.npy`, and `galaxy_final_sample_mask.npy`.

These files store the matched records and sample selections used in the later calculations.

---

## Main Results

The results show that satellite galaxies have higher quenched fractions than central galaxies across all five stellar-mass bins examined.

### Comparing Central and Satellite Galaxies

| Stellar-mass bin | Central quenched fraction | Satellite quenched fraction | Quenching efficiency |
|---|---|---|---|
| 10.0–10.1 | 0.200 | 0.594 | 0.493 |
| 10.1–10.2 | 0.239 | 0.604 | 0.480 |
| 10.2–10.3 | 0.294 | 0.642 | 0.493 |
| 10.3–10.4 | 0.358 | 0.654 | 0.460 |
| 10.4–10.5 | 0.475 | 0.698 | 0.425 |

Two trends stand out.

First, both central and satellite galaxies become more likely to be quenched as stellar mass increases.

Second, satellites consistently have higher quenched fractions than centrals, even when their stellar masses are similar.

### Dependence on Host-Halo Mass

I also found that satellite galaxies generally become more likely to be quenched as their host-halo mass increases, particularly in more massive group and cluster environments.

This suggests that the environment surrounding a satellite galaxy plays an important role in its star-formation activity.

These trends are consistent with environmental quenching mechanisms such as ram-pressure stripping and starvation, although the analysis does not directly identify which mechanism is responsible.

---

## Final Figures

The three main figures from the analysis are included in this repository.

**Figure 1 — Quenched Fraction vs. Stellar Mass**

`FINAL_FIGURE_1.pdf`

Compares the quenched fractions of central and satellite galaxies across five stellar-mass bins, including 95% group-bootstrap confidence intervals.

**Figure 2 — Satellite Quenching vs. Host-Halo Mass**

`FINAL_FIGURE_2.pdf`

Shows how the satellite quenched fraction changes with host-halo mass for galaxies within the stellar-mass range:

$$
10.2 \leq \log_{10}(M_*/M_\odot) < 10.3
$$

A central-galaxy quenched fraction is included as a reference.

**Figure 3 — Satellite-Quenching Efficiency**

`FINAL_FIGURE_3.pdf`

Shows how satellite-quenching efficiency depends on both stellar mass and host-halo mass.

Blank regions indicate bins with fewer than 30 satellites, where measurements were not reported.

### Numerical Results

The numerical values used in the figures are saved in:

- `final_quenching_results.csv`
- `final_environmental_results.csv`

---

## Checking the Reliability of the Results

An important part of the project was checking whether the conclusions depended strongly on the analysis choices.

I carried out several tests to investigate this.

**Statistical uncertainty:** I used bootstrap resampling of galaxy groups to calculate 95% confidence intervals. This helps account for the fact that galaxies in the same group may not be statistically independent.

**Quenching threshold:** I repeated the quenching calculations using thresholds of −11.2, −11.0, and −10.8. The measured efficiencies changed, but the satellite-quenching excess remained positive in every tested stellar-mass bin.

**Redshift selection:** I compared different redshift intervals to investigate whether the observed trends were sensitive to the distance range of the sample.

**Redshift adjustment:** I adjusted the central-galaxy reference fractions to better match the redshift distributions of the satellite populations.

**Survey-volume weighting:** I also compared unweighted measurements with exploratory inverse-$V_{\max}$ weighted estimates.

These checks helped me identify which trends were relatively stable and where observational selection could influence the measurements.

---

## Limitations

Although the results support an environmental contribution to galaxy quenching, there are several limitations.

The SDSS sample is magnitude-limited, meaning some galaxies may be underrepresented at greater distances. The stellar-mass completeness selection reduces this concern but does not eliminate it entirely.

Central/satellite classifications and host-halo masses are also estimates obtained from a galaxy group finder, rather than direct measurements without uncertainty.

The measured quenching efficiencies depend partly on the chosen sSFR threshold, and the reported statistical confidence intervals do not include every possible systematic error.

Finally, the quenched fraction tells us whether star formation has been suppressed, but not exactly how it happened. Distinguishing mechanisms such as ram-pressure stripping and starvation would require additional measurements or analysis.

---

## Scientific Background and References

This project builds on previous observational and theoretical research into galaxy quenching.

Some of the main studies used to guide the analysis include:

- Gunn & Gott (1972) — Ram-pressure stripping
- Peng et al. (2012) — Environmental quenching and satellite-quenching efficiency
- Wetzel et al. (2013) — Satellite-quenching timescales
- Zu & Mandelbaum (2015) — Galaxy stellar-mass selection and completeness
- Wright et al. (2022) — Environmental quenching in simulations
- Goubert et al. (2024) — Central and satellite quenching in observations and simulations

These studies provide the scientific background needed to interpret the results and compare them with previous work.

---

## Acknowledgements

This project uses publicly available data from the Sloan Digital Sky Survey, the MPA-JHU galaxy-property catalogues, and the Tinker galaxy group catalogue.

I acknowledge the survey teams and researchers responsible for collecting, processing, and publishing these data.

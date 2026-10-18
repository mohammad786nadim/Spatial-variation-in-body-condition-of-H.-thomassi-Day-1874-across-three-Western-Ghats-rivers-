# Spatial variation in body condition of *Hypselobarbus thomassi* (Day, 1874) across three Western Ghats rivers

## Overview

This repository contains the R script used to examine length–weight relationships (LWR), length–length relationships (LLR), and body condition of *Hypselobarbus thomassi* (Day, 1874) from three rivers in the Western Ghats, India: Netravati, Periyar, and Kallada.

The analyses focus on differences in morphometric relationships and body condition among the three river populations. They also provide a statistical framework for comparing growth patterns and assessing variation in fish condition across rivers.

## Statistical analyses

The R script includes the following analyses:

1. **Length–weight relationships (LWR):** Length and weight measurements were used to estimate the regression coefficients (a and b) and coefficients of determination (R²). The relationships were fitted using log-transformed data.

2. **Length–length relationships (LLR):** Relationships among total length (TL), standard length (SL), and fork length (FL) were examined using linear regression. Outliers were identified and removed based on studentized residuals before fitting the final regression models.

3. **Body condition indices:** Fulton's condition factor (K), relative condition factor (Kn) and scaled mass index (SMI) were calculated to assess and compare fish condition among the three rivers.

4. **Analysis of covariance (ANCOVA):** ANCOVA was used to examine differences in regression slopes among river populations, followed by pairwise comparisons.

5. **Scaled Mass Index (SMI):** The SMI was calculated using standardised major axis (SMA) regression to assess body condition while accounting for differences in body length.

6. **Statistical comparisons:** Shapiro–Wilk tests were used to assess normality, while Kruskal–Wallis tests and Dunn's post hoc tests with Benjamini–Hochberg adjustment were used to compare body condition indices among rivers.

7. **Data visualisation:** The script generates individual and combined regression plots, along with boxplots illustrating variation in body condition indices among river populations.

## Software and R packages

All analyses were carried out using R.

The following R packages are required:

- `dplyr`
- `ggplot2`
- `tidyr`
- `emmeans`
- `FSA`
- `smatr`

The required packages can be installed using:

```r
install.packages(c("dplyr", "ggplot2", "tidyr", "emmeans", "FSA", "smatr"))
```

## Repository contents

The repository contains the following files:

- `LWR, LLR, SMI.R` – R script used for morphometric analyses, calculation of body condition indices, statistical comparisons, and visualisation.
- `README.md` – Overview of the analyses and instructions for using the script.
- `LICENSE` – Terms and conditions governing the use and redistribution of the code.

## Input data

The analyses are based on morphometric measurements of *H. thomassi* collected from the Netravati, Periyar, and Kallada rivers.

The script requires the following CSV files:

- `Netravati.csv`
- `Periyar.csv`
- `Kallada.csv`
- `Netravati_LL.csv`
- `Periyar_LL.csv`
- `Kallada_LL.csv`

These datasets are not included in the repository.

## Running the analysis

To reproduce the analyses using the required datasets:

1. Download or clone the GitHub repository.
2. Install the required R packages.
3. Place the six input CSV files in the working directory.
4. Open `LWR, LLR, SMI.R` in RStudio.
5. Run the script to obtain the statistical results and figures.

The input files must contain the length and weight measurements in the format expected by the script.

## Citation

Researchers who use or adapt this R script are requested to acknowledge the source by citing the associated research publication or the archived software release.

The repository will be archived in Zenodo, and the DOI and recommended citation will be added once they are available.

## License

This repository is made available under the MIT License, which permits the use, modification, and redistribution of the code, provided that the original copyright notice and licence terms are retained. See the `LICENSE` file for details.

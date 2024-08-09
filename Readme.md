
tlaepple@awi
2024_08_08

This directory contains the code to reproduce the results and figures
for Laepple et al., Matters Arising in response to Jones et al., Nature 2024

The code requires some libraries that can be installed via CRAN:

zoo; tidyr; ggplot2 and dplyr

if further requires PaleoSpec (an R package of the AWI Earth System Diagnostics lab to assist in the spectral analysis of timeseries) that can be installed via github
https://github.com/EarthSystemDiagnostics/paleospec

The data folder contains all the dataset used for the analysis; if they are not the original
datasets supplied by the original authors, a header is added to explain their origin.

The code to do the calculations and produce the figures are in the src folder
src/Figure1andExtendedFigure1.R
src/Figure2andExtendedDataFigure2.R
src/Figure3andExtendedDataFigure3_DiffusionLength.R

In these standalone files, please adapt the path ('basedrive') to the path where you extracted directories



This directory contains `R` code to reproduce the results and figures for the
Laepple et al. Matters Arising in response to the Jones et al., Nature 2023
publication (doi:
[10.1038/s41586-022-05411-8](https://doi.org/10.1038/s41586-022-05411-8).

The analysis code was written by Thomas Laepple, with contributions to the
library functions by Thomas Münch; both Alfred Wegener Institute, Helmholtz
Centre for Polar and Marine Research (AWI).

The code requires some `R` libraries that can be installed via CRAN:
```
install.packages(zoo)
install.packages(tidyr)
install.packages(ggplot2)
install.packages(dplyr)
```

It further requires `PaleoSpec` (an `R` package of the AWI Earth System
Diagnostics lab to assist in the spectral analysis of timeseries) that can be
installed from its [GitHub
repository](https://github.com/EarthSystemDiagnostics/paleospec).

The data folder contains all the datasets used for the analysis; if they are not
the original datasets supplied by the original authors, a header was added to
explain the data origin.

The code to do the calculations and produce the figures are in the `src` folder:
```
src/Figure1andExtendedFigure1.R
src/Figure2andExtendedDataFigure2.R
src/Figure3andExtendedDataFigure3_DiffusionLength.R
```

In these standalone files, please adapt the path (`'basedrive'`) to the path
of the directory into which you downloaded or cloned this project.

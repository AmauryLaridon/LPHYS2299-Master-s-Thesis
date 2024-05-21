---
# SURFER v3.2
---

<p align="center">
<img src="https://github.com/AmauryLaridon/LPHYS2299-Master-s-Thesis/assets/58213378/5ced037d-3522-4b51-a871-87df8b400603" width="600" />
</p>

SURFER: a simple climate model linking CO2 emissions and Solar Radiation Management to sea level rise and ocean acidification.

For a full explanation of the model please check: *LPHYS2299___MastersThesis.pdf* in the main repository

We present the SURFER code and examples to reproduce all figures in the paper using a **jupyter notebook**: https://jupyter.org/.

### What is new in the version 3.2

- AMOC-GIS Cascade Calibration Module
- Representation of the impacts on the carbon cycle of an AMOC collapse

### Requisites to running the jupyter notebook

The jupyter notebook containing the code and examples has been written in Python

### Data from other sources

The repository contains data from other sources. The data is used to force the model (emission scenarios data) or to compare SURFER's output with other models (ZECMIP data).

#### SSP emission scenarios

The data contained in folder `scenarios/SSP/` is based on the SSP database hosted by the IIASA Energy Program at https://tntcat.iiasa.ac.at/SspDb.

#### RCP emission scenarios

The data contained in folder `scenarios/RCP/` is based on the RCP databased hosted by IIASA Energy Program at https://tntcat.iiasa.ac.at/RcpDb

#### ZECMIP results

The data contained in folder `ZECMIP_data/` is freely available and has been downloaded from http://terra.seos.uvic.ca/ZECMIP/index.html

#### Data not in the repository

In the paper, SURFER's output is contrasted to results from other references:

[1] A. M. R. Bakker, T. E. Wong, K. L. Ruckert, and K. Keller, “Sea-level projections representing the deeply uncertain contribution of the West Antarctic ice sheet,” Sci. Rep., vol. 7, no. 1, pp. 1–7, 2017.

[2] T. Frederikse et al., “The causes of sea-level rise since 1900,” Nature, vol. 584, no. 7821, pp. 393–397, 2020.

[3] A. Robinson, R. Calov, and A. Ganopolski, “Multistability and critical thresholds of the Greenland ice sheet,” Nat. Clim. Chang., vol. 2, no. 6, pp. 429–432, 2012.

[4] J. Van Breedam, H. Goelzer, and P. Huybrechts, “Semi-equilibrated global sea-level change projections for the next 10 000 years,” Earth Syst. Dyn., vol. 11, no. 4, pp. 953–976, 2020.

Data in Refs. [1] and [2] can be found by following the references and data in Refs. [3] and [4] was obtained by private communication with the authors.

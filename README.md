# Introduction

Modeling photovoltaic (PV) output using models with hourly resolution can overestimate output by a few percent; this is the Average-then-Clip (AtC) Bias [1]. 
Use the code here to correct hourly modeling output for the AtC Bias as described in [2] and [3].

This project consists of Python scripts and associated data in CSV files. Users are expected be familiar with Python and PV modeling. The users specifies file locations and other inputs in the User Input section of each script and runs these in a Python console. This 
project uses pvlib [4] for PV calculations.

The AtC Bias correction method in [2], also described as "Classification Method 2" in [3], uses a correction matrix based on clipping potential, CP, and a clearness index, gamma_DNI or simply gamma. The matrix used in [2] and [3], as well as NREL SAM 2023.12.17 [5], is provided here in the file [method2.csv](method2.csv). Sections below provide short desciptions of how to apply this or a similar matrix, and how to create a new one, for example, using additional observation data beyond what was included in [2] and [3].

# Validation

This AtC Bias correction method was demonstrated and compared with a model based on minute-scale data in [6].

# Apply Matrix
Estimate AtC Bias by applying the correction matrix to hourly PV modeling output using the method described in [2] and [3]. The steps in the process are:

1. Run tmy_refsky.py to generate a reference sky radiation time series in the TMY3 format.
2. Run PV modeling software like SAM or PVLib using hourly inputs for these cases:
    - Plant specifications with TMY meteorological data
    - Plant with additional modules so dc:ac ratio is 1.0 with TMY meteorological data
    - Plant specifications with reference sky meteorological data
3. Export time series of the 3 model runs as CSV files. These contain:
    - SAM
      - Time stamp
      - Irradiance DNI from weather file | (W/m2)
      - Inverter DC input power | (kW)
      - Inverter AC output power | (kW)
    - PVSyst
      - date
      - H_sol
      - GlobHor
      - DiffHor
      - EOutArray
      - EOutInv
4. Run atcbias_apply_matrix.py to generate estimates AtC Bias

Note that the hourly estimates of average scaled AtC Bias are expected to 
match the annual (and perhaps month-hour) average AtC Bias, but are not
expected to match the hourly values. See [2] and [3] for a detailed explanation.

# Create Matrix

Run atcbias_create_matrix.py to calculate a correction matrix for the AtC 
bias based on minute-scale Pdc inputs using the method described in [2] and [3]. 

The user provides inverter specifications and multiple CSV files with minute PV 
data. The CSV files are named [site]_[pv_system].csv. These contain the columns:

- time: ISO datetime string in format %Y%m%d%H%M%S, e.g. 20090101144800
- zenith: Solar zenith angle in degrees
- dni: Direct Normal Irradiance, DNI, (W/m2)
- vdc: voltage from pv installation in V
- pdc: dc power from pv installation in W
- dni_refsky: DNI for reference sky, i.e. turbidity = 1.0 (W/m2)
- pdc_refsky: dc power from pv installation under reference sky

These files should contain whole years of data to avoid seasonal weighting. 
Times are in increments of one minute. Time zones do not affect the matrix
creation. pdc_refsky must have the same units as pdc.

The AtC Bias correction matrix will be saved to a CSV file which can be used
to estimate annual AtC Bias, for example using atcbias_apply_matrix.py.

# References
[1] J. O. Allen and W. B. Hobbs, "The effect of short-term inverter saturation on modeled
hourly PV output using minute DC power measurements
", Journal of Renewable and Sustainable Energy 14, 063503, 2022.
https://doi.org/10.1063/5.0130265

[2] J. O. Allen, W. B. Hobbs, and M. Bolen, "Classification Method to Predict the Effect of Short-Term Inverter Saturation on PV Performance Modeling", PV Performance Modeling Workshop Salt Lake City, 2022.

[3] J. O. Allen, “Improved PV Plant Energy Production
(Phases 1 and 2): The Effect of Short-term Inverter Saturation on PV Performance Modeling”. EPRI Technical Report
3002018708. EPRI, Palo Alto, CA. 2022.
https://www.epri.com/research/products/000000003002018708. 

[4] Anderson, K., Hansen, C., Holmgren, W., Jensen, A., Mikofski, M., and Driesse, A. “pvlib python: 2023 project update.” Journal of Open Source Software, 8(92), 5994, (2023). [DOI: 10.21105/joss.05994](http://dx.doi.org/10.21105/joss.05994).

[5] System Advisor Model Version 2023.12.17, Revision 1 (SAM 2023.12.17r1). National Renewable Energy Laboratory. Golden, CO. https://https://sam.nrel.gov .

[6] Prilliman, M. J., J. M. F. Keith, and W. B. Hobbs. 2024. Side-by-Side
Comparison of Subhourly Clipping Models: Preprint. Golden, CO: National Renewable
Energy Laboratory. NREL/CP-7A40-90054. https://www.nrel.gov/docs/fy24osti/90054.pdf
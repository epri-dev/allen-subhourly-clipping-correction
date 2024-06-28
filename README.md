# Introduction

Modeling photovoltaic (PV) output using models with hourly resolution can overestimate output by a few percent; this is the Average-then-Clip (AtC) Bias [1]. 
Use the code here to correct hourly modeling output for the AtC Bias as described in [2].

This project consists of Python scripts and associated data in CSV files. Users are expected be familiar with Python and PV modeling. The users specifies file locations and other inputs in the User Input section of each script and runs these in a Python console. This 
project uses pvlib [https://pvlib-python.readthedocs.io/en/stable/#] for PV calculations.

# Apply Matrix
Estimate AtC Bias by applying the correction matrix to hourly PV modeling output using the method described in [2]. The steps in the process are:

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
expected to match the hourly values. See [2] for a detailed explanation.

# Create Matrix

Run atcbias_create_matrix.py to calculate a correction matrix for the AtC 
bias based on minute-scale Pdc inputs using the method described in [2]. 

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

[2] J. O. Allen, W. B. Hobbs, and M. Bolen, "Classification Method to
   Predict the Effect of Short-Term Inverter Saturation on PV Performance
   Modeling", PV Performance Modeling Workshop Salt Lake City, 2022.




#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to create correction matrix for Average then Clip (AtC) bias
"""

import os
import errno
import re
import copy
import numpy as np
import pandas as pd

"""
BSD 3-Clause License

Copyright (c) 2024 Allen Analytics LLC

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

  Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

__author__ = "Jonathan Allen"
__copyright__ = "Copyright 2024, Allen Analytics LLC"
__license__ = "BSD"
__email__ = "jon@allen-analytics.com"

"""
Project funded by EPRI, see ***
"""

"""
Use this script to calculate a correction matrix for the Average then Clip (AtC) 
bias using the method described in [1]. 

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

References
----------
.. [1] J. O. Allen, W. B. Hobbs, and M. Bolen, "Classification Method to
   Predict the Effect of Short-Term Inverter Saturation on PV Performance
   Modeling", PV Performance Modeling Workshop Salt Lake City, 2022.

"""

# ------------------------------------------------------------------------------
# User Inputs
# ------------------------------------------------------------------------------
# Sites and PV systems to use for matrix calculation
sites = ["PSU", "SXF", "TBL"]
pv_systems = ["PV1Axis", "PVS25"]

# Location of input csv files, and output csv file
project_dir = r"C:\a2\epri_atcbias3\data\rabin"

# Nominal dc capacity of modules connected to inverter
# module Pdc0 * number of modules * number of strings
# units must match pdc values in csv files
dc_capacity = 320 * 384 * 18  # W

# dc:ac ratios to include in matrix calculations
dcac_ratios = np.arange(1.2, 1.8, 0.1)

# Inverter specifications using Sandia model
inverter = {
    "Model": "Sandia",
    "Name": "SMA America: SC-2200-US 665V [CEC 2016]",
    "Vac": 665,
    "Pac0": 2079000,
    "Pdc0": 2130260,
    "Vdc0": 665,
    "Ps0": 4991.09,
    "C0": -7.83848E-09,
    "C1": 1.41278E-05,
    "C2": 0.00210083,
    "C3": 0.000112255,
}

# One minus - Combined ac losses
ac_derate = 0.99

# AtC Bias Correction Matrix output file
matrix_fn = os.path.join(project_dir, "matrix.csv")

cp_bins = [-2.0, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4,
           0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6]

gamma_bins = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
              0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.10]


# ------------------------------------------------------------------------------
# Utility Functions
# ------------------------------------------------------------------------------
def sandia_inverter(v_dc, p_dc, inv_param, inv_scale=1.0):
    """
    Calculate the inverter AC power using the Sandia model

    Parameters
    ----------
    v_dc : numeric
        DC voltage input to the inverter. [V]

    p_dc : numeric
        DC power input to the inverter. [W]

    inv_param : dict-like
        Parameters for the inverter model in [1].

    inv_scale : float
        Inverter capacity scale factor from [2].

    Returns
    -------
    p_ac : numeric
        AC power output. [W]

    Notes
    -----

    This function differs from the pv_lib version in that:
        - inverter sized can be scaled
        - unclipped power is returned
        - parameters are named with tailing zeros, not 'oh's.

    References
    ----------
    [1] D. King, S. Gonzalez, G. Galbraith, W. Boyson, "Performance Model
       for Grid-Connected Photovoltaic Inverters", SAND2007-5036, Sandia
       National Laboratories.
    [2] J. O. Allen and W. B. Hobbs, "The effect of short-term inverter
        saturation on modeled hourly PV output using minute DC power
        measurements", Journal of Renewable and Sustainable Energy
        14, 063503 (2022). https://doi.org/10.1063/5.0130265
    """

    # Scale constants in Sandia inverter model
    pac0 = inv_param["Pac0"] * inv_scale
    pdc0 = inv_param["Pdc0"] * inv_scale
    ps0 = inv_param["Ps0"] * inv_scale
    c0 = inv_param["C0"] / inv_scale

    a = pdc0 * (1 + inv_param["C1"] * (v_dc - inv_param["Vdc0"]))
    b = ps0 * (1 + inv_param["C2"] * (v_dc - inv_param["Vdc0"]))
    c = c0 * (1 + inv_param["C3"] * (v_dc - inv_param["Vdc0"]))

    pac_unclipped = (pac0 / (a - b) - c * (a - b)) * (p_dc - b) + c * (p_dc - b) ** 2

    # Clip Pac
    pac_clipped = np.minimum(pac0, pac_unclipped)

    return pac_clipped, pac_unclipped


# ------------------------------------------------------------------------------
# Create 3-Level Dictionary to Store Results
# ------------------------------------------------------------------------------
# Results are stored in x, a 3-level dictionary with these levels:
#   1 - site
#   2 - pv_system
#   3 - content:
#       - TMYFileName
#       - MinuteFileName
#       - NominalPac [scalar]
#       - ValidHour [Hourly Pandas Series of bools]
#       - Gamma [Hourly Pandas Series of floats]
#       - CP (Clipping Potential) [NumPy Array, Num DCAC Ratios x Num Hours]
#       - PacCtA (AC Power, Clip then Average) [NumPy Array, Num DCAC Ratios x Num Hours]
#       - PacAtC (AC Power, Average then Clip) [NumPy Array, Num DCAC Ratios x Num Hours]
#       - Delta (Scaled AtC Bias) [NumPy Array, Num DCAC Ratios x Num Hours]
x = {}

# ------------------------------------------------------------------------------
# Add Data Source File Names and Verify That These Exist
# ------------------------------------------------------------------------------
for site in sites:
    x[site] = {}
    for pv_system in pv_systems:
        x[site][pv_system] = {}

        fn = f'{site}_{pv_system}_TMY.csv'
        fn = os.path.join(project_dir, fn)
        if os.path.isfile(fn):
            x[site][pv_system]['TMYFileName'] = fn
        else:
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), fn)

        fn = f'{site}_{pv_system}_Minute.csv'
        fn = os.path.join(project_dir, fn)
        if os.path.isfile(fn):
            x[site][pv_system]['MinuteFileName'] = fn
        else:
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), fn)

# ------------------------------------------------------------------------------
# Calculate Nominal PV Output - Full Year TMY with DC:AC = 1.0
# ------------------------------------------------------------------------------
for site in sites:
    for pv_system in pv_systems:
        tmy_df = pd.read_csv(x[site][pv_system]['TMYFileName'])

        # Run calculations at dcac ratio of 1.0
        dcac_ratio = 1.0
        inv_scale = dc_capacity / inverter["Pac0"] / dcac_ratio
        pac_clipped, pac_unclipped = sandia_inverter(
            v_dc=tmy_df["Vdc"],
            p_dc=tmy_df["Pdc"],
            inv_param=inverter,
            inv_scale=inv_scale)

        x[site][pv_system]['NominalPac'] = pac_clipped.sum() * ac_derate

# ------------------------------------------------------------------------------
# Hourly Calculations of Pac from Minute Data
# ------------------------------------------------------------------------------
for site in sites:
    for pv_system in pv_systems:
        # Read minute data from CSV file
        minute_df = pd.read_csv(x[site][pv_system]['MinuteFileName'], index_col='Time')
        minute_df.index = pd.to_datetime(minute_df.index)
        num_hours = len(minute_df.index) // 60

        # ------------------------------------------------------------------------------
        # Find Valid Hours of Data
        # ------------------------------------------------------------------------------
        # Exclude hours with less than 50% data coverage for any variable
        nan_df = minute_df.isna()
        nan_count_df = nan_df.resample('H').sum()
        valid_df = nan_count_df < 30

        # Exclude night and sunrise/sunset hours
        max_zenith = minute_df['Zenith'].resample('H').max()
        valid_df['DayTime'] = max_zenith < 90

        x[site][pv_system]['ValidHour'] = valid_df.all(axis=1)

        # ------------------------------------------------------------------------------
        # Calculate Gamma DNI
        # ------------------------------------------------------------------------------
        # This is the same for all pv_systems at site, but duplicate calculation
        # and data so that comparable data are stored in the same node.
        dni_hour = minute_df['DNI'].resample('H').mean()
        dni_refsky_hour = minute_df['DNI_RefSky'].resample('H').mean()
        x[site][pv_system]['Gamma'] = dni_hour / dni_refsky_hour

        # ------------------------------------------------------------------------------
        # Calculate Clipping Potential
        # ------------------------------------------------------------------------------
        cp = np.empty((len(dcac_ratios), num_hours))
        cp[:] = np.nan
        pdc_refsky_hour = minute_df['Pdc_RefSky'].resample('H').mean()
        for i, dcac_ratio in enumerate(dcac_ratios):
            pac0 = dc_capacity / dcac_ratio
            cp[i, :] = (pdc_refsky_hour - pac0) / pac0

        x[site][pv_system]['CP'] = cp

        # ------------------------------------------------------------------------------
        # Calculate Pac Clip then Average
        # ------------------------------------------------------------------------------
        pac_cta = np.empty((len(dcac_ratios), num_hours))
        pac_cta[:] = np.nan
        for i, dcac_ratio in enumerate(dcac_ratios):
            inv_scale = dc_capacity / inverter["Pac0"] / dcac_ratio
            minute_df['Pac'] = sandia_inverter(
                v_dc=minute_df["Vdc"],
                p_dc=minute_df["Pdc"],
                inv_param=inverter,
                inv_scale=inv_scale)[0]
            pac_cta[i, :] = minute_df['Pac'].resample('H').mean()

        x[site][pv_system]['PacCtA'] = pac_cta

        # ------------------------------------------------------------------------------
        # Calculate Pac Average then Clip
        # ------------------------------------------------------------------------------
        pac_atc = np.empty((len(dcac_ratios), num_hours))
        pac_atc[:] = np.nan
        vdc_hour = minute_df['Vdc'].resample('H').mean()
        pdc_hour = minute_df['Pdc'].resample('H').mean()

        for i, dcac_ratio in enumerate(dcac_ratios):
            inv_scale = dc_capacity / inverter["Pac0"] / dcac_ratio
            pac_atc[i, :] = sandia_inverter(
                v_dc=vdc_hour,
                p_dc=pdc_hour,
                inv_param=inverter,
                inv_scale=inv_scale)[0]

        x[site][pv_system]['PacAtC'] = pac_atc

        # ------------------------------------------------------------------------------
        # Calculate Scaled AtC Bias
        # ------------------------------------------------------------------------------
        x[site][pv_system]['Delta'] = (pac_atc - pac_cta) / x[site][pv_system]['NominalPac']

# ------------------------------------------------------------------------------
# Print AtC Bias by Site - PV System - Year
# ------------------------------------------------------------------------------
for site in sites:
    for pv_system in pv_systems:
        valid_hour = x[site][pv_system]['ValidHour']
        delta = x[site][pv_system]['Delta'][:, valid_hour.to_numpy()]
        num_years = len(valid_hour) / 8760  # misses leap days

        print(f"{site}, {pv_system}")
        for i, dcac_ratio in enumerate(dcac_ratios):
            print(f"    dcac ratio = {dcac_ratio:4.2f}, sum delta = {np.sum(delta[i, :]) / num_years:6.4f}")

# ------------------------------------------------------------------------------
# Collect Hourly Data from All Sites and PV Systems
# ------------------------------------------------------------------------------
gamma = np.empty(0)
cp = np.empty(0)
delta = np.empty(0)
for site in sites:
    for pv_system in pv_systems:
        valid_hour = x[site][pv_system]['ValidHour'].to_numpy()

        y = x[site][pv_system]['Gamma'].loc[valid_hour].to_list()
        gamma = np.append(gamma, np.tile(y, len(dcac_ratios)))
        y = x[site][pv_system]['CP'][:, valid_hour]
        cp = np.append(cp, np.reshape(y, -1))
        y = x[site][pv_system]['Delta'][:, valid_hour]
        delta = np.append(delta, np.reshape(y, -1))

# ------------------------------------------------------------------------------
# Create Matrix
# ------------------------------------------------------------------------------
num_gamma_bins = len(gamma_bins) - 1
num_cp_bins = len(cp_bins) - 1
delta_matrix = np.zeros([num_gamma_bins, num_cp_bins])

gamma_bin_index = np.digitize(gamma, gamma_bins)
cp_bin_index = np.digitize(cp, cp_bins)

for gamma_bin in range(num_gamma_bins):
    for cp_bin in range(num_cp_bins):
        index = np.all(np.array([gamma_bin_index == gamma_bin,
                                 cp_bin_index == cp_bin]), axis=0)
        delta_matrix[gamma_bin, cp_bin] = delta[index].mean()

        if gamma_bin == 13 and cp_bin == 8:
            index2 = copy.copy(index)

# ------------------------------------------------------------------------------
# Export Matrix
# ------------------------------------------------------------------------------
matrix_df = pd.DataFrame(delta_matrix)
matrix_df = matrix_df.rename(index=lambda k: f"{gamma_bins[k]} - {gamma_bins[k+1]}")
matrix_df = matrix_df.rename(columns=lambda k: f"{cp_bins[k]} - {cp_bins[k+1]}")
matrix_df.to_csv(matrix_fn)


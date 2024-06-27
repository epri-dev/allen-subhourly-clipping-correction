#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to apply matrix correction for Average then Clip (AtC) bias
"""

import os
import re
import numpy as np
import pandas as pd
import pvlib

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
Use this script to estimate hourly corrections for the Average then Clip (AtC) 
bias using the method described in [1]. 

The user provides three CSV files containing hourly PV modeling inputs and 
outputs. The user enters the file names in the script below.

The CSV files are typically created using PV modeling tools, e.g. SAM or PVSyst;
they must contain: 
    - Hourly modeling using plant specifications and realistic meteorological data
        o time stamp
        o DNI (W/m2) (or GHI and DHI for PVSyst)
        o Pac
    - Hourly modeling using dc:ac ratio = 1.0 and realistic meteorological data
        o time stamp
        o Pac
    - Hourly modeling using plant specifications and reference sky radiation 
        o time stamp
        o DNI (W/m2) (or GHI and DHI for PVSyst)
        o Pdc

The CSV files will be joined based on the row number, so these must be 
consistent across the files. Units of power must also be consistent across 
the CSV files. The hourly modeling with dc:ac ratio = 1.0 is generated by 
increasing the number of modules.
 
The user also provides a file with an AtC Bias correction matrix [1] which 
contains average scaled AtC Bias (delta) for combinations of radiation index
(gamma) and clipping potential.

The script outputs an estimates of annual AtC Bias and a CSV file with the
time series average scaled AtC Bias. The rows of the output CSV file match
those of the input files.

Note that the hourly estimates of average scaled AtC Bias are expected to 
match the annual (and perhaps month-hour) average AtC Bias, but are not
expected to match the hourly values. See [1] for a detailed explanation.

References
----------
.. [1] J. O. Allen, W. B. Hobbs, and M. Bolen, "Classification Method to
   Predict the Effect of Short-Term Inverter Saturation on PV Performance
   Modeling", PV Performance Modeling Workshop Salt Lake City, 2022.

"""

# ------------------------------------------------------------------------------
# User Inputs
# ------------------------------------------------------------------------------
project_dir = r"C:\a2\epri_atcbias3\data\rabin"

# dc:ac ratio of plant specifications; see note below.
dc_ac_ratio = 1.4

# Format of PV modeling files
# For PVSyst modify the col_name lists in read_pvsyst_csv()
pv_model = "PVSyst"

# Hourly modeling using plant specifications and realistic meteorological data
pac_csv = os.path.join(project_dir, "PVModelPVSyst1.4.csv")
# Hourly modeling using dc:ac ratio = 1.0 and realistic meteorological data
pac_dcac1_csv = os.path.join(project_dir, "PVModel_PVSyst_1.csv")
# Hourly modeling using plant specifications and reference sky radiation
pac_refsky_csv = os.path.join(project_dir, "PVModel_PVSyst_refsky.csv")

# AtC Bias Correction Matrix
matrix_fn = os.path.join(project_dir, "method2.csv")

# Output CSV file name
output_csv = os.path.join(project_dir, "atc_bias_output.csv")


# ------------------------------------------------------------------------------
# Utility Functions
# ------------------------------------------------------------------------------
def read_limits(limit_lines):
    # Convert text list of bin limits to list
    # So for
    #   limit_lines = ['0 - 1', '1 - 5', '5 - 20']
    #   limit_list = [0, 1, 5, 20]

    limit_list = []
    n = []
    for limit_line in limit_lines:
        n = re.findall(r"[+-]?\d*\.*\d+", limit_line)
        limit_list.append(float(n[0]))
    limit_list.append(float(n[1]))

    return limit_list


# ------------------------------------------------------------------------------
# Read Input Files from SAM
# ------------------------------------------------------------------------------
def read_sam_one_csv(csv_fn, col_dict):
    csv_df = pd.read_csv(csv_fn)
    csv_df = csv_df.rename(columns=col_dict)
    csv_df = csv_df.drop(columns=csv_df.columns.difference(col_dict.values()))
    return csv_df


def read_sam_csv(pac_csv, pac_dcac1_csv, pac_refsky_csv):
    # Pac and DNI for plant specifications
    col_dict = {
        "Irradiance DNI from weather file | (W/m2)": "dni",
        "Inverter AC output power | (kW)": "pac",
    }
    df = read_sam_one_csv(pac_csv, col_dict)

    # Pac for dc:ac ratio = 1.0
    col_dict = {
        "Inverter AC output power | (kW)": "pac_dcac1",
    }
    df2 = read_sam_one_csv(pac_dcac1_csv, col_dict)
    df = df.join(df2)

    # Pdc and DNI for reference sky
    col_dict = {
        "Irradiance DNI from weather file | (W/m2)": "dni_refsky",
        "Inverter DC input power | (kW)": "pdc_refsky",
    }
    df2 = read_sam_one_csv(pac_refsky_csv, col_dict)
    df = df.join(df2)

    return df


# ------------------------------------------------------------------------------
# Read Input Files from PVSyst
# ------------------------------------------------------------------------------
def read_pvsyst_one_csv(csv_fn):
    # Read column names on row 11
    skip_rows = 10
    encoding = "latin-1"
    csv_df1 = pd.read_csv(csv_fn, skiprows=skip_rows, encoding=encoding)

    # Read data starting on row 13
    skip_rows = 13
    encoding = "latin-1"
    csv_df2 = pd.read_csv(csv_fn, skiprows=skip_rows, encoding=encoding,
                          header=None, names=csv_df1.columns)
    csv_df2.index = pd.to_datetime(csv_df2['date'], dayfirst=True, format='mixed')

    return csv_df2


def read_pvsyst_csv(pac_csv, pac_dcac1_csv, pac_refsky_csv):
    # Expected columns of data (not all are needed for each of the 3 csv files)
    col_dict = {'GlobHor': 'ghi', 'DiffHor': 'dhi', 'HSol': 'solar_elevation',
                'EOutArray': 'pdc', 'EOutInv': 'pac'}

    # Pac and DNI for plant specifications
    # Required columns - 'date', 'zenith', 'ghi', 'dhi', 'pac'
    df = read_pvsyst_one_csv(pac_csv)
    df = df.rename(columns=col_dict)
    df = df.sort_index()

    # Calculate DNI as (GHI - DHI) / cos(zenith)
    df["dni"] = (df["ghi"] - df["dhi"]) / pvlib.tools.cosd(90 - df["solar_elevation"])

    # Pac for dc:ac ratio = 1.0
    # Required columns - 'date', 'pac'
    df2 = read_pvsyst_one_csv(pac_dcac1_csv)
    df2 = df2.rename(columns=col_dict)
    df = df.join(df2, rsuffix='_dcac1')

    # Pdc and DNI for reference sky
    # Required columns - 'date', 'ghi', 'dhi', 'pdc'
    df2 = read_pvsyst_one_csv(pac_refsky_csv)
    df2 = df2.rename(columns=col_dict)
    df = df.join(df2, rsuffix='_refsky')

    # Calculate DNI as (GHI - DHI) / cos(zenith)
    df["dni_refsky"] = (df["ghi_refsky"] - df["dhi_refsky"]) / pvlib.tools.cosd(90 - df["solar_elevation"])

    return df


# ------------------------------------------------------------------------------
# Read Input Data from CSV Files
# ------------------------------------------------------------------------------
match pv_model:
    case "SAM":
        pac_df = read_sam_csv(pac_csv, pac_dcac1_csv, pac_refsky_csv)

    case "PVSyst":
        pac_df = read_pvsyst_csv(pac_csv, pac_dcac1_csv, pac_refsky_csv)

    case _:
        pac_df = pd.DataFrame(["ghi", "dhi", "pac", "pac_dcac1", "ghi_refsky",
                               "dhi_refsky", "pac_refsky"])

# ------------------------------------------------------------------------------
# Read Correction Matrix
# ------------------------------------------------------------------------------
matrix_df = pd.read_csv(matrix_fn)
cp_limit_headers = matrix_df.iloc[:, 1:].columns
cp_limits = read_limits(cp_limit_headers)

gamma_limit_lines = matrix_df.iloc[:, 0]
gamma_limits = read_limits(gamma_limit_lines)

corr_matrix = matrix_df.iloc[:, 1:]

# ------------------------------------------------------------------------------
# Calculate gamma and CP
# ------------------------------------------------------------------------------
pac_df["gamma"] = pac_df["dni"] / pac_df["dni_refsky"]

# Get inverter capacity estimated from the maximum observed hourly output.
# This will only work for installations with dc:ac ratio > 1.0.
# The user can also enter the actual value here.
pac0 = pac_df["pac"].max()

pac_df["cp"] = (pac_df["pdc_refsky"] - pac0) / pac0

# ------------------------------------------------------------------------------
# Create AtC Bias Correction Time Series
# ------------------------------------------------------------------------------
# Find the gamma and CP bins for each hour
gamma_bin = pd.cut(pac_df.gamma, gamma_limits, labels=False).fillna(-1).astype("int")
cp_bin = pd.cut(pac_df.cp, cp_limits, labels=False).fillna(-1).astype("int")

# Fill in the estimated correction for each hour
correction = np.zeros(len(pac_df))
# Use for loop to fill correction values from matrix with these complications:
#   - bin indices may be invalid (-1)
#   - no corrections for first and last hour of PV output (hours of sunrise and sunset)
for i in range(len(pac_df)):
    # Fill corrections from matrix
    if gamma_bin[i] >= 0 and cp_bin[i] >= 0:
        correction[i] = corr_matrix.iloc[gamma_bin[i], cp_bin[i]]

    # Exclude hours of sunrise and sunset
    if 0 < i < (len(pac_df) - 1):
        if np.isnan(pac_df["pac"][i - 1]) or np.isnan(pac_df["pac"][i + 1]):
            correction[i] = 0

# Scale correction by dc:ac ratio because correction matrix was developed by
# scaling inverter _down_ whereas dc output using PV models is calculated by
# scaling modules _up_. Remove '* dc_ac_ratio' in next line if dc output
# generated by scaling inverter down.
pac_df["bias_correction"] = correction * dc_ac_ratio
pac_df["pac_bias_Correction"] = correction * pac_df["pac_dcac1"].sum()

# ------------------------------------------------------------------------------
# Create AtC Bias Correction Time Series
# ------------------------------------------------------------------------------
n_years = len(pac_df) / 8760
print(f'Annual AtC Bias = {pac_df["bias_correction"].sum() / n_years * 100}%')
pac_df.to_csv(output_csv)

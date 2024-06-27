#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to calculate reference sky radiation for a TMY site
"""

import os
import numpy as np
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

# ------------------------------------------------------------------------------
# User Input of TMY File
# ------------------------------------------------------------------------------
tmy_dir = r"C:\a2\tmy3\alltmy3a"
# Full path to TMY input and output files
tmy_id = "725128"  # State College, PA
tmy_fn = os.path.join(tmy_dir, f"{tmy_id}TYA.CSV")


# ------------------------------------------------------------------------------
# Read TMY File
# ------------------------------------------------------------------------------
# pvlib returns TZ aware tmy_df.index
# TMY has end-of-period time stamp
tmy_df, tmy_meta = pvlib.iotools.read_tmy3(tmy_fn, map_variables=False)

# Parse meta data
latitude = tmy_meta["latitude"]
longitude = tmy_meta["longitude"]
altitude = tmy_meta["altitude"]

# ------------------------------------------------------------------------------
# Calculate Reference Sky Radiation
# ------------------------------------------------------------------------------
# Calculate solar position
# Make calculations at middle of hour, i.e. 30 min prior to time index
delta_t = np.timedelta64(1, "h") / 2.0
t = tmy_df.index - delta_t
solar_df = pvlib.solarposition.get_solarposition(
    t,
    latitude=latitude,
    longitude=longitude,
    altitude=altitude,
    method="nrel_numpy",
)

# Calculate air mass and extra terrestrial DNI needed for reference sky calcs
solar_df["airmass_relative"] = pvlib.atmosphere.get_relative_airmass(
    zenith=solar_df["zenith"]
)
solar_df["airmass_absolute"] = pvlib.atmosphere.get_absolute_airmass(
    airmass_relative=solar_df["airmass_relative"], pressure=tmy_df["Pressure (mbar)"]
)
solar_df["dni_extra"] = pvlib.irradiance.get_extra_radiation(solar_df.index)

# Calculate reference sky radiation
linke_turbidity = 1.0
clear_sky_df = pvlib.clearsky.ineichen(
    apparent_zenith=solar_df["apparent_zenith"],
    airmass_absolute=solar_df["airmass_absolute"],
    linke_turbidity=linke_turbidity,
    altitude=altitude,
    dni_extra=solar_df["dni_extra"],
)

# ------------------------------------------------------------------------------
# Write Reference Sky Radiation in TMY Format
# ------------------------------------------------------------------------------
# Replace radiation in tmy_df with reference sky values
tmy_df["DHI (W/m^2)"] = clear_sky_df["dhi"]
tmy_df["DNI (W/m^2)"] = clear_sky_df["dni"]
tmy_df["GHI (W/m^2)"] = clear_sky_df["ghi"]

# Write new file in TMY format with metadata header
refsky_fn = f"{tmy_id}_RefSky.CSV"
f = open(refsky_fn, "w")
f.write(
    f'{tmy_meta["USAF"]}, '
    f'{tmy_meta["Name"]}, '
    f'{tmy_meta["State"]}, '
    f'{tmy_meta["TZ"]}, '
    f'{tmy_meta["latitude"]}, '
    f'{tmy_meta["longitude"]}, '
    f'{tmy_meta["altitude"]} '
    "\n"
)
tmy_df.to_csv(f, index=False)
f.close()

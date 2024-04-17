#!/usr/bin/env python3

import astropy.units as u
from astropy.coordinates import get_body, SkyCoord
from astropy.time import Time
from sunpy.coordinates import Helioprojective

import hi_c_targeter

"""
Targets
A. AR 13644 (Previously East limb): N12E78 (-914”,  218”)
B. AR 13639: N30E35                        (-476",  541")
C-1. AR 13641: N11W15                      ( 243",  268")
C-2. AR 13635: N22W28                      ( 417",  431")
D. AR 13643: S10E49                        (-712", -107")
E-1. AR 13637: S12E25                      (-396", -118")
E-2. AR 13638: S17E30                      (-458", -203")
E-3. (region that flared today):           (-537”, -69”)
F. AR 13634: N25W56                        ( 719",  448")
G. AR 13636: S18E06                        (- 95", -209")
The region positions are valid on 16-Apr-2024 23:30 UT
"""

# Define target coordinates.
obstime = Time("2024-04-16 23:30", scale="utc")
frame = Helioprojective(observer=get_body("Earth", obstime))
coord_units = [u.arcsec, u.arcsec]
planning_targets = {"Alpha":     SkyCoord(Tx= -914, Ty=  218, unit=coord_units, frame=frame),
                    "Bravo":     SkyCoord(Tx= -476, Ty=  541, unit=coord_units, frame=frame),
                    "Charlie-1": SkyCoord(Tx=  243, Ty=  268, unit=coord_units, frame=frame),
                    "Charlie-2": SkyCoord(Tx=  417, Ty=  431, unit=coord_units, frame=frame),
                    "Delta":     SkyCoord(Tx= -712, Ty= -107, unit=coord_units, frame=frame),
                    "Echo-1":    SkyCoord(Tx= -396, Ty= -118, unit=coord_units, frame=frame),
                    "Echo-2":    SkyCoord(Tx= -458, Ty= -203, unit=coord_units, frame=frame),
                    "Echo-3":    SkyCoord(Tx= -537, Ty= -69, unit=coord_units, frame=frame),
                    "Foxtrot":   SkyCoord(Tx=  719, Ty=  448, unit=coord_units, frame=frame),
                    "Golf":      SkyCoord(Tx= -95,  Ty= -209, unit=coord_units, frame=frame)}

# Define payload orientations.
payload_orientations = ["east"] * len(planning_targets)

# Define median time of launch window.
median_launch_time = Time("2024-04-17 21:53", scale="utc")

# Predict target positions during window and 5" shifts as a function of time.
x_shifts, predicted_positions = hi_c_targeter.calculate_sparcs_x_offsets_for_targets(
    planning_targets, median_launch_time, payload_orientations=payload_orientations)

print(x_shifts)

x_shifts.write(f"HiC-Science-Target_Tables_{median_launch_time.value[:4]}_{median_launch_time.value[5:7]}_{median_launch_time.value[8:10]}.csv")

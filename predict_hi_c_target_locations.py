#!/usr/bin/env python3

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from sunpy.coordinates import HeliographicStonyhurst

import hi_c_targeter

"""
Targets
A. NE limb: N20E84 (-895",335")
B. AR 13632: N26E05 (-75",509")
C. AR 13629: N05W67 (879",122")
D. AR 13633: S08E47 (-695",-64")
E. AR 13628: N08E02 (-33",233")
F. AR 13630: S11W65 (853",-140")
G. AR 13634: N26E57 (-720", 470")
H. AR 13631: N11W49 (711",247")
"""

# Define target coordinates.
obstime = Time("2024-04-08 19:30", scale="utc")
frame = HeliographicStonyhurst
coord_units = [u.deg, u.deg]
planning_targets = {"A": SkyCoord(lat=20, lon=-84, unit=coord_units, obstime=obstime, frame=frame),
                    "B": SkyCoord(lat=26, lon=-5, unit=coord_units, obstime=obstime, frame=frame),
                    "C": SkyCoord(lat=5, lon=67, unit=coord_units, obstime=obstime, frame=frame),
                    "D": SkyCoord(lat=-8, lon=-47, unit=coord_units, obstime=obstime, frame=frame),
                    "E": SkyCoord(lat=8, lon=-2, unit=coord_units, obstime=obstime, frame=frame),
                    "F": SkyCoord(lat=-11, lon=65, unit=coord_units, obstime=obstime, frame=frame),
                    "G": SkyCoord(lat=26, lon=-57, unit=coord_units, obstime=obstime, frame=frame),
                    "H": SkyCoord(lat=11, lon=49, unit=coord_units, obstime=obstime, frame=frame)}
planning_targets = list(planning_targets.values())

# Define payload orientations
payload_orientations = ["north"] * len(planning_targets)

# Define median time of launch window.
median_launch_time = Time("2024-04-09 21:53", scale="utc")

# Predict target positions during window and 5" shifts as a function of time.
x_shifts, predicted_positions = hi_c_targeter.calculate_sparcs_x_offsets_for_targets(
    planning_targets, median_launch_time, payload_orientations=payload_orientations)

print(x_shifts)

x_shifts.write("hi_c_sparcs_targets.csv")

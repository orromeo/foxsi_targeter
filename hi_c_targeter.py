import warnings

import astropy.coordinates
import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.table import Table
from sunpy.coordinates import Helioprojective, RotatedSunFrame
from sunpy.util import SunpyUserWarning

TIME_EDGES = [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2] * u.hour
TIME_CENTERS = (TIME_EDGES[1:] - TIME_EDGES[:-1]) / 2 + TIME_EDGES[:-1]
UTC_TO_LOCAL = -8 * u.hour
SPARCS_X_INCREMENT = 5 * u.arcsec


def calculate_sparcs_x_offsets_for_targets(planning_targets, median_launch_time, payload_orientations=None, diff_rotate=None, **kwargs):
    # Sanitize inputs
    n_targets = len(planning_targets)
    if not payload_orientations:
        payload_orientations = ["north"] * n_targets
    if not diff_rotate:
        diff_rotate = [True] * n_targets
    
    # Build times for which target positions are needed
    date_edges = TIME_EDGES + median_launch_time
    date_centers = TIME_CENTERS + median_launch_time
    date_edges_local = date_edges + UTC_TO_LOCAL
    date_centers_local = date_centers + UTC_TO_LOCAL

    # Build table of results.
    n_intervals = len(TIME_CENTERS)
    idx_median = int(n_intervals / 2)
    sparcs_buttons = {"north": ["Bottom Right"] * idx_median + ["N/A"] + ["Top Left"] * idx_median,
                      "west": ["Bottom Left"] * idx_median + ["N/A"] + ["Top Right"] * idx_median,
                      "south": ["Top Left"] * idx_median + ["N/A"] + ["Bottom Right"] * idx_median,
                      "east": ["Top Right"] * idx_median + ["N/A"] + ["Bottom Left"] * idx_median}
    Delta_str = "\N{GREEK CAPITAL LETTER DELTA}"
    initial_rows = [tuple(["Interval End [UTC]"] + [d[11:16] for d in date_edges[1:].value] + [""]),
                    tuple(["Interval End [AKDT]"] + [d[11:16] for d in date_edges_local[1:].value] + [""]),
                    tuple([""] * (n_intervals + 2)),
                    tuple(["Target Name/SPARCS button"] + sparcs_buttons[payload_orientations[0]] + [""])]
    idxs = list(range(len(date_centers)))
    colnames = tuple([date_edges[0].utc.value[:10]]
                     + [f"Interval {i} {Delta_str}W" for i in idxs[:idx_median]]
                     + ["Median (N, W)"]
                     + [f"Interval {i} {Delta_str}W" for i in idxs[idx_median+1:]]
                     + ["Payload Orientation"])

    results = Table(rows=initial_rows, names=colnames)
    target_pointings = {}

    # Get position of each target, calculate shift from median position, and enter results into table.
    hpc_unit = u.arcsec
    warnings.filterwarnings("error")
    factor = 1
    for ((target_name, planning_target), payload_orientation) in zip(planning_targets.items(), payload_orientations):
        coord_off_disk = True
        while coord_off_disk:
            try:
                sparcs_x_offsets, flight_pointings = calculate_sparcs_x_offsets(
                    planning_target, date_centers, payload_orientation=payload_orientation, **kwargs)
                coord_off_disk = False
            except SunpyUserWarning as warning:
                if "because off-disk" in warning.args[0]:
                    factor -= 0.001
                    cu = u.arcsec
                    new_Tx, new_Ty = contract_to_disk(planning_target.Tx.to_value(cu),
                                                      planning_target.Ty.to_value(cu),
                                                      factor=factor)
                    planning_target = SkyCoord(Tx=new_Tx, Ty=new_Ty, unit=[cu]*2, frame=planning_target.frame)
                else:
                    coord_off_disk = False
        dx = (np.round(sparcs_x_offsets.to_value(SPARCS_X_INCREMENT.unit) / SPARCS_X_INCREMENT.value, decimals=0) * SPARCS_X_INCREMENT).astype(int)
        if payload_orientation in {"west", "south"}:
            dx *= -1
        dx_str = [str(x.value) + '"' for x in dx]
        new_row = ([target_name]
                   + dx_str[:idx_median]
                   + [f'({round(flight_pointings.Ty[idx_median].to_value(hpc_unit))}", {round(flight_pointings.Tx[idx_median].to_value(hpc_unit))}")']
                   + dx_str[idx_median+1:]
                   + [payload_orientation])
        results.add_row(new_row)
        target_pointings[target_name] = flight_pointings
    warnings.resetwarnings()

    return results, target_pointings


def calculate_sparcs_x_offsets(planning_target, times, payload_orientation=None, **kwargs):
    # Get predicted positions of target at key times throughout launch window.
    idx_median = int(len(times) / 2)
    target_positions = predict_target_position(planning_target, times, **kwargs)
    sparcs_x_offsets = target_positions.Tx - target_positions[idx_median].Tx
    if payload_orientation is not None and payload_orientation.lower() == "south":
        sparcs_x_offsets *= -1
    return sparcs_x_offsets, target_positions 


def predict_target_position(planning_target, times, **kwargs):
    # Get coordinate of FOXSI at apogee.
    poker_lon, poker_lat, poker_altitude = -147.47 * u.deg, 65.12 * u.deg, 197 * u.m
    foxsi_apogee = 400 * u.km
    foxsi_loc = EarthLocation(lon=poker_lon, lat=poker_lat, height=poker_altitude + foxsi_apogee)
    foxsi_skycoord = SkyCoord(foxsi_loc.get_gcrs(times))

    # Differentially rotate target to launch time.
    target_diffrot = SkyCoord(RotatedSunFrame(base=planning_target.frame, rotated_time=foxsi_skycoord.obstime,
                                              rotation_model=kwargs.pop("rotation_model", None), **kwargs))
    # Transform to FOXSI observer view.
    target_prediction = target_diffrot.transform_to(Helioprojective(observer=foxsi_skycoord))

    return target_prediction


def contract_to_disk(Tx, Ty, factor=1):
        rsun_arcsec = 959.63
        Tx_unit = Tx / np.sqrt(Tx**2 + Ty**2)
        Ty_unit = Ty / np.sqrt(Tx**2 + Ty**2)
        Tx_onlimb = Tx_unit * rsun_arcsec * factor # Factor of solar radii to consider
        Ty_onlimb = Ty_unit * rsun_arcsec * factor
        return Tx_onlimb, Ty_onlimb

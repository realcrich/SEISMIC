# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 22:12:51 2025

@author: datco
"""

import os
import pandas as pd

def parse_my_vs_folder(folder, num_geophones, fill_value=-999.031250):
    """
    Parse xxx_picks.vs files into a dictionary of structure:

    {
        shot_key: [
            {"geophone": int, "x": float, "z": float,
             "t_obs": float, "t_calc": float},
            ...
        ]  # always length = num_geophones
    }
    """

    my_shot_dict = {}

    for fname in os.listdir(folder):
        if fname.endswith("_picks.vs"):
            shot_key = fname.split("_")[0]   # e.g., "30000-30007"
            filepath = os.path.join(folder, fname)

            with open(filepath, "r") as f:
                lines = f.readlines()

            # Actual data lines (skip header/footer)
            data_lines = lines[2:-3]

            picks = []

            # ---- Parse real picks first ----
            for geoph_idx, line in enumerate(data_lines, start=1):
                if geoph_idx > num_geophones:
                    break  # avoid overflow

                parts = line.split()
                if len(parts) < 2:
                    continue

                x     = float(parts[0])
                t_obs = float(parts[1])/1000    # ms -> s
                z     = 0.0

                picks.append({
                    "#": geoph_idx,
                    "x": x,
                    "z": z,
                    "t_obs": t_obs,
                    "t_calc": fill_value
                })

            # ---- Pad missing receivers ----
            for geoph_idx in range(len(picks)+1, num_geophones+1):
                picks.append({
                    "#": geoph_idx,
                    "x": fill_value,      # or set to expected X positions if known
                    "z": 0.0,
                    "t_obs": fill_value,
                    "t_calc": fill_value
                })

            my_shot_dict[shot_key] = picks

    return my_shot_dict

def write_shots_to_txt(my_shot_dict, output_path, elev_file, data_summary_file):
    """
    Writes the full dictionary of shots to a formatted text file.

    my_shot_dict = {
        "100": [ {geophone, x, z, t_obs, t_calc}, ... ],
        "101": [...],
        ...
    }

    shot_locations (optional) = {
        "100": (x, z),
        ...
    }
    If not provided, defaults to (0.0, 0.0)
    """

    shot_keys = list(my_shot_dict.keys())
    total_shots = len(shot_keys)
    
    dx = pd.read_csv(dsum_f)['Distance (m)'][:-3]

    with open(output_path, "w") as f:

        # ---- (1) Elevation block at top ----
        if elev_file is not None:
            f.write("# ELEVATION:\n")
            with open(elev_file, "r") as ef:
                for line in ef:
                    f.write(line.rstrip() + "\n")
            f.write("\n")  # blank line after elevation block

        for idx,shot_key in enumerate(shot_keys):
            picks = my_shot_dict[shot_key]
            num_receivers = len(picks)

            sx,sz = dx[idx],0.0
            
            # ---- Shot header ----
        
            f.write(f"TOTAL_NUMBER_OF_SHOTS: {total_shots}\n")
            f.write(f"SHOT_NUMBER: {shot_key}\n")
            f.write(f"SHOT_LOCATION[X,Z]: {sx:.3f} {sz:.3f}\n")
            f.write(f"NUMBER_OF_RECEIVERS: {num_receivers}\n")
            f.write("     #                 X                 Z             T_OBS             T_CAL\n")

            # ---- Receiver lines ----
            for p in picks:
                f.write(
                    f"{p['#']:6d}"
                    f"{p['x']:20.3f}"
                    f"{p['z']:20.3f}"
                    f"{p['t_obs']:16.6f}"
                    f"{p['t_calc']:16.6f}\n"
                )

            # blank line between shots
            f.write("\n")

# ---- Read in all directories and filepaths to create/write out all shots ---- #
    ''' 
    home_dir -> directory where all seismic data is housed
    line -> defines which geophone line we are processing
    ddir -> data directory for line within home_dir
    pdir -> directory where all pick files are stored
    elev_f -> elevation data file for line X 
    dsum_f -> data summary file as .csv to read in shot distance along transect
    out_f -> output filename 
    '''
    
home_dir = 'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/'
line = 'T2' 
ddir = home_dir+line+'/'
pdir = ddir+'picks/'
elev_f = ddir+'SS_P304_070922_L2-topo.txt'
dsum_f = ddir+'BCZCN_Seismic_RES_DataSheet_Summer2022_P304_Line2All.csv'
out_f = ddir+'/{}_all_shots_firstbreaks.txt'.format(line)

my_shot_dict = parse_my_vs_folder(pdir, num_geophones=240)
write_shots_to_txt(my_shot_dict, out_f, elev_f, dsum_f)
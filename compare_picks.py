# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 12:51:11 2025

@author: Collin
"""
import os
import re
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
'''
shot_dict = defaultdict(list)

with open("C:/Users/Collin/Downloads/BW_080714_2/BW_080714_2/Seismic/comp_picks.txt", "r") as f:
    lines = f.readlines()

current_shot = None
for line in lines:
    line = line.strip()

    # detect start of shot block
    if line.startswith("SHOT_NUMBER:"):
        current_shot = int(line.split(":")[1].strip())

    # detect receiver rows (they start with an index number)
    elif re.match(r"^\d+", line):
        parts = line.split()
        # parts = [#, X, Z, T_OBS, T_CAL]
        x = float(parts[1])
        t_obs = float(parts[3])
        if current_shot is not None:
            shot_dict[current_shot].append((x, t_obs))

# convert defaultdict to dict
shot_dict = dict(shot_dict)

# quick test: show first 5 receivers of first shot
first_shot = list(shot_dict.keys())[0]
print(f"Shot {first_shot} → {shot_dict[first_shot][:5]}")
'''
def parse_test_picks(filename):
    """
    Parse test picks into a dictionary like:
    {
        shot_number: [
            {"x": float, "t_obs": float},
            ...
        ]
    }
    """
    shot_dict = {}
    with open(filename, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("SHOT_NUMBER:"):
            shot_number = int(line.split(":")[1].strip())
            i += 2  # skip SHOT_LOCATION line
            n_receivers = int(lines[i].split(":")[1].strip())
            i += 1  # move to receiver header line
            picks = []
            for _ in range(n_receivers):
                i += 1
                parts = lines[i].split()
                x = float(parts[1])
                t_obs = float(parts[3])
                picks.append({"x": x, "t_obs": t_obs})
            shot_dict[shot_number] = picks
        i += 1

    return shot_dict

def parse_my_vs_folder(folder):
    """
    Parse xxx_picks.vs files into a dictionary:
    {
        shot_number: [
            {"x": float, "t_obs": float},
            ...
        ]
    }
    """
    my_shot_dict = {}

    for fname in os.listdir(folder):
        if fname.endswith("_picks.vs"):
            shot_number = fname.split("_")[0]
            #shot_number = int(fname.split("_")[0])  # e.g. "600_picks.vs" → 600
            filepath = os.path.join(folder, fname)

            with open(filepath, "r") as f:
                lines = f.readlines()

            # Data starts from line 3 (skip headers)
            picks = []
            for line in lines[3:-3]:
                parts = line.split()
                if len(parts) < 2:
                    continue
                x = float(parts[0])
                t_obs = float(parts[1])
                picks.append({"x": x, "t_obs": t_obs})

            my_shot_dict[shot_number] = picks

    return my_shot_dict

def compare_shots(shot_dict, my_shot_dict, shot_number, ms=True):
    """
    Compare test vs. my picks for a given shot.
    shot_dict: reference picks
    my_shot_dict: your picks
    shot_number: int
    ms: if True, convert seconds → milliseconds
    """
    if shot_number not in shot_dict or shot_number not in my_shot_dict:
        print(f"Shot {shot_number} not found in one of the dictionaries.")
        return

    # Extract test picks
    test_x = [rec["x"] for rec in shot_dict[shot_number]]
    test_t = [rec["t_obs"] for rec in shot_dict[shot_number]]

    # Extract your picks
    my_x = [rec["x"] for rec in my_shot_dict[shot_number]]
    my_t = [rec["t_obs"] for rec in my_shot_dict[shot_number]]

    # Convert to ms if needed
    if ms:
        test_t = [t*1000 for t in test_t]
        my_t = [t for t in my_t]
    
    test_x = np.array(test_x)    
    test_t = np.array(test_t)
    
    my_x = np.array(my_x)
    my_t = np.array(my_t)

    # Plot
    plt.figure(figsize=(8,5))
    cond = test_t > 0
    my_cond = my_t > 0
    plt.plot(test_x[cond], test_t[cond], "ko-", label="Test Picks")   # black circles
    plt.plot(my_x[my_cond], my_t[my_cond], "rx-", label="My Picks")         # red x's
    plt.xlabel("Receiver X (m)")
    plt.ylabel("First Break Time (ms)" if ms else "First Break Time (s)")
    plt.title(f"Shot {shot_number} Comparison")
    plt.legend()
    plt.grid(True)
    plt.savefig('C:/Users/Collin/SEISMIC/compare_picks/shot_{}_comparison.png'.format(shot_number))
    #plt.show()
    plt.close()
    
    #return test_x, test_t, my_x, my_t
    
#test_x, test_t, my_x, my_t = compare_shots(shot_dict, my_shot_dict, shot_number=600)

def read_out_shots(my_shot_dict, shot, outdir):
    
    '''
    if shot not in my_shot_dict:
        print(f"Shot {shot} not found in my_shot_dict.")
        return
    '''
    # Define receiver x positions (96 receivers, 0–237.5 m in 2.5 m steps)
    #x_positions = np.arange(0, 240, 2.5)
    # Define receiver x positions (X receivers, 0–237.5 m in 1 0 m steps)
    x_positions = np.arange(0, 240, 1.0)

    # Build lookup dict for fast access
    picks = {rec["x"]: rec["t_obs"] for rec in my_shot_dict[shot]}
    #picks = {rec["x"]: rec["t_obs"]/1000 for rec in my_shot_dict[shot]}  # convert ms → s
    fill_value = -999.031250

    # Construct output rows
    data = [["x", "T_OBS"]]
    for x in x_positions:
        t_obs = picks.get(x, fill_value)
        data.append([x, t_obs])

    # Write to file
    outfile = f"{outdir}{shot}_firstbreaks.txt"
    with open(outfile, "w") as f:
        for row in data:
            f.write(f"{row[0]}\t{row[1]}\n")

    print(f"Saved: {outfile}")
    return data

def write_all_shots(my_shot_dict, outfile):
    # Define receiver x positions (96 receivers, 0–237.5 m in 2.5 m steps)
    #x_positions = np.arange(0, 240, 2.5)
    # Define receiver x positions (96 receivers, 0–237.5 m in 2.5 m steps)
    x_positions = np.arange(0, 240, 1.0)
    fill_value = -999.031250

    with open(outfile, "w") as f:
        for shot in sorted(my_shot_dict.keys()):  # loop through shots in order
            f.write(f"SHOT_NUMBER = {shot}\n")
            f.write("x\tT_OBS\n")

            # Build lookup dict for fast access
            picks = {rec["x"]: rec["t_obs"]for rec in my_shot_dict[shot]}
            #picks = {rec["x"]: rec["t_obs"]/1000 for rec in my_shot_dict[shot]}  # ms → s

            # Write rows for all receivers
            for x in x_positions:
                t_obs = picks.get(x, fill_value)
                f.write(f"{x:.3f}\t{t_obs:.6f}\n")

            f.write("\n")  # blank line between shots

    print(f"Saved all shots to {outfile}")

#shot_dict = parse_test_picks("C:/Users/Collin/Downloads/BW_080714_2/BW_080714_2/Seismic/comp_picks.txt")
#folder = "C:/Users/Collin/Downloads/BW_080714_2/BW_080714_2/Seismic/Deployment1"
folder = 'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T2/picks'
my_shot_dict = parse_my_vs_folder(folder)

#for shot in range(600,631):
for shot in my_shot_dict.keys():
    #compare_shots(shot_dict, my_shot_dict, shot_number=shot)
    read_out_shots(my_shot_dict, shot, outdir=folder+'/all_picks/')
    
write_all_shots(my_shot_dict, outfile=folder+'/all_firstbreaks.ttx')
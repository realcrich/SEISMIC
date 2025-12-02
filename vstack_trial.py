# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 11:34:19 2025

@author: datco
"""

# -*- coding: utf-8 -*-
"""
Clean version of the stacking script.
Stacks files 8 at a time, computes mean, and writes out SEG-Y.
"""

import os
import numpy as np
import obspy as op

# ----- CONFIG -----
'''
ddirs = [
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T7/Raw1/',
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T7/Raw2/'
]

ddirs = [
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T2/Raw1/',
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T2/Raw2/',
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T2/Raw3/'
]
'''
ddirs = [
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T1/Raw1/',
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T1/Raw2/',
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T1/Raw3/',
    'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T1/Raw4/',
]

stack_size = 8
skip_vals = ['32999', '32256', '32569'], 
             #'38993', '38994', '38995', '38996', '38997',
             #'38998', '38999', '39000']#, 
             #'32930', '32931', '32932', '32933', '32934',
             #'32935', '32936']

outdir = 'C:/Users/datco/Downloads/southern_sierra/data/contents/P304/Seismic/T1/stacks/'
os.makedirs(outdir, exist_ok=True)

# ----- Derive trace/sampling structure from a sample file -----
f_samp = os.path.join(ddirs[0], os.listdir(ddirs[0])[1])
sample_st = op.read(f_samp)
n_tr, n_pts = sample_st[0].stats.npts * 0 + len(sample_st), sample_st[0].stats.npts


# =====================================================================
#                           MAIN LOOP
# =====================================================================

for ddir in ddirs:

    flist = [
        os.path.join(ddir, f) for f in os.listdir(ddir)
        if f.endswith(".dat") and f[-9:-4] not in skip_vals
    ]
    flist = sorted(flist)

    print(f"\nDirectory: {ddir}")
    print(f"Found {len(flist)} usable files.")

    i = 0
    nfiles = len(flist)

    while i < nfiles:

        # Proposed 8-file chunk
        cstack = flist[i:i + stack_size]

        # --- CASE 1: Special 7-file stack starting with 32930 ---
        if '32930' in cstack[0]:
            cstack = cstack[:-1]      # force 7 files
            step = 7                  # advance by 7
            print(f"Special 7-file stack triggered at index {i}")

        # --- CASE 2: Normal full 8 files ---
        elif len(cstack) == stack_size:
            step = stack_size         # advance by 8

        # --- CASE 3: Incomplete leftover chunk (skip it) ---
        else:
            print(f"Skipping last partial chunk ({len(cstack)} files).")
            break

        # Extract shot numbers from first/last file in cstack
        shot0 = cstack[0][-9:-4]
        shotf = cstack[-1][-9:-4]
        print(f"Stacking {len(cstack)} files: {shot0}–{shotf}")

        # Allocate stack based on actual number of files
        d_stack = np.empty((len(cstack), n_tr, n_pts))

        # Load all traces into d_stack
        for j, fname in enumerate(cstack):
            st = op.read(fname)
            d = np.array([tr.data for tr in st])

            if len(d) == n_tr:
                d_stack[j] = d
            else:
                d_stack[j] = np.zeros((n_tr, n_pts))
                d_stack[j][:len(d)] = d

        # Compute the average
        d_out = np.mean(d_stack, axis=0)

        # Use template to create output SEGY
        st_template = op.read(f_samp)
        for k, tr in enumerate(st_template):
            tr.data = d_out[k].astype(np.float32)

        outname = f"{shot0}-{shotf}_stack.segy"
        outpath = os.path.join(outdir, outname)
        st_template.write(outpath, format="SEGY")

        print(f" → Wrote {outname}")

        # Advance index correctly
        i += step

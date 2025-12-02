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
skip_vals = ['32999', '32256', '32569', 
             '38993', '38994', '38995', '38996', '38997',
             '38998', '38999', '39000']#, 
             #'32930', '32931', '32932', '32933', '32934'
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

    #skips_T1: = ['32930','32931','32932','32933','32934','32935','32936']
    
    # Collect & filter files
    flist = [
        os.path.join(ddir, f) for f in os.listdir(ddir)
        if f.endswith('.dat') and f[-9:-4] not in skip_vals
        #and f[-9:-4] not in skips_T1
        ]

    # Ensure correct order (critical!)
    flist = sorted(flist)

    print(f"\nDirectory: {ddir}")
    print(f"Found {len(flist)} usable files.")

    # Loop through files in 8-file chunks
    for i in range(0, len(flist), stack_size): 

        cstack = flist[i:i + stack_size]

        # Skip last incomplete chunk
        if len(cstack) < stack_size:
            print(f"Skipping last partial chunk ({len(cstack)} files).")
            continue

        # Extract shot numbers from filenames
        shot0 = cstack[0][-9:-4]
        shotf = cstack[-1][-9:-4]

        print(f"Stacking {stack_size} files: {shot0}–{shotf}")

        # Allocate stack array
        d_stack = np.empty((stack_size, n_tr, n_pts))

        # Read all 8 files into stack array
        for j, fname in enumerate(cstack):
            st = op.read(fname)
            d = np.array([tr.data for tr in st])
            if len(d) == n_tr:
                d_stack[j] = d
            else:
                d_stack[j] = np.zeros((n_tr, n_pts))
                d_stack[j][:len(d)] = d

        # Compute average
        d_out = np.mean(d_stack, axis=0)

        # Build output stream from template
        st_template = op.read(f_samp)
        for k, tr in enumerate(st_template):
            tr.data = d_out[k].astype(np.float32)

        # Write SEGY
        outname = f"{shot0}-{shotf}_stack.segy"
        outpath = os.path.join(outdir, outname)
        st_template.write(outpath, format='SEGY')

        print(f" → Wrote {outname}")
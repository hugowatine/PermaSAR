#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine (CRPG)
################################################################################

"""flatsim_quality_listdate.py
-----------------------------

Add a new quality column to list_date using values from list_date and aps_corr2date:
new_col = (Norm(exp(-abs(bp)/100)) + Norm(1/aps)) / 2

Usage: add_column_from_two_txt.py --list_date=<file> --aps_date=<file> --output=<file>
                                  [-h]

Options:
  --file1=<file>           Input main file
  --file2=<file>           Second file providing extra columns
  --output=<file>          Output file with new appended column
  -h --help                Show this screen and exit.
"""

import docopt
import pandas as pd
import numpy as np

# ---- Parse arguments ----
arguments = docopt.docopt(__doc__)

list_date = arguments["--list_date"]
aps_date = arguments["--aps_date"]
output_file = arguments["--output"]

# ---- Load files ----
list_date = pd.read_csv(list_date, delim_whitespace=True, header=None)
aps_date  = pd.read_csv(aps_date,  delim_whitespace=True, header=None)

# Rename columns for clarity
list_date.columns = ["col0", "date", "other", "bp"]
aps_date.columns  = ["date", "aps"]

# ---- Merge based on date ----
merged = list_date.merge(aps_date, on="date", how="left")

# Check for missing matches
if merged["aps"].isna().any():
    print("⚠️ Warning: some dates from list_date have no match in aps_date")

# Extract needed columns AFTER merge
bp  = merged["bp"].values.astype(float)
aps = merged["aps"].values.astype(float)

# ---- Normalization functions ----
def normalize(x):
    return (x / np.max(x)) 

# Norm(exp(-abs(col2)/100))
term1 = normalize(np.exp(-np.abs(bp) / 100.0))

# Norm(col3)
term2 = np.zeros_like(aps)
aps[aps == 0] = 0.2
term2 = normalize(1.0 / aps)

# New column
new_col = ((1+term1) * (1+term2))

# ---- Append to df1 ----
list_date["new_col"] = new_col

# ---- Save ----
list_date.to_csv(output_file, sep=" ", index=False, header=False)

print(f"✅ New column added and saved to {output_file}")

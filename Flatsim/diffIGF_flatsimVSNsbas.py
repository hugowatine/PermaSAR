#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine (CRPG)
################################################################################

"""compare_pairs.py
-------------------

Compare two lists of interferometric pairs and list differences.

Usage: compare_pairs.py [--RMSinterfero=<file>] [--inter_pairs=<file>] [--output=<file>] [-h]

Options:
  --RMSinterfero=<file>    Path to RMSinterfero file (4 columns) [default: RMSinterfero.txt]
  --inter_pairs=<file>     Path to inter_pair.rsc file (2 columns) [default: inter_pair.rsc]
  --output=<file>          Output file to write differences [default: comparaisonIFGs_Flatsim_NSBAS.txt]
  -h --help                Show this screen and exit.

The script compares interferogram date pairs from two files and prints:
- Pairs present in both files.
- Pairs missing from either file.
"""

import docopt
import numpy as np

arguments = docopt.docopt(__doc__)

rms_file = arguments["--RMSinterfero"]
inter_file = arguments["--inter_pairs"]
if arguments["--output"] == None:
    output_file = 'comparaisonIFGs_Flatsim_NSBAS.txt'
else:
    output_file = arguments["--output"]

# Load pairs from RMSinterfero file (columns 2 and 3)
rms_data = np.loadtxt(rms_file, usecols=(1, 2), dtype=str)
rms_pairs = {f"{start}_{end}" for start, end in rms_data}

# Load pairs from inter_pair.rsc file (columns 1 and 2)
inter_data = np.loadtxt(inter_file, usecols=(0, 1), dtype=str)
inter_pairs = {f"{start}_{end}" for start, end in inter_data}

# Compare pairs
common = rms_pairs & inter_pairs
only_in_rms = rms_pairs - inter_pairs
only_in_inter = inter_pairs - rms_pairs

# Console output
print(f"\n✅ Number of common pairs       : {len(common)}")
print(f"❌ Pairs only in FLATSIM ({rms_file})         : {len(only_in_rms)}")
print(f"❌ Pairs only in NSBAS ({inter_file})         : {len(only_in_inter)}\n")

if only_in_rms:
    print("🔸 Pairs only in RMSinterfero (FLATSIM):")
    for p in sorted(only_in_rms):
        print(f"  {p}")

if only_in_inter:
    print("\n🔹 Pairs only in inter_pair.rsc (NSBAS):")
    for p in sorted(only_in_inter):
        print(f"  {p}")

# Write to output file
with open(output_file, 'w') as f:
    f.write(f"Number of common pairs: {len(common)}\n")
    f.write(f"Pairs only in {rms_file}: {len(only_in_rms)}\n")
    f.write(f"Pairs only in {inter_file}: {len(only_in_inter)}\n\n")

    if only_in_rms:
        f.write("Pairs only in RMSinterfero (FLATSIM):\n")
        for p in sorted(only_in_rms):
            f.write(f"  {p}\n")
        f.write("\n")

    if only_in_inter:
        f.write("Pairs only in inter_pair.rsc (NSBAS):\n")
        for p in sorted(only_in_inter):
            f.write(f"  {p}\n")

print(f"\n📝 Summary written to {output_file}")


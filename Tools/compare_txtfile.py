#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine (CRPG)
################################################################################

"""compare_txtfile.py
-------------------

Compare two lists of interferometric pairs, compute correlation, and optionally plot.

Usage: compare_pairs.py [--file1=<file>] [--file2=<file>] [--output=<file>]
                        [--file1_date1_col=<int>] [--file1_date2_col=<int>] [--file1_val_col=<int>]
                        [--file2_date1_col=<int>] [--file2_date2_col=<int>] [--file2_val_col=<int>]
                        [--plot=<yesno>]
                        [-h]

Options:
  --file1=<file>            Path to first file (default: RMSinterfero.txt)
  --file2=<file>            Path to second file (default: inter_pairs.rsc)
  --output=<file>           Output file (default: corelation_txtfile.txt)
  --file1_date1_col=<int>   Column index for date1 in file1 [default: 0]
  --file1_date2_col=<int>   Column index for date2 in file1 [default: 1]
  --file1_val_col=<int>     Column index for value in file1 [default: 2]
  --file2_date1_col=<int>   Column index for date1 in file2 [default: 1]
  --file2_date2_col=<int>   Column index for date2 in file2 [default: 2]
  --file2_val_col=<int>     Column index for value in file2 [default: 3]
  --plot=<yesno>            Show scatter plot of val1 vs val2 [default: no]
  -h --help                 Show this screen and exit.
"""

import docopt
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import numpy as np

# --- Parse arguments ---
arguments = docopt.docopt(__doc__)

file1 = arguments["--file1"]
file2 = arguments["--file2"]
output_file = arguments["--output"]
plot_enabled = arguments["--plot"].lower() == 'yes'

# Convert column indices to int
f1_d1 = int(arguments["--file1_date1_col"])
f1_d2 = int(arguments["--file1_date2_col"])
f1_val = int(arguments["--file1_val_col"])

f2_d1 = int(arguments["--file2_date1_col"])
f2_d2 = int(arguments["--file2_date2_col"])
f2_val = int(arguments["--file2_val_col"])

# --- Load files ---
df1 = pd.read_csv(file1, delim_whitespace=True, comment='#', header=None)
df2 = pd.read_csv(file2, delim_whitespace=True, comment='#', header=None)

# Rename columns for merge
df1 = df1.rename(columns={f1_d1:'date1', f1_d2:'date2', f1_val:'val1'})
df2 = df2.rename(columns={f2_d1:'date1', f2_d2:'date2', f2_val:'val2'})

# --- Merge on date1 and date2 ---
df_merged = pd.merge(df1[['date1','date2','val1']], df2[['date1','date2','val2']], on=['date1','date2'])

corr=None
pval=None
# --- Compute correlation ---
#if len(df_merged) > 0:
#    corr, pval = pearsonr(df_merged['val1'], df_merged['val2'])
#else:
#    corr, pval = None, None
corr=None

# --- Write results ---
if output_file != None:
    with open(output_file, 'w') as f:
        f.write("# Merged pairs (common to both files)\n")
        df_merged.to_csv(f, sep='\t', index=False)
        f.write("\n# Pearson correlation between val1 and val2:\n")
        if corr is not None:
            f.write(f"Correlation: {corr:.4f}, p-value: {pval:.4e}\n")
        else:
            f.write("No common pairs found.\n")
        print(f"Results written to {output_file}")

print(f"Processed {len(df_merged)} common pairs.")
if corr is not None:
    print(f"Pearson correlation: {corr:.4f}, p-value: {pval:.4e}")
else:
    print("No common pairs found.")


# --- Plotting ---
if plot_enabled and len(df_merged) > 0:
    plt.figure(figsize=(6,6))
    plt.scatter(abs(df_merged['val1']), df_merged['val2'], color='dodgerblue', edgecolor='k', alpha=0.7)
    plt.xlabel(f'{file1}')
    plt.ylabel(f'{file2}')
    #plt.title(f'Correlation: {corr:.2f} (p={pval:.2e})')

    #a, b = np.polyfit(df_merged['val1'], df_merged['val2'], 1)
    #x_vals = np.array([df_merged['val1'].min(), df_merged['val1'].max()])
    #y_vals = a * x_vals + b
    #plt.plot(x_vals, y_vals, 'r-', label=f'Regression: y={a:.2f}x+{b:.2f}')
    plt.legend()
    plt.tight_layout()
    plt.show()

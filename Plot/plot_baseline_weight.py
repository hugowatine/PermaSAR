#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine (CRPG)
################################################################################

"""
plot_baseline_hist.py
---------------------

Plot the temporal baseline distribution from a text file, and add a second Y axis
showing the weight as a function of baseline:
      weight(x) = exp(-0.2 / x)

Usage:
  plot_baseline_hist.py --file=<file> --col=<int> [--binsize=<int>] [--show=<yesno>]
  plot_baseline_hist.py -h

Options:
  --file=<file>         Path to the input text file containing baseline values.
  --col=<int>           Column index (0-based OR 1-based) of the baseline.
  --binsize=<int>       Width of histogram bins in days [default: 15].
  --show=<yesno>        Display plot on screen? (yes/no) [default: yes].
  -h --help             Show this help message.
"""

import docopt
import numpy as np
import matplotlib.pyplot as plt


def main():

    # -------------------------------------------------------------------------
    # Parse arguments
    # -------------------------------------------------------------------------
    args = docopt.docopt(__doc__)
    filepath   = args["--file"]
    col_index  = int(args["--col"])
    binsize    = int(args["--binsize"])
    show_plot  = args["--show"].lower() == "yes"

    # Convert 1-based indexing to Python indexing
    if col_index >= 1:
        col_index -= 1

    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------
    try:
        data = np.loadtxt(filepath)
    except Exception as e:
        print(f"❌ Error reading file '{filepath}': {e}")
        return

    baseline = data[:, col_index]

    # Keep baseline between 0 and 365
    baseline = baseline[(baseline > 0) & (baseline <= 365)]

    # -------------------------------------------------------------------------
    # Create bins
    # -------------------------------------------------------------------------
    bins = np.arange(0, 365 + binsize, binsize)

    # -------------------------------------------------------------------------
    # Plot histogram
    # -------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(10, 6))

    n, bins, patches = ax1.hist(
        baseline,
        bins=bins,
        edgecolor="black",
        alpha=0.7
    )

    ax1.set_xlabel("Temporal baseline (days)")
    ax1.set_ylabel("Number of interferograms")
    ax1.set_title("Temporal Baseline Distribution with Weight Function")
    ax1.grid(axis="y", linestyle="--", alpha=0.4)

    # -------------------------------------------------------------------------
    # Second axis for weight function
    # -------------------------------------------------------------------------
    ax2 = ax1.twinx()

    # Dense x for smooth curve
    x_vals = np.linspace(1, 365, 500)  # avoid division by zero

    weight = np.exp(-72 / x_vals) + 0.1

    ax2.plot(x_vals, weight, color="red", linewidth=2, label="exp(-0.2/x)")
    ax2.set_ylabel("Weight exp(-0.2/x)", color="red")
    ax2.tick_params(axis="y", labelcolor="red")

    # Legend on second axis
    ax2.legend(loc="upper right")

    # -------------------------------------------------------------------------
    # Show or save
    # -------------------------------------------------------------------------
    if show_plot:
        plt.show()
    else:
        out_png = "baseline_hist_with_weight.png"
        plt.savefig(out_png, dpi=300)
        print(f"✅ Histogram saved to {out_png}")


if __name__ == "__main__":
    main()

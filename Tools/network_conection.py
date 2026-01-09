#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figure_pairs.py
------------

Usage: figure_pairs.py <table>

Options:
-h | --help         Show this screen
"""

import docopt
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from cmcrameri import cm


def read_pairs_file(file):
    try:
        # Try to load d1 / d2 / bp / bt
        try:
            pairs = np.loadtxt(file, skiprows=2, usecols=(0, 1, 2, 3), delimiter='\t', dtype=str)
        except:
            pairs = np.loadtxt(file, skiprows=2, usecols=(0, 1, 2, 3), delimiter=' ', dtype=str)
    except:
        # Else load d1 / d2
        try:
            pairs = np.loadtxt(file, skiprows=2, usecols=(0, 1), delimiter='\t', dtype=str)
        except:
            pairs = np.loadtxt(file, skiprows=2, usecols=(0, 1), delimiter=' ', dtype=str)
    return pairs


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    table = arguments["<table>"]

    pairs = read_pairs_file(table)
    print("Number of pairs", len(pairs))
    dates1 = pairs[:, 0]
    dates2 = pairs[:, 1]
    all_dates = list(set(list(dates1) + list(dates2)))
    print("Number of dates", len(all_dates))

    dates1 = [datetime.strptime(d, "%Y%m%d") for d in dates1]
    dates2 = [datetime.strptime(d, "%Y%m%d") for d in dates2]
    all_dates = [datetime.strptime(d, "%Y%m%d") for d in all_dates]
    all_dates.sort()

    step_count_1month = np.zeros(len(all_dates) - 1)
    step_count_3month = np.zeros(len(all_dates) - 1)
    step_count_6month = np.zeros(len(all_dates) - 1)
    step_count_more = np.zeros(len(all_dates) - 1)
    for d in range(len(all_dates) - 1):
        step1 = all_dates[d]
        step2 = all_dates[d + 1]

        for p in range(len(pairs)):
            if dates1[p] < dates2[p]:
                d1 = dates1[p]
                d2 = dates2[p]
            else:
                d1 = dates2[p]
                d2 = dates1[p]
            if not (d2 < step1 or d1 > step2):
                baseline = (dates2[p] - dates1[p]).days
                if baseline < 30:
                    step_count_1month[d] += 1
                elif baseline < 90:
                    step_count_3month[d] += 1
                elif baseline < 180:
                    step_count_6month[d] += 1
                else:
                    step_count_more[d] += 1

    all_dates = np.array(all_dates)
    inter = all_dates[:-1] + (all_dates[:-1] - all_dates[1:]) / 2
    cum1 = step_count_1month
    cum2 = cum1 + step_count_3month
    cum3 = cum2 + step_count_6month
    cumtot = cum3 + step_count_more
    cmap = cm.navia
    plt.fill_between(inter, cum1, color=cmap(0), alpha=0.6, label='<1month')
    plt.fill_between(inter, cum2, cum1, color=cmap(1/6), alpha=0.6, label='<3months')
    plt.fill_between(inter, cum3, cum2, color=cmap(2/6), alpha=0.6, label='<6months')
    plt.fill_between(inter, cumtot, cum3, color=cmap(1/2), alpha=0.6, label='>6months')
    plt.plot(inter, cumtot, color='k')
    plt.title("Pair Coverage")
    plt.xlabel("Date Interval")
    plt.legend()
    plt.tight_layout()
    plt.show()
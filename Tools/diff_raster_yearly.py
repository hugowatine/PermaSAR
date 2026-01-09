#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package
# written in Python-Gdal
#
############################################
# Author        : Hugo Watine
############################################

"""\
diff_raster_yearly.py
-------------
Compute yearly differences between two dates from a multi-band raster cube

Usage:
    diff_raster_yearly.py --infile=<path> --date1=<MMDD> --date2=<MMDD> [--vmin=<value>] [--vmax=<value>]

Options:
    -h --help        Show this help message
    --infile PATH    Multi-band raster cube
    --date1 DATE     Start date (YYYYMMDD)
    --date2 DATE     End date (YYYYMMDD)
    --vmin VALUE     Min value for color scale [default: auto]
    --vmax VALUE     Max value for color scale [default: auto]
"""

print()
print("Author: Hugo Watine")
print()

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from osgeo import gdal
import docopt
import math

arguments = docopt.docopt(__doc__ if "__doc__" in globals() else "")
cube_file = arguments["--infile"]
month1 = int(arguments["--date1"])
month2 = int(arguments["--date2"])
vmin_arg = arguments.get("--vmin")
vmax_arg = arguments.get("--vmax")

# -----------------------------
# Functions
# -----------------------------
def read_cube_with_dates(path):
    ds = gdal.Open(path)
    nb = ds.RasterCount
    dates = []
    for i in range(1, nb+1):
        desc = ds.GetRasterBand(i).GetDescription()
        if "tserie" in desc:
            date_str = desc.split(":")[1].strip()
            dates.append((i, datetime.strptime(date_str, "%Y%m%d")))
    return ds, sorted(dates, key=lambda x: x[1])

def compute_yearly_differences(ds, dates_sorted, month1, month2):
    """
    Pour chaque année du dataset :
        - trouve la bande la plus proche du mois1 de l'année N
        - trouve la bande la plus proche du mois2 de l'année N+1
        - calcule la différence (date2 - date1)
    """
    diffs = {}
    # obtenir toutes les années disponibles
    years = sorted({d.year for _, d in dates_sorted})
    for y in years:
        # bande la plus proche de mois1 dans l'année y
        candidates1 = [(idx, d) for idx, d in dates_sorted if d.year == y]
        if not candidates1:
            continue
        b1_idx, b1_date = min(candidates1, key=lambda x: abs(x[1].month - month1))

        # bande la plus proche de mois2 dans l'année y+1
        candidates2 = [(idx, d) for idx, d in dates_sorted if d.year == y+1]
        if not candidates2:
            continue
        b2_idx, b2_date = min(candidates2, key=lambda x: abs(x[1].month - month2))

        arr1 = ds.GetRasterBand(b1_idx).ReadAsArray().astype(float)
        arr2 = ds.GetRasterBand(b2_idx).ReadAsArray().astype(float)
        diffs[y] = (arr2 - arr1, b1_date, b2_date)

    return diffs

def plot_yearly_differences(diffs, exclude_year=None):
    if exclude_year is not None:
        diffs = {year: v for year, v in diffs.items() if year != exclude_year}

    n = len(diffs)
    if n == 0:
        print("Aucune différence trouvée !")
        return

    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 4*rows))
    axes = np.array(axes).reshape(rows, cols)

    # ---- 1) Nettoyage + collecte des données pour les percentiles globaux ----
    all_arrays = []

    for year, (arr, d1, d2) in diffs.items():
        arr = arr.astype(float).copy()
        arr[arr == 0.0] = np.nan
        all_arrays.append(arr.ravel())
        diffs[year] = (arr, d1, d2)  # mise à jour dans la structure

    # ---- 2) Calcul global des percentiles ----
    all_data = np.concatenate(all_arrays)
    vmin = np.nanpercentile(all_data, 2)
    vmax = np.nanpercentile(all_data, 98)

    # ---- 3) Affichage ----
    for ax, (year, (arr, d1, d2)) in zip(axes.flatten(), sorted(diffs.items())):
        im = ax.imshow(arr, vmin=vmin, vmax=vmax, cmap="RdBu")
        ax.set_title(f"{year}: {d1.strftime('%Y%m%d')} – {d2.strftime('%Y%m%d')}")
        ax.axis("off")

        cbar = fig.colorbar(im, ax=ax, orientation="vertical",
                            fraction=0.046, pad=0.04)
        cbar.set_label("Déformation (m)")

    # Désactiver les axes vides
    for ax in axes.flatten()[n:]:
        ax.axis("off")

    plt.tight_layout()
    plt.show()

# -----------------------------
# Main
# -----------------------------
ds, dates_sorted = read_cube_with_dates(cube_file)
diffs = compute_yearly_differences(ds, dates_sorted, month1, month2)
plot_yearly_differences(diffs, exclude_year=2014)

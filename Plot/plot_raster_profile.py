#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot_cube_profile.py
--------------
Plot a profile line of a given cube with the option to select a specific band and visualize the band with the profile overlay.

Usage:
    plot_cube_profile.py --cube=<path> --lstart=<value> --lend=<value> --list_dates=<path> [--scale=<value>] [--band=<value>] [--aspect=<path>] [--amplitude=<path>] [--cbarmarks=<values>] [--dem=<path>] [--window=<value>] [--ylim=<value>] [--crop=<values>] [--ymin_max=<values>]
    plot_cube_profile.py --help

Options:
-h | --help             Show this screen
--cube                  Path to cube file
--lstart                Start coordinate of line (x,y)
--lend                  End coordinate of line (x,y)
--list_dates            Path to images_retenues file
--band                  Band to read from the cube (1-based index)
--aspect                Path to aspect map
--scale                 Scaling factor to apply to the band [Default: 1]
--amplitude             Path to displacement amplitude
--cbarmarks             Dates to mark in colorbar given date1,date2,date3,... in YYYYMMDD [Default: None]
--dem                   Path to DEM file to plot profile
--window                Window size to smooth the profile [Default: 1]
--ylim                  Set ylim extent as min,max [Default: automatic]
--crop                  Define crop area for cumulative displacement map as xmin,xmax,ymin,ymax [Default: full image]
--ymin_max              Constrain the y axis for the profile line (min,max)
"""

##########
# IMPORT #
##########

import os
import numpy as np
from osgeo import gdal
from pathlib import Path
import docopt
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

#############
# FUNCTIONS #
#############

def read_tif(input_file):
    ds = gdal.OpenEx(input_file, allowed_drivers=['GTiff'])
    ds_band = ds.GetRasterBand(1)
    values = ds_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
    return values

def get_cube_dimension(cube_file):
    ds = gdal.Open(cube_file)
    ncols, nlines = ds.RasterXSize, ds.RasterYSize
    n_img = ds.RasterCount
    return nlines, ncols, n_img

def read_cube_band(cube_file, band_number):
    """
    Reads a specific band from the cube file.
    """
    print(f'Reading band {band_number} from cube...')
    ds = gdal.Open(cube_file)
    if band_number < 1 or band_number > ds.RasterCount:
        raise ValueError(f"Band {band_number} is out of range. Cube has {ds.RasterCount} bands.")
    band = ds.GetRasterBand(band_number)
    nodata = band.GetNoDataValue()
    print(f'Band {band_number} read successfully.')
    return band.ReadAsArray(), nodata

def compute_percentile_limits(data, low_percentile=2, high_percentile=98):
    """
    Compute the 2nd and 98th percentiles of the data for visualization.
    """
    vmin = np.percentile(data, low_percentile)
    vmax = np.percentile(data, high_percentile)
    return vmin, vmax

def plot_band_with_profile(band_data, lstart, lend, profile, band_number, window_size=20):
    """
    Plots the band with the profile overlay and the extracted profile with statistics.
    """
    rows = np.linspace(lstart[1], lend[1], 100).astype(int)
    cols = np.linspace(lstart[0], lend[0], 100).astype(int)

    # Calcul des statistiques de fenêtre
    medians, stds = calculate_profile_stats(profile, window_size)

    # Plot the band with profile overlay
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

    try:
        from matplotlib.colors import LinearSegmentedColormap
        cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
        cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs + "roma.txt"))
        cmap = cmap.reversed()
    except:
        cmap = plt.cm.rainbow

    # Plot the band
    im = ax1.imshow(band_data, cmap=cmap, vmin=np.nanpercentile(band_data, 2), vmax=np.nanpercentile(band_data, 98), interpolation='nearest')
    fig.colorbar(im, ax=ax1, label='Value')  # Ajout de la colorbar au subplot

    ax1.set_title(f'Band {band_number} with Profile Overlay\n(2nd-98th Percentile Scaling)')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')

    # Overlay the profile line
    ax1.plot(cols, rows, color='red', label='Profile')
    
    # Ajouter des marqueurs (grey dots) tous les n_dots
    num_steps = max(abs(lend[1] - lstart[1]), abs(lend[0] - lstart[0]))
    n_dots = int(num_steps / 10)
    for j in range(0, len(cols), n_dots):
        print("coucou", j)
        ax1.plot(cols[j], rows[j], marker='o', color='grey', markersize=3)

    # Adjusting ticks and grid on profile plot
    ax2.set_xticks(np.arange(0, num_steps, n_dots))  # x-axis ticks every n pixels
    ax2.grid(which='both', color='gray', linestyle='--', linewidth=0.5)

    # Affichage du profil avec les statistiques
    ax2.plot(medians, label='Median (Window)', color='green')
    ax2.fill_between(range(len(profile)), medians - stds, medians + stds, color='orange', alpha=0.4, label='± 1 Std. Dev.')
    ax2.set_xlabel('Distance along the profile')
    ax2.set_ylabel('Value')
    ax2.set_title(f'Profile plot from band {band_number}')
    ax2.legend()
    ax2.grid()
    
    ymin, ymax = None, None
    if "--ymin_max" in arguments and arguments["--ymin_max"]:
        ymin, ymax = tuple(map(int, arguments["--ymin_max"].replace(',', ' ').split()))

    ax2.set_ylim([ymin,ymax])
    # Affichage
    plt.tight_layout()
    plt.show()

def calculate_profile_stats(profile, window_size=20):
    """
    Calcule la médiane et l'écart-type pour une fenêtre glissante le long d'un profil.
    Args:
        profile (numpy array): Les valeurs le long du profil.
        window_size (int): Taille de la fenêtre en pixels.
    Returns:
        medians (numpy array): Médiane pour chaque point du profil.
        stds (numpy array): Écart-type pour chaque point du profil.
    """
    half_window = window_size // 2
    medians = np.zeros_like(profile)
    stds = np.zeros_like(profile)

    for i in range(len(profile)):
        start = max(0, i - half_window)
        end = min(len(profile), i + half_window + 1)
        window = profile[start:end]
        medians[i] = np.nanmedian(window)
        stds[i] = np.nanstd(window)

    return medians, stds


########
# MAIN #
########

arguments = docopt.docopt(__doc__)

cube_file = arguments['--cube']
list_dates_file = arguments['--list_dates']
lstart = tuple(map(int, arguments['--lstart'].split(',')))
lend = tuple(map(int, arguments['--lend'].split(',')))
scale_factor = float(arguments['--scale']) if arguments['--scale'] else 1

# Get band to read, default is None
band_number = int(arguments['--band']) if arguments['--band'] else None

# Get cube dimensions
nlines, ncols, n_img = get_cube_dimension(cube_file)

# Read the specific band if specified, otherwise read the entire cube
if band_number:
    band_data, nodata = read_cube_band(cube_file, band_number)
    band_data = band_data*scale_factor
    band_data[band_data == nodata] = np.nan
    print('Nodata = ', nodata)

else:
    print("No band specified. Please use --band=<value> to select a specific band.")
    exit(1)

# Extract the profile

window_size = int(arguments['--window']) if arguments['--window'] else 20

rows = np.linspace(lstart[1], lend[1], 100).astype(int)
cols = np.linspace(lstart[0], lend[0], 100).astype(int)
profile = band_data[rows, cols]

# Plot the band with profile overlay and the profile
plot_band_with_profile(band_data, lstart, lend, profile, band_number, window_size)


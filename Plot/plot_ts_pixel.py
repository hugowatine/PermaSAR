#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine (CRPG) 
################################################################################

"""\
invers_disp_pixel.py
-------------
Temporal decomposition of the time series delays of selected pixels (used depl_cumule (BIP format) and images_retenues, output of invers_pixel). 

Usage: plot_ts__pixel.py --cols=<values> --ligns=<values> [--cube=<path>] [--list_images=<path>] [--lectfile=<path>] [--aps=<path>] [--vmax_vmin=<values>] [<iref>] [<jref>]

-h --help               Show this screen
--ncols VALUE           Pixel column numbers (eg. 200,400,450)
--nligns VALUE          Pixel lines numbers  (eg. 1200,1200,3000)
--cube PATH             Path to displacement file [default: depl_cumul_flat]
--list_images PATH      Path to list images file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues]
--aps PATH              Path to the APS file giving the error associated to each dates [default: No weigthing]
--vmax_vmin             Set the max,min values of the TS for all pixels
iref                  colum numbers of the reference pixel [default: None] 
jref                  lign number of the reference pixel [default: None]
"""

from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
import sys
import docopt
import struct
import pandas as pd
from datetime import datetime
from matplotlib.gridspec import GridSpec
import os

arguments = docopt.docopt(__doc__)

if arguments["--list_images"] ==  None:
    listim = "images_retenues"
else:
    listim = arguments["--list_images"]

if arguments["--cube"] ==  None:
    cubef = "depl_cumule"
else:
    cubef = arguments["--cube"]
if arguments["--lectfile"] ==  None:
    infile = "lect.in"
else:
    infile = arguments["--lectfile"]
if arguments["--aps"] ==  None:
    apsf = None
else:
    apsf = arguments["--aps"]
if arguments["<iref>"] ==  None:
    iref = None
else:
    iref = int(arguments["<iref>"])
if arguments["<jref>"] ==  None:
    jref = None
else:
    jref = int(arguments["<jref>"])

## Lecture des lignes et colonnes
ipix = list(map(int,arguments["--cols"].replace(',',' ').split()))
jpix = list(map(int,arguments["--ligns"].replace(',',' ').split()))
if len(ipix) != len(jpix):
    print("ncols and nligns lists are not the same size")
    sys.exit()

# number of pixels
Npix = len(ipix)

print(f"Colonne : {ipix}, Ligne : {jpix}")

## Lecture du fichier liste_image
idates,rms=np.loadtxt(listim, comments='#', usecols=(3,4), unpack=True, dtype='f,f')
formatted_dates = idates
#formatted_dates = [datetime.strptime(str(date), "%Y%m%d").strftime("%Y/%m/%d") for date in idates]

ds = gdal.Open(cubef)
if ds is None:
    print(f"Impossible d'ouvrir le fichier cube {cubef}")
    sys.exit()

# Récupérer le nombre de bandes
n_bands = ds.RasterCount
print(f"Le cube contient {n_bands} bandes (séries temporelles)")

# Initialiser une liste pour stocker les séries temporelles de chaque pixel
series_temporales = []

# Extraction des valeurs pour chaque pixel spécifié
for i in range(Npix):
    col = ipix[i]
    row = jpix[i]
    print(f"Extraction des valeurs pour la colonne {col}, ligne {row}...")

    # Initialiser une liste pour stocker les valeurs temporelles de ce pixel
    ts_values = []

    # Lire directement la valeur du pixel pour chaque bande en une seule fois avec ReadRaster
    for band_idx in range(1, n_bands + 1):
        band = ds.GetRasterBand(band_idx)
        data = band.ReadRaster(col, row, 1, 1, buf_type=gdal.GDT_Float32)
        value = struct.unpack('f', data)[0]  # Décompresser les données
        ts_values.append(value)

    series_temporales.append(ts_values)

ref = np.copy(series_temporales[0])*0
if iref != None and jref != None:
    ref = []
    for band_idx in range(1, n_bands + 1):
        band = ds.GetRasterBand(band_idx)
        data = band.ReadRaster(iref, jref, 1, 1, buf_type=gdal.GDT_Float32)
        value = struct.unpack('f', data)[0]  # Décompresser les données
        ref.append(value)

fig = plt.figure(figsize=(12, 8))
gs = GridSpec(Npix, 3, width_ratios=[1, 1, 2])

for i in range(Npix):
    ax = fig.add_subplot(gs[i, :2])
    ax.errorbar(formatted_dates, np.array(series_temporales[i])-np.array(ref), yerr=rms, fmt='o', alpha=0.5, label=f"Pixel {str(i+1)}  ({ipix[i]}, {jpix[i]})")
    ax.set_ylabel(f'Cumulative displacement (mm)')
    ax.set_xlabel('Time (Year/month/day)')
    ax.legend(loc='upper left')
    if arguments["--vmax_vmin"] != None:
        ymax, ymin = map(int,arguments["--vmax_vmin"].replace(',',' ').split())
        ax.set_ylim(ymin, ymax)

# Affichage de la carte de la dernière bande
last_band = ds.GetRasterBand(n_bands-1).ReadAsArray()
nodata_value = ds.GetRasterBand(n_bands-1).GetNoDataValue()
ax_map = fig.add_subplot(gs[:, 2])

try:
    last_band[last_band == nodata_value] = float('NaN')
except:
    pass
masked_band_data = np.ma.array(last_band, mask=np.isnan(last_band))

try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
    cmap = cmap.reversed()
except:
    cmap=cm.rainbow

cax = ax_map.imshow(masked_band_data, cmap=cmap, vmin=np.nanpercentile(last_band, 2), vmax=np.nanpercentile(last_band, 98), interpolation='nearest')
fig.colorbar(cax, ax=ax_map, orientation='vertical')
ax_map.set_title(f'Map of band : {n_bands-1}')
ax_map.set_xlim(np.min(ipix)-500, np.max(ipix)+500)
ax_map.set_ylim(np.max(jpix)+500, np.min(jpix)-500)

# Ajouter les points des pixels sur la carte
for i in range(Npix):
    ax_map.plot(ipix[i], jpix[i], marker='o', color='k', markersize=8)
    ax_map.text(ipix[i] + 5, jpix[i] + 5, str(i+1), color='k', fontsize=10)
if iref != None and jref != None:
    ax_map.plot(iref, jref, marker='X', color='red', markersize=8)
    ax_map.text(iref - 5, jref - 5, 'ref', color='red', fontsize=10)

# Ajuster la mise en page
plt.tight_layout()
plt.show()

# Fermer le dataset GDAL
ds = None








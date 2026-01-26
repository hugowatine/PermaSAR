#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
col2cor.py
-------------
create .cor and .rsc files from colinear interfero.

Usage: col2cor.py --col=<path> [--plot]\

Options:
-h --help           Show this screen.
--coh PATH          path of col .int file
--plot          Display the input coh map and the generated .cor file.
"""

print()
print()
print('Author: Hugo WATINE')
print()
print()

import numpy as np
from osgeo import gdal, osr, gdalconst
import matplotlib.pyplot as plt

gdal.UseExceptions()

try:
    from nsbas import docopt
except:
    import docopt

import sys
import re
import os

def open_raster(raster_path):
    ds = gdal.Open(raster_path, gdalconst.GA_ReadOnly)
    data = ds.GetRasterBand(1).ReadAsArray()
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize
    geotrans = ds.GetGeoTransform()
    proj = ds.GetProjection()

    return data, Xsize, Ysize, geotrans, proj


# read arguments
arguments = docopt.docopt(__doc__)

coh_path = arguments["--col"]
coh, Xifg, Yifg, geotrans, proj = open_raster(coh_path)

basename = os.path.basename(coh_path)

# Extraction des dates YYYYMMDD-YYYYMMDD
match_dates = re.search(r'(\d{8}-\d{8})', basename)
if not match_dates:
    sys.exit("Erreur: impossible de trouver les dates dans le nom du fichier")

dates = match_dates.group(1)

# Extraction du facteur multilook (ex: 8rlks)
match_rlks = re.search(r'(\d+rlks)', basename)
if not match_rlks:
    sys.exit("Erreur: impossible de trouver le facteur rlks dans le nom du fichier")

rlks = match_rlks.group(1)

# Construction du nom final
output_name = f"{dates}_{rlks}.cor"
output_path = os.path.join(os.path.dirname(coh_path), output_name)



ds = gdal.Open(coh_path)
drv = gdal.GetDriverByName("roi_pac")
dst_ds = drv.Create(output_path, Xifg, Yifg, 2, gdal.GDT_Float32)
band1 = dst_ds.GetRasterBand(1)
band2 = dst_ds.GetRasterBand(2)
band = ds.GetRasterBand(1)
data = band.ReadAsArray()

amp = np.absolute(data)
phi = np.angle(data)
data = [amp, phi]

band1.WriteArray(amp)
band2.WriteArray(amp)

dst_ds.SetGeoTransform(geotrans)
dst_ds.SetProjection(proj)

band1.FlushCache()
band2.FlushCache()
del band1
del band2
del dst_ds

print(f".cor saved: {output_path}")

## Modification of the .rsc file

#RLOOK = re.search(r'_(\d+)rlks', output_path).group(1)
with open(output_path + '.rsc', 'r') as f:
    lines = f.readlines()
#lines.append(f"{'RLOOK'.ljust(41)}{RLOOK}\n")
with open(output_path + '.rsc', 'w') as f:
    f.writelines(lines)

print(f".rsc saved: {output_path+'.rsc'}")

do_plot = arguments["--plot"]
if do_plot:
    print("\n Displaying input and output images...\n")

    # Reload output .cor file (band 1)
    ds_cor = gdal.Open(output_path, gdalconst.GA_ReadOnly)
    cor = ds_cor.GetRasterBand(1).ReadAsArray()
    del ds_cor

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    im0 = axes[0].imshow(coh, cmap="viridis")
    axes[0].set_title("Input COH (tif)")
    plt.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(cor, cmap="viridis")
    axes[1].set_title(".cor output")
    plt.colorbar(im1, ax=axes[1])

    for ax in axes:
        ax.axis("off")

    plt.suptitle(f"flatsim2cor preview: {output_name}")
    plt.tight_layout()
    plt.show()
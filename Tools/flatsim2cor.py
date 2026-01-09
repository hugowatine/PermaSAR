#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
flatsim2cor.py
-------------
create .cor and .rsc files from flatsim data.

Usage: flatsim2int.py --coh=<path> [--plot]\

Options:
-h --help           Show this screen.
--coh PATH          path of Coh .tif file with, in band (1), the spatial coherence. Need to be in the same geometry and resolution than InW
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

coh_path = arguments["--coh"]
coh, Xifg, Yifg, geotrans, proj = open_raster(coh_path)

basename = os.path.splitext(os.path.basename(coh_path))[0]
match = re.search(r'\d{8}_\d{8}.*', basename)
if match:
    basename = match.group()
match = re.search(r'\d{8}_\d{8}', basename)
if match:
    new_dates = match.group().replace('_', '-')
    basename_modified = basename.replace(match.group(), new_dates, 1)
output_path = os.path.join(os.path.dirname(coh_path), basename_modified + ".cor")

drv = gdal.GetDriverByName("roi_pac")
dst_ds = drv.Create(output_path, Xifg, Yifg, 2, gdal.GDT_Float32)

band1 = dst_ds.GetRasterBand(1)
band2 = dst_ds.GetRasterBand(2)

band1.WriteArray(coh)
band2.WriteArray(coh)

dst_ds.SetGeoTransform(geotrans)
dst_ds.SetProjection(proj)

band1.FlushCache()
band2.FlushCache()
del band1
del band2
del dst_ds

print(f".cor saved: {output_path}")

## Modification of the .rsc file

RLOOK = re.search(r'_(\d+)rlks', output_path).group(1)
with open(output_path + '.rsc', 'r') as f:
    lines = f.readlines()
lines.append(f"{'RLOOK'.ljust(41)}{RLOOK}\n")
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

    plt.suptitle(f"flatsim2cor preview: {basename_modified}")
    plt.tight_layout()
    plt.show()
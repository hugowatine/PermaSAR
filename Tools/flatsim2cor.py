#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
flatsim2cor.py
-------------
create .cor and .rsc files from flatsim data.

Usage: flatsim2int.py --coh=<path>\

Options:
-h --help           Show this screen.
--coh PATH          path of Coh .tif file with, in band (1), the spatial coherence. Need to be in the same geometry and resolution than InW
"""

print()
print()
print('Author: Hugo WATINE')
print()
print()

import numpy as np
from osgeo import gdal, osr, gdalconst

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


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
unw_nodata0.py
-------------
Change phase value to 0. por pixel with 0. coherence 

Usage: unw_nodata0.py --unw=<path> \

Options:
-h --help           Show this screen.
--unw PATH          path of .unw file with.
"""
print()
print()
print('Author: Hugo WATINE')
print()
print()

import numpy as np
from osgeo import gdal, osr, gdalconst
import sys
import re
import os

gdal.UseExceptions()

try:
    from nsbas import docopt
except:
    import docopt

def open_raster(raster_path, band_idx=1):
    ds = gdal.Open(raster_path, gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError("Impossible d'ouvrir le fichier")
    band = ds.GetRasterBand(int(band_idx))
    arr = band.ReadAsArray().astype(np.float32)
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize
    geotrans = ds.GetGeoTransform()
    proj = ds.GetProjection()

    return arr, Xsize, Ysize, geotrans, proj

# read arguments
arguments = docopt.docopt(__doc__)

unw_path = arguments["--unw"]

phi, Xifg, Yifg, geotrans, proj = open_raster(unw_path, band_idx=2)
coh, Xcoh, Ycoh, _, _ = open_raster(unw_path, band_idx=1)

if Xifg != Xcoh or Yifg != Ycoh:
    print('ERROR : not the same size')
    sys.exit

phi[coh == 0.] = 0.

output_path = os.path.join(os.path.dirname(unw_path), unw_path)

drv = gdal.GetDriverByName("roi_pac")
dst_ds = drv.Create(output_path, Xifg, Yifg, 2, gdal.GDT_Float32)

band1 = dst_ds.GetRasterBand(1)
band2 = dst_ds.GetRasterBand(2)

band1.WriteArray(coh)
band2.WriteArray(phi)

dst_ds.SetGeoTransform(geotrans)
dst_ds.SetProjection(proj)

band1.FlushCache()
band2.FlushCache()
del band1
del band2
del dst_ds

print(f".unw saved: {output_path}")

## Modification of the .rsc file

RLOOK = re.search(r'_(\d+)rlks', output_path).group(1)
with open(output_path + '.rsc', 'r') as f:
    lines = f.readlines()
lines.append(f"{'RLOOK'.ljust(41)}{RLOOK}\n")
with open(output_path + '.rsc', 'w') as f:
    f.writelines(lines)

print(f".rsc saved: {output_path+'.rsc'}")















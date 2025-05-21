#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
flatsim2int.py
-------------
create .int and .rsc files from flatsim data.

Usage: flatsim2int.py --ifg=<path> --coh=<path>\

Options:
-h --help           Show this screen.
--ifg PATH          path of InW .tif file with, in band (1), the wrapped phase
--coh PATH          path of Coh .tif file with, in band (1), the spatial coherence. Need to be in the same geometry and resolution than InW
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

gdal.UseExceptions()

try:
    from nsbas import docopt
except:
    import docopt

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

ifg_path = arguments["--ifg"]
coh_path = arguments["--coh"]

phi, Xifg, Yifg, geotrans, proj = open_raster(ifg_path)
coh, Xcoh, Ycoh, _, _ = open_raster(coh_path)

if Xifg != Xcoh or Yifg != Ycoh:
    print('ERROR : not the same size')
    sys.exit

complex_data = coh * np.exp(1j*phi)
output_path = ifg_path.split('.')[0] + '.int'

drv = gdal.GetDriverByName("roi_pac")
dst_ds = drv.Create(output_path, Xifg, Yifg, 1, gdal.GDT_CFloat32)

band = dst_ds.GetRasterBand(1)
band.WriteArray(complex_data.astype(np.complex64))
dst_ds.SetGeoTransform(geotrans)
dst_ds.SetProjection(proj)

band.FlushCache()
del band
del dst_ds

print(f".int saved: {output_path}")

## Modification of the .rsc file

RLOOK = re.search(r'_(\d+)rlks', output_path).group(1)
with open(output_path + '.rsc', 'r') as f:
    lines = f.readlines()
lines.append(f"{'RLOOK'.ljust(41)}{RLOOK}\n")
with open(output_path + '.rsc', 'w') as f:
    f.writelines(lines)

print(f".rsc saved: {output_path+'.rsc'}")















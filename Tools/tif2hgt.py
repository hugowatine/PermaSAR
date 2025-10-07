#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
tif2hgt.py
----------
Convert a single-band GeoTIFF to ROI_PAC .hgt format with .rsc file.

Usage: tif2hgt.py --tif=<path>\

Options:
-h --help           Show this screen.
--tif PATH          Path to input single-band GeoTIFF file.
"""

import numpy as np
from osgeo import gdal, gdalconst
import sys
import os
import re

try:
    from nsbas import docopt
except:
    import docopt

def open_raster(raster_path):
    ds = gdal.Open(raster_path, gdalconst.GA_ReadOnly)
    if ds is None:
        print(f"ERROR: Cannot open {raster_path}")
        sys.exit(1)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize
    geotrans = ds.GetGeoTransform()
    proj = ds.GetProjection()
    return data, Xsize, Ysize, geotrans, proj

def write_rsc(rsc_path, Xsize, Ysize, geotrans):
    """
    Write minimal ROI_PAC .rsc file.
    geotrans = (top left x, pixel width, 0, top left y, 0, pixel height)
    Note: pixel height usually negative.
    """
    x_first = geotrans[0]
    y_first = geotrans[3]
    x_step = geotrans[1]
    y_step = geotrans[5]

    with open(rsc_path, 'w') as f:
        f.write(f"WIDTH           {Xsize}\n")
        f.write(f"FILE_LENGTH     {Ysize}\n")
        f.write(f"X_FIRST         {x_first}\n")
        f.write(f"Y_FIRST         {y_first}\n")
        f.write(f"X_STEP          {x_step}\n")
        f.write(f"Y_STEP          {y_step}\n")

def main():
    arguments = docopt.docopt(__doc__)
    tif_path = arguments["--tif"]
    data, Xsize, Ysize, geotrans, proj = open_raster(tif_path)

    basename = os.path.splitext(os.path.basename(tif_path))[0]
    out_hgt_path = os.path.join(os.path.dirname(tif_path), basename + ".hgt")
    out_rsc_path = out_hgt_path + ".rsc"

    # ROI_PAC requires float32, 1 band
    drv = gdal.GetDriverByName("roi_pac")

    dst_ds = drv.Create(out_hgt_path, Xsize, Ysize, 2, gdal.GDT_Float32)
    band1 = dst_ds.GetRasterBand(1)
    band2 = dst_ds.GetRasterBand(2)

    band1.WriteArray(data)
    band2.WriteArray(data)

    dst_ds.SetGeoTransform(geotrans)
    dst_ds.SetProjection(proj)

    print(f"ROI_PAC .hgt saved: {out_hgt_path}")

    # Write .rsc file
    write_rsc(out_rsc_path, Xsize, Ysize, geotrans)
    print(f".rsc saved: {out_rsc_path}")

if __name__ == "__main__":
    main()


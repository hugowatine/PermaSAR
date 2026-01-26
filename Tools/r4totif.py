#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Simon Daout (CRPG)
# Modification by Hugo Watine in September 2025
############################################

"""\
r4totiff.py
-------------
Change real4 file to a Tiff file

Usage: r4totiff.py --infile=<path> --outfile=<path> --ref_file=<path> [--lectfile=<path>]
r4totiff.py -h | --help

Options:
-h --help           Show this screen
--infile PATH       r4 file
--outfile PATH      output file
--lectfile PATH     Path of the lect.in file [default: lect.in]
--ref_file PATH     Path of the ref file
"""

print()
print()
print('Author: Hugo Watine')
print()
print()

import numpy as np
from osgeo import gdal
try:
    from nsbas import docopt
except:
    import docopt
arguments = docopt.docopt(__doc__)

def pixel_to_latlon(dataset, pixel_x, pixel_y):
    # Get the geotransform of the dataset
    geotransform = dataset.GetGeoTransform()

    # Geotransform format:
    # [origin_x, pixel_width, 0, origin_y, 0, pixel_height]

    origin_x = geotransform[0]  # Top-left x (longitude for geographic projection)
    pixel_width = geotransform[1]  # Pixel size in x direction
    origin_y = geotransform[3]  # Top-left y (latitude for geographic projection)
    pixel_height = geotransform[5]  # Pixel size in y direction (usually negative)

    # Calculate the geographic coordinates (lat/lon) of the pixel
    lon = origin_x + pixel_x * pixel_width
    lat = origin_y + pixel_y * pixel_height
    print(geotransform)
    return (lat, lon)

if arguments["--lectfile"] ==  None:
   lecfile = "lect.in"
else:
   lecfile = arguments["--lectfile"]

infile = arguments["--infile"]
outfile = arguments["--outfile"]
ref_file= arguments["--ref_file"]

# read lect.in
ncols, nlines = list(map(int, open(lecfile).readline().split(None, 2)[0:2]))
fid = open(infile, 'r')
m = np.fromfile(fid,dtype=np.float32).reshape((nlines,ncols))
#m[m==0.0] = float('NaN')

#Take metadata from ref_file
ds = gdal.Open(ref_file)
ds_geo = ds.GetGeoTransform()
proj = ds.GetProjection()

drv = gdal.GetDriverByName('GTiff')
dst_ds = drv.Create(outfile, ncols, nlines, 1, gdal.GDT_Float32)
dst_ds.SetGeoTransform(ds_geo)
dst_ds.SetProjection(proj)

dst_band = dst_ds.GetRasterBand(1)
dst_band.WriteArray(m)
dst_band.FlushCache()

print(f'The file "{outfile}" has been successfully saved.')


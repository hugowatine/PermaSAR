#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo Watine (CRPG)
############################################

"""\
int2tif.py
-------------
Change .int file to a Tiff file

Usage: r4totiff.py --infile=<path> --outfile=<path> --ref_file=<path> [--lectfile=<path>]
r4totiff.py -h | --help

Options:
-h --help           Show this screen
--infile PATH       .int file
--outfile PATH      output file
--ref_file PATH     Path of the ref file
"""

print()
print()
print('Author: Hugo Watine')
print()
print()

import numpy as np
from osgeo import gdal
import sys
import os
try:
    from nsbas import docopt
except:
    import docopt

arguments = docopt.docopt(__doc__)


def open_gdal(file, band=1, supp_ndv=None, complex=False):
    """
    Use GDAL to open band as real value
    """
    print('-----')
    print(file)
    if not os.path.isfile(file):
        raise FileNotFoundError('File does not exists: {}'.format(file))
    ds = gdal.Open(file)

    print('dims', ds.RasterXSize, ds.RasterYSize, ds.RasterCount)
    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    data = band.ReadAsArray()
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize
    if complex:
        amp = np.absolute(data)
        phi = np.angle(data)
        data = [amp, phi]
        ndv = 0.0
        if ndv is not None and ndv != np.nan:
            amp[amp==ndv] = np.nan
            phi[phi==ndv] = np.nan
        if supp_ndv is not None and supp_ndv != np.nan:
            amp[amp==supp_ndv] = np.nan
            phi[phi==supp_ndv] = np.nan
    else:
        if ndv is not None and ndv != np.nan:
            data[data==ndv] = np.nan
        if supp_ndv is not None and supp_ndv != np.nan:
            data[data==supp_ndv] = np.nan
    return data, Xsize, Ysize

infile = arguments["--infile"]
outfile = arguments["--outfile"]
ref_file= arguments["--ref_file"]

complex, ncols, nlines = open_gdal(infile, complex=True)

# Read georeferencing from reference file
ds_ref = gdal.Open(ref_file)
geo = ds_ref.GetGeoTransform()
proj = ds_ref.GetProjection()

# Create output 2-band GeoTIFF
driver = gdal.GetDriverByName("GTiff")
dst = driver.Create(outfile, ncols, nlines, 2, gdal.GDT_Float32)

dst.SetGeoTransform(geo)
dst.SetProjection(proj)

dst.GetRasterBand(1).WriteArray(complex[0])
dst.GetRasterBand(1).SetDescription("Amplitude")
dst.GetRasterBand(2).WriteArray(complex[1])
dst.GetRasterBand(2).SetDescription("Phase")

dst.FlushCache()

print(f"\nThe file '{outfile}' has been successfully saved.\n")



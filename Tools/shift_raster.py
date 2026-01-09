#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shift raster pixels in X and Y directions

Usage:
    shift_pixels.py <infile> [--shift_x=<sx>] [--shift_y=<sy>] [--outfile=<outfile>]

Options:
    -h --help           Display this message
    infile              Input raster (.tif)
    shift_x             Shift in X (columns), + right, - left [default: 0]
    shift_y             Shift in Y (rows), + down, - up [default: 0]
    outfile             Output raster (.tif), default: infile_shifted.tif
"""

import numpy as np
from osgeo import gdal
import docopt
import sys
import os

gdal.UseExceptions()

def open_gdal(file, band=1):
    """
    Open a raster band with GDAL and return array, geotransform and projection
    """
    ds = gdal.Open(file, gdal.GA_ReadOnly)
    if ds is None:
        raise FileNotFoundError(f"Cannot open {file}")
    arr = ds.GetRasterBand(band).ReadAsArray()
    geo = ds.GetGeoTransform()
    proj = ds.GetProjection()
    return arr, geo, proj

def shift_array(arr, shift_x=0, shift_y=0):
    """
    Shift array in X and Y
    Positive shift_x -> right
    Positive shift_y -> down
    New positions filled with 0
    """
    out = np.zeros_like(arr)

    # Y (rows)
    if shift_y >= 0:
        y_src_start = 0
        y_src_end = arr.shape[0] - shift_y
        y_dst_start = shift_y
        y_dst_end = arr.shape[0]
    else:
        y_src_start = -shift_y
        y_src_end = arr.shape[0]
        y_dst_start = 0
        y_dst_end = arr.shape[0] + shift_y

    # X (cols)
    if shift_x >= 0:
        x_src_start = 0
        x_src_end = arr.shape[1] - shift_x
        x_dst_start = shift_x
        x_dst_end = arr.shape[1]
    else:
        x_src_start = -shift_x
        x_src_end = arr.shape[1]
        x_dst_start = 0
        x_dst_end = arr.shape[1] + shift_x

    out[y_dst_start:y_dst_end, x_dst_start:x_dst_end] = arr[y_src_start:y_src_end, x_src_start:x_src_end]
    return out

def save_raster(outfile, arr, geo, proj):
    """
    Save array as GeoTIFF with GDAL
    """
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(outfile, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)
    out_ds.SetGeoTransform(geo)
    out_ds.SetProjection(proj)
    out_ds.GetRasterBand(1).WriteArray(arr)
    out_ds.FlushCache()
    out_ds = None
    print(f"✅ Raster saved: {outfile}")

def main():
    args = docopt.docopt(__doc__)
    infile = args["<infile>"]
    shift_x = int(args["--shift_x"]) if args["--shift_x"] is not None else 0
    shift_y = int(args["--shift_y"]) if args["--shift_y"] is not None else 0

    outfile = args["--outfile"] if args["--outfile"] else infile.replace(".tif", "_shifted.tif")

    arr, geo, proj = open_gdal(infile)
    shifted = shift_array(arr, shift_x, shift_y)
    save_raster(outfile, shifted, geo, proj)

if __name__ == "__main__":
    main()


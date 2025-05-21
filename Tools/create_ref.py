#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine (CRPG)
################################################################################


"""\
create_ref.py
-------------
Create a .tif ref file from a cropped tif file.

Usage: create_ref.py --infile=<path> --crop=<value> --output=<value>

Options:

-h | --help             Show this screen
--infile                Path to file that will be cropped and save as a ref file
--crop                  Crop dimensions ymin,ymax,xmin,xmax
--output                Name of the output file
"""
import os, sys
import numpy as np
#from correct_ts_from_gacos import geotransform
from numpy.lib.stride_tricks import as_strided
from osgeo import gdal, osr
import pandas as pd
from pathlib import Path
import shutil
import scipy.ndimage
from dateutil import parser
from matplotlib import pyplot as plt
import time
import multiprocessing as mp

from osgeo.ogr import wkbTINM

#from theano.tensor.random.basic import TruncExponentialRV

try:
    from nsbas import docopt
except:
    import docopt

arguments = docopt.docopt(__doc__)

ds_ref = gdal.OpenEx(arguments["--infile"], allowed_drivers=['GTiff'])
translated_raster = arguments["--outfile"]
ymin,ymax,xmin,xmax = map(int, arguments['--crop'].split(','))

gdal.Translate(translated_raster, ds_ref, srcWin=[xmin,ymin, xmax-xmin, ymax-ymin ])
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
depl_cumul2tif.py
-------------
Convert an output from `invers_disp2coef.py` of Pygdalsar into a .tif cube

Usage: depl_cumul2tif.py --infile=<path> --outfile=<path> [--lectfile=<path>]

Options:
-h --help           Show this screen.
--infile PATH       depl_cumul file to convert [default: depl_cumul_flat]
--outfile PATH      output file
--lectfile PATH     Path of the lect.in file [default: lect.in]
"""

print()
print()
print('Author: Hugo WATINE')
print()
print()

import sys
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
try:
    from nsbas import docopt
except:
    import docopt

# read arguments
arguments = docopt.docopt(__doc__)

if arguments["--infile"] == None:
    infile = 'depl_cumule_flat'
else :
    infile = arguments["--infile"]
if arguments["--infile"] == None:
    lectfile = 'lect_ts.in'
else :
    lectfile = arguments["--infile"]

#ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
#print(ds.RasterXSize)
fid = open(infile, 'r')
print(fid)
ncols, nlines = list(map(int, open(lectfile).readline().split(None, 2)[0:2]))
phi = np.fromfile(fid,dtype=np.float32)[:nlines*ncols].reshape((nlines,ncols))
print("> Driver:   REAL4  band file")
print("> Size:     ", ncols,'x',nlines,'x')
print("> Datatype: FLOAT32")


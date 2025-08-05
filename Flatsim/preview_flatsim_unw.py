#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""
preview_flatsim_unw.py
========================
This script plot jpeg file for a list of interferograms .unw

Usage:
  preview_unw.py  --outputdir=<path> [--homedir=<path>] [--int_list=<path>] [--int_path=<path>] [--prefix=<value>] [--suffix=<value>] [--rlook=<value>] [--nproc=<value>] [--ext=<vaue>] [--wrap=<vaue>] [--dpi=<vaue>]

Options:
  --outputdir PATH    Output directory where .jpeg are saved
  --homedir PATH      Path to home directory  [default: ./]
  --int_list PATH     Text file containing list of interferograms dates in two colums, $data1 $date2 [default: interf_pair.rsc]
  --int_path PATH     Path to interfeorgrams [default: int]
  --prefix=<value>    Prefix name ${prefix}$date1-$date2${suffix}_${rlook}.{ext} [default: '']
  --suffix=<vaue>     Suffix name ${prefix}$date1-$date2${suffix}_${rlook}.{ext}  [default: '']
  --rlook VALUE       Multilook number ${prefix}$date1-$date2${suffix}_${rlook}.unw [default: 2]
  --ext=<vaue>        Extension of the file [default: unw]
  --nproc Value       Number of processor [default: 4]
  --wrap=<value>      Wrapped phase between value for unwrapped files 
  --dpi=<value>       dpi for jpg file [default: 300]
  -h --help           Show this screen
"""

import docopt
import numpy as np
import os, subprocess, glob, sys, shutil
import multiprocessing
import matplotlib.pyplot as plt
from osgeo import gdal
gdal.UseExceptions()
    
def generate_preview(infile, outjpeg, dpi, wrap=None):
    try:
        ds = gdal.OpenEx(infile, allowed_drivers=["ROI_PAC"])
        if ds is None:
            raise ValueError(f"Impossible d’ouvrir le fichier: {infile}")

        band1 = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
        band2 = ds.GetRasterBand(2).ReadAsArray().astype(np.float32)

        # Remplacer les zéros par NaN
        band1[band1 == 0.0] = np.nan
        band2[band2 == 0.0] = np.nan

        # Normalisation bande 1
        band1_norm = (band1 - np.nanmin(band1)) / (np.nanmax(band1) - np.nanmin(band1) + 1e-6)

        vmin = np.nanpercentile(band2, 2)
        vmax = np.nanpercentile(band2, 98)

        if wrap is not None:
            wrap = float(wrap)
            band2 = np.mod(band2 + wrap, 2 * wrap) - wrap
            vmin = -wrap
            vmax = wrap
            band2_norm = (band2 - vmin) / (vmax - vmin + 1e-6)
            cmap = plt.get_cmap('hsv')
        else:
            band2_norm = np.clip((band2 - vmin) / (vmax - vmin + 1e-6), 0, 1)
            cmap = plt.get_cmap('Spectral')
        
        rgba = cmap(band2_norm)

        rgb_mod = np.empty_like(rgba)
        rgb_mod[..., :3] = rgba[..., :3] * band1_norm[..., np.newaxis]
        rgb_mod[..., 3] = 1.0

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(rgb_mod, interpolation='none')
        ax.axis('off')

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)

        plt.savefig(outjpeg, bbox_inches='tight', pad_inches=0, dpi=dpi)
        plt.close()
        return True

    except Exception as e:
        print(f"[Erreur] {infile} -> {e}")
        return False
    

# read arguments
arguments = docopt.docopt(__doc__)
if arguments["--homedir"] == None:
  homedir = os.path.abspath('.')
else:
  homedir=os.path.abspath(arguments["--homedir"])

outputdir=os.path.join(homedir, arguments["--outputdir"]) + '/'

if arguments["--int_path"] == None:
    int_path= os.path.join(homedir, 'int')
else:
    int_path=os.path.join(homedir, arguments["--int_path"]) + '/'

if arguments["--int_list"] == None:
  int_list = os.path.join(homedir,'interf_pair.rsc')
else:
  int_list = os.path.join(homedir,arguments["--int_list"])
if arguments["--prefix"] == None:
    prefix = ''
else:
    prefix=arguments["--prefix"]
if arguments["--suffix"] == None:
    suffix = ''
else:
    suffix=arguments["--suffix"]
if arguments["--rlook"] == None:
    rlook = 2
else:
    rlook = '_' + arguments["--rlook"] + 'rlks'
if arguments["--nproc"] == None:
  nproc = 4
else:
  nproc = int(arguments["--nproc"])

if arguments["--ext"] == None:
   ext = 'unw'
else:
   ext = arguments["--ext"]
if arguments["--dpi"] == None:
   dpi = 300
else:
   dpi = int(arguments["--dpi"])

print('')

date_1,date_2=np.loadtxt(int_list,comments="#",unpack=True,usecols=(0,1),dtype='i,i')
kmax=len(date_1)
# cleanif 
if os.path.exists(outputdir):
  # print('{0}/{1}*{2}{3}.jpeg'.format(outputdir, prefix, suffix, rlook))
  jpeg_files = glob.glob('{0}/{1}*{2}{3}.jpeg'.format(outputdir, prefix, suffix, rlook))
  for f in jpeg_files:
    os.remove(f)
else:
  os.makedirs(outputdir)

if os.path.exists(os.path.join(outputdir, "interf_pair_problems.txt")):
    os.remove(os.path.join(outputdir, "interf_pair_problems.txt"))
    os.remove(os.path.join(outputdir, "interf_pair_success.txt"))
def preview(kk):
    successf = open(os.path.join(outputdir, "interf_pair_success.txt"), "a")
    failf =  open(os.path.join(outputdir, "interf_pair_problems.txt"), "a")
    date1, date2 = date_1[kk], date_2[kk]
    idate = str(date1) + '-' + str(date2) 
    folder =  'int_'+ str(date1) + '_' + str(date2) + '/'
    infile = int_path + folder +  prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.' + ext
    jpeg = int_path + folder +  prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.jpeg'
    outjpeg = outputdir + prefix + str(date1) + '-' + str(date2) + suffix + rlook + '.jpeg'

    subprocess.call("length.pl "+str(infile),shell=True)

    try:
      success = generate_preview(infile, outjpeg, dpi, arguments["--wrap"])
      if not success:
        raise Exception("preview failed")
      print('Create:', outjpeg)
      successf.write("%s %s\n" % (str(date1), str(date2)))
    except:
      failf.write("%s %s\n" % (str(date1), str(date2)))

    successf.close()
    failf.close()

pool = multiprocessing.Pool(nproc)
work = [(kk) for kk in range(kmax)]
pool.map(preview, work)

print() 
print('Write successed IFG in: interf_pair_success.txt')
print('Write failed IFG in: interf_pair_problems.txt')


# print('{0}/int_*/{1}*{2}{3}.jpeg'.format(int_path, prefix, suffix, rlook))
jpeg_files = glob.glob('{0}/int_*/{1}*{2}{3}.jpeg'.format(int_path, prefix, suffix, rlook))
print()
print('Move files into:', outputdir)
for f in jpeg_files:
    shutil.move(f,outputdir)
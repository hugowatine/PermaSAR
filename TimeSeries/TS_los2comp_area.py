#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
#
# PermaSAR: An InSAR post-processing package
# written in Python-Gdal
#
############################################
# Author        : Hugo Watine / Simon Daout
############################################

"""\
TS_los2comp_area.py
-------------
Compute the mean time series of an area and the mod of the seasonality. Read cube and geometrical data in geocoded.

Usage: TS_los2comp_area.py --input_file=<path> [--crop=<values>] [--crop_shp=<path>]

Options:
-h --help               Show this screen.
--input_file<file>      Path to file with all the data for the decomposition  
--crop=<values>         Crop data [default: 0,nlines,0,ncol]
--crop_shp=<file>       Crop data with a shapefile [default: None]


Ajouter option masque

"""

import numpy as np
from osgeo import gdal, ogr
from matplotlib import pyplot as plt
import subprocess
import sys

try:
    from nsbas import docopt
except:
    import docopt

print()
print()
print('Author: Hugo Watine')
print()
print()


### DEFINITIION DES ARGUMENTS ###

arguments = docopt.docopt(__doc__)
if arguments["--input_file"] == '':
    print("No input_file, exit")
    sys.exit()

if arguments["--crop_shp"] ==  None:
    shp = None
else:
    shp = arguments["--crop_shp"]

#if arguments["--crop"] ==  None:
#    crop = [0,nlines,0,ncols]
#else:
#    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
#ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])


### DEFINITION DES FONCTIONS ###

def extract_data_in_shapefile_zone(raster_path, shapefile_path):
    # Ouvrir le raster
    raster_ds = gdal.Open(raster_path)
    if raster_ds is None:
        print(f"Erreur: Impossible d'ouvrir le raster {raster_path}")
        return None

    # Ouvrir le shapefile et obtenir la couche
    shapefile_ds = ogr.Open(shapefile_path)
    if shapefile_ds is None:
        print(f"Erreur: Impossible d'ouvrir le shapefile {shapefile_path}")
        return None

    layer = shapefile_ds.GetLayer()
    feature = layer.GetNextFeature()

    # Obtenir la géométrie du shapefile
    geom = feature.GetGeometryRef()

    # Définir le chemin de sortie pour le nouveau raster découpé
    new_name = f"{raster_path.split('.')[0]}_{shapefile_path.split('/')[-1].split('.')[0]}.tif"

    # Paramètres de découpe pour gdal.Warp
    warp_options = gdal.WarpOptions(
        format='GTiff',
        cutlineDSName=shapefile_path,
        cutlineWhere=None,
        cropToCutline=True,
        dstNodata=np.nan,
        outputBounds=None,
        width=0,
        height=0,
        srcSRS=None,
        dstSRS=None,
        multithread=False,
        resampleAlg='nearest',
        srcAlpha=False,
        dstAlpha=False,
        warpOptions=None,
        errorThreshold=None,
        creationOptions=None,
        callback=None,
        callback_data=None
    )

    # Découper le raster en utilisant gdal.Warp
    gdal.Warp(new_name, raster_ds, options=warp_options)

    return new_name

### MAIN ###


# Ouverture du fichier d'entrée comportant la liste insar et le path pour le DEM

exec(open(arguments["--input_file"]).read())

# Pour chacun des fichier, crop via shapefile

if shp != None:
    
    print("Crop data based on shapefile")
    crop_insar = []

    for i in range(len(insar)):
        crop_insar.append([])
        for j in range(len(insar[i])):
            if insar[i][j].split('.')[-1] in ['tif', 'tiff']:
                crop_insar[i].append(extract_data_in_shapefile_zone(insar[i][j], shp))
            else:
                crop_insar[i].append(insar[i][j])

            print(crop_insar[i][j])
else:
    crop_insar = insar

# Ouverture des fichier raster
cubes = [] #liste des cubes ouverts 3D
dates = [] #liste des dates de cubes 1D
looks = [] #liste des look 2D
heads = [] #liste des head 2D

for i in range(len(crop_insar)):
    print(i)
    cubes.append([]), dates.append([]), looks.append([]), heads.append([])
    
    cubes[i].append(gdal.Open(crop_insar[i][0]).ReadAsArray())
    looks[i].append(gdal.Open(crop_insar[i][2]).ReadAsArray())
    heads[i].append(gdal.Open(crop_insar[i][3]).ReadAsArray())
    looks[i].append(np.loadtxt(crop_insar[i][1], comments='#', usecols=(3), unpack=True,dtype='f'))

print(looks[0])
    

#amp_map=gdal.Open(arguments["--ampfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]











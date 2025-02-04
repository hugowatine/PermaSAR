#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
ENUp2LookHeading.py
-------------
Create a look and head angle file from CosENU FLATSIM AUX file

Usage: ENUp2LookHeading.py --CosENUfile=<path> --outpath=<path> --name=<path>\

Options:
-h --help           Show this screen.
--CosENUfile PATH   tif file with 3 bands : Est, north and Up coeficient.
--outpath PATH      output path
--name value
"""
print()
print()
print('Author: Hugo WATINE')
print()
print()

import numpy as np
from osgeo import gdal, osr
import matplotlib.pyplot as plt
import math

try:
    from nsbas import docopt
except:
    import docopt

# read arguments
arguments = docopt.docopt(__doc__)

def calculate_theta_phi(input_raster, output_theta, output_phi):
    """
    Calcule les rasters theta et phi à partir d'un raster d'entrée composé de 3 bandes (E, N, Up).

    Parameters:
        input_raster (str): Chemin du fichier raster d'entrée.
        output_theta (str): Chemin du fichier raster de sortie pour theta.
        output_phi (str): Chemin du fichier raster de sortie pour phi.
    """
    # Ouvrir le raster d'entrée
    dataset = gdal.Open(input_raster)
    if dataset is None:
        raise FileNotFoundError(f"Impossible d'ouvrir le fichier raster : {input_raster}")

    # Lire les bandes
    band_e = dataset.GetRasterBand(1).ReadAsArray()
    band_n = dataset.GetRasterBand(2).ReadAsArray()
    band_up = dataset.GetRasterBand(3).ReadAsArray()

    # Vérifier que toutes les bandes ont la même taille
    if not (band_e.shape == band_n.shape == band_up.shape):
        raise ValueError("Les dimensions des bandes ne correspondent pas.")

    # Calculer theta et phi
    #theta = np.degrees(np.arctan2(band_e, band_up))  # theta = tan-1(E/Up)
    #phi = np.degrees(np.arctan2(band_e, band_n))     # phi = tan-1(E/N)

    #theta = np.arctan2(band_e, band_up)  # theta = tan-1(E/Up)
    #phi = np.arctan2(band_e, band_n)     # phi = tan-1(E/N)

    r_t = band_e/band_up
    r_p = band_e/band_n

    theta = np.degrees(np.arctan(r_t))
    phi = np.degrees(np.arctan(r_p))


    # Obtenir les métadonnées du raster d'entrée
    cols, rows = dataset.RasterXSize, dataset.RasterYSize
    geotransform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()

    # Créer les rasters de sortie
    create_output_raster(output_theta, theta, cols, rows, geotransform, projection)
    create_output_raster(output_phi, phi, cols, rows, geotransform, projection)

    print(f"Les fichiers raster ont été créés :\n  Theta : {output_theta}\n  Phi : {output_phi}")

# Fonction pour créer un raster de sortie
def create_output_raster(output_path, data, cols, rows, geotransform, projection):
    driver = gdal.GetDriverByName("GTiff")
    out_raster = driver.Create(output_path, cols, rows, 1, gdal.GDT_Float32)
    out_raster.SetGeoTransform(geotransform)
    out_raster.SetProjection(projection)
    out_band = out_raster.GetRasterBand(1)
    out_band.WriteArray(data)
    out_band.SetNoDataValue(np.nan)
    out_band.FlushCache()
    out_raster = None

calculate_theta_phi(arguments["--CosENUfile"], arguments["--outpath"] +'/' + arguments["--name"] + "_Look.tif", arguments["--outpath"] +'/' + arguments["--name"] + "_Head.tif")


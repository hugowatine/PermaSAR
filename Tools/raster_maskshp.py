#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# Raster Masking: Mask a raster file using a shapefile
#
############################################
# Author        : Hugo Watine
############################################

"""\
raster_mask.py
--------------
Mask raster image using a shapefile polygon without changing resolution or extent.

Usage: raster_mask.py --raster=<path> --shapefile=<path> --output=<path> 

Options:
-h --help             Show this screen.
--raster=<file>       Raster to be masked 
--shapefile=<file>    Path to the shapefile containing polygon mask
--output=<file>       Path to the output masked raster
"""

import os
from osgeo import gdal, ogr
import sys
import numpy as np

# Read arguments
import docopt
gdal.UseExceptions()

arguments = docopt.docopt(__doc__)

raster_path = arguments["--raster"]
shapefile_path = arguments["--shapefile"]
output_path = arguments["--output"]


# Ouvrir les fichiers
raster = gdal.Open(raster_path)
shapefile = ogr.Open(shapefile_path)
layer = shapefile.GetLayer()

band = raster.GetRasterBand(1)
nodata_value = band.GetNoDataValue()
if nodata_value is None:
    nodata_value = 0.0  # Valeur fallback par défaut
print(f"NoData value: {nodata_value}")


# Créer le raster de sortie avec les mêmes dimensions
driver = gdal.GetDriverByName("GTiff")
out_raster = driver.Create(output_path, raster.RasterXSize, raster.RasterYSize, raster.RasterCount, gdal.GDT_Float32)
out_raster.SetGeoTransform(raster.GetGeoTransform())
out_raster.SetProjection(raster.GetProjection())

# Créer le raster temporaire de masque dans le même dossier que le shapefile
shapefile_dir = os.path.dirname(shapefile_path)
temp_mask_path = os.path.join(shapefile_dir, "temp_mask.tif")

mask_ds = driver.Create(temp_mask_path, raster.RasterXSize, raster.RasterYSize, 1, gdal.GDT_Byte)
mask_ds.SetGeoTransform(raster.GetGeoTransform())
mask_ds.SetProjection(raster.GetProjection())
gdal.RasterizeLayer(mask_ds, [1], layer, burn_values=[1])

# Lire le masque en array numpy
mask_array = mask_ds.GetRasterBand(1).ReadAsArray()

# Appliquer le masque à chaque bande
for i in range(1, raster.RasterCount + 1):
    band = raster.GetRasterBand(i)
    data = band.ReadAsArray()

    data_masked = np.where(mask_array == 1, data, nodata_value)

    out_band = out_raster.GetRasterBand(i)
    out_band.WriteArray(data_masked)
    out_band.SetNoDataValue(nodata_value)

# Fermer les fichiers
out_raster.FlushCache()
raster = None
shapefile = None
mask_ds = None
out_raster = None

# Supprimer le raster temporaire
if os.path.exists(temp_mask_path):
    os.remove(temp_mask_path)

print(f"Masked raster saved successfully to: {output_path}")

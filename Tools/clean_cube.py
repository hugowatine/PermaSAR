#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
clean_cube.py
-------------
Clean a cube file using the fonction clean_raster.py from PygdalSAR

Usage: clean_raster.py --infile=<path> --outfile=<path> [--plot=<yes/no>] [--filter=<HP/LP>] [--fwindsize=<value>] \
[--mask=<path>] [--threshold=<value>]

Options:
-h --help           Show this screen.
--infile PATH       tif file to clean
--outfile PATH      output file
--filter=<HP/LP>    Apply a high pass (HP) or a low pass (LP) gaussian filter to the image
--mask PATH         .tif file used as mask, ex amplitude data (default: None)
--threshold=<value> threshold value on mask file (Keep pixel with mask > threshold)
--fwindsize=<value> Filter window size (default: 16)
--plot=<yes/no>     plot intermediate result (default:no)
"""

print()
print()
print('Author: Hugo WATINE /  Simon DAOUT')
print()
print()

import time
import sys
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal
try:
    from nsbas import docopt
except:
    import docopt

# read arguments
arguments = docopt.docopt(__doc__)

def clean_raster(band, arguments, nodata):
    if arguments["--filter"] == 'HP':
        m_filter = np.copy(band)
        sum_coef = 0*np.copy(band) +1 # The goal is to take into acount nodata value in the filter

        index = (band == nodata) | np.isnan(band)
        
        m_filter[index] = 0.
        sum_coef[index] = 0.

        band = band - ndimage.gaussian_filter(m_filter, int(arguments["--fwindsize"]))/ndimage.gaussian_filter(sum_coef, int(arguments["--fwindsize"]))

        band[index] = float('nan')
    
    elif arguments["--filter"] == 'LP':
        m_filter = np.copy(band)
        index = np.isnan(band)
        m_filter[index] = 0.
        band = ndimage.gaussian_filter(m_filter, int(arguments["--fwindsize"]))
        band[index] = float('nan')
    
    if arguments["--plot"] == 'yes':
        fig, ax = plt.subplots()
        
        # Définir les limites de la colorbar (vmin, vmax) selon tes préférences
        vmax, vmin = np.nanpercentile(band,98), np.nanpercentile(band,2)
        cmap = 'viridis'  # Ou ou autre colormap
        
        # Afficher la bande filtrée
        cax = ax.imshow(band, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_title('Filtered Band')

        # Masquer les labels des ticks x
        plt.setp(ax.get_xticklabels(), visible=False)

        # Ajouter une colorbar
        divider = make_axes_locatable(ax)
        c = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax, cax=c)

        # Afficher la figure
        plt.show()
    return band

# Ouvrir le fichier raster cube
dataset = gdal.Open(arguments["--infile"], gdal.GA_ReadOnly)
if dataset is None:
    print(f"Erreur lors de l'ouverture du fichier {arguments['--infile']}")

print("> Driver:   ", dataset.GetDriver().ShortName)
print("> Size:     ", dataset.RasterXSize,'x',dataset.RasterYSize,'x',dataset.RasterCount)

# Extraire les dimensions du raster
ncols = dataset.RasterXSize
nrows = dataset.RasterYSize
nbands = dataset.RasterCount
gt = dataset.GetGeoTransform()
proj = dataset.GetProjectionRef()
driver = gdal.GetDriverByName('GTiff')

# Créer un nouveau cube vide pour stocker les bandes traitées
processed_cube = np.zeros((nbands, nrows, ncols))

# Traiter chaque bande individuellement

l_rms = []
l_max = []
l_min = []
for band_index in range(1, nbands + 1):
    print(f"Traitement de la bande {band_index}/{nbands}...")
    start_time= time.time()
    # Lire la bande actuelle
    
    band = dataset.GetRasterBand(band_index)
    nodata_value = band.GetNoDataValue()
    band = band.ReadAsArray(0, 0, ncols, nrows, ncols, nrows)
    
    # Appliquer un traitement à la bande (par exemple, nettoyer ou filtrer)
    min_b = np.nanmin(band.flatten())
    max_b = np.nanmax(band.flatten())
    print('min = ', min_b, ' ; max', max_b)
    print('rms*2 = ', np.nanmean(band.flatten()**2))
    print('rms = ', np.sqrt(np.nanmean(band.flatten()**2)))         
    
    l_rms.append(np.sqrt(np.nanmean(band.flatten()**2)))
    l_max.append(max_b)
    l_min.append(min_b)

    #processed_band = clean_raster(band, arguments, nodata_value)
    
    # Ajouter la bande traitée au cube final
    #processed_cube[band_index - 1, :, :] = processed_band
    
    elapsed_time = time.time() - start_time
    print(f"Temps de traitement pour la bande : {elapsed_time:.2f} secondes")

print('------')
print('data ; min,max :',min(l_min), max(l_max)) 
print('rms ; min,max,mean :', min(l_rms), max(l_rms), np.mean(l_rms))
print('------')

sys.exit()

# Réassembler le cube et l'enregistrer
driver = gdal.GetDriverByName('GTiff')  # Choisir le format de sortie (ex: GeoTIFF)
out_dataset = driver.Create(arguments["--outfile"], ncols, nrows, nbands, gdal.GDT_Float32)
out_dataset.SetGeoTransform(gt)
out_dataset.SetProjection(proj)

# Sauvegarder chaque bande traitée dans le fichier de sortie
for band_index in range(1, nbands + 1):
    out_band = out_dataset.GetRasterBand(band_index)
    out_band.WriteArray(processed_cube[band_index - 1, :, :])
    band_metadata = dataset.GetRasterBand(band_index).GetMetadata()
    if band_metadata:
        out_band.SetMetadata(band_metadata)
# Fermer les fichiers
out_dataset.FlushCache()
dataset = None
out_dataset = None
print(f"Cube traité et enregistré dans {arguments['--outfile']}")

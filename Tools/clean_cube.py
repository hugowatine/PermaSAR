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
import matplotlib.cm as cm

# read arguments
arguments = docopt.docopt(__doc__)

def clean_raster(band, arguments, nodata, mask):
    if arguments["--filter"] == 'HP' and arguments["--mask"] is None:
        save_band = np.copy(band)
        m_filter = np.copy(band)
        sum_coef = 0*np.copy(band) +1 # The goal is to take into acount nodata value in the filter

        index = (band == nodata) | np.isnan(band) 
        
        m_filter[index] = 0.
        sum_coef[index] = 0.

        filtre_sum_coef = ndimage.gaussian_filter(sum_coef, int(arguments["--fwindsize"]))
        filtre_m_filter = ndimage.gaussian_filter(m_filter, int(arguments["--fwindsize"]))

        LP = filtre_m_filter / filtre_sum_coef
        band = band - LP

        band[index] = float('nan')

    if arguments["--filter"] == 'HP' and arguments["--mask"] is not None:
        save_band = np.copy(band)
        m_filter = np.copy(band)
        sum_coef = np.ones_like(band)  # The goal is to take into acount nodata value in the filter

        index = (band == nodata) | np.isnan(band) 
        index_mask = (mask > float(arguments["--threshold"]))
        
        m_filter[index] = 0.
        m_filter[index_mask] = 0.
        sum_coef[index] = 0.
        sum_coef[index_mask] = 0.

        filtre_sum_coef = ndimage.gaussian_filter(sum_coef, int(arguments["--fwindsize"]))
        filtre_m_filter = ndimage.gaussian_filter(m_filter, int(arguments["--fwindsize"]))

        LP = filtre_m_filter / filtre_sum_coef
        band = band - LP

        band[index] = float('nan')
        #band[index_mask] = save_band[index_mask] + band[index_mask]
    

    elif arguments["--filter"] == 'LP':
        m_filter = np.copy(band)
        index = np.isnan(band)
        m_filter[index] = 0.
        band = ndimage.gaussian_filter(m_filter, int(arguments["--fwindsize"]))
        band[index] = float('nan')
    
    if arguments["--plot"] == 'yes':

        try:
            from matplotlib.colors import LinearSegmentedColormap
            cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
            cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
            cmap = cmap.reversed()
        except:
            cmap=cm.rainbow

        vmin, vmax = -8, 8
        # Définir les limites de la colorbar (vmin, vmax) en excluant les valeurs extrêmes
        #vmax, vmin = np.nanpercentile(save_band, 98), np.nanpercentile(save_band, 2)
        #cmap = 'viridis'  # Colormap utilisée

        # --- 1. Image originale (save_band) ---
        plt.figure(figsize=(6, 5))
        plt.title('Original (save_band)')
        cax1 = plt.imshow(save_band, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        plt.colorbar(cax1, fraction=0.046, pad=0.04)
        plt.show()

        plt.figure(figsize=(6, 5))
        plt.title('filtre masqué (1)')
        cax1 = plt.imshow(filtre_m_filter, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        plt.colorbar(cax1, fraction=0.046, pad=0.04)
        plt.show()

        plt.figure(figsize=(6, 5))
        plt.title('somme pixel coef (2)')
        cax1 = plt.imshow(filtre_sum_coef, cmap=cmap, vmin=0.4, vmax=1, interpolation='nearest')
        plt.colorbar(cax1, fraction=0.046, pad=0.04)
        plt.show()

        # --- 2. Image filtrée (band) ---
        plt.figure(figsize=(6, 5))
        plt.title('Filtered (band)')
        cax2 = plt.imshow(band, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        plt.colorbar(cax2, fraction=0.046, pad=0.04)
        plt.show()

        plt.figure(figsize=(6, 5))
        plt.title('Low pass (band)')
        cax2 = plt.imshow(LP, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        plt.colorbar(cax2, fraction=0.046, pad=0.04)
        plt.show()

        if arguments["--mask"] is not None:
            # --- 3. Masque utilisé (index_mask) ---
            plt.figure(figsize=(6, 5))
            plt.title('Mask (index_mask)')
            cax3 = plt.imshow(index_mask, cmap='gray', interpolation='nearest')
            plt.colorbar(cax3, fraction=0.046, pad=0.04)
            plt.show()

    return band

def open_cube(path):
    dataset = gdal.Open(path, gdal.GA_ReadOnly)
    if dataset is None:
        print(f"Erreur lors de l'ouverture du fichier {arguments['--infile']}")
        sys.exit()

    print('Opening raster file : ', path)
    print("> Driver:   ", dataset.GetDriver().ShortName)
    print("> Size:     ", dataset.RasterXSize,'x',dataset.RasterYSize,'x',dataset.RasterCount)
    print("")
    return dataset

# Ouvrir le fichier raster cube
dataset = open_cube(arguments["--infile"])

if arguments["--mask"] is not None: 
    maskdata =  open_cube(arguments["--mask"])
    if dataset.RasterXSize != maskdata.RasterXSize:
        print("Erreur, mask and cube file don't have the same X Size")
        sys.exit()
    if dataset.RasterYSize != maskdata.RasterYSize:
        print("Erreur, mask and cube file don't have the same Y Size")
        sys.exit()

# Extraire les dimensions du raster
ncols = dataset.RasterXSize
nrows = dataset.RasterYSize
nbands = dataset.RasterCount
gt = dataset.GetGeoTransform()
proj = dataset.GetProjectionRef()
driver = gdal.GetDriverByName('GTiff')

# Créer un nouveau cube vide pour stocker les bandes traitées
processed_cube = np.zeros((nbands, nrows, ncols))
mask =  np.zeros((nbands, nrows, ncols))

if arguments["--mask"] is not None: 
    mask = maskdata.GetRasterBand(1).ReadAsArray()

# Traiter chaque bande individuellement

#l_rms = []
#l_max = []
#l_min = []
for band_index in range(1, nbands + 1):
    print(f"Traitement de la bande {band_index}/{nbands}...")
    start_time= time.time()
    # Lire la bande actuelle
    #band_index = 170
    band = dataset.GetRasterBand(band_index)
    nodata_value = band.GetNoDataValue()
    band = band.ReadAsArray(0, 0, ncols, nrows, ncols, nrows)
    
    # Appliquer un traitement à la bande (par exemple, nettoyer ou filtrer)
    #min_b = np.nanmin(band.flatten())
    #max_b = np.nanmax(band.flatten())
    #print('min = ', min_b, ' ; max', max_b)
    #print('rms*2 = ', np.nanmean(band.flatten()**2))
    #print('rms = ', np.sqrt(np.nanmean(band.flatten()**2)))         
    
    #l_rms.append(np.sqrt(np.nanmean(band.flatten()**2)))
    #l_max.append(max_b)
    #l_min.append(min_b)

    processed_band = clean_raster(band, arguments, nodata_value, mask)
    
    # Ajouter la bande traitée au cube final
    processed_cube[band_index - 1, :, :] = processed_band
    
    elapsed_time = time.time() - start_time
    print(f"Temps de traitement pour la bande : {elapsed_time:.2f} secondes")

#print('------')
#print('data ; min,max :',min(l_min), max(l_max)) 
#print('rms ; min,max,mean :', min(l_rms), max(l_rms), np.mean(l_rms))
#print('------')

#sys.exit()

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

#!/usr/bin/env python
# -*- coding: utf-8 -*-

from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

# liste harcodé des fichiers pour proxy
fichiers_entree = [
    "RMS_geo_8rlks.tiff",
    "NDSI_reprojected.tif",
    "NDWI_reprojected.tif",
    "inverted_TCoh_geo_8rlks.tiff",
    "CNES_DEM_geo_8rlks_slope_clean.tif"
]

bands_normalized = []
for fichier in fichiers_entree:
    dataset = gdal.Open(fichier)
    band = dataset.GetRasterBand(1).ReadAsArray()
    p2, p98 = np.nanpercentile(band, (2, 98))
    band_stretch = np.clip((band - p2) / (p98 - p2), 0, 1) + 0.1
    bands_normalized.append(band_stretch)


# Multiplication des bands normalisées
resultat = np.ones(np.shape(band_stretch))
for band in bands_normalized:
    resultat = band * resultat

# normalisation
cut = resultat/np.nanpercentile(resultat,98) # higher the value is, higher the cohrence is decreased for unwrapping
coh = np.clip(1. - cut, 0, 1)   # new weight for filtering

# sauvegarde du cutfile et de la coherence
dataset_ref = gdal.Open(fichiers_entree[0])
geotransform = dataset_ref.GetGeoTransform()
projection = dataset_ref.GetProjection()
driver = gdal.GetDriverByName("roi_pac")

dataset_cut = driver.Create("proxy.hgt", dataset_ref.RasterXSize, dataset_ref.RasterYSize, 2, gdal.GDT_Float32)
dataset_cut.SetGeoTransform(geotransform)
dataset_cut.SetProjection(projection)
band_cut = dataset_cut.GetRasterBand(1)
band_cut.WriteArray(cut.astype(np.float32))
band_cut = dataset_cut.GetRasterBand(2)
band_cut.WriteArray(cut.astype(np.float32))

dataset_coh = driver.Create("coherence.cor", dataset_ref.RasterXSize, dataset_ref.RasterYSize, 2, gdal.GDT_Float32)
dataset_coh.SetGeoTransform(geotransform)
dataset_coh.SetProjection(projection)
band_coh = dataset_coh.GetRasterBand(1)
band_coh.WriteArray(coh.astype(np.float32))
band_coh = dataset_coh.GetRasterBand(2)
band_coh.WriteArray(coh.astype(np.float32))

# Fermeture des fichiers
dataset_cut = None
dataset_coh = None
dataset_ref = None

# Affichage des bandes et du résultat avec matplotlib
fig, axs = plt.subplots(3, 2, figsize=(12, 8))
fig.suptitle("Visualisation des bandes et du résultat")

# Affichage des bandes d'entrée
for i, band in enumerate(bands_normalized):
    cax = axs[i//2, i%2].imshow(band, cmap='jet')
    axs[i//2, i%2].set_title(fichiers_entree[i])
    axs[i//2, i%2].axis('off')
    divider = make_axes_locatable(axs[i//2, i%2])
    cax_colorbar = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax, cax=cax_colorbar)

plt.tight_layout()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))
cax = ax1.imshow(cut, cmap="jet", vmax=0.5, vmin=0)
ax1.set_title('Proxy')
ax1.axis('off')
divider = make_axes_locatable(ax1)
cax_colorbar = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=cax_colorbar)

cax = ax2.imshow(coh, cmap="jet", vmin=0, vmax=1)
ax2.set_title('Coherence')
ax2.axis('off')
divider = make_axes_locatable(ax2)
cax_colorbar = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=cax_colorbar)

plt.tight_layout()
plt.show()

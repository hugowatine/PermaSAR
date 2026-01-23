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
resultat = resultat/np.nanmax(resultat)

# sauvegarde resultats
dataset_ref = gdal.Open(fichiers_entree[0])
geotransform = dataset_ref.GetGeoTransform()
projection = dataset_ref.GetProjection()
driver = gdal.GetDriverByName("GTiff")
dataset_sortie = driver.Create("proxy.tif", dataset_ref.RasterXSize, dataset_ref.RasterYSize, 1, gdal.GDT_Float32)
dataset_sortie.SetGeoTransform(geotransform)
dataset_sortie.SetProjection(projection)
band_sortie = dataset_sortie.GetRasterBand(1)
band_sortie.WriteArray(resultat.astype(np.float32))

# Fermeture des fichiers
dataset_sortie = None
dataset_ref = None

# Affichage des bandes et du résultat avec matplotlib
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
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

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
fig.suptitle("Proxy")
# Affichage du résultat
cax = ax.imshow(resultat, cmap="jet", vmax=0.5, vmin=0)
ax.axis('off')
divider = make_axes_locatable(ax)
cax_colorbar = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=cax_colorbar)

plt.tight_layout()
plt.show()

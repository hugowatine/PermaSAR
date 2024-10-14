#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from osgeo import gdal
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import math
from datetime import datetime
from matplotlib.colors import LinearSegmentedColormap

# Chemin du fichier cube
cube_file = '/data2/scratch/A129_V3/CNES_DTs_geo_8rlks.tiff'
dates, idates,rms=np.loadtxt('/data2/scratch/A129_V3/list_images.txt', comments='#', usecols=(1,3,4), unpack=True, dtype='f,f,f')
formatted_dates = idates
#formatted_dates = [datetime.strptime(str(int(date)), "%Y%m%d").strftime("%Y/%m/%d") for date in dates]

# Ouverture du cube avec GDAL
ds = gdal.Open(cube_file, gdal.GA_ReadOnly)
n_bands = ds.RasterCount  # Nombre de bandes (dates)
x_size = ds.RasterXSize   # Largeur du cube (colonnes)
y_size = ds.RasterYSize   # Hauteur du cube (lignes)

# Définition de la zone (col, row) à analyser
# Par exemple, une zone de 10x10 pixels à partir de (50, 50)
zone_start_col = 2300
zone_start_row = 3540
zone_width = 3
zone_height = 400
theta = np.radians(90)

s = [math.sin(theta),- math.cos(theta), 0]
n = [math.cos(theta), math.sin(theta), 0]


# Liste pour stocker les séries temporelles de chaque pixel
all_series_temporales = []
all_velocity = []
analysed_pixel = []
# Extraction des valeurs pour chaque pixel dans la zone
for row_offset in range(zone_height):
    all_series_temporales.append([])
    all_velocity.append([])

    for col_offset in range(zone_width):
        # Calcul des coordonnées de la zone
        orig_col = col_offset  # Colonne d'origine dans la zone
        orig_row = row_offset   # Ligne d'origine dans la zone

        # Rotation des coordonnées par rapport à l'origine (x0, y0)
        #insarx = (orig_col - (zone_width / 2))  # Centrer autour de (0, 0)
        #insary = (orig_row - (zone_height / 2))   # Centrer autour de (0, 0)
        insarx = orig_col
        insary = orig_row
        rotated_xp = (insarx * s[0]) + (insary * n[0])
        rotated_yp = (insarx * s[1]) + (insary * n[1])

        # Ajustement pour l'origine
        rotated_col = int(round(rotated_xp + zone_start_col))
        rotated_row = int(round(rotated_yp + zone_start_row))
        # Vérification que les coordonnées sont dans les limites du raster
        if 0 <= rotated_col < x_size and 0 <= rotated_row < y_size:
            if (rotated_col, rotated_row) in analysed_pixel:
                print(f"Pixel col {rotated_col}, line {rotated_row} already analyzed, not taken into account")
            else :
                print(f"Extraction des valeurs pour la colonne {rotated_col}, ligne {rotated_row}...")
                analysed_pixel.append((rotated_col, rotated_row))
                # Initialiser une liste pour stocker les valeurs temporelles de ce pixel
                ts_values = []
                # Lire directement la valeur du pixel pour chaque bande
                for band_idx in range(1, n_bands + 1):
                    band = ds.GetRasterBand(band_idx)
                    nodata_value = band.GetNoDataValue()
                    data = band.ReadRaster(rotated_col, rotated_row, 1, 1, buf_type=gdal.GDT_Float32)
                    value = struct.unpack('f', data)[0]# Décompresser les données
                    if value == nodata_value:
                        value = np.nan
                    ts_values.append(value)
                if not np.all(np.isnan(ts_values)):
                    all_series_temporales[-1].append(ts_values)
                    valid_mask = ~np.isnan(ts_values)
                    if np.sum(valid_mask) > 1:  # Vérifier qu'il y a au moins deux valeurs valides
                        slope, intercept = np.polyfit(np.array(idates)[valid_mask], np.array(ts_values)[valid_mask], 1)
                        all_velocity[-1].append(slope)
                    else:
                        all_velocity[-1].append(np.nan)

                else:
                    print(f"All NaN values for pixel ({rotated_col}, {rotated_row})")
                    all_series_temporales[-1].append([np.nan] * n_bands)
                    all_velocity[-1].append(np.nan)

#print(all_series_temporales)
#print(all_velocity)
# Calcul de la moyenne des séries temporelles sur la zone
mean_ts_values = np.nanmedian(all_series_temporales, axis=1)
mean_vel_values = np.nanmedian(all_velocity, axis=1)
std_vel_values = np.nanstd(all_velocity, axis=1)

print('Ploting....')
# Affichage 

fig = plt.figure(figsize=(15, 10))  # Taille de la figure
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])  # Grille 2x2

ax1 = fig.add_subplot(gs[0, 0])

ax1.plot(np.arange(1, zone_height + 1), mean_vel_values - 2 * std_vel_values, linestyle='--', color='k', alpha=0.5)
ax1.plot(np.arange(1, zone_height + 1), mean_vel_values, linestyle='-', color='k', alpha=1, label='median')
ax1.plot(np.arange(1, zone_height + 1), mean_vel_values + 2 * std_vel_values, linestyle='--', color='k', alpha=0.5, label='+- 2 sigma')
ax1.set_title('Velocity profile')
ax1.set_xlabel('Distance (pixels)')
ax1.set_ylabel('Velocity (?/years)')
ax1.set_ylim([-2,2])
ax1.grid(True)
ax1.legend()


try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
    cmap = cmap.reversed()
except:
    cmap=cm.rainbow

ax2 = fig.add_subplot(gs[1, 0])
norm = plt.Normalize(vmin=1, vmax=mean_ts_values.shape[1])
for i in range(mean_ts_values.shape[1]):
    color = cmap(norm(i + 1))  # Obtenir la couleur correspondante à la ligne
    ax2.plot(np.arange(1, zone_height + 1), mean_ts_values[:, i], color=color, alpha=0.5)


#axs[1].plot(np.arange(1, zone_height + 1), mean_ts_values, linestyle='-', color='k', alpha=0.5)
ax2.set_title(f'Cube profile ({zone_width}x{zone_height})')
ax2.set_xlabel('Distance (pixels)')
ax2.set_ylabel('Median displacement (?)')
ax2.grid(True)
ax2.set_ylim([-20,30])

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax2, orientation='horizontal', pad=0.3)
tick_locs = np.linspace(1, len(formatted_dates), len(formatted_dates))[::20]  # Un tick sur 10
tick_labels = formatted_dates[::20]  # Correspondance des labels (une date sur 10)

# Appliquer les ticks et labels
cbar.set_ticks(tick_locs)
cbar.set_ticklabels(tick_labels, rotation=45)
cbar.set_label('Dates')


last_band = ds.GetRasterBand(n_bands-1).ReadAsArray()
nodata_value = ds.GetRasterBand(n_bands-1).GetNoDataValue()

try:
    last_band[last_band == nodata_value] = float('NaN')
except:
    pass
masked_band_data = np.ma.array(last_band, mask=np.isnan(last_band))

ax3 = fig.add_subplot(gs[:, 1])

cax = ax3.imshow(masked_band_data, cmap=cmap, vmin=-6, vmax=6, interpolation='nearest')
#cax = ax3.imshow(masked_band_data, cmap=cmap, vmin=np.nanpercentile(last_band, 2), vmax=np.nanpercentile(last_band, 98), interpolation='nearest')
ax3.set_title(f'Map of Band: {n_bands}')
ax3.set_xlim(zone_start_col - 500, zone_start_col + 500)
ax3.set_ylim(zone_start_row + 500, zone_start_row - 500)  # Inverser l'axe Y pour correspondre aux indices de l'image

line_x = [zone_start_col, zone_start_col + 0* s[0] + (zone_height) * n[0]]
line_y = [zone_start_row, zone_start_row + 0* s[1] + ( zone_height) * n[1]] 
ax3.plot(line_x, line_y, color='black', lw=2, label='Zone analysée')

# Ajout de la barre de couleur pour la carte
fig.colorbar(cax, ax=ax3, orientation='vertical')
plt.tight_layout()
print('Showing...')
plt.show()

# Fermeture du dataset
ds = None


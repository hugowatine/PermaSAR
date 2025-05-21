#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine / Simon DAOUT (CRPG-ENSG)
################################################################################

"""\
invers_los2comp_pixel.py
-------------
Decomposition of the ascending and descending time series delays of selecte pixel into slope parallel and slope normal direction (used depl_cumule (BIP format) and images_retenues, output of invers_pixel).

Usage: invers_los2comp_pixel.py --lat=<values> --lon=<values> --cubeAsc=<path> --cubeDsc=<path> --list_imagesAsc=<path> --list_imagesDsc=<path> --lookAsc=<path> --lookDsc=<path> --headAsc=<path> --headDsc=<path> --slope=<path> --aspect=<path> [--apsAsc=<path>] [--apsDsc=<path>] [--windowsize=<value>]

invers_los2comp_pixel.py -h | --help
Options:

-h --help               Show this screen
...
"""

print()
print()
print('Authors: Hugo Watine / Simon Daout')
print()
print()

import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.transform import rowcol
from pyproj import Transformer
from osgeo import gdal, osr
import docopt

def get_pixel_value(raster_path, longitude, latitude, w_size=1):
    """
    Extrait la valeur du pixel à une longitude et latitude données à partir d'un fichier raster.

    Parameters:
        raster_path (str): Chemin vers le fichier raster.
        longitude (float): Longitude de la position du pixel.
        latitude (float): Latitude de la position du pixel.

    Returns:
        np.ndarray: Valeurs des pixels (une valeur si une seule bande, ou un tableau si plusieurs bandes).
    """

    # Ouvrir le fichier raster
    dataset = gdal.Open(raster_path)
    proj = dataset.GetProjection()
    transform = dataset.GetGeoTransform()

    # Convertir les coordonnées géographiques en coordonnées du raster
    source = osr.SpatialReference(wkt=proj)
    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)  # WGS84
    coord_transform = osr.CoordinateTransformation(target, source)
    x, y, _ = coord_transform.TransformPoint(longitude, latitude)  # Nous avons seulement besoin des coordonnées x et y

    # Convertir les coordonnées géographiques en indices de pixel
    col = int((x - transform[0]) / transform[1])
    row = int((y - transform[3]) / transform[5])

    # Lire la valeur du pixel pour chaque bande
    #values = [dataset.GetRasterBand(b).ReadAsArray(col, row, 1, 1)[0, 0] for b in range(1, dataset.RasterCount + 1)]

    # Déterminer les dimensions de la fenêtre
    half_w_size = w_size // 2

    # Liste pour stocker les médianes pour chaque bande
    median_values = []
    # Lire les données pour chaque bande
    for b in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(b)
        # Lire la fenêtre de w_size x w_size autour du pixel central
        array = band.ReadAsArray(col - half_w_size, row - half_w_size, w_size, w_size) #Ne prend le pixel central que si w_size est impair (pour avec 3x3 faut mettre w_size à 3)
        # Calculer la médiane des valeurs dans la fenêtre
        median_value = np.median(array)
        median_values.append(median_value)

    values = median_values

    return values

# read arguments
arguments = docopt.docopt(__doc__)
if None in arguments:
    print("One of the argument is None, run syntetic test")

    t_asc = [-0.5, 0, 1.5, 2, 4, 6, 8]
    asc = [2, 2, 4, 1, 5, 3, 4]
    phi_a = np.radians(-166.7)
    theta_a = np.radians(50)

    t_dsc = [-3, -1, 1, 3, 5, 7, 9]
    dsc = [2, 3, 0, 4, 6, 2, 3]
    phi_d = np.radians(-13.3)
    theta_d = np.radians(50)

    slope = 20
    rot = 40

else:
     if arguments["--windowsize"] == None:
         w_size = 1
     else:
         w_size = int(arguments["--windowsize"])
     
     lon = float(arguments["--lon"])
     lat = float(arguments["--lat"])
     print("Run data for pixel : lon =",lon, " ; lat =",lat)

     t_asc = np.loadtxt(arguments["--list_imagesAsc"], comments='#', usecols=(3), unpack=True,dtype='f')
     asc = get_pixel_value(arguments["--cubeAsc"], lon, lat, w_size=w_size)
     phi_a = np.radians(- get_pixel_value(arguments["--headAsc"], lon, lat, w_size=w_size)[0] - 90) #Phi = -90 -Head
     theta_a = np.radians(- get_pixel_value(arguments["--lookAsc"], lon, lat, w_size=w_size)[0] + 90) #theta = 90 -Look 
     if arguments["--apsAsc"] == None:
        apsAsc = None
     else:
         apsAsc = np.loadtxt(arguments["--apsAsc"], comments='#', unpack=True,dtype='f')

     t_dsc = np.loadtxt(arguments["--list_imagesDsc"], comments='#', usecols=(3), unpack=True,dtype='f')
     dsc = get_pixel_value(arguments["--cubeDsc"], lon, lat, w_size=w_size)
     phi_d = np.radians(- get_pixel_value(arguments["--headDsc"], lon, lat, w_size=w_size)[0] - 90) #Phi = -90 -H
     theta_d = np.radians(- get_pixel_value(arguments["--lookDsc"], lon, lat, w_size=w_size)[0] + 90) #theta = 90 -Look
     if arguments["--apsDsc"] == None:
        apsDsc = None
     else:
         apsDsc = np.loadtxt(arguments["--apsDsc"], comments='#', unpack=True,dtype='f')

     slope = np.radians(get_pixel_value(arguments["--slope"], lon, lat, w_size=w_size)[0])
     rot = - np.radians(get_pixel_value(arguments["--aspect"], lon, lat, w_size=w_size)[0]) # -Aspect (-90 south, +90 north)
     
     print()
     print("Ascending : nb_data=", len(asc), " ; phi=", np.rad2deg(phi_a), " ; theta=", np.rad2deg(theta_a))
     print("Descending : nb_data=", len(dsc), " ; phi=", np.rad2deg(phi_d), " ; theta=", np.rad2deg(theta_d))
     print("Slope=", np.rad2deg(slope))
     print("Rot=", np.rad2deg(rot))
     print()

## 0:Fonction ##
def interpolation_lineaire(x, y, x_interpolation, error=None):
    """
    Calcule les valeurs y_interpolation pour chaque point dans x_interpolation
    en utilisant l'interpolation linéaire seulement dans la plage des valeurs de x.
    Si une valeur de x_interpolation est hors de cette plage, elle est ignorée (renvoie None).
    """
    y_interpolation = []
    y_error=[]
    
    for xi in x_interpolation:
        if x[0] <= xi <= x[-1]:
            # Interpolation linéaire avec numpy
            yi = np.interp(xi, x, y)
            y_interpolation.append(yi)
            if error is None:
                y_error.append(np.nan)
            else:
                yi_error = np.interp(xi, x, error)
                y_error.append(yi_error)
        else:
            # Si xi est en dehors de la plage de x, on met None ou ignore
            y_interpolation.append(None)
    
    return y_interpolation, y_error
    
def projection(slope, theta, phi, rot):
    return [
        np.cos(slope)*np.cos(theta)*(np.cos(rot)*np.cos(phi) - np.sin(rot)*np.sin(phi)) - np.sin(slope)*np.sin(theta),
        np.cos(theta)*(np.sin(rot)*np.cos(phi) + np.cos(rot)*np.sin(phi)),
        np.sin(slope)*np.cos(theta)*(np.cos(rot)*np.cos(phi) - np.sin(rot)*np.sin(phi)) + np.cos(slope)*np.sin(theta)
    ]

def invertion(asc, dsc, phi_a, phi_d, theta_a, theta_d, slope, rot, error_asc, error_dsc):
    proj_asc = projection(slope, theta_a, phi_a, rot)
    proj_dsc = projection(slope, theta_d, phi_d, rot)
    para = []
    perp = []
    sig_para = []
    sig_perp = []
    for i in range(len(asc)):
        data = np.array([asc[i], dsc[i]])
        G = np.array([[proj_asc[0], proj_asc[2]], [proj_dsc[0], proj_dsc[2]]])
        Cd = np.diag([error_asc[i]**2, error_dsc[i]**2])# Diagonal matrix
        Cd[np.isnan(Cd)] = 1.
        
        #Coef
        pars = np.dot(np.linalg.inv(np.dot(np.dot(G.T,np.linalg.inv(Cd)),G)),np.dot(np.dot(G.T,np.linalg.inv(Cd)),data))
        #pars = np.dot(np.linalg.inv(G),data)
        para.append(pars[0])
        perp.append(pars[1])
        
        # Error
        sigma = np.sqrt(np.linalg.inv(np.dot(np.dot(G.T,np.linalg.inv(Cd)),G)))
        sig_para.append(sigma[0][0])
        sig_perp.append(sigma[1][1])

    return para, perp, sig_para, sig_perp
    
## 1:Interpolation ##
interpolation_asc, error_asc_valid = interpolation_lineaire(t_asc, asc, t_dsc, error=apsAsc)

# Filtrer les points où l'interpolation a donné None (hors de portée)
t_dsc_valid = [t for t, y in zip(t_dsc, interpolation_asc) if y is not None]
interpolation_asc_valid = [y for y in interpolation_asc if y is not None]

# Créer la liste dsc_valid avec la même taille et abscisses que interpolation_asc_valid
dsc_valid = [d for d, y in zip(dsc, interpolation_asc) if y is not None]
error_dsc_valid = [d for d, y in zip(apsDsc, interpolation_asc) if y is not None]

## 2:Invertion ##
para, perp, sig_para, sig_perp = invertion(interpolation_asc_valid, dsc_valid, phi_a, phi_d, theta_a, theta_d, slope, rot, error_asc_valid, error_dsc_valid)

## 3:Plot ##
fig, axs = plt.subplots(4, figsize=(8, 10))

axs[0].plot(t_asc, asc, 'ro-', label='asc')
axs[0].errorbar(t_asc, asc, yerr=apsAsc, fmt='none', color='r', alpha=0.5)
axs[0].plot(t_dsc_valid, interpolation_asc_valid, 'g+', label='interpolation asc')
axs[0].errorbar(t_dsc_valid, interpolation_asc_valid, yerr=error_asc_valid, fmt='none', color='g', alpha=0.5)

axs[0].plot(t_dsc, dsc, 'bo-', label='dsc')
axs[0].errorbar(t_dsc, dsc, yerr=apsDsc, fmt='none', color='b', alpha=0.5)
axs[0].plot(t_dsc_valid, dsc_valid, 'y+', label='interpolation dsc')
axs[0].errorbar(t_dsc_valid, dsc_valid, yerr=error_dsc_valid, fmt='none', color='y', alpha=0.5)

axs[1].plot(t_dsc_valid, para, 'mo-', label="para")
axs[1].errorbar(t_dsc_valid, para, yerr=sig_para, fmt='none', color='m', alpha=0.5)

axs[2].plot(t_dsc_valid, perp, 'co-', label="perp")
axs[2].errorbar(t_dsc_valid, perp, yerr=sig_perp, fmt='none', color='c', alpha=0.5)

t = np.arange(min(t_dsc_valid), max(t_dsc_valid), 100)
para = interpolation_lineaire(t_dsc_valid, para, t)
perp = interpolation_lineaire(t_dsc_valid, perp, t)
para_cumulee = np.cumsum(np.nan_to_num(para, nan=0.0))
axs[3].plot(para_cumulee,perp, 'ko')

axs[0].set_title('Données et Interpolation')
axs[1].set_title('Composante Parallèle')
axs[2].set_title('Composante Perpendiculaire')
axs[0].legend()
axs[1].legend()
axs[2].legend()

plt.subplots_adjust(hspace=0.4)
plt.show()



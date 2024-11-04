#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
invert_cube2slope.py
--------------
Project each map of LOS cube in line of steepest slope and saves it in cube. Possible to crop the result to a given extent

Usage:
    invert_cube2slope.py --input_data=<path> --dem=<path> --dest=<path> --ext=<value> [--crop=<value>] [--plot=<yes,no>]
    invert_cube2slope.py --test
    invert_cube2slope.py -h | --help

Options:
-h | --help             Show this screen
--input_data            Path to file with input data (same crop and same size) : geocoded ascending and descending cube, incidence map and heading angle for both.
--dem                   Path to dem file to derive the slope and aspect (same crop and same size)
--dest                  Path to destination 
--ext                   Naming extension of output file 
--crop                  Crop dimensions ymin,ymax,xmin,xmax
--test                  Run synthetic test
--plot                  plot process (yes)

"""

##########
# IMPORT #
##########

import os, sys
import numpy as np
from numpy.lib.stride_tricks import as_strided
from osgeo import gdal
import pandas as pd
from pathlib import Path
import shutil
import scipy.ndimage
from dateutil import parser
from matplotlib import pyplot as plt
import time
try:
    from nsbas import docopt
except:
    import docopt

#############
# FUNCTIONS #
#############

t1 = time.time()

## TEST FUNCTIONS ##
# See flo code

## MAIN FUNCTIONS #

def interpolate(dates, pixel_ts, master_index, plot=False):
    """
    Interpole les séries temporelles de pixels par rapport à la série maîtresse
    spécifiée par master_index, en utilisant l'interpolation linéaire.

    :param dates: Liste de listes contenant les dates pour chaque série temporelle.
    :param pixel_ts: Liste de listes contenant les valeurs des pixels pour chaque série temporelle.
    :param master_index: Index de la série maîtresse à utiliser pour l'interpolation.
    :return: (dates_interpolated, pixel_ts_interpolated)
    """
    # Dates de la série maîtresse
    master_dates = dates[master_index]
    master_ts = pixel_ts[master_index]

    # Prépare les listes pour les valeurs interpolées
    dates_interpolated = []
    pixel_ts_interpolated = []

    for i in range(len(dates)):
        # Ignore la série maîtresse, car elle est déjà à la bonne échelle
        if i == master_index:
            dates_interpolated.append(master_dates)
            pixel_ts_interpolated.append(master_ts)
            continue

        # Dates et valeurs pour la série actuelle
        current_dates = dates[i]
        current_ts = pixel_ts[i]

        # Interpolation linéaire pour les valeurs de la série actuelle
        interpolated_ts = np.interp(master_dates, current_dates, current_ts)

        # Stocke les résultats
        dates_interpolated.append(master_dates)
        pixel_ts_interpolated.append(interpolated_ts)

    if plot:
        for i in range(len(dates)):
            plt.plot(dates[i], pixel_ts[i], 'o-', alpha=0.5)
            plt.plot(dates_interpolated[i], pixel_ts_interpolated[i], '+--', alpha=0.5)
        plt.show()

    return dates_interpolated, pixel_ts_interpolated

def projection(slope, theta, phi, rot):
    return [
        np.cos(slope)*np.cos(theta)*(np.cos(rot)*np.cos(phi) - np.sin(rot)*np.sin(phi)) - np.sin(slope)*np.sin(theta),
        np.cos(theta)*(np.sin(rot)*np.cos(phi) + np.cos(rot)*np.sin(phi)),
        np.sin(slope)*np.cos(theta)*(np.cos(rot)*np.cos(phi) - np.sin(rot)*np.sin(phi)) + np.cos(slope)*np.sin(theta)
    ]

def invertion(dates, pixel_ts, phi, theta, slope, rot ,plot=False):
    proj = [projection(slope, theta[i], phi[i], rot) for i in range(len(dates))]
    TS_U0 = []
    TS_U2 = []

    for day in range(len(dates[0])):
        data = np.array([pixel_ts[i][day] for i in range(len(dates))])
        G = np.array([[proj[i][0], proj[i][2]] for i in range(len(dates))])

        try:
            pars = np.linalg.lstsq(G, data, rcond=None)[0]  # Utilisation de least squares
            TS_U0.append(pars[0])
            TS_U2.append(pars[1])
        except np.linalg.LinAlgError:
            TS_U0.append(np.nan)
            TS_U2.append(np.nan)

    if plot:
        # Affichage des résultats de décomposition
        plt.plot(dates[0], TS_U0, 's-', color='red', label="U0 (composante parallèle)")
        plt.plot(dates[0], TS_U2, 'd-', color='blue', label="U2 (composante perpendiculaire)")

        plt.xlabel("Temps (yr)")
        plt.ylabel("Déplacement (mm)")
        plt.legend()
        plt.title(f"Décomposition spatiale des séries temporelles : slope : {np.rad2deg(slope)}° - aspect {np.rad2deg(rot)}°")
        plt.show()
    return TS_U0, TS_U2

def compute_slope_aspect(path, crop, min_slope=0):

    global filter, Px, Py, slope, aspect

    #### LOAD DEM
    ds = gdal.Open(path,gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    topo = band.ReadAsArray()
    ncols, nlines = ds.RasterYSize, ds.RasterXSize

    fwindsize = 2. ## Filtrage de la topo -> posibilité de mettre une variable ici
    
    filter = scipy.ndimage.gaussian_filter(topo,fwindsize)
    #filter = scipy.ndimage.gaussian_filter(topo,2.)
    gt = ds.GetGeoTransform()
    projref = ds.GetProjectionRef()
    drv = gdal.GetDriverByName('GTiff')

    # Get middle latitude
    data1 = ds.GetGeoTransform()
    lats = data1[3] + (np.arange(nlines) * data1[5])
    lat_moy = np.mean(lats)

    # Get resolution depending on whether the
    res = data1[1]*40075e3/360
    print("Spatial resolution in deg: {}, in meter: {}".format(data1[1],res))
    if res<1 or res>500:
        print("Spatial resolution seems unrealistic. Exit!")
        exit()

    # Calcul gradient
    Py, Px = np.gradient(filter, res, res*np.cos(np.deg2rad(lat_moy)))
    Px = Px.astype(float); Py = Py.astype(float)
    # Slope
    slope = np.arctan(np.sqrt(Py**2+Px**2))
    # smooth slope to avoid crazy values
    slope[np.rad2deg(slope)>80]=0.
    slope[np.rad2deg(slope)<min_slope]=0.
    # Line of max slope
    aspect = np.arctan2(Py,-Px)

    # Create aspect and slope files
    dst = drv.Create('slope.tif', nlines, ncols, 1, gdal.GDT_Float32)
    bandt = dst.GetRasterBand(1)
    bandt.WriteArray(np.rad2deg(slope))
    dst.SetGeoTransform(gt)
    dst.SetProjection(projref)
    bandt.FlushCache()

    dst = drv.Create('aspect.tif', nlines, ncols, 1, gdal.GDT_Float32)
    bandt = dst.GetRasterBand(1)
    bandt.WriteArray(np.rad2deg(aspect))
    dst.SetGeoTransform(gt)
    dst.SetProjection(projref)
    bandt.FlushCache()

    if(crop):
        return -aspect[crop[0]:crop[1], crop[2]:crop[3]], slope[crop[0]:crop[1], crop[2]:crop[3]]
    else:
        return -aspect, slope

def get_cube_dimension(cube_file):
    ds = gdal.Open(cube_file)
    ncols, nlines = ds.RasterXSize, ds.RasterYSize
    n_img = ds.RasterCount
    nodata_value = ds.GetRasterBand(1).GetNoDataValue()

    img_data = (nlines, ncols, n_img, nodata_value)

    return img_data

def read_tif(input_file, crop):

    name = input_file.split('/')[-1]
    
    ds = gdal.OpenEx(input_file, allowed_drivers=['GTiff'])
    print()  
    print(f'Start reading raster file {name}')
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',1)
    ds_band = ds.GetRasterBand(1)
    values = ds_band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)

    if(crop):
        return values[crop[0]:crop[1], crop[2]:crop[3]]
    else:
        return values


def read_cube(cube_file, crop):
    name = cube_file.split('/')[-1]
    print()
    print(f'Start reading cube file {name}')

    ds = gdal.Open(cube_file, gdal.GA_ReadOnly)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(ds.GetRasterBand(1).DataType))
    print("> NoData: ", ds.GetRasterBand(1).GetNoDataValue())
    
    cube = ds.ReadAsArray()
    cube = np.transpose(cube, (1, 2, 0))
    
    cube = np.where(np.isin(cube, [ds.GetRasterBand(1).GetNoDataValue()]), np.nan, cube)

    if crop:
        maps = cube[crop[0]:crop[1], crop[2]:crop[3], :]
    else:
        maps = cube

    return maps

def save_cube(dest_path, out_filename, maps):
    print('Writing cube')

    fid = open(os.path.join(dest_path, out_filename), 'wb')
    maps[:,:,:].flatten().astype('float32').tofile(fid)
    fid.close()

def save_cube_metadata(dest_path, out_filename, img_data):
    nrow, ncol, nimg = img_data[0], img_data[1], img_data[2]

    # be careful here with ncol, nrow (normal: lines=nrow; samples=ncol)
    with open(os.path.join(dest_path, '{}.hdr'.format(out_filename)), 'w') as hdr_file:
        hdr_file.write("ENVI\n")
        hdr_file.write("samples = {}\n".format(ncol))
        hdr_file.write("lines = {}\n".format(nrow))
        hdr_file.write("bands = {}\n".format(nimg))
        hdr_file.write("header offset = 0\n")
        hdr_file.write("file type = ENVI Standard\n")
        hdr_file.write("data type = 4\n")  # 4 represents float32, adjust if needed
   
    # save also .in file to plot pixel
    with open(os.path.join(dest_path, 'lect_{}.in'.format(out_filename)), 'w') as lect_file:
        lect_file.write('\t{}\t{}\t{}'.format(ncol, nrow, nimg))

########
# MAIN #
########

arguments = docopt.docopt(__doc__)
if arguments["--plot"] == None:
    plot=False
else:
    plot=True


### 0: Set input data

##  Get input parameter from file 
if arguments["--input_data"] == '':
    print("No input_data, exit")
    sys.exit()

exec(open(arguments["--input_data"]).read()) # Contient la variable insar

## Get others parameters
dem_path = arguments['--dem']
dest_path = arguments['--dest']
ext = arguments['--ext']

###
## Read and prepare 2D and 1D  input data
###

## Prepare ##
# get cube dimensions from one of the cubes 
# full extent to initially read the cube, keep this if no crop option
img_data = get_cube_dimension(insar[0][0])
nlines, ncols, n_img, nodata = img_data[0], img_data[1], img_data[2], img_data[3]

# set crop parameters if crop is input
if(arguments['--crop']):
    # (xmin, xmax, ymin, ymax)
    crop = tuple(map(int, arguments['--crop'].split(',')))
    
    nlines_crop = crop[1] - crop[0]
    ncols_crop = crop[3] - crop[2]

    print('Crop option is set. New cube dimensions are ({},{},{})'.format(nlines_crop, ncols_crop, n_img))
else:
    crop = None
    print('No crop option is set. Process the full extent of cube')

## Read ##
# read topo files

dem = read_tif(dem_path, crop)
rot, slope = compute_slope_aspect(dem_path, crop, min_slope=0)

if plot:
    ## PLOT TOPO ##
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 1 ligne, 3 colonnes

    im0 = axes[0].imshow(dem, cmap='terrain')
    axes[0].set_title("DEM")
    fig.colorbar(im0, ax=axes[0], orientation="horizontal", fraction=0.046, pad=0.1)

    im1 = axes[1].imshow(rot, cmap='hsv')  # Choisir un colormap adapté
    axes[1].set_title("Aspect (0°Est - rad)")
    fig.colorbar(im1, ax=axes[1], orientation="horizontal", fraction=0.046, pad=0.1)

    im2 = axes[2].imshow(slope, cmap='viridis')  # Colormap pour la pente
    axes[2].set_title("Slope (rad)")
    fig.colorbar(im2, ax=axes[2], orientation="horizontal", fraction=0.046, pad=0.1)

    plt.tight_layout()
    plt.show()
    ### PLOT TOPO ###

###
## read insar data
###

crop_insar = [] 
# Will be organised as folow : [data1, data2...]
# with data1 = [cube (3D), dates (1D), look (2D), head (2D)]

for i in range(len(insar)):
    crop_insar.append([])
    for j in range(len(insar[0])):
        if insar[i][j].split('.')[-1] in ['tif', 'tiff']:
            if j == 0: #Read an asc or dsc 3D cube
                crop_insar[i].append(read_cube(insar[i][j], crop))
            else :
                crop_insar[i].append(read_tif(insar[i][j], crop))
        else :
            crop_insar[i].append(np.loadtxt(insar[i][j], comments='#', usecols=(3), unpack=True,dtype='f'))
    print('----------------------------------------')

if plot:
    ## PLOT ASC and DSC data ##
    fig, axes = plt.subplots(len(crop_insar), 3, figsize=(15, 10))  # 1 ligne, 3 colonnes

    for i in range(len(insar)):
        last_band_index = crop_insar[i][0].shape[2] - 1
        im0 = axes[i][0].imshow(crop_insar[i][0][:, :, last_band_index], vmin=np.nanpercentile(crop_insar[i][0][:, :, last_band_index],2), vmax=np.nanpercentile(crop_insar[i][0][:, :, last_band_index],98), cmap='rainbow', interpolation='nearest')
        axes[i][0].set_title(insar[i][0].split('/')[-1], pad=1)  # Ajoute de l'espace sous le titre
        fig.colorbar(im0, ax=axes[i][0], orientation="horizontal", fraction=0.046, pad=0.1)

        im1 = axes[i][1].imshow(crop_insar[i][2], cmap='YlGnBu')
        axes[i][1].set_title("Look angle", pad=1)  # Ajoute de l'espace sous le titre
        fig.colorbar(im1, ax=axes[i][1], orientation="horizontal", fraction=0.046, pad=0.1)

        im2 = axes[i][2].imshow(crop_insar[i][3], cmap='YlGnBu')
        axes[i][2].set_title("Heading angle", pad=1)  # Ajoute de l'espace sous le titre
        fig.colorbar(im2, ax=axes[i][2], orientation="horizontal", fraction=0.046, pad=0.1)

    plt.subplots_adjust(wspace=0.3, hspace=0.4)  # Ajuste wspace et hspace selon tes besoins
    plt.tight_layout()
    plt.show()
    ## PLOT ASC and DSC data ##

###
## Define angles for projection
###

# convert head, look to angle phi, theta in rad
print('Average LOOK:{0:.5f}, HEADING:{1:.5f} angles'.format(np.nanmean(crop_insar[i][2]),np.nanmean(crop_insar[i][3])))

for i in range(len(insar)):
    crop_insar[i][2] = np.deg2rad(90. - crop_insar[i][2]) # theta: angle between LOS and horizontal (positive in anti-clockwise direction)
    crop_insar[i][3] = np.deg2rad(-90. - crop_insar[i][3]) #phi: horizontal angle between LOS and East (positive in anti-clockwise direction)

print('Average THETA:{0:.5f}, PHI:{1:.5f} angles'.format(np.nanmean(np.rad2deg(crop_insar[i][2])), np.nanmean(np.rad2deg(crop_insar[i][3]))))

###
##### 3:  Interpolate and inverse all pixel of all bands 
###

# Define the number of bands of the output cubes = cube with lower amount of bands
nimg_final = float('inf')  # Utilise l'infini pour trouver le minimum
index = -1

for i in range(len(crop_insar)):
    nimg = crop_insar[i][0].shape[2]  # Obtenir le nombre de bandes
    if nimg < nimg_final:
        nimg_final = nimg
        index = i

nlines_final = crop_insar[index][0].shape[0]
ncols_final = crop_insar[index][0].shape[1]
img_data_final = (nlines_final, ncols_final, nimg_final)

print()
print(f"Each time series of each cube will be interpolated based on the master {insar[index][0].split('/')[-1]} with {nimg_final} dates")

# Create output cube
U0 = np.zeros((nlines_final, ncols_final, nimg_final))
U2 = np.zeros((nlines_final, ncols_final, nimg_final))

# Create the list the date of each bands for each cube
dates = [crop_insar[i][1] for i in range(len(crop_insar))]

# Interpolate and decompse each pixel
## IL FAUT EMPECHER DE FAIRE CES OP2RATIONS SI Y'A PAS DEUX DONNEES
## OU SI ON MASQUE LES ZONES

# Mask data 
mask_ignore = np.zeros((nlines_final, ncols_final), dtype=bool)
for j in range(nlines_final):
    for k in range(ncols_final):
        if slope[j, k] < np.deg2rad(2.):
            mask_ignore[j, k] = True
            continue
        
        if (np.deg2rad(-120.) <= rot[j, k] <= np.deg2rad(-60.)) or (np.deg2rad(60.) <= rot[j, k] <= np.deg2rad(120.)):
            mask_ignore[j, k] = True
            continue

for j in range(nlines_final):
    if j%20 ==0 :
        print(f'Start inverting line {j}')
    for k in range(ncols_final):
        if mask_ignore[j, k]:
            U0[j, k, :] = np.nan
            U2[j, k, :] = np.nan
            continue

        pixel_ts = []
        phi = []
        theta = []
        # Interpolate TS
        for i in range(len(crop_insar)):
            pixel_ts.append(crop_insar[i][0][j, k])
            phi.append(crop_insar[i][3][j, k])
            theta.append(crop_insar[i][2][j, k])

        if j == 0 and k == 0 and plot:
            dates_interpolated, pixel_ts_interpolated = interpolate(dates, pixel_ts, index, plot=plot)
        else:
            dates_interpolated, pixel_ts_interpolated = interpolate(dates, pixel_ts, index, plot=False)
        
        # Decompose TS
        pixel_U0, pixel_U2 = invertion(dates_interpolated, pixel_ts_interpolated, phi, theta, slope[j,k], rot[j,k] ,plot)
        U0[j, k, :] = pixel_U0
        U2[j, k, :] = pixel_U2

band_index = -1  # dernière bande (ajustez selon votre besoin)

plt.figure(figsize=(10, 5))

# Affichage de la bande de U0
plt.subplot(1, 2, 1)
plt.imshow(U0[:, :, band_index], cmap='rainbow', interpolation='nearest', vmin=np.nanpercentile(U0[:, :, band_index],2), vmax=np.nanpercentile(U0[:, :, band_index],98))
plt.colorbar(label="Déplacement U0")
plt.title(f"Bande U0 index {band_index}")

# Affichage de la bande de U2
plt.subplot(1, 2, 2)
plt.imshow(U2[:, :, band_index], cmap='rainbow', interpolation='nearest', vmin=np.nanpercentile(U2[:, :, band_index],2), vmax=np.nanpercentile(U2[:, :, band_index],98))
plt.colorbar(label="Déplacement U2")
plt.title(f"Bande U2 index {band_index}")

plt.show()

## 2: Save the inverted cubes as file

if(crop):
    U0_filename = 'depl_cumule_U0_{}_crop_{}_{}_{}_{}'.format(ext, crop[0], crop[1], crop[2], crop[3])
    U2_filename = 'depl_cumule_U2_{}_crop_{}_{}_{}_{}'.format(ext, crop[0], crop[1], crop[2], crop[3])
else:
    uslope_filename = 'depl_cumule_uslope_{}'.format(ext)
    uz_filename = 'depl_cumule_uz_{}'.format(ext)

save_cube(dest_path, U0_filename, U0)
save_cube_metadata(dest_path, U0_filename, img_data_final)

save_cube(dest_path, U2_filename, U2)
save_cube_metadata(dest_path, U2_filename, img_data_final)

print(f"Time : {time.time() - t1}") 

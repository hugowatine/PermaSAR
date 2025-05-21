#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
invert_cube2slope.py
--------------
Project each map of LOS cube in line of steepest slope and saves it in cube. Possible to crop the result to a given extent

Usage:
    invert_cube2slope.py --input_data=<path> --dem=<path> --dest=<path> --ext=<value> [--crop=<value>] [--plot=<yes,no>] [--proc=<value>] [--block_size=<value]
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
--proc                  number of processes to run at the same time [default: 4 proc]
--block_size            Size of block annalysed by one proc in pixels [default: 200 lines]

"""

##########
# IMPORT #
##########

import os, sys
import numpy as np
#from correct_ts_from_gacos import geotransform
from numpy.lib.stride_tricks import as_strided
from osgeo import gdal, osr
import pandas as pd
from pathlib import Path
import shutil
import scipy.ndimage
from dateutil import parser
from matplotlib import pyplot as plt
import time
import multiprocessing as mp

from osgeo.ogr import wkbTINM

#from theano.tensor.random.basic import TruncExponentialRV

try:
    from nsbas import docopt
except:
    import docopt


t1 = time.time()
#############
# FUNCTIONS #
#############

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

        if np.all(np.isnan(current_ts)):
            interpolated_ts = np.full(len(master_dates), np.nan)
            continue
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
            if np.isnan(G).any() or np.isnan(data).any():
                np.nan_to_num(G, nan=0.0)
                np.nan_to_num(data, nan=0.0)
            pars = np.linalg.lstsq(G, data, rcond=None)[0]  # Utilisation de least squares
            TS_U0.append(pars[0])
            TS_U2.append(pars[1])
        except np.linalg.LinAlgError:
            print('Warning LinAlgError')
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
        return -aspect[crop[0]:crop[1]+1, crop[2]:crop[3]+1], slope[crop[0]:crop[1]+1, crop[2]:crop[3]+1]
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
        return values[crop[0]:crop[1]+1, crop[2]:crop[3]+1]
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
        maps = cube[crop[0]:crop[1]+1, crop[2]:crop[3]+1, :]
    else:
        maps = cube

    return maps

########
# MAIN #
########

arguments = docopt.docopt(__doc__)
if arguments["--plot"] == None or arguments["--plot"] =='no':
    plot=False
else:
    plot=True

####
## 0: Set input data
####

##  Get input parameter from file 
if arguments["--input_data"] == '':
    print("No input_data, exit")
    sys.exit()

exec(open(arguments["--input_data"]).read()) # Contient la variable insar

## Get others parameters
dem_path = arguments['--dem']
dest_path = arguments['--dest']
ext = arguments['--ext']

####
## 1: Read and prepare 2D and 1D input data
####

## Prepare ##
# get cube dimensions from one of the cubes 
# full extent to initially read the cube, keep this if no crop option
img_data = get_cube_dimension(insar[0][0])
nlines, ncols, n_img, nodata = img_data[0], img_data[1], img_data[2], img_data[3]

# set crop parameters if crop is input
if(arguments['--crop']):
    # (xmin, xmax, ymin, ymax)
    crop = tuple(map(int, arguments['--crop'].split(',')))

    nlines_crop = crop[1] - crop[0] + 1
    ncols_crop = crop[3] - crop[2] + 1

    print('Crop option is set. New cube dimensions are ({},{},{})'.format(nlines_crop, ncols_crop, n_img))
else:
    crop = None
    print('No crop option is set. Process the full extent of cube')

## Read ##
# read topo files

dem = read_tif(dem_path, crop)
rot, slope = compute_slope_aspect(dem_path, crop, min_slope=0)

# read satelite angle and dates

cubes = []
head = []
look = []
dates = []
for i in range(len(insar)):
    cubes.append(gdal.Open(insar[i][0]))
    dates.append(np.loadtxt(insar[i][1], comments='#', usecols=(3), unpack=True,dtype='f'))
    look.append(read_tif(insar[i][2], crop))
    head.append(read_tif(insar[i][3], crop))

print()

####
## 2: Plot input data
####

if plot:
    print('--')
    print('Plot is enabled: reading 3D rasters and plotting 2D and 3D rasters')
    print('--')
    
    ## Read cubes to plot them
    cubes_array = []
    for i in range(len(insar)):
      cubes_array.append(read_cube(insar[i][0], crop))  # Lire uniquement la dernière bande + utiliser le ds

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
    
    ## PLOT ASC and DSC data ##
    fig, axes = plt.subplots(len(insar), 3, figsize=(15, 10))  # 1 ligne, 3 colonnes

    for i in range(len(insar)):
        last_band_index = cubes_array[i].shape[2] - 1
        im0 = axes[i][0].imshow(cubes_array[i][:, :, last_band_index], vmin=np.nanpercentile(cubes_array[i][:, :, last_band_index],2), vmax=np.nanpercentile(cubes_array[i][:, :, last_band_index],98), cmap='rainbow', interpolation='nearest')
        axes[i][0].set_title(insar[i][0].split('/')[-1], pad=1)  # Ajoute de l'espace sous le titre
        fig.colorbar(im0, ax=axes[i][0], orientation="horizontal", fraction=0.046, pad=0.1)

        im1 = axes[i][1].imshow(look[i], cmap='YlGnBu')
        axes[i][1].set_title("Look angle", pad=1)  # Ajoute de l'espace sous le titre
        fig.colorbar(im1, ax=axes[i][1], orientation="horizontal", fraction=0.046, pad=0.1)

        im2 = axes[i][2].imshow(head[i], cmap='YlGnBu')
        axes[i][2].set_title("Heading angle", pad=1)  # Ajoute de l'espace sous le titre
        fig.colorbar(im2, ax=axes[i][2], orientation="horizontal", fraction=0.046, pad=0.1)

    plt.subplots_adjust(wspace=0.3, hspace=0.4)  # Ajuste wspace et hspace selon tes besoins
    plt.tight_layout()
    plt.show()
    
    del cubes_array
    print()

####
## 3: Define angles for projection
####

theta = []
phi = []

# convert head, look to angle phi, theta in rad
for i in range(len(insar)):
    print(insar[i][0].split('/')[-1])
    print('Average LOOK:{0:.5f}, HEADING:{1:.5f} angles'.format(np.nanmean(look[i]),np.nanmean(head[i])))
    theta.append(np.deg2rad(90. - look[i])) # theta: angle between LOS and horizontal (positive in anti-clockwise direction)
    phi.append(np.deg2rad(-90. - head[i])) #phi: horizontal angle between LOS and East (positive in anti-clockwise direction)
    print('Average THETA:{0:.5f}, PHI:{1:.5f} angles'.format(np.nanmean(np.rad2deg(theta[i])), np.nanmean(np.rad2deg(phi[i]))))
    print()

####
## 4:  Interpolate and inverse all pixel of all bands 
####

# Define the number of bands of the output cubes = cube with lower amount of bands
nimg_final = float('inf')  # Utilise l'infini pour trouver le minimum
index = -1

for i in range(len(cubes)):
    nimg = cubes[i].RasterCount  # Obtenir le nombre de bandes
    if nimg < nimg_final:
        nimg_final = nimg
        index = i

print(f"Each time series of each cube will be interpolated based on the master {insar[index][0].split('/')[-1]} with {nimg_final} dates")
print()


def process_chunk(slope, rot, theta, phi, data, dates, master):
    """
    Applies a method of linear interpolation on the time series based on a master,
    then applies a method of spatial decomposition on the data, date by date.

    :param slope, rot: topographic angle - 1D array (pixels)
    :param theta, phi: satellite angle - 2D array (track x pixels)
    :param data: contains the time series of the cubes - 3D array (track x pixels x time series)
    :param dates: contains the dates of the cubes - 2D array (tracks x dates)
    :param master: contains the index of the master dates
    :return: list of decomposed pixels to insert into the final products
    """
    plot_yes = False

    U0 = np.full((len(slope), len(data[master])), np.nan)
    U2 = np.full((len(slope), len(data[master])), np.nan)
    print(U0.shape)
    print(len(slope))
    print(len(data[master]))
    
    for j in range(len(slope)):
        if np.all(np.isnan(data[master][:, j])): # Explications
            continue

        all_nan = 0
        for k in range(master):
            if not np.all(np.isnan(data[k][:, j])): # Explications
                all_nan += 1
                break
        if not all_nan:
            for k in range(master +1, len(data)):
                if not np.all(np.isnan(data[k][:, j])):  # Explications
                    all_nan += 1
                    break

        if all_nan == 0:
            continue

        if slope[j] < np.deg2rad(2.):
            continue

        if (np.deg2rad(-120.) <= rot[j] <= np.deg2rad(-60.)) or (np.deg2rad(60.) <= rot[j] <= np.deg2rad(120.)):
            continue
       
        dates_interpolated, pixel_ts_interpolated = interpolate(dates, [data[i][:, j] for i in range(len(data))], master, plot=plot_yes)

        pixel_U0, pixel_U2 = invertion(dates_interpolated, pixel_ts_interpolated, [phi[i][j] for i in range(len(data))], [theta[i][j] for i in range(len(data))], slope[j], rot[j], plot=plot_yes)
       
        
        U0[j, :] = pixel_U0
        U2[j, :] = pixel_U2
   
    return U0, U2

# Define multiprocess parameters
if arguments["--proc"] == None:
    num_processes = 4
else :
    num_processes = int(arguments["--proc"])

if arguments["--block_size"] == None:
    block_size = 200
else :
    block_size = int(arguments["--block_size"])

# Define area
start_line = crop[0] if crop else 0
end_line = crop[1] if crop else nlines
block_height = end_line - start_line + 1

start_col = crop[2] if crop else 0
end_col = crop[3] if crop else ncols
block_width = end_col - start_col + 1

# Create output cube
U0 = np.full((block_height, block_width, nimg_final), np.nan)
U2 = np.full((block_height, block_width, nimg_final), np.nan)

with mp.Pool(processes=num_processes) as pool:
    for line in range(0, block_height, block_size):
        end_block_line = min(line + block_size, block_height)
        block_size_local = end_block_line - line  #taille (en lignes) du bloc analysée dans cette boucle
        if line%20 == 0:
            print(f'Start inverting line {line}')

        slope_block = np.array(slope[line:end_block_line, 0:end_col])
        #plt.imshow(slope_block)
        #plt.show()
        slope_block = np.reshape(slope_block, (block_width*block_size_local))
        slope_block = np.array_split(slope_block, num_processes, axis=0) 
        
        rot_block = np.array(rot[line:end_block_line, 0:end_col])
        #plt.imshow(rot_block)
        #plt.show()
        rot_block = np.reshape(rot_block, (block_width*block_size_local))
        rot_block = np.array_split(rot_block, num_processes, axis=0) 
        
        theta_block = []
        phi_block = []
        ts = []

        for i in range(len(cubes)):
            theta_i = np.array(theta[i][line:end_block_line, 0:end_col])
            phi_i = np.array(phi[i][line:end_block_line, 0:end_col])
            
            theta_i = np.reshape(theta_i, (block_width*block_size_local))
            phi_i = np.reshape(phi_i, (block_width*block_size_local))

            theta_block.append(np.array_split(theta_i, num_processes, axis=0))
            phi_block.append(np.array_split(phi_i, num_processes, axis=0))
            
            data_block = np.array([cubes[i].GetRasterBand(b+1).ReadAsArray(start_col, start_line + line, block_width, block_size_local) for b in range(cubes[i].RasterCount)]) 
            #plt.imshow(cubes[i].GetRasterBand(10).ReadAsArray(start_col, start_line+line, block_width, block_size_local))
            #plt.show()
            data_block = np.where(np.isin(data_block, [cubes[i].GetRasterBand(1).GetNoDataValue()]), np.nan, data_block)

            data_block = np.reshape(data_block, (cubes[i].RasterCount, block_width*block_size_local)) # flaten each band
            data_block = np.squeeze(data_block) # remove dimentions of size 1 if there is some
            data_block = np.array_split(data_block, num_processes, axis=1) # divide data into x 'num_processes' array
            ts.append(data_block)
        
        # Création des paramètres pour le processus
        parameters = []
        for i in range(num_processes):
            theta_input = []
            phi_input = []
            ts_input = []
            for j in range(len(cubes)):
                theta_input.append(theta_block[j][i])
                phi_input.append(phi_block[j][i])
                ts_input.append(ts[j][i])
            parameters.append((
                slope_block[i], 
                rot_block[i], 
                np.array(theta_input), 
                np.array(phi_input), 
                ts_input,
                dates,
                index
            ))

        print(
            f"slope_block[i].shape = {slope_block[i].shape if isinstance(slope_block[i], np.ndarray) else 'Not a NumPy array'}")
        print(
            f"rot_block[i].shape = {rot_block[i].shape if isinstance(rot_block[i], np.ndarray) else 'Not a NumPy array'}")
        print(f"theta_input shape = {np.array(theta_input).shape}")
        print(f"phi_input shape = {np.array(phi_input).shape}")
        print(f"ts_input shape = {len(ts_input)}, {len(ts_input[0])}, {len(ts_input[0][0])}")
        print(f"dates shape = {len(dates)}, {len(dates[0])} ")


        results = pool.starmap(process_chunk, parameters)

        U0[line:end_block_line, :, :] = np.concatenate([res[0] for res in results], axis=0).reshape((block_size_local, block_width, nimg_final)) 
        U2[line:end_block_line, :, :] = np.concatenate([res[1] for res in results], axis=0).reshape((block_size_local, block_width, nimg_final))


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

print(f"Time : {time.time() - t1}")

if(crop):
    U0_filename = 'depl_cumule_U0_{}_crop_{}_{}_{}_{}'.format(ext, crop[0], crop[1], crop[2], crop[3])
    U2_filename = 'depl_cumule_U2_{}_crop_{}_{}_{}_{}'.format(ext, crop[0], crop[1], crop[2], crop[3])
else:
    uslope_filename = 'depl_cumule_uslope_{}'.format(ext)
    uz_filename = 'depl_cumule_uz_{}'.format(ext)

if(crop):
    ds_ref = gdal.OpenEx(insar[0][2], allowed_drivers=['GTiff'])
    translated_raster = "ref_file.tif"
    gdal.Translate(translated_raster, ds_ref, srcWin=[start_col, start_line, block_width, block_height])

    ds_ref = gdal.OpenEx(translated_raster, allowed_drivers=['GTiff'])
    geotransform = ds_ref.GetGeoTransform()
    projection = ds_ref.GetProjection()

else:
    geotransform = cubes[0].GetGeoTransform()
    projection = cubes[0].GetProjection()

def save_cube_tif(dest_path, out_filename, maps, geotransform=None, projection=None):
    nrow, ncol, nimg = maps.shape

    print(maps.shape)
    driver = gdal.GetDriverByName("GTiff")
    out_path = os.path.join(dest_path, f"{out_filename}.tif")
    dataset = driver.Create(out_path, ncol, nrow, nimg, gdal.GDT_Float32)

    # Ajouter les métadonnées géospatiales si fournies
    if geotransform:
        dataset.SetGeoTransform(geotransform)
    if projection:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(projection)
        dataset.SetProjection(srs.ExportToWkt())

    # Écrire chaque bande dans le fichier
    for i in range(nimg):
        band = dataset.GetRasterBand(i + 1)
        band.WriteArray(maps[:, :, i])
        band.SetNoDataValue(np.nan)

    dataset.FlushCache()
    dataset = None  # Fermer le fichier
    print(f"Cube sauvegardé sous {out_path}")

save_cube_tif(dest_path, U0_filename, U0, geotransform=geotransform, projection=projection)
save_cube_tif(dest_path, U2_filename, U2, geotransform=geotransform, projection=projection)





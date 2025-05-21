#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
# Author        : Simon Daout / Hugo Watine
############################################

"""\
compute_modseas_geo.py
-------------
Compute the mod of the seasonality. Read output files from invers_dips2coef.py but georef. 

Usage: clean_phi-amp.py [--cube=<path> ] [--ampfile=<path>] \
[--crop=<values>] [--crop_shp=<path>] [--images_retenues=<path>]\
[--sigampfile=<path>] [--minamp=<value>] [--maxamp=<value>] [--name=<value>]


Options:
-h --help             Show this screen.
--images_retrenues    Path to dates file
--cube=<file>         Path to time series [default: depl_cumule_flat]
--ampfile=<file>      Amplitude map file [default: ampwt_coeff.tif]
--crop=<values>       Crop data [default: 0,nlines,0,ncol]
--crop_shp=<file>     Crop data with a shapefile [default: None]  
--sigampfile=<file>   Uncertainty amplitude map file [default: ampwt_sigcoeff.tif]
--minamp=<value>      Mask on minimum Amplitude in mm [default: 3.5]
--maxamp=<value>      Threshold Amplitude limit in mm [default: 10.]
--name=<value>        Output file name
"""

import numpy as np
from numpy.lib.stride_tricks import as_strided

import scipy.stats as stat
import scipy.optimize as opt
import scipy.linalg as lst
from scipy import ndimage
from osgeo import gdal, ogr
import sys

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import subprocess

#from Plot.plot_profile_tif import nodata_value

#np.warnings.filterwarnings('ignore')
try:
    from nsbas import docopt
except:
    import docopt

print()
print()
print('Author: Simon Daout')
print('Please cite:')
print('Daout, S., Dini, B., Haeberli, W., Doin, M. P., & Parsons, B. (2020). Ice loss in the Northeastern Tibetan Plateau permafrost as seen by 16 yr of ESA SAR missions. Earth and Planetary Science Letters, 545, 116404')
print()
print()

arguments = docopt.docopt(__doc__)
if arguments["--minamp"] ==  None:
    minamp = 2.5
else:
    minamp = float(arguments["--minamp"])
if arguments["--maxamp"] ==  None:
    maxamp = 10
else:
    maxamp = float(arguments["--maxamp"])
if arguments["--ampfile"] ==  None:
    arguments["--ampfile"] = 'ampwt_coeff.tif'
if arguments["--sigampfile"] ==  None:
    arguments["--sigampfile"] = 'ampwt_sigcoeff.tif'
if arguments["--cube"] ==  None:
    arguments["--cube"] = 'depl_cumule_flat'
if arguments["--crop_shp"] ==  None:
    shp = None
else:
    shp = arguments["--crop_shp"]
if arguments["--images_retenues"] == None:
    fimages = 'images_retenues'
else:
    fimages = arguments["--images_retenues"]

imref = 0
rad2mm = -4.4563

##################################
def supr_file(name):
    subprocess.call(f"rm {name}", shell=True)


def extract_data_in_shapefile_zone1(raster_path, shapefile_path):
    # Ouvrir le raster
    raster_ds = gdal.Open(raster_path)

    if raster_ds is None:
        print(f"Erreur: Impossible d'ouvrir le raster {raster_path}")
        return None

    nodata = raster_ds.GetRasterBand(1).GetNoDataValue()
    if nodata is None:
        print(f"Aucune valeur NoData définie pour {raster_path}. L'utilisation de 0 par défaut.")
        nodata = 0

    # Découper le raster en utilisant gdal.Warp
    new_name = f"{raster_path.split('.')[0]}_{shapefile_path.split('/')[-1].split('.')[0]}.tif"

    print(nodata)
    # Appliquer le découpage avec un shapefile
    gdal.Warp(new_name, raster_ds, cutlineDSName=shapefile_path, cropToCutline=True, dstNodata=np.nan,srcNodata=np.nan, resampleAlg='nearest')

    return new_name


def extract_data_in_shapefile_zone(raster_path, shapefile_path):
    # Ouvrir le raster
    raster_ds = gdal.Open(raster_path)
    if raster_ds is None:
        print(f"Erreur: Impossible d'ouvrir le raster {raster_path}")
        return None
    
    # Ouvrir le shapefile et obtenir la couche
    shapefile_ds = ogr.Open(shapefile_path)
    if shapefile_ds is None:
        print(f"Erreur: Impossible d'ouvrir le shapefile {shapefile_path}")
        return None
    
    layer = shapefile_ds.GetLayer()
    feature = layer.GetNextFeature()

    # Obtenir la géométrie du shapefile
    geom = feature.GetGeometryRef()

    # Définir le chemin de sortie pour le nouveau raster découpé
    new_name = f"{raster_path.split('.')[0]}_{shapefile_path.split('/')[-1].split('.')[0]}.tif"

    # Paramètres de découpe pour gdal.Warp
    warp_options = gdal.WarpOptions(
        format='GTiff',
        cutlineDSName=shapefile_path,
        cutlineWhere=None,
        cropToCutline=True,
        dstNodata=np.nan,
        outputBounds=None,
        width=0,
        height=0,
        srcSRS=None,
        dstSRS=None,
        multithread=False,
        resampleAlg='nearest',
        srcAlpha=False,
        dstAlpha=False,
        warpOptions=None,
        errorThreshold=None,
        creationOptions=None,
        callback=None,
        callback_data=None
    )

    # Découper le raster en utilisant gdal.Warp
    gdal.Warp(new_name, raster_ds, options=warp_options)

    return new_name

# lect cube
nb,idates,dates,base=np.loadtxt(fimages, comments='#', usecols=(0,1,3,5), unpack=True,dtype='i,i,f,f')
ds = gdal.Open(arguments["--cube"])

ncols = ds.RasterXSize
nlines = ds.RasterYSize
N = ds.RasterCount
print("Nombre de colonnes (ncols):", ncols)
print("Nombre de lignes (nlines):", nlines)
print("Nombre de bandes (N):", N)

if arguments["--crop"] ==  None:
    crop = [0,nlines,0,ncols]
else:
    crop = list(map(float,arguments["--crop"].replace(',',' ').split()))
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])
# print(ibeg,iend,jbeg,jend)

if shp != None:
    arguments["--ampfile"] = extract_data_in_shapefile_zone(arguments["--ampfile"], shp)
    #arguments["--sigampfile"]= extract_data_in_shapefile_zone(arguments["--sigampfile"], shp)##

# Open mapsi"
amp_map=gdal.Open(arguments["--ampfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]

#sigamp_map=gdal.Open(arguments["--sigampfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]

# if shp != None:
#     ibeg,iend,jbeg,jend = 0, len(amp_map), 0, len(amp_map[0])
#
#     print("Type of files that will be deleted:")
#     print(arguments["--ampfile"])
#     # Asking for confirmation
#     #confirm_delete = input("Do you want to delete these files? (Type 'Yes' to confirm) ")
#     confirm_delete = 'yes'
#     if confirm_delete.strip().lower() == 'yes':
#         supr_file(arguments["--ampfile"])
#         supr_file(arguments["--sigampfile"])
#
#         if arguments["--amp_increase_file"] != None:
#             supr_file(arguments["--amp_increase_file"])
#     else:
#         print("Operation cancelled. No files were deleted.")



# slope in %
#slope_map = slope_map*100

# Initialisation
dmod=[]
flat_disp_pos = []
flat_disp_neg = []
time_pos,time_neg = [],[]
amp_pos,amp_neg = [],[]

##################################

## function invers seasonality
def invers_seas(x,y,std):
    G=np.zeros((len(x),3))
    G[:,0]=np.cos(2*np.pi*(x))
    G[:,1]=np.sin(2*np.pi*(x))
    G[:,2]=np.ones(len(x))

    x0 = lst.lstsq(G,y)[0]
    _func = lambda x: np.sum(((np.dot(G,x)-y)/std)**2)
    _fprime = lambda x: 2*np.dot(G.T/std, (np.dot(G,x)-y)/std)
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=20000,full_output=True,iprint=0)[0] 

    Cd = np.diag(std_neg**2,k=0)

    Cov = np.linalg.inv(Cd)

    a,b = pars[0],pars[1]
    phi=np.arctan2(b,a)
    amp = np.sqrt(a**2+b**2)
    cova = np.linalg.inv(np.dot(G.T,np.dot(Cov,G)))
    siga,sigb,sigc = np.sqrt(np.diag(cova))
    sigamp = np.sqrt(siga**2+sigb**2)
    sigphi = (a*siga+b*sigb)/(a**2+b**2)

def invers_linseas(x, y, std):
    G = np.zeros((len(x), 3))
    G[:, 0] = np.cos(2 * np.pi * (x))
    G[:, 1] = np.sin(2 * np.pi * (x))
    G[:, 2] = np.ones(len(x))

    x0 = lst.lstsq(G, y)[0]
    _func = lambda x: np.sum(((np.dot(G, x) - y) / std) ** 2)
    _fprime = lambda x: 2 * np.dot(G.T / std, (np.dot(G, x) - y) / std)
    pars = opt.fmin_slsqp(_func, x0, fprime=_fprime, iter=20000, full_output=True, iprint=0)[0]

    Cd = np.diag(std_neg ** 2, k=0)

    Cov = np.linalg.inv(Cd)

    a, b = pars[0], pars[1]
    phi = np.arctan2(b, a)
    amp = np.sqrt(a ** 2 + b ** 2)
    cova = np.linalg.inv(np.dot(G.T, np.dot(Cov, G)))
    siga, sigb, sigc = np.sqrt(np.diag(cova))
    sigamp = np.sqrt(siga ** 2 + sigb ** 2)
    sigphi = (a * siga + b * sigb) / (a ** 2 + b ** 2)

    return pars, amp, phi, sigamp, sigphi

def seasonal(time,a,b,c):
    return a*np.cos(2*np.pi*time) + b*np.sin(2*np.pi*time) + c

def coef_determination(y, pred):
    u = ((y- pred)**2).sum()
    v = ((y- y.mean())**2).sum()
    return 1- u/v

##################################

# # plot
# fig_ampphi = plt.figure(0,figsize=(14,4))
#
# ax = fig_ampphi.add_subplot(1,3,1)
# im = ax.imshow(slope_map,interpolation='nearest', cmap=cm.Greys,vmax=np.nanpercentile(slope_map,90),vmin=np.nanpercentile(slope_map,10),zorder=1)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# plt.colorbar(im, cax=cax)
# ax.set_title('Slope',fontsize=12)
#
# ax = fig_ampphi.add_subplot(1,3,2)
# cmap = cm.rainbow
# cax = ax.imshow(slope_map,cmap=cm.Greys,vmax=np.nanpercentile(slope_map,92),vmin=np.nanpercentile(slope_map,8),interpolation='nearest', zorder=1)
# # no point to take negatif amplitude
# print(len(amp_map))
# print(len(np.ma.array(amp_map, mask=np.isnan(amp_map))))
# im = ax.imshow(np.ma.array(amp_map, mask=np.isnan(amp_map)),interpolation='nearest', cmap=cmap,vmin=minamp,\
#     vmax=np.nanpercentile(amp_map,98), alpha=0.5, zorder=2)
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# plt.colorbar(im, cax=cax)
# ax.set_title('Amplitude (mm)',fontsize=12)
#
# ax = fig_ampphi.add_subplot(1,3,3)
# cmap_phi = cm.tab20b
# cax = ax.imshow(slope_map,cmap=cm.Greys,zorder=1,vmax=np.nanpercentile(slope_map,92),interpolation='nearest',vmin=np.nanpercentile(slope_map,8))
# im = ax.imshow(np.ma.array(phi_map, mask=np.isnan(phi_map)),cmap=cmap_phi,vmin=0,vmax=2*np.pi,alpha=0.8,zorder=2, interpolation='nearest')
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# plt.colorbar(im, cax=cax)
# ax.set_title('Timing (rad)',fontsize=12)
# fig_ampphi.tight_layout()
# fig_ampphi.savefig('{}-modmaps.pdf'.format(arguments["--name"]), format='PDF',dpi=80)
#
# plt.show()
# sys.exit(0)

# Open cube

if shp != None:
    arguments["--cube"] = extract_data_in_shapefile_zone1(arguments["--cube"], shp)


maps=gdal.Open(arguments["--cube"])
maps = maps.ReadAsArray()[:, ibeg:iend, jbeg:jend]
maps= np.transpose(maps, (1, 2, 0))

print(maps.shape)
print(amp_map.shape)

last_band = maps[:, :, -1]  # Sélection de la dernière bande
plt.imshow(last_band, cmap='viridis',  vmin=np.nanpercentile(last_band, 2), vmax=np.nanpercentile(last_band, 98), interpolation='nearest')  # Choix du colormap 'viridis', ajustez si nécessaire
plt.colorbar(label="Valeurs")  # Ajouter une échelle de couleur
plt.title("Dernière bande du cube")
plt.xlabel("Colonnes")
plt.ylabel("Lignes")
plt.show()

plt.imshow(amp_map, cmap='viridis',  vmin=np.nanpercentile(amp_map, 2), vmax=np.nanpercentile(amp_map, 98), interpolation='nearest')  # Choix du colormap 'viridis', ajustez si nécessaire
plt.colorbar(label="Valeurs")  # Ajouter une échelle de couleur
plt.title("Amplitude")
plt.xlabel("Colonnes")
plt.ylabel("Lignes")
plt.show()




ibeg, iend, jbeg, jend = 0, maps.shape[0], 0, maps.shape[1]

#cube = np.fromfile(arguments["--cube"],dtype=np.float32)[:nlines*ncols*N]*rad2mm
#maps = cube.reshape((nlines,ncols,N))[ibeg:iend,jbeg:jend]

dmodt = np.fmod(dates,1)
for t in range((N)):
    dmod.append(dmodt[t])

print("iend-ibeg = ", iend-ibeg)
print("jend-jbeg = ", jend-jbeg)

from scipy.optimize import curve_fit
def model(t, a, b, c, phi):
    return a * t + b + c * np.cos(2 * np.pi * t + phi)

for i in range(iend-ibeg):
    for j in range(jend-jbeg):
        if not np.all(np.isnan(maps[i, j, :])) and amp_map[i, j] < maxamp and amp_map[i, j] > minamp:
            time_series = maps[i, j, :]
            dates_series = dates[:]

            # Enlever les NaN pour le fitting
            mask = ~np.isnan(time_series)
            filtered_dates = dates_series[mask]
            filtered_time_series = time_series[mask]

            # Ajustement du modèle complet
            try:
                # Estimations initiales pour les paramètres (a, b, c, phi)
                initial_guess = [0, np.mean(filtered_time_series),
                                 (np.max(filtered_time_series) - np.min(filtered_time_series)) / 2, 0]
                params, _ = curve_fit(model, filtered_dates, filtered_time_series, p0=initial_guess)
                a, b, c, phi = params

                fitted_model = model(filtered_dates, a, b, c, phi)

                # Calcul du terme linéaire (a * t + b) seulement
                linear_trend = a * filtered_dates + b

                # Soustraction du terme linéaire
                detrended_series = filtered_time_series - linear_trend

                # Normalisation par l'amplitude c du terme cosinus
                normalized_series = detrended_series / c

                # Affichage combiné des données brutes, du modèle ajusté et de la série normalisée
                plot = False
                if plot:
                    plt.figure(figsize=(12, 6))

                    # Plot des données brutes et du modèle ajusté
                    plt.subplot(1, 2, 1)
                    plt.plot(filtered_dates, filtered_time_series, 'bo', label="Données brutes")
                    plt.plot(filtered_dates, fitted_model, '-', label="Modèle ajusté")
                    plt.xlabel("Dates")
                    plt.ylabel("Valeur")
                    plt.title(f"Série temporelle et modèle pour le pixel ({i}, {j})")
                    plt.legend()

                    # Plot de la série normalisée
                    plt.subplot(1, 2, 2)
                    plt.plot(filtered_dates, normalized_series, 'bo', label="Série normalisée")
                    plt.xlabel("Dates")
                    plt.ylabel("Valeur normalisée")
                    plt.title(f"Série temporelle normalisée pour le pixel ({i}, {j})")
                    plt.legend()

                    plt.tight_layout()
                    plt.show()

                for t in range(len(normalized_series)):
                    flat_disp_pos.append(normalized_series[t])
                    time_pos.append(dmodt[t])
                    amp_pos.append(amp_map[i, j])

            except RuntimeError:
                print(f"L'ajustement du modèle a échoué pour le pixel ({i}, {j}).")

        #temp_pos = (maps[i,j,:] - lin_map[i,j]*(dates[:]-dates[imref]) - ref_map[i,j] - bperp_map[i,j]*(base[:]-base[imref]))/amp_map[i,j]


# wrap everything
time_pos = np.array(time_pos).flatten()
amp_pos = np.array(amp_pos).flatten()
#time_neg = np.array(time_neg).flatten()
#amp_neg = np.array(amp_neg).flatten()
flat_disp_pos = np.array(flat_disp_pos).flatten()
#flat_disp_neg = np.array(flat_disp_neg).flatten()
dmod=np.unique(np.array(dmod).flatten())

# check
#print(len(amp_pos), amp_pos)
#print(len(amp_neg), amp_neg)
#sys.exit()

mean_pos,std_pos, dmod_pos=[],[],[]
#mean_neg,std_neg, dmod_neg=[],[],[]

if flat_disp_pos.size == 0 or time_pos.size == 0 or dmod.size == 0:
    print("Un ou plusieurs tableaux sont vides pour les données positives")
else:
    for d in dmod:
        uu = np.flatnonzero(time_pos==d)
        if len(uu) > 0:
            if not np.isnan(flat_disp_pos[uu]).all():
                mean_pos.append(np.nanmedian(flat_disp_pos[uu]))
                std_pos.append(np.nanstd(flat_disp_pos[uu]))
                dmod_pos.append(d)
            else:
                print(f"All values are NaN for d={d}. Skipping.")
        else:
            print(f"No matches found for d={d} in time_neg. Skipping.")

# if flat_disp_neg.size == 0 or time_neg.size == 0 or dmod.size == 0:
#     print("Un ou plusieurs tableaux sont vides pour les données negatives")
# else:
#     for d in dmod:
#         uu = np.flatnonzero(time_neg==d)
#         if len(uu) > 0:
#             if not np.isnan(flat_disp_neg[uu]).all():
#                 mean_neg.append(np.nanmedian(flat_disp_neg[uu]))
#                 std_neg.append(np.nanstd(flat_disp_neg[uu]))
#                 dmod_neg.append(d)
#             else:
#                 print(f"All values are NaN for d={d}. Skipping.")
#         else:
#             print(f"No matches found for d={d} in time_neg. Skipping.")
mean_pos,std_pos,dmod_pos  = np.array(mean_pos),np.array(std_pos),np.array(dmod_pos)
#mean_neg,std_neg,dmod_neg  = np.array(mean_neg),np.array(std_neg),np.array(dmod_neg)

# plot slope postive
fig=plt.figure(1,figsize=(14,5))
ax=fig.add_subplot(1,2,1)
ax.plot(dmod_pos,mean_pos,'o',c='blue',ms=6.)
ax.errorbar(dmod_pos,mean_pos,yerr=std_pos,fmt='none',ecolor='blue',alpha=0.1)

#print(len(mean_pos),len(std_pos),len(dmod_pos))
#print(len(mean_neg),len(std_neg),len(dmod_neg))

try:
    pars, amp, phi, sigamp, sigphi = invers_seas(dmod_pos,mean_pos,std_pos) 
    t = np.arange(1,100)/100.
    r2 =  coef_determination(mean_pos, seasonal(dmod_pos,pars[0],pars[1],pars[2]))
    ax.plot(t,seasonal(t,pars[0],pars[1],pars[2]),'red',\
        lw=2,label='{:0.1f}+-{:0.1f} * cos(wt - {:0.1f}+-{:0.1f}) - R2 = {:0.1f}'.format(amp,sigamp,phi,sigphi,r2))
except Exception as e:
    print(f"Erreur lors du calcul avec invers_seas : {e}")
    pass




ax.set_ylim([-3,3])
ax.set_xlim([0,1])
ax.set_xticks(np.arange(0,1, 1./12))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.legend(loc='best')

# Adding vertical dashed lines at 0.360 and 0.690
ax.axvline(x=0.360, color='black', linestyle='--')
ax.axvline(x=0.690, color='black', linestyle='--')

# # plot slope negative
ax2=fig.add_subplot(1,2,2)
# ax2.plot(dmod_neg,mean_neg,'o',c='blue',ms=6.,label='{} < {}'.format(threshold,maxthreshold))
# ax2.errorbar(dmod_neg,mean_neg,yerr=std_neg,fmt='none',ecolor='blue',alpha=0.1)
#
# try:
#     pars, amp, phi, sigamp, sigphi = invers_seas(dmod_neg,mean_neg,std_neg)
#     t = np.arange(1,100)/100.
#     r2 =  coef_determination(mean_neg, seasonal(dmod_neg,pars[0],pars[1],pars[2]))
#     ax2.plot(t,seasonal(t,pars[0],pars[1],pars[2]),'red',\
#         lw=2,label='{:0.1f}+-{:0.1f} * cos(wt - {:0.1f}+-{:0.1f}) - R2 = {:0.1f}'.format(amp,sigamp,phi,sigphi,r2))
# except Exception as e:
#     print(f"Erreur lors du calcul avec invers_seas : {e}")
#     pass

ax2.set_ylim([-3,3])
ax2.set_xlim([0,1])
ax2.set_xticks(np.arange(0,1, 1./12))
ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
plt.legend(loc='best')

# Adding vertical dashed lines at 0.360 and 0.690
ax2.axvline(x=0.360, color='black', linestyle='--')
ax2.axvline(x=0.690, color='black', linestyle='--')

fig.tight_layout()
fig.savefig('{}-modseas.pdf'.format(arguments["--name"]), format='PDF',dpi=80)
plt.show()

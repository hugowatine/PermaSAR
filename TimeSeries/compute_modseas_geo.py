#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon Daout
############################################

"""\
compute_modseas_geo.py
-------------
Compute the mod of the seasonality. Read output files from invers_dips2coef.py but georef. 

Usage: clean_phi-amp.py [--cube=<path> ] [--ampfile=<path>] [--amp_increase_file=<path>] [--phifile=<path>] [--linfile=<path>] [--demerrfile=<path>] \
[--ref_file=<path>] [--slopefile=<path>] [--outampfile=<path>] [--outphifile=<path>] [--crop=<values>] [--crop_shp=<path>] [--lectfile=<path>] \
[--sigampfile=<path>] [--sigphifile=<path>] [--minamp=<value>] [--maxamp=<value>] [--maxslope=<value>] [--perc_sig=<value>] [--name=<value>] [--topofile=<path>] [--threshold=<amp/slope>] [--minelev=<value>] [--aspectfile=<path>] [--mask_aspect=<value>]


Options:
-h --help             Show this screen.
--cube=<file>         Path to time series [default: depl_cumule_flat]
--ampfile=<file>      Amplitude map file [default: ampwt_coeff.tif]
--amp_increase_file=<file>   Increase of amplitude map file, if you add it, the Amplitude used will be amp+3*ampi [default: None]
--phifile=<file>      phi map file [default: phiwt_coeff.tif]
--linfile=<file>      Linear map file [default: lin_coeff.tif]
--demerrfile=<file>   Path to the dem error file [default: corrdem_coeff.tif]
--ref_file=<file>     Path to the reference file [default: ref_coeff.tif]
--slopefile=<file>    Elevation slope [default: Slope.tif]
--lectfile=<file>     Path of the lect.in file for r4 format [default: lect_ts.in]
--crop=<values>       Crop data [default: 0,nlines,0,ncol]
--crop_shp=<file>     Crop data with a shapefile [default: None]  
--sigampfile=<file>   Uncertainty amplitude map file [default: ampwt_sigcoeff.tif]
--sigphifile=<file>   Uncertainty phi map file [default: phiwt_sigcoeff.tif]
--minamp=<value>      Mask on minimum Amplitude in mm [default: 3.5]
--maxamp=<value>      Threshold Amplitude limit in mm [default: 10.]
--minelev=<value>     Mask based on elevation [default:3500]
--maxslope=<value>    Threshold on Slope in % [default: 3.]
--perc_sig=<value>    Percentile uncertainty for map cleaning [default: 99.]
--name=<value>        Output file name
--topofile=<file>     Topographic file [default: dem]
--threshold=<amp/slope> Choose if you want to apply the threshold on the amplitude (ampmax) or the slope (ampslope) [default: amp]
--aspectfile=<file>   aspect map file [default: None]
--mask_aspect=<value> Mask on aspect range [default:0,360]
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
if arguments["--ampfile"] ==  None:
    arguments["--ampfile"] = 'ampwt_coeff.tif'
if arguments["--phifile"] ==  None:
    arguments["--phifile"] = 'phiwt_coeff.tif'
if arguments["--ref_file"] ==  None:
    arguments["--ref_file"] = 'ref_coeff.tif'
if arguments["--linfile"] ==  None:
    arguments["--linfile"] = 'lin_coeff.tif'
if arguments["--lectfile"] ==  None:
    arguments["--lectfile"] = "lect_ts.in"
if arguments["--sigampfile"] ==  None:
    arguments["--sigampfile"] = 'ampwt_sigcoeff.tif'
if arguments["--sigphifile"] ==  None:
    arguments["--sigphifile"] = 'phiwt_sigcoeff.tif'
if arguments["--cube"] ==  None:
    arguments["--cube"] = 'depl_cumule_flat'
if arguments["--minamp"] ==  None:
    minamp = 2.5
else:
    minamp = float(arguments["--minamp"])
if arguments["--maxamp"] ==  None:
    maxamp = 10
else:
    maxamp = float(arguments["--maxamp"])
if arguments["--maxslope"] ==  None:
    maxslope = 3
else:
    maxslope = float(arguments["--maxslope"])
if arguments["--minelev"] ==  None:
    minelev = -1
else:
    minelev = float(arguments["--minelev"])
if arguments["--topofile"] == None:
    demf = None
else:
    demf = arguments["--topofile"]
if arguments["--slopefile"] ==  None:
    slopef = None
else:
    slopef = arguments["--slopefile"]
if arguments["--threshold"] == None:
    threshold = "Amp"
    maxthreshold = maxamp
elif arguments["--threshold"] == "amp":
    threshold = "Amp"
    maxthreshold = maxamp
elif arguments["--threshold"] == "slope":
    threshold = "Slope"
    maxthreshold = maxslope
else:
    "threshold argument not recognise. Exit!"
    sys.exit()
if arguments["--crop_shp"] ==  None:
    shp = None
else:
    shp = arguments["--crop_shp"]

if arguments["--mask_aspect"] ==  None:
    mask_aspect = [0, 360]
else:
    mask_aspect = list(map(float,arguments["--mask_aspect"].replace(',',' ').split()))
    print(mask_aspect)
fimages='images_retenues'
imref = 0
rad2mm = -4.4563

##################################
def supr_file(name):
    subprocess.call(f"rm {name}", shell=True)

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
if not ds:
  print ('.hdr file time series cube {0}, not found, open {1}'.format(arguments["--cube"],arguments["--lectfile"]))
  ncols, nlines, N = list(map(int, open(arguments["--lectfile"]).readline().split(None, 3)[0:3]))
else:
  hdr = arguments["--cube"] + '.hdr'
  print ('Read ', hdr)
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
    arguments["--phifile"] = extract_data_in_shapefile_zone(arguments["--phifile"], shp)
    arguments["--linfile"]= extract_data_in_shapefile_zone(arguments["--linfile"], shp)
    arguments["--sigampfile"]= extract_data_in_shapefile_zone(arguments["--sigampfile"], shp)
    arguments["--sigphifile"]= extract_data_in_shapefile_zone(arguments["--sigphifile"], shp)
    arguments["--ref_file"]= extract_data_in_shapefile_zone(arguments["--ref_file"], shp)
    if arguments["--amp_increase_file"] != None:
        arguments["--amp_increase_file"] = extract_data_in_shapefile_zone(arguments["--amp_increase_file"], shp)

# Open mapsi"
amp_map=gdal.Open(arguments["--ampfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]
phi_map=gdal.Open(arguments["--phifile"]).ReadAsArray()[ibeg:iend, jbeg:jend]
lin_map=gdal.Open(arguments["--linfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]
sigamp_map=gdal.Open(arguments["--sigampfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]
sigphi_map=gdal.Open(arguments["--sigphifile"]).ReadAsArray()[ibeg:iend, jbeg:jend]
ref_map=gdal.Open(arguments["--ref_file"]).ReadAsArray()[ibeg:iend, jbeg:jend]
if arguments["--amp_increase_file"] != None:
    amp_increase_map=gdal.Open(arguments["--amp_increase_file"]).ReadAsArray()[ibeg:iend, jbeg:jend]
    amp_map = abs(amp_map + 3*amp_increase_map)

if shp != None:
    ibeg,iend,jbeg,jend = 0, len(amp_map), 0, len(amp_map[0])
    
    print("Type of files that will be deleted:")
    print(arguments["--ampfile"])
    # Asking for confirmation
    #confirm_delete = input("Do you want to delete these files? (Type 'Yes' to confirm) ")
    confirm_delete = 'yes'
    if confirm_delete.strip().lower() == 'yes':
        supr_file(arguments["--ampfile"])
        supr_file(arguments["--phifile"])
        supr_file(arguments["--linfile"])
        supr_file(arguments["--sigampfile"])
        supr_file(arguments["--sigphifile"])
        supr_file(arguments["--ref_file"])

        if arguments["--amp_increase_file"] != None:
            supr_file(arguments["--amp_increase_file"])
    else:
        print("Operation cancelled. No files were deleted.")

if arguments["--demerrfile"] is not None:
    if shp != None:
        arguments["--demerrfile"] = extract_data_in_shapefile_zone(arguments["--demerrfile"], shp)
        bperp_map = gdal.Open(arguments["--demerrfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]
        supr_file(arguments["--demerrfile"])
    else:
        bperp_map = gdal.Open(arguments["--demerrfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]
else:
    bperp_map = np.zeros((nlines,ncols))[ibeg:iend,jbeg:jend]

if demf ==  None:
    try:
        demf = 'dem'
        if shp != None:
            dem = extract_data_in_shapefile_zone(arguments["--dem"], shp)
            dem_map = gdal.Open(dem).ReadAsArray()[ibeg:iend, jbeg:jend]
            supr_file(dem)
        else:
            dem_map = gdal.Open('dem').ReadAsArray()[ibeg:iend, jbeg:jend]
    except:
        print('DEM file is not readible. Set elelvation to zeros.')
        dem_map = np.zeros((nlines,ncols))[ibeg:iend,jbeg:jend]
else:
    if shp != None:
        dem = extract_data_in_shapefile_zone(demf, shp)
        dem_map = gdal.Open(dem).ReadAsArray()[ibeg:iend, jbeg:jend]
        supr_file(dem)
    else:
        dem_map = gdal.Open(demf).ReadAsArray()[ibeg:iend, jbeg:jend]
print('Min Elev: {}, Max Elve: {}'.format(np.nanmin(dem_map), np.nanmax(dem_map)))

if slopef ==  None:
    print('Slope file is not readible. Set relief to one.')
    slope_map = np.ones((nlines,ncols))[ibeg:iend,jbeg:jend]
else:
    if shp != None:
        slopef_shp = extract_data_in_shapefile_zone(slopef, shp)
        slope_map = gdal.Open(slopef_shp).ReadAsArray()[ibeg:iend, jbeg:jend]
        supr_file(slopef_shp)
    else:
        slope_map = gdal.Open(slopef).ReadAsArray()[ibeg:iend, jbeg:jend]

# # clean
if arguments["--perc_sig"] ==  None:
    perc_sig = 100.
else:
    perc_sig = float(arguments["--perc_sig"])

# cleaning 

print(len(amp_map))
print(len(amp_map[0]))
print(len(slope_map))
print(len(slope_map[0]))
print(np.nanmax(slope_map))
print(np.nanmin(slope_map))

nan_mask_combined = np.isnan(phi_map) & np.isnan(amp_map) & np.isnan(lin_map)

sigamp_map[nan_mask_combined] = float('NaN')
sigphi_map[nan_mask_combined] = float('NaN')
phi_map[nan_mask_combined] = float('NaN')
amp_map[nan_mask_combined] = float('NaN')
lin_map[nan_mask_combined] = float('NaN')

phi_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = float('NaN')
amp_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = float('NaN')
lin_map[sigamp_map>np.nanpercentile(sigamp_map,perc_sig)] = float('NaN')

phi_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = float('NaN')
amp_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = float('NaN')
lin_map[sigphi_map>np.nanpercentile(sigphi_map,perc_sig)] = float('NaN')

# conversion
# Il n'est pas nécéssaire d'appliquer ces modifications car les données d'entrée ont déja été traitées
#amp_map,sigamp_map,sigphi_map = amp_map*abs(rad2mm),sigamp_map*abs(rad2mm),sigphi_map*abs(rad2mm)
#lin_map = lin_map*rad2mm
#ref_map = ref_map*rad2mm

phi_map[amp_map<float(minamp)] = float('NaN')
lin_map[amp_map<float(minamp)] = float('NaN')

phi_map[dem_map<minelev] = float('NaN')
amp_map[dem_map<minelev] = float('NaN')
lin_map[dem_map<minelev] = float('NaN')


# slope in %
#slope_map = slope_map*100

phi_map[slope_map<maxslope] = float('NaN')
amp_map[slope_map<maxslope] = float('NaN')
lin_map[slope_map<maxslope] = float('NaN')

if arguments["--aspectfile"] != None:
    if shp != None:
        aspectf_shp = extract_data_in_shapefile_zone(arguments["--aspectfile"], shp)
        aspect_map = gdal.Open(aspectf_shp).ReadAsArray()[ibeg:iend, jbeg:jend]
        supr_file(aspectf_shp)
    else:
        aspect_map = gdal.Open(arguments["--aspectfile"]).ReadAsArray()[ibeg:iend, jbeg:jend]
    
    mask = np.logical_or(aspect_map < mask_aspect[0], aspect_map > mask_aspect[1])
    phi_map[mask] = float('NaN')
    amp_map[mask] = float('NaN')
    lin_map[mask] = float('NaN')

# convert phi between 0 and 2pi
phi_map[phi_map<0] = phi_map[phi_map<0] + 2*np.pi

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

    return pars, amp, phi, sigamp, sigphi

def seasonal(time,a,b,c):
    return a*np.cos(2*np.pi*time) + b*np.sin(2*np.pi*time) + c

def coef_determination(y, pred):
    u = ((y- pred)**2).sum()
    v = ((y- y.mean())**2).sum()
    return 1- u/v

##################################

# plot
fig_ampphi = plt.figure(0,figsize=(14,4))

ax = fig_ampphi.add_subplot(1,3,1)
im = ax.imshow(slope_map,interpolation='nearest', cmap=cm.Greys,vmax=np.nanpercentile(slope_map,90),vmin=np.nanpercentile(slope_map,10),zorder=1)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Slope',fontsize=12)

ax = fig_ampphi.add_subplot(1,3,2)
cmap = cm.rainbow
cax = ax.imshow(slope_map,cmap=cm.Greys,vmax=np.nanpercentile(slope_map,92),vmin=np.nanpercentile(slope_map,8),interpolation='nearest', zorder=1)
# no point to take negatif amplitude
print(len(amp_map))
print(len(np.ma.array(amp_map, mask=np.isnan(amp_map))))
im = ax.imshow(np.ma.array(amp_map, mask=np.isnan(amp_map)),interpolation='nearest', cmap=cmap,vmin=minamp,\
    vmax=np.nanpercentile(amp_map,98), alpha=0.5, zorder=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Amplitude (mm)',fontsize=12)

ax = fig_ampphi.add_subplot(1,3,3)
cmap_phi = cm.tab20b
cax = ax.imshow(slope_map,cmap=cm.Greys,zorder=1,vmax=np.nanpercentile(slope_map,92),interpolation='nearest',vmin=np.nanpercentile(slope_map,8))
im = ax.imshow(np.ma.array(phi_map, mask=np.isnan(phi_map)),cmap=cmap_phi,vmin=0,vmax=2*np.pi,alpha=0.8,zorder=2, interpolation='nearest')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
ax.set_title('Timing (rad)',fontsize=12)
fig_ampphi.tight_layout()
fig_ampphi.savefig('{}-modmaps.pdf'.format(arguments["--name"]), format='PDF',dpi=80)

plt.show()
# sys.exit(0)

# Open cube

if shp != None:
    arguments["--cube"] = extract_data_in_shapefile_zone(arguments["--cube"], shp)

maps=gdal.Open(arguments["--cube"])
maps=maps.ReadAsArray()[:, ibeg:iend, jbeg:jend]
maps = np.transpose(maps, (1, 2, 0))

#cube = np.fromfile(arguments["--cube"],dtype=np.float32)[:nlines*ncols*N]*rad2mm
#maps = cube.reshape((nlines,ncols,N))[ibeg:iend,jbeg:jend]

dmodt = np.fmod(dates,1)
for t in range((N)):
    dmod.append(dmodt[t])

for i in range(iend-ibeg):
    for j in range(jend-jbeg):
        if threshold == "Amp":  
            if (~np.isnan(amp_map[i,j]) and ~np.isnan(phi_map[i,j]) and amp_map[i,j]>float(maxamp)):
                temp_pos = (maps[i,j,:] - lin_map[i,j]*(dates[:]-dates[imref]) - ref_map[i,j] \
                - bperp_map[i,j]*(base[:]-base[imref]))/amp_map[i,j]
                for t in range((N)):
                    flat_disp_pos.append(temp_pos[t])
                    time_pos.append(dmodt[t])
                    amp_pos.append(amp_map[i,j])
            elif (~np.isnan(amp_map[i,j]) and ~np.isnan(phi_map[i,j]) and amp_map[i,j]<float(maxamp)):
                temp_neg = (maps[i,j,:] - lin_map[i,j]*(dates[:]-dates[imref]) - ref_map[i,j] \
                - bperp_map[i,j]*(base[:]-base[imref]))/amp_map[i,j]
                for t in range((N)):
                    flat_disp_neg.append(temp_neg[t])
                    time_neg.append(dmodt[t])
                    amp_neg.append(amp_map[i,j])
      
        elif threshold == "Slope":  
            if (~np.isnan(amp_map[i,j]) and ~np.isnan(phi_map[i,j]) and slope_map[i,j]>maxslope):
                temp_pos = (maps[i,j,:] - lin_map[i,j]*(dates[:]-dates[imref]) - ref_map[i,j] \
                - bperp_map[i,j]*(base[:]-base[imref]))/amp_map[i,j]
                for t in range((N)):
                    flat_disp_pos.append(temp_pos[t])
                    time_pos.append(dmodt[t])
                    amp_pos.append(amp_map[i,j])
            elif (~np.isnan(amp_map[i,j]) and ~np.isnan(phi_map[i,j]) and slope_map[i,j]<maxslope):
                temp_neg = (maps[i,j,:] - lin_map[i,j]*(dates[:]-dates[imref]) - ref_map[i,j] \
                - bperp_map[i,j]*(base[:]-base[imref]))/amp_map[i,j]
                for t in range((N)):
                    flat_disp_neg.append(temp_neg[t])
                    time_neg.append(dmodt[t])
                    amp_neg.append(amp_map[i,j])

# wrap everything
time_pos = np.array(time_pos).flatten()
amp_pos = np.array(amp_pos).flatten()
time_neg = np.array(time_neg).flatten()
amp_neg = np.array(amp_neg).flatten()
flat_disp_pos = np.array(flat_disp_pos).flatten()
flat_disp_neg = np.array(flat_disp_neg).flatten()
dmod=np.unique(np.array(dmod).flatten())

# check
#print(len(amp_pos), amp_pos)
#print(len(amp_neg), amp_neg)
#sys.exit()

mean_pos,std_pos, dmod_pos=[],[],[]
mean_neg,std_neg, dmod_neg=[],[],[]

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

if flat_disp_neg.size == 0 or time_neg.size == 0 or dmod.size == 0:
    print("Un ou plusieurs tableaux sont vides pour les données negatives")
else:
    for d in dmod:
        uu = np.flatnonzero(time_neg==d)
        if len(uu) > 0:
            if not np.isnan(flat_disp_neg[uu]).all():
                mean_neg.append(np.nanmedian(flat_disp_neg[uu]))
                std_neg.append(np.nanstd(flat_disp_neg[uu]))
                dmod_neg.append(d)
            else:
                print(f"All values are NaN for d={d}. Skipping.")
        else:
            print(f"No matches found for d={d} in time_neg. Skipping.")
mean_pos,std_pos,dmod_pos  = np.array(mean_pos),np.array(std_pos),np.array(dmod_pos)
mean_neg,std_neg,dmod_neg  = np.array(mean_neg),np.array(std_neg),np.array(dmod_neg)

# plot slope postive
fig=plt.figure(1,figsize=(14,5))
ax=fig.add_subplot(1,2,1)
ax.plot(dmod_pos,mean_pos,'o',c='blue',ms=6.,label='{} > {}'.format(threshold,maxthreshold)) 
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

# plot slope negative
ax2=fig.add_subplot(1,2,2)
ax2.plot(dmod_neg,mean_neg,'o',c='blue',ms=6.,label='{} < {}'.format(threshold,maxthreshold))
ax2.errorbar(dmod_neg,mean_neg,yerr=std_neg,fmt='none',ecolor='blue',alpha=0.1)

try:
    pars, amp, phi, sigamp, sigphi = invers_seas(dmod_neg,mean_neg,std_neg) 
    t = np.arange(1,100)/100.
    r2 =  coef_determination(mean_neg, seasonal(dmod_neg,pars[0],pars[1],pars[2]))
    ax2.plot(t,seasonal(t,pars[0],pars[1],pars[2]),'red',\
        lw=2,label='{:0.1f}+-{:0.1f} * cos(wt - {:0.1f}+-{:0.1f}) - R2 = {:0.1f}'.format(amp,sigamp,phi,sigphi,r2))
except Exception as e:
    print(f"Erreur lors du calcul avec invers_seas : {e}")
    pass

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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
flatism_compute_mean_miscolsure.py
-------------
Compute misclosure of interferograms (in rad) averaged by bins of temporal baseline.

Usage: flatism_compute_mean_miscolsure.py --output=<file> [--list_pair=<path>]  [--weight] [--plot]
flatism_compute_mean_miscolsure.py -h | --help

Options:
-h --help           Show this screen
--list_pair PATH      Path to the list of interferogram col1=date1 col2=date2 and, if it exist, col3=weights [default:list_pair]
--output PATH         Path for the output figure .pdf and the raster .tiff
--weight              Use weighted mean using column 3 of list_pair
--plot                Do plot
"""

import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import os
import docopt
from datetime import datetime
gdal.UseExceptions()

arguments = docopt.docopt(__doc__)
use_weight = arguments["--weight"]

def open_gdal(file, band=1, supp_ndv=None, complex=False):
    """
    Use GDAL to open band as real value or complex interferogram.
    Returns (data, Xsize, Ysize)
    If complex=True: data = [amplitude, phase]
    """
 
    if not os.path.isfile(file):
        raise FileNotFoundError('File does not exists: {}'.format(file))
    ds = gdal.Open(file)

    #print('dims', ds.RasterXSize, ds.RasterYSize, ds.RasterCount)
    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    data = band.ReadAsArray()
    #Xsize = ds.RasterXSize
    #Ysize = ds.RasterYSize
    if ndv is not None and ndv != np.nan:
        data[data==ndv] = np.nan
    if supp_ndv is not None and supp_ndv != np.nan:
        data[data==supp_ndv] = np.nan
    return data

if arguments["--list_pair"] == None:
    list_pair = "list_pair"
else:
    list_pair=arguments["--list_pair"]

list_pair=np.loadtxt(list_pair, dtype=str)
dates1 = list_pair[:, 0]
dates2 = list_pair[:, 1]
weights = list_pair[:, 2].astype(float) if list_pair.shape[1] > 2 else None

def yyyymmdd_to_datetime(d):
    return datetime.strptime(d, "%Y%m%d")

dt1 = np.array([yyyymmdd_to_datetime(d) for d in dates1])
dt2 = np.array([yyyymmdd_to_datetime(d) for d in dates2])

baseline_days = np.array([(d2 - d1).days for d1, d2 in zip(dt1, dt2)])
unique_baselines = np.unique(baseline_days)

print("Baselines temporelles uniques (jours) :")
print(unique_baselines)

groups = {"G1_12_24":  [12, 24], "G2_36_60":  [36, 48, 60], "G3_72_120": [72, 84, 96, 120], "G4_324_360": [324, 336, 348, 360]}

group_results = {}
group_count = {}
global_sum = None
global_weight = 0.0

for group_name, baseline_list in groups.items():
    print(f"\n=== Processing {group_name} ===")
    rms_sum = None
    count = 0
    for b in baseline_list:
        indices = np.where(baseline_days == b)[0]

        for idx in indices:
            d1 = dates1[idx]
            d2 = dates2[idx]
            fname = f"RMSpixel_{d1}_{d2}"
            if not os.path.exists(fname):
                print(f"⚠️  Fichier manquant : {fname}")
                continue
            #print(f"   Loading {fname}")
            rms = open_gdal(fname)
            w = weights[idx] if (use_weight and weights is not None) else 1.0

            if rms_sum is None:
                rms_sum = rms.astype(np.float64)*w
            else:
                rms_sum += rms * w
            count += w

            if global_sum is None:
                global_sum = rms.astype(np.float64)*w
            else:
                global_sum += rms * w
            global_weight += w


    if count == 0:
        print(f"⚠️ Aucun fichier trouvé pour {group_name}")
        continue
    print(f"   Loading {count} data")
    rms_mean = rms_sum / count
    group_results[group_name] = rms_mean
    group_count[group_name] = count

if global_weight > 0:
    global_mean = global_sum / global_weight

group_results["GLOBAL"] = global_mean
group_count["GLOBAL"] = global_weight

#all_data = np.array([arr for arr in group_results.values()])
#vmin = min(arr.min() for arr in group_results.values())
#vmax = max(arr.max() for arr in group_results.values())
#max_abs = max(abs(vmin), abs(vmax))

if arguments["--plot"]:
    n = len(group_results)
    fig, axes = plt.subplots(2, n, figsize=(4*n, 8), gridspec_kw={'height_ratios':[3, 1]})

    for i, (group_name, arr) in enumerate(group_results.items()):
        # Carte RMS mean
        vmin, vmax = np.percentile(arr, [2, 98])
        max_abs = max(abs(vmin), abs(vmax))
        max_abs = 0.1
        im = axes[0, i].imshow(arr, cmap='RdBu', vmin=-max_abs, vmax=max_abs)
        axes[0, i].set_title(f"{group_name.replace('_',' ')} : {group_count[group_name]}", fontsize=12)
        axes[0, i].axis('off')
        
        # Colorbar pour chaque image
        cbar = fig.colorbar(im, ax=axes[0, i], fraction=0.046, pad=0.04)
        cbar.set_label("RMS mean")
        
        # Histogramme sous la carte
        arr_nonzero = arr[arr != 0.0]
        mean_val = arr_nonzero.mean()
        axes[1, i].hist(arr_nonzero, range=(-0.2, 0.2), bins=100, color='k', alpha=0.5, density=True)
        axes[1, i].axvline(mean_val, color='k', linestyle='--', linewidth=2)
        axes[1, i].text(0.95, 0.95, f"{mean_val:.3f}", color='k', fontweight='bold',
                        ha='right', va='top', transform=axes[1, i].transAxes)
        axes[1, i].set_xlabel("RMS mean")
        axes[1, i].set_ylabel("Normalized count")
        axes[1, i].grid(True, linestyle='--', alpha=0.5)

        plt.tight_layout()
        plt.savefig(str(arguments["--output"]+'.pdf'), dpi=300)
        plt.show()

ds_template = gdal.Open(fname)
driver = gdal.GetDriverByName("GTiff")
ds = driver.CreateCopy(str(arguments["--output"]+'.tiff'), ds_template)

band_obj = ds.GetRasterBand(1)
ndv = band_obj.GetNoDataValue()
global_mean = np.where(global_mean == 0.0, np.nan, global_mean)

if ndv is not None:
    data = np.where(np.isnan(global_mean), ndv, global_mean)
    band_obj.SetNoDataValue(ndv)
else:
    data = global_mean

band_obj.WriteArray(data)
band_obj.FlushCache()
ds = None
print("save:", str(arguments["--output"]+'.tiff'))







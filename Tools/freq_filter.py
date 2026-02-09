#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aapply a filter in the spectral domain to an input image 

Usage:
    freq_filter.py <input> <filter> [--filter_size=<value>] [--outfile=<outfile>] [--plot]
    freq_filter.py -h | --help

Options:
    -h --help           Display command details
    input               Input raster to be filter in gdal format
    filter              Low pass or hight pass filtering
    filter_size         size of the filter in pixel
"""

import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import docopt
import sys
import os

gdal.UseExceptions()

def arg2value(value, conversion=None, default=None):
    """Convert a string argument if exists otherwise use default value"""
    if value is None:
        return default
    elif conversion is not None:
        return conversion(value)
    return value

def open_gdal(file, band=1, supp_ndv=None):
    """
    Use GDAL to open band as real value or complex interferogram.
    Returns (data, Xsize, Ysize)
    If complex=True: data = [amplitude, phase]
    """
    print('-----')
    print(file)
    if not os.path.isfile(file):
        raise FileNotFoundError('File does not exists: {}'.format(file))
    ds = gdal.Open(file)

    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    data = band.ReadAsArray()
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize

    print('dims', Xsize, Ysize, ds.RasterCount, '/ nodata : ', ndv)
    print('-----\n')

    if ndv is not None and ndv != np.nan:
        data[data==ndv] = np.nan
    if supp_ndv is not None and supp_ndv != np.nan:
        data[data==supp_ndv] = np.nan

    return data, Xsize, Ysize

def plot_comparison(original, filtered, title_filter=""):
    """
    Plot original image, filtered image and residual
    """

    residual = original - filtered

    # Common color scale for original & filtered
    vmin = np.nanpercentile(original, 2)
    vmax = np.nanpercentile(original, 98)

    # Residual scale centered on 0
    rmax = np.nanpercentile(np.abs(residual), 98)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    im0 = axes[0].imshow(original, cmap='RdBu_r', vmin=vmin, vmax=vmax)
    axes[0].set_title("Original")
    plt.colorbar(im0, ax=axes[0], shrink=0.8)

    vmin = np.nanpercentile(filtered, 2)
    vmax = np.nanpercentile(filtered, 98)
    im1 = axes[1].imshow(filtered, cmap='RdBu_r', vmin=vmin, vmax=vmax)
    axes[1].set_title(f"Filtered ({title_filter})")
    plt.colorbar(im1, ax=axes[1], shrink=0.8)

    #vmin = np.nanpercentile(residual, 2)
    #vmax = np.nanpercentile(residual, 98)
    im2 = axes[2].imshow(residual, cmap='RdBu_r', vmin=-rmax, vmax=rmax)
    #im2 = axes[2].imshow(residual, cmap='RdBu_r', vmin=vmin, vmax=vmax)
    axes[2].set_title("Original − Filtered")
    plt.colorbar(im2, ax=axes[2], shrink=0.8)

    for ax in axes:
        ax.axis("off")

    plt.show()

def save_gdal(outfile, data, template=None, band=1, alt_ndv=None, roi_pac=True):
    """
    Create a new raster with same dims as the template or overwrite if exists.
    If data contains complex values (complex dtype), will write magnitude/phase splitted is not implemented.
    We assume template is a complex ROI_PAC-like stack when writing complex arrays.
    """

    final_outfile = outfile

    # Create dataset
    ds_template = gdal.Open(template)
    driver = gdal.GetDriverByName("ROI_PAC" if roi_pac else "GTiff")
    ds = driver.CreateCopy(final_outfile, ds_template)

    band_obj = ds.GetRasterBand(band)
    ndv = alt_ndv if alt_ndv is not None else band_obj.GetNoDataValue()
    if ndv is not None:
        band_obj.SetNoDataValue(ndv)
        data = np.where(np.isnan(data), ndv, data)

    band_obj.WriteArray(data)
    ds.FlushCache()

    print("save:", final_outfile)

def filter_freq(array, filter_type, filter_size, plot='no'):
    """
    Apply a low-pass or high-pass frequency filter to a 2D array

    Parameters
    ----------
    array : 2D numpy array Input image
    filter_type : str 'LP' or 'HP'
    filter_size : int Radius of the filter in pixels (frequency domain)

    Returns
    -------
    filt_array : 2D numpy array  Filtered image
    """

    if filter_type not in ("LP", "HP"):
        raise ValueError("filter_type must be 'LP' or 'HP'")
    if filter_size is None:
        raise ValueError("filter_size must be provided")

    nan_mask = np.isnan(array)
    # Replace NaNs by 0 for FFT
    arr = np.copy(array)
    arr[nan_mask] = 0.0
    # FFT
    fft = np.fft.fft2(arr) # On passe en freq
    fft_shift = np.fft.fftshift(fft) # On met les basses frequences au centre 

    ny, nx = arr.shape
    cy, cx = ny // 2, nx // 2 # Centre de l'image en frequence

    # Frequency grid
    y, x = np.ogrid[:ny, :nx]
    r2 = (y - cy)**2 + (x - cx)**2 # Distance au carrée du centre de l'image

    sigma = filter_size
    gauss = np.exp(-r2 / (2 * sigma**2)) #Filtre gaussien si on est au centre (r2 = 0) on à 1 et plus on s'éloigne moins c'est fort
    # Sigma donnc l'écart-type du filtre gaussien, 1 ecart type couvre 68.2% de la données / 2 sigma = 95.4% / 3 sigma = 99.7%
    # Moi après ~70 pixel je n'ai plus de relation de corélation (semis-vario) donc je souhaite avoir dans la cloche toutes mes valeurs
    # => 3 sigma = 70 donc 1 sigma = 70/3 = ~23
    # => 2 sigma = 70 donc 1 sigma = 70/2 = 35

    if filter_type == "HP":
        gauss = 1.0 - gauss

    fft_filt = fft_shift * gauss

    filt = np.fft.ifft2(np.fft.ifftshift(fft_filt)).real
    filt[nan_mask] = np.nan

    if plot=='yes':
        spec_before = np.log10(np.abs(fft_shift) + 1)
        spec_after  = np.log10(np.abs(fft_filt) + 1)

        vmax = np.nanpercentile(spec_before, 99)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

        im0_vmax= np.nanpercentile(spec_before, 98)
        im0_vmin= np.nanpercentile(spec_before, 50)
        im0 = axes[0].imshow(spec_before, cmap='gray', vmin=im0_vmin, vmax=np.max(spec_before))
        axes[0].set_title("Spectrum before filtering")
        plt.colorbar(im0, ax=axes[0], shrink=0.8)

        im1 = axes[1].imshow(gauss, cmap='gray')
        axes[1].set_title("Frequency mask")
        plt.colorbar(im1, ax=axes[1], shrink=0.8)

        im2 = axes[2].imshow(spec_after, cmap='gray', vmin=0, vmax=vmax)
        axes[2].set_title("Spectrum after filtering")
        plt.colorbar(im2, ax=axes[2], shrink=0.8)

        for ax in axes:
            ax.axis("off")

        plt.show()

    return filt

if __name__ == '__main__':
    arguments = docopt.docopt(__doc__)
    #print(arguments)


    plot = 'yes' if arguments["--plot"] else 'no'
    filter_size = arg2value(arguments["--filter_size"], float, 20)
    outfile = arg2value(arguments["--outfile"], str, 'filt_raster.tif')
    filter_type = arguments["<filter>"]

    if filter_size <= 0:
        raise ValueError("filter_size must be > 0")

    array, Xsize, Ysize = open_gdal(arguments["<input>"], band=1, supp_ndv=0.0)

    filt_array = filter_freq(array, filter_type, filter_size, plot)

    save_gdal(outfile, filt_array, template=arguments["<input>"], band=1, roi_pac=False)
    residual_path = outfile.split('.')[0] + '_residual.' + outfile.split('.')[1]
    save_gdal(residual_path, array - filt_array, template=arguments["<input>"], band=1, roi_pac=False)

    if plot=='yes':
        plot_comparison(array, filt_array, title_filter=f"{filter_type} / size={filter_size}px")


# r4totif.py --infile=slope_corr_var --outfile=slope_corr_var.tif --lectfile=lect.in --ref_file=../NSBAS_TS-PKG_S1_TIBET-HIM-A158CENTRE-VV-2014-2022_IW123_2014-10-16_2022-05-31/CNES_MV-LOS_geo_8rlks.tiff
# ~/PermaSAR/Tools/freq_filter.py slope_corr_var.tif LP --plot
    






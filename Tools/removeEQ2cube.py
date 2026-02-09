#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove earthquake pattern from cumulative displacement data.

Usage:
    removeEQ2cube.py <cube> <YYYYMMDD> <x> <y> [--HPfilter_size=<value>] [--LPfilter_size=<value>] [--gaussian_mask=<outfile>] [--outfile=<outfile>] [--plot] [-f]
    removeEQ2cube.py -h | --help

Options:
    -h --help               Display command details
    cube                    Input cube to be cleaned
    YYYYMMDD                Date of the earthquake
    x                       position of the earthquake (col)
    y                       posiition of the earthquake (line)
    --LPfilter_size=<value>   Size of the filter in pixels (default: 100)
    --HPfilter_size=<value>   Size of the filter in pixels (default: 9)
    --gaussian_mask=<value> Sigma of the gaussian mask apply in the process (default: 300)
    --outfile=<outfile>     Output file name [default: cube_eqclean.tif]
    --plot                  Display diagnostic plots
    -f                      Overright
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from scipy import ndimage
from osgeo import gdal
import docopt
import math

import os

gdal.UseExceptions()

def save_gdal(outfile, data, template=None, band=1, alt_ndv=None):
    """
    Create a new raster with same dims as the template or overwrite if exists.
    If data contains complex values (complex dtype), will write magnitude/phase splitted is not implemented.
    We assume template is a complex ROI_PAC-like stack when writing complex arrays.
    """

    final_outfile = outfile

    # Create dataset
    ds_template = gdal.Open(template)
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.CreateCopy(final_outfile, ds_template)

    band_obj = ds.GetRasterBand(band)
    ndv = alt_ndv if alt_ndv is not None else band_obj.GetNoDataValue()
    if ndv is not None:
        band_obj.SetNoDataValue(ndv)
        data = np.where(np.isnan(data), ndv, data)

    band_obj.WriteArray(data)
    ds.FlushCache()

    print("save:", final_outfile)

def arg2value(value, conversion=None, default=None):
    """Convert a string argument if it exists, otherwise use default value."""
    if value is None:
        return default
    if conversion is not None:
        return conversion(value)
    return value


def open_gdal(file, band=1, supp_ndv=None):
    """
    Open a raster band using GDAL.

    Returns
    -------
    data : ndarray
    Xsize, Ysize : int
    """
    print('-----')
    print(file)

    if not os.path.isfile(file):
        raise FileNotFoundError(f'File does not exist: {file}')

    ds = gdal.Open(file)
    band = ds.GetRasterBand(band)

    ndv = band.GetNoDataValue()
    data = band.ReadAsArray().astype(float)

    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize

    print('dims:', Xsize, Ysize, ds.RasterCount, '/ nodata:', ndv)
    print('-----\n')

    if ndv is not None:
        data[data == ndv] = np.nan
    if supp_ndv is not None:
        data[data == supp_ndv] = np.nan

    return data, Xsize, Ysize

def plot_raster(raster, plot=False, title='', show='yes', crop='yes', freq='no', sigma=None, h=400):
    """Construct the raster display"""
    global x
    global y 
    if plot: 
        try:
            from matplotlib.colors import LinearSegmentedColormap
            cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
            cpt = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
            cpt = cpt.reversed()
        except:
            cpt=cm.rainbow

        if crop == 'yes':
            ny, nx = raster.shape
            x0, x1 = max(0, x-h), min(nx, x+h)
            y0, y1 = max(0, y-h), min(ny, y+h)
            raster = raster[y0:y1, x0:x1]
            xp, yp = x-x0, y-y0
        else:
            xp, yp = x, y

        if freq == 'yes':
            nan_mask = np.isnan(raster)
            raster_fft = raster.copy()
            raster_fft[nan_mask] = 0.0
            fft = np.fft.fft2(raster_fft)
            fft_shift = np.fft.fftshift(fft)
            raster = np.log10(np.abs(fft_shift) + 1)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1,1,1)
        vmin = np.nanpercentile(raster, 2)
        vmax = np.nanpercentile(raster, 98)
        if freq == 'yes':
            hax = ax.imshow(raster, 'grey', interpolation='nearest', vmin=np.min(raster), vmax=np.max(raster))
        else:
            vmin = -10
            vmax = 10
            hax = ax.imshow(raster, cpt, interpolation='nearest', vmin=vmin, vmax=vmax)
            ax.plot(xp, yp, marker='+', color='k', markersize=15, mew=2)
        ax.set_title(title)
        divider = make_axes_locatable(ax)
        c = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(hax, cax=c)

        if sigma!= None:
            circle_plot = plt.Circle((xp, yp), sigma, color='r', fill=False, linewidth=2)
            ax.add_patch(circle_plot)
        
        plt.tight_layout()
        if show=='yes':
            plt.show()
    else:
        pass

def extract_bands(output, ds, bands):
    if os.path.isfile(output):
        print(f"→ {output} already exists, skipping extraction")
        return
    
    print(f"→ extracting band {bands} → {output}")
    gdal.Translate(output, ds, bandList=bands, creationOptions=["TILED=YES", "COMPRESS=DEFLATE"])


def generate_eq_model(cube_path, date_eq, plot=False):
    """Placeholder for earthquake model generation."""
    
    ds = gdal.Open(cube_path)
    md = ds.GetMetadata()

    # {date: band_number}
    date_to_band = {
        int(v): int(k.split("_")[1])
        for k, v in md.items()
        if k.startswith("Band_")
    }

    sorted_dates = sorted(date_to_band.keys())
    after_date_eq = next((d for d in sorted_dates if d > date_eq), None)


    if after_date_eq is None:
        print(f"No date found after {date_eq}")
    else:
        print("Next date :", after_date_eq)
        print("Corresponding band :", date_to_band[after_date_eq])
        band_after = date_to_band[after_date_eq]
        band_before = band_after -1

    intermediate_output_after = cube_path.split('.')[0] + f'_band{band_after}_eq{date_eq}.' + cube_path.split('.')[1]
    intermediate_output_before = cube_path.split('.')[0] + f'_band{band_before}_eq{date_eq}.' + cube_path.split('.')[1]
    extract_bands(intermediate_output_before, ds, [band_before])
    extract_bands(intermediate_output_after, ds, [band_after])

    data_b, _, _ = open_gdal(intermediate_output_before)
    data_a, _, _ = open_gdal(intermediate_output_after)

    diff = data_a - data_b

    plot_raster(data_a, title=f'after eq - {date_eq}', plot=plot, show='no')
    plot_raster(data_b, title=f'before eq - {date_eq}', plot=plot, show='no')
    plot_raster(diff, title=f'Diff - {date_eq}', plot=plot, show='no')

    return diff, data_b, data_a

def compute_sigma_map(shape, ref, sigma_min, sigma_max, dmax):
    ny, nx = shape
    y_px, x_px = np.indices((ny, nx))
    y0, x0 = ref

    dist = np.sqrt((y_px - y0)**2 + (x_px - x0)**2)

    # Exemple : sigma croît linéairement avec la distance
    sigma = sigma_min + (sigma_max - sigma_min) * np.clip(dist / dmax, 0, 1)

    return sigma

def gaussian_lp_variable_sigma(mf, sigma_map, n_bins=10):
    mf_out = np.zeros_like(mf)
    nan_mask = np.isnan(mf)

    # Préparation NaNs
    mf_filled = mf.copy()
    mf_filled[nan_mask] = 0.0
    ones = np.ones_like(mf)
    ones[nan_mask] = 0.0

    # Binning des sigmas
    sigmas = np.linspace(sigma_map.min(), sigma_map.max(), n_bins)

    for i in range(len(sigmas) - 1):
        s0, s1 = sigmas[i], sigmas[i+1]
        mask = (sigma_map >= s0) & (sigma_map < s1)

        if not np.any(mask):
            continue

        sigma_eff = 0.5 * (s0 + s1)

        val_f = ndimage.gaussian_filter(mf_filled, sigma_eff)
        one_f = ndimage.gaussian_filter(ones, sigma_eff)

        local_lp = val_f / one_f
        mf_out[mask] = local_lp[mask]
        print(i, s0, s1, np.sum(mask))


    mf_out[nan_mask] = np.nan
    return mf_out


def filt_eq_model(mf, filter_size, plot, method='freqLP'):
    global x
    global y 
    if method == 'HP' or method == 'LP':
        m_filter_vals = np.copy(mf)
        m_filter_vals[np.isnan(mf)] = 0.
        m_lp_vals = ndimage.gaussian_filter(m_filter_vals, int(filter_size))

        m_filter_ones = 0*np.copy(mf)+1
        m_filter_ones[np.isnan(mf)] = 0.
        m_lp_ones = ndimage.gaussian_filter(m_filter_ones, int(filter_size))

        mf_lp = m_lp_vals/m_lp_ones

        if method == 'HP':
            mf_final = mf - mf_lp
        else:
            mf_final = mf_lp

        mf_final[np.isnan(mf)] = np.nan

    if method == 'freqLP' or method== 'freqHP' :
        nan_mask = np.isnan(mf)
    
        arr = np.copy(mf)
        arr[nan_mask] = 0.0 # Replace NaNs by 0 for FFT

        fft = np.fft.fft2(arr) # On passe en freq
        fft_shift = np.fft.fftshift(fft) # On met les basses frequences au centre 
        ny, nx = arr.shape
        cy, cx = ny // 2, nx // 2 # Centre de l'image en frequence
        # Frequency grid
        y_px, x_px = np.ogrid[:ny, :nx]
        r2 = (y_px - cy)**2 + (x_px - cx)**2 # Distance au carrée du centre de l'image
        sigma = filter_size
        gauss = np.exp(-r2 / (2 * sigma**2))
        
        if method == "freqHP":
            gauss = 1.0 - gauss

        fft_filt = fft_shift * gauss

        mf_final = np.fft.ifft2(np.fft.ifftshift(fft_filt)).real
        mf_final[nan_mask] = np.nan

    if method == 'LP_var':
        ref = (x, y)
        sigma_min=3
        sigma_max=50
        dmax=200
        sigma_map = compute_sigma_map(mf.shape, ref=ref, sigma_min=sigma_min, sigma_max=sigma_max, dmax=dmax)
        print("sigma_map:", sigma_map.min(), sigma_map.max())
        mf_final = gaussian_lp_variable_sigma(mf, sigma_map, n_bins=12)

    #plot_raster(mf_final, title=f'filter {method} - size {filter_size}', plot=plot, show='no')
    return mf_final

def apply_gaussiankernel(raster, sigma, order=4):
    global x
    global y

    ny, nx = raster.shape
    y_px, x_px = np.indices((ny, nx))

    r2 = (y_px - y)**2 + (x_px - x)**2

    # Super-gaussien
    gaussian_mask = np.exp(-(r2 / (2 * sigma**2))**order)

    # Normalisation
    gaussian_mask /= gaussian_mask.max()

    result = raster * gaussian_mask

    # Rayon à 0.9 (approximation utile)
    r_09 = sigma * (2 * np.log(1 / 0.9))**(1 / (2 * order))

    return result, gaussian_mask, r_09

def remove_model_from_cube(cube, model, band):
    """Placeholder for model removal from cube."""
    ds = gdal.Open(file)
    band = ds.GetRasterBand(band)

    pass

if __name__ == '__main__':

    arguments = docopt.docopt(__doc__)

    plot = arguments["--plot"]
    x = int(arguments["<x>"])
    y = int(arguments["<y>"])
    cube_path = arguments["<cube>"]
    date_eq = int(arguments["<YYYYMMDD>"])
    LPfilter_size = arg2value(arguments["--LPfilter_size"], int, 0)
    HPfilter_size = arg2value(arguments["--HPfilter_size"], int, 0)
    gaussian_mask = arg2value(arguments["--gaussian_mask"], int, 200)
    outfile = arg2value(arguments["--outfile"], str, 'cube_eqclean.tif')

    row_model, before_model, after_model = generate_eq_model(cube_path, date_eq, plot)
    #plot_raster(row_model, title=f'freq', plot=plot, show='no', freq='yes', crop='no')

    print('--- High Pass ---')

    if HPfilter_size == 0:
        print("choose the good size")
        list_HP_size = [1, 2, 4, 6, 8, 10, 12, 14, 16]
        result = []
        std = []
        for elt in list_HP_size:
            raster_output = filt_eq_model(row_model, elt, plot=plot, method='freqHP')
            plot_raster(raster_output, title=f'model HP filt {elt}', plot=plot, show='no', crop='yes', h=500)
            result.append(raster_output)
        plt.show()
        Hp_size = int(input("Best HP filter size: "))
        i = list_HP_size.index(Hp_size)
        filt_HP_model = result[i]
    else:
        print(f"HP filtering wiht sigma = {HPfilter_size}")
        filt_HP_model = filt_eq_model(row_model, HPfilter_size, plot=plot, method='freqHP')
    
    

    print('')
    print('--- Low Pass ---')

    if LPfilter_size == 0:
        list_LP_size = [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300]
        result = []
        std = []
        for elt in list_LP_size:
            raster_output = filt_eq_model(filt_HP_model, elt, plot=plot, method='freqLP')
            plot_raster(raster_output, title=f'model LP filt {elt}', plot=plot, show='no', crop='yes', h=500)
            result.append(raster_output)
        plt.show()
        Lp_size = int(input("Best LP filter size: "))
        i = list_LP_size.index(Lp_size)
        filt_LP_model = result[i]
    else:
        print(f"LP filtering wiht sigma = {LPfilter_size}")
        filt_LP_model = filt_eq_model(filt_HP_model, LPfilter_size, plot=plot, method='freqLP')

    print('')
    print('--- Mask with gaussian kernel ---')
    print('')

    happy_gaussian_kernel = 'no'
    while happy_gaussian_kernel == 'no':
        sigma = gaussian_mask
        final_model, kernel, r_90 = apply_gaussiankernel(filt_LP_model, sigma=sigma)
        plot_raster(filt_LP_model, title=f'model LP filt {LPfilter_size}', plot=plot, show='no', crop='yes', sigma=r_90, h=700)
        plot_raster(final_model, title=f'final model', plot=plot, show='yes', crop='yes', sigma=r_90, h=700)
        plot_raster(kernel, title=f'kernel', plot=plot, show='yes', crop='yes', sigma=r_90, h=700)
        answer = input("Happy with the mask ? sigma={sigma} (y/n): ")
        if answer == 'y':
            happy_gaussian_kernel = 'yes'
        else : 
            gaussian_mask = int(input("New value for sigma gaussian mask (int): "))

    # Save model, ne pas faire tout ça quand on a le model, tej le model si force 
    # Ensuite retirer le model de toutes les bandes après le eq
    # Faudra check ce qu'on fait en regardant une ts
            






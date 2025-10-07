#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model coefficient estimation and removal for wrapped interferograms

Usage:
    flatten_stack.py <phase> <phase_filt> <model> [--outfile=<outfile>] [--nreg=<float>] [--thresh_amp=<float>] [--thresh_model=<float>] [--thresh_std_model=<float>] [--thresh_min_pixel=<int>] [--cyclmax=<float>]
    flatten_stack.py -h | --help
    flatten_stack.py estimate <phase_filt> <model> --outfile=<outfile> [--outfile=<outfile>] [--nreg=<nreg>] [--thresh_amp=<float>] [--thresh_model=<float>] [--thresh_std_model=<float>] [--thresh_min_pixel=<int>] [--cyclmax=<float>]
    flatten_stack.py add <unwrapped> <model> --outfile=<outfile> --coeff=<coeff>
    flatten_stack.py remove <phase> <model> --outfile=<outfile> --coeff=<coeff>

Options:
    -h --help           Display command details
    phase               None filtered interferogram (.int)       
    phase_filt          Filtered interferogram (.int)
    model               Model for the propotionality estimation (.r4)
    outfile             Outfile name
    nreg                Number of region in range, put 1 for 1 global region (default: 28)
    thresh_amp          Minimal value of amplitude for pixel selection (default: 0.1)
    thresh_model        Minimal value of model for pixel selection (default: 0.4)
    thresh_std_model    Minimal standard deviation within a window for estimation (default: 0.3)
    thresh_min_pixel    Minimal percent of pixel per windows for estimation (default: 10)
    cyclmax             Maximum number of phase cycle (default: )
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from osgeo import gdal
from numba import jit
import docopt
import sys
import os
from scipy.optimize import minimize


def open_gdal(file, band=1, supp_ndv=None, complex=False):
    """
    Use GDAL to open band as real value
    """
    print('-----')
    print(file)
    if not os.path.isfile(file):
        raise FileNotFoundError('File does not exists: {}'.format(file))
    ds = gdal.Open(file)
    # ds = gdal.OpenEx(file, allowed_drivers = ["ROI_PAC"])
    print('ds', ds)
    print('dims', ds.RasterXSize, ds.RasterYSize, ds.RasterCount)
    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    data = band.ReadAsArray()
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize
    if complex:
        amp = np.absolute(data)
        phi = np.angle(data)
        data = [amp, phi]
        ndv = 0.0
        if ndv is not None and ndv != np.nan:
            amp[amp==ndv] = np.nan
            phi[phi==ndv] = np.nan
        if supp_ndv is not None and supp_ndv != np.nan:
            amp[amp==supp_ndv] = np.nan
            phi[phi==supp_ndv] = np.nan
    else:
        if ndv is not None and ndv != np.nan:
            data[data==ndv] = np.nan
        if supp_ndv is not None and supp_ndv != np.nan:
            data[data==supp_ndv] = np.nan
    return data, Xsize, Ysize


def open_complex(file, band=1, supp_ndv=None):
    """
    Use GDAL to open band as complex value
    """
    if not os.path.isfile(file):
        raise FileNotFoundError('File does not exists: {}'.format(file))
    ds = gdal.Open(file)
    band = ds.GetRasterBand(band)
    data = band.ReadAsArray()
    amp = np.absolute(data)
    phi = np.angle(data)
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize
    ndv = 0.0
    if ndv is not None and ndv != np.nan:
        amp[amp==ndv] = np.nan
        phi[phi==ndv] = np.nan
    if supp_ndv is not None and supp_ndv != np.nan:
        amp[amp==supp_ndv] = np.nan
        phi[phi==supp_ndv] = np.nan
    return [amp, phi], Xsize, Ysize

def save_gdal(outfile, template, data, band=1, alt_ndv=None):
    """
    Create a new raster with same dims as the template
    """
    ds_template = gdal.Open(template)
    ds = gdal.GetDriverByName('ROI_PAC').CreateCopy(outfile, ds_template)
    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    if alt_ndv is not None:
        ndv = alt_ndv
        band.SetNoDataValue(ndv)
    if ndv is not None:
        data[data == np.nan] = ndv
    band.WriteArray(data)
    ds.FlushCache()


def former_flatten_stack(phase, phase_filt, model, outfile, **kwargs):
    """
    Remove a pattern from a wrapped signal, as in flatten_stack.f from NSBAS
    """
    coeff = estimate_coeff(phase_filt, model, **kwargs)
    model_removed = remove_model(phase, model, coeff)
    save_gdal(outfile, phase, model_removed)


def estimate_coeff(phase_filt, model, **kwargs) -> float:
    """
    Estimate the proportionnality coefficient between the model and the wrapped phase values
    minimizing the complex product (maximizing the coherence)
    """
    amp_filt = np.copy(phase_filt[0])
    phi_filt = np.copy(phase_filt[1])
    
    index = np.isnan(amp_filt)| np.isnan(phi_filt) | np.isnan(model) | (amp_filt<kwargs['thresh_amp']) | (model<kwargs['thresh_model'])

    amp_filt[index] = np.nan
    phi_filt[index] = np.nan
    model[index] = np.nan

    # fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    # im0 = axes[0].imshow(phi_filt, cmap='rainbow')
    # axes[0].set_title('Phase filtrée')
    # fig.colorbar(im0, ax=axes[0])
    # im1 = axes[1].imshow(amp_filt, cmap='gray_r')
    # axes[1].set_title('Amplitude filtrée')
    # fig.colorbar(im1, ax=axes[1])
    # im2 = axes[2].imshow(model, cmap='rainbow')
    # axes[2].set_title('Modèle')
    # fig.colorbar(im2, ax=axes[2])

    # plt.show()

    ## Creation of sub-area
    nreg = kwargs["nreg"]

    if nreg != 1:
        amp_filt = reg_pad(amp_filt, nreg)
        phi_filt = reg_pad(phi_filt, nreg)
        model = reg_pad(model, nreg)

    # Initialisation
    cyclmax = 2
    npascycl = 30
    #Dstack = 10
    ismax = int(cyclmax * npascycl)
    #fac = 2 * np.pi * cyclmax / (Dstack * ismax)
    #amax = -np.inf
    #best_slope = 0
    list_coeff = np.full(shape=(model.shape[0], model.shape[1]), fill_value=np.nan)


    print(model.shape)
    for i in range(model.shape[0]):     
        for j in range(model.shape[1]):
            # for each region
            phi_filt_region = phi_filt[i,j]
            model_region = model[i,j]

            if np.count_nonzero(~np.isnan(phi_filt_region)) < kwargs['thresh_min_pixel']:
                # print("skip not enough px", i, j)
                continue
            elif np.nanstd(model_region) < kwargs['thresh_std_model']:
                # print("skip not enough std", i, j)
                continue
            # print("launch lsqtsq ", i, j)
            #coeff = least_squares(minimize_complex_product, [0.0], bounds=[[-ismax], [ismax]], args=[phi_filt_region, model_region])
            #print(coeff.x)
            #list_coeff[i,j] = coeff.x

            #res = minimize(lambda c: -np.abs(np.nanmean(np.exp(1j*phi_filt_region) * np.exp(-1j*(np.mod(c*model_region + np.pi, 2*np.pi) - np.pi)))),
            #   x0=0.0, bounds=[(-ismax, ismax)])
            res = minimize(lambda c: 1/(np.abs(np.nanmean(np.exp(1j*phi_filt_region) * np.exp(-1j*(c*model_region))))),
               x0=0.0, bounds=[(-ismax, ismax)])
            coeff = res.x[0]
            print(coeff)
            list_coeff[i,j] = coeff


            #res = minimize(lambda params: -np.abs(np.nanmean(np.exp(1j * phi_filt_region) * np.exp(-1j*(np.mod(params[0]*model_region + params[1] + np.pi, 2*np.pi) - np.pi)))),
            #    x0=[0.0, 0.0],  bounds=[(-ismax, ismax), (-np.pi, np.pi)])
            #res = minimize(lambda params: -np.abs(np.nanmean(np.exp(1j * (phi_filt_region - (params[0]*model_region + params[1]))))),
            #    x0=[0.0, 0.0],  bounds=[(-ismax, ismax), (-np.pi, np.pi)])
            #coeff, offset = res.x
            #print(coeff, offset)


            if True:
                plt.figure(figsize=(5, 3))
                plt.scatter(np.angle(np.exp(1j*model_region)), np.angle(np.exp(1j*phi_filt_region)), s=10, c='k', alpha=0.6, label="Data")
                #plt.scatter(np.angle(np.exp(1j*phi_filt_region)) - np.angle(np.exp(1j*coeff*model_region)), np.angle(np.exp(1j*phi_filt_region)), s=10, c='r', alpha=0.6, label="Data")
                x_fit = np.linspace(np.nanmin(np.angle(np.exp(1j*model_region))), np.nanmax(np.angle(np.exp(1j*model_region))), 100)
                y_fit = np.angle(np.exp(1j*coeff*x_fit))

                x = np.angle(np.exp(1j*model_region))
                y = np.angle(np.exp(1j*phi_filt_region))

                y_fit2 = np.angle(np.exp(1j*coeff*model_region))
                cst =  np.nanmean(y - y_fit2)
                print(cst)

                #y_fit = np.mod(coeff*x_fit + offset + np.pi, 2*np.pi) - np.pi
                #plt.plot(x_fit, y_fit, 'r-', lw=2, label=f"Fit: c={coeff:.3f}, off={offset:.2f} rad")
                #plt.plot(x, y_fit2, 'r-', lw=2, label=f"Fit: coeff={list_coeff[i,j]:.3f}")
                plt.plot(x_fit, y_fit + cst, 'r-', lw=2, label=f"Fit: coeff={list_coeff[i,j]:.3f}")
                #plt.xlim(-4, 4)
                #plt.ylim(-4, 4)
                plt.xlabel("Model")
                plt.ylabel("Phase filtered")
                plt.legend()
                plt.title(f"Window ({i}, {j})")
                plt.grid(True, ls='--', alpha=0.4)
                plt.tight_layout()
                plt.show()
    
    print(np.nanmedian(list_coeff))
    print("all")
    
    plt.imshow(list_coeff, cmap='RdBu_r')
    plt.colorbar()
    plt.show()
    return np.median(list_coeff)
    pass


#@jit
def minimize_complex_product(x, *args):
    coeff = x[0]
    wrapped = args[0]
    model = args[1]
    # need to be between -pi and pi TODO    
    #return np.absolute(np.nanmean(np.exp(1j * wrapped) * np.exp(-1j * coeff * model)))
    return - np.absolute(np.nanmean(np.exp(1j * wrapped) * np.exp(-1j * coeff * (np.mod(model + np.pi, 2 * np.pi) - np.pi))))
    return np.absolute(np.nanmean(np.exp(1j * wrapped) * np.exp(-1j * coeff * np.mod(model, 2 * np.pi))))

    #phase_model = np.mod(coeff * model, 2 * np.pi)
    #phase_model[phase_model > np.pi] -= 2 * np.pi
    #return np.nanmean(wrapped-coeff * np.mod(model, 2 * np.pi))
    #return np.absolute(np.nanmean(np.exp(1j * wrapped) * np.exp(-1j * coeff * model)))
    #return np.nanmean(np.exp(1j * wrapped) * np.exp(-1j * coeff * np.mod(model, 2 * np.pi)))
    # return np.nanmean(np.exp(1j * wrapped) * np.exp(-1j * phase_model))


def reg_pad(array, block_nb_x):
    """
    Regionalize an array from a given block number in the x dimension.
    The squared block size is determined from the x dimension, by padding the array
    to ensure nx % block_size = 0.
    The number of blocks in the y dimension is determined by conserving the block size 
    constant between x and y axis.
    """
    ny, nx = array.shape
    # taille de bloc à nombre de blocs fixé
    dxy = nx // block_nb_x + 1
    pad_x = dxy * block_nb_x - nx

    # Nombre de blocs à taille de bloc fixée
    nby = ny // dxy + 1
    pad_y = dxy * nby - ny
    
    # Pad to ensure that it can perfectly be divided into blocks
    pad_array = np.pad(array, ((0, pad_y), (0, pad_x)), constant_values=np.nan)
    # Cut the different blocks into different sub-arrays
    y_blocks = np.array(np.array_split(pad_array, nby, axis=0))
    blocks = np.array([np.array_split(b, nreg, axis=1) for b in y_blocks])
    
    return blocks


@jit
def remove_model(phase, model, coeff):
    """
    Remove a model from wrapped values with a given proportionnality coefficient
    """
    return phase - np.mod(coeff* model, 2 * np.pi)


@jit
def add_model(amplitude, model, coeff):
    """
    Add a model to values given a proportionnality coefficient
    """
    return amplitude + coeff * model


def arg2value(value, conversion=None, default=None):
    """Convert a string argument if exists otherwise use default value"""
    if value is None:
        return default
    elif conversion is not None:
        return conversion(value)
    return value


if __name__ == '__main__':
    arguments = docopt.docopt(__doc__)
    print(arguments)
    
    # Additionnal parameters
    nreg = arg2value(arguments["--nreg"], int, 28)
    thresh_amp = arg2value(arguments["--thresh_amp"], float, 0.1)
    thresh_model = arg2value(arguments["--thresh_model"], float, 0.4)
    thresh_std_model = arg2value(arguments["--thresh_std_model"], float, 0.3)
    thresh_min_pixel = arg2value(arguments["--thresh_min_pixel"], float, 10) # nb px (TODO en %)
    cyclmax = arg2value(arguments["--cyclmax"], float, 0)

    if arguments["estimate"]:
        # Find the optimal coefficient of proportionality between the phase and the model
        phase, pXsize, pYsize = open_gdal(arguments["<phase>"], band=1, complex=True)
        phase_filt, pfXsize, pfYsize = open_gdal(arguments["<phase_filt>"], band=1, complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])
        if not (pXsize == pfXsize == mXsize and pYsize == pfYsize == mYsize):
            sys.exit("Error: input rasters (phase, phase_filt, model) do not have the same dimensions.")
        coeff = estimate_coeff(phase_filt, model, nreg=nreg, thresh_amp=thresh_amp, thresh_model=thresh_model, 
                       thresh_std_model=thresh_std_model, thresh_min_pixel=thresh_min_pixel, cyclmax=cyclmax)
        print(coeff)
    elif arguments["add"]:
        # Add the model to an unwrapped file
        unwrapped, unwXsize, unwYsize = open_gdal(arguments["<unwrapped>"], complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])
        if not (unwXsize == mXsize and unwYsize == mYsize):
            sys.exit("Error: input rasters (unw, model) do not have the same dimensions.")
        coeff = float(arguments["<coeff>"])
        outfile = arguments["--outfile"]
        unw_plus_model = add_model(unwrapped, model, coeff)
        save_gdal(outfile, unwrapped, unw_plus_model)
    elif arguments["remove"]:
        # Remove the model from wrapped signal
        phase, pXsize, pYsize = open_gdal(arguments["<phase>"], band=1, complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])
        if not (pXsize == mXsize and pYsize == mYsize):
            sys.exit("Error: input rasters (phase, model) do not have the same dimensions.")
        
        coeff = float(arguments["<coeff>"])
        outfile = arguments["--outfile"]
        phase_minus_model = remove_model(phase, model, coeff)
        save_gdal(outfile, phase, phase_minus_model, band=2)
    else:
        # Remove a pattern from a wrapped signal, as in flatten_stack.f from NSBAS
        phase, pXsize, pYsize = open_gdal(arguments["<phase>"], complex=True)
        # print(phase[1].dtype)
        # print('MIN/MAX', np.nanmin(phase[1]), np.nanmax(phase[1]))
        # plt.imshow(phase[1], cmap='rainbow')
        # plt.show()
        phase_filt, pfXsize, pfYsize = open_gdal(arguments["<phase_filt>"], complex=True)
        model, mXsize, mYsize = open_gdal(arguments["<model>"])
        # plt.imshow(model, cmap='Greys_r')
        # plt.show()
        if not (pXsize == pfXsize == mXsize and pYsize == pfYsize == mYsize):
            sys.exit("Error: input rasters (phase, phase_filt, model) do not have the same dimensions.")
        outfile = arguments["--outfile"]
        phase_minus_model = former_flatten_stack(phase, phase_filt, model, 
                                outfile, nreg=nreg, thresh_amp=thresh_amp, thresh_model=thresh_model, 
                                thresh_std_model=thresh_std_model, thresh_min_pixel=thresh_min_pixel, cyclmax=cyclmax)


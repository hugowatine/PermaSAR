#!/usr/bin/env python3

"""
vert-vs-topo.py
-------------------
Reads in raster files - vertical and topo
Plots topo on x axis vs vertical vel on y axis, pixel by pixel
Inputs must be the same size and the same number of pixels

.. Author:
    Roberta Wilkinson, University of Oxford; April 2023.
    Repurposes code from Simon Daout's plot_raster.py.
.. Last Modified:
    10th May 2023
Usage:
    vert-vs-topo.py --vertical=<path> --topo=<path> --outfile=<path> [--ymin=<value>] [--ymax=<value>] [--tp=<value>] [<ibeg>] [<iend>] [<jbeg>] [<jend>]

Options:
    -h --help             Show this screen
    --vertical=PATH       Path to vertical velocity file.
    --topo=PATH           Path to topography file
    --ibeg=<value>        Line number (smallest y-value) bounded the cutting zone [default: 0]
    --iend=<value>        Line number (biggest y-value) bounded the cutting zone
                          [default: nlines]
    --jbeg=<value>        Column number (smallest x-value) bounded the cutting zone
                          [default: 0]
    --jend=<value>        Column number (biggest x-value) bounded the cutting zone 
                          [default: ncols]
    --outfile=PATH        Output file name or path
    --ymin=<value>        Minimum value to display on vert vs topo plot (must also specify ymax)
    --ymax=<value>        Maximum value to display on vert vs topo plot (must also specify ymin)
    --tp=<value>          Set transparency for vert vs topo plot [default: 0.01]
"""

import os, sys, logging, glob
import docopt
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from osgeo.gdal import GA_ReadOnly
from numpy.lib.stride_tricks import as_strided
import matplotlib.cm as cm
from pylab import setp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle


#####################################################################################
# READ INPUT PARAM
#####################################################################################

# read input parameters and arguments
arguments = docopt.docopt(__doc__)

if arguments["--vertical"] != None:
    vertical = arguments["--vertical"]

if arguments["--topo"] != None:
    topo = arguments["--topo"]

if arguments["--outfile"] != None:
    outfile = arguments["--outfile"]

if arguments["--tp"] != None:
    tp = float(arguments["--tp"])
else:
    tp = 0.01

yflag = False
if arguments["--ymin"] != None and arguments["--ymax"] != None:
    ymin = float(arguments["--ymin"])
    ymax = float(arguments["--ymax"])
    yflag = True
elif arguments["--ymin"] != None and arguments["--ymax"] == None:
    print('WARNING: Specified ymix will not be used as ymax was not specified')
elif arguments["--ymin"] == None and arguments["--ymax"] != None:
    print('WARNING: Specified ymax will not be used as ymin was not specified')

sformat = 'GTIFF'
band = 1
if sformat == "GTIFF":
    #### Vertical velocity
    ds = gdal.Open(vertical, gdal.GA_ReadOnly)
    phase_band = ds.GetRasterBand(band)
    # Attributes
    print("Vertical velocity file")
    print("> Driver:   ", ds.GetDriver().ShortName)
    print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(phase_band.DataType))
    nlines, ncols = ds.RasterYSize, ds.RasterXSize

    #### Topography
    dstopo = gdal.Open(topo, gdal.GA_ReadOnly)
    phase_bandtopo = dstopo.GetRasterBand(band)
    # Attributes
    print("Topography file")
    print("> Driver:   ", dstopo.GetDriver().ShortName)
    print("> Size:     ", dstopo.RasterXSize,'x',ds.RasterYSize,'x',dstopo.RasterCount)
    print("> Datatype: ", gdal.GetDataTypeName(phase_bandtopo.DataType))
    nlinestopo, ncolstopo = dstopo.RasterYSize, dstopo.RasterXSize

if arguments["<ibeg>"] ==  None:
     ibeg = 0
else:
    ibeg = int(arguments["<ibeg>"])
if arguments["<iend>"] ==  None:
    iend = nlines
else:
    iend = int(arguments["<iend>"])
if arguments["<jbeg>"] ==  None:
    jbeg = 0
else:
    jbeg = int(arguments["<jbeg>"])
if arguments["<jend>"] ==  None:
    jend = ncols
else:
    jend = int(arguments["<jend>"])
crop = [0,jend,0,iend]

print('Using vertical velocity file:', vertical)
print('Using topography file:', topo)



## fonction ##

def fit_linear_model(x, data, RMS):

    G = np.vstack([np.ones_like(x)]).T
    Cd = np.diag(RMS ** 2)  # Matrice de covariance basée sur RMS

    # Calcul du paramètre de régression a et b en minimisant les moindres carrés pondérés
    params = np.dot(np.linalg.inv(np.dot(np.dot(G.T, np.linalg.inv(Cd)), G)),
                    np.dot(np.dot(G.T, np.linalg.inv(Cd)), data))

    return params

#####################################################################################
# CHECK THE FILES EXIST
#####################################################################################
try:
    topoopen = gdal.Open(topo)
except:
    print('The topo file does not exist.')
    sys.exit(0)

try:
    vertopen = gdal.Open(vertical)
except:
    print('The vertical file does not exist.')
    sys.exit(0)

#####################################################################################
# PLOT TIFF MAPS
#####################################################################################
try:
    from matplotlib.colors import LinearSegmentedColormap
    cm_locs = os.environ["PYGDALSAR"] + '/contrib/python/colormaps/'
    cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
    cmap = cmap.reversed()
    cmaptopo = cm.gray
except:
    cmap=cm.rainbow

# Initialize a matplotlib figure
fig = plt.figure(1,figsize=(15,8))

rad2mm = 1

####### Plot the vertical velocity - Big context map
axcontext = fig.add_subplot(1,3,1)
# Put the values of the vertical velocity array into phi
phi = phase_band.ReadAsArray(0, 0,
       ds.RasterXSize, ds.RasterYSize,
       ds.RasterXSize, ds.RasterYSize)

vmaxcontext = np.nanpercentile(phi,98)
vmincontext = np.nanpercentile(phi,2)

# replace 0 by nan
try:
    phi[phi==0] = float('NaN')
except:
    pass
masked_arraycontext = np.ma.array(phi, mask=np.isnan(phi))

caxcontext = axcontext.imshow(masked_arraycontext, cmap, interpolation='nearest',vmax=vmaxcontext,vmin=vmincontext)
#cax = ax.imshow(masked_array, cmap, interpolation='none',vmax=vmax,vmin=vmin)
divider = make_axes_locatable(axcontext)
ccontext = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(caxcontext, cax=ccontext)

axcontext.set_title('Vertical velocity', fontsize=14)

# Add the subplot outline
deli=abs(iend-ibeg)
delj=abs(jend-jbeg)
axcontext.add_patch(Rectangle((jbeg, ibeg), delj, deli, fc='none',color="black",linewidth=2))

####### Plot the vertical velocity - cropped map
ax = fig.add_subplot(1,3,2)
# Trim the maps
cutphi = as_strided(phi[ibeg:iend,jbeg:jend])*rad2mm

vmax = np.nanpercentile(cutphi,98)
vmin = np.nanpercentile(cutphi,2)

# replace 0 by nan
try:
    cutphi[cutphi==0] = float('NaN')
except:
    pass
masked_array = np.ma.array(cutphi, mask=np.isnan(cutphi))

cax = ax.imshow(masked_array, cmap, interpolation='nearest',vmax=vmax,vmin=vmin)
#cax = ax.imshow(masked_array, cmap, interpolation='none',vmax=vmax,vmin=vmin)
divider = make_axes_locatable(ax)
c = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cax, cax=c)

ax.set_title('Vertical velocity (cropped)', fontsize=14)

###### Plot the topography
axtopo = fig.add_subplot(1,3,3)

# Put the values of the vertical velocity array into phi
phitopo = phase_bandtopo.ReadAsArray(0, 0,
       dstopo.RasterXSize, dstopo.RasterYSize,
       dstopo.RasterXSize, dstopo.RasterYSize)

cutphitopo = as_strided(phitopo[ibeg:iend,jbeg:jend])*rad2mm

vmaxtopo = np.nanpercentile(cutphitopo,98)
vmintopo = np.nanpercentile(cutphitopo,2)

# replace 0 by nan
try:
    cutphitopo[cutphitopo==0] = float('NaN')
except:
    pass
masked_arraytopo = np.ma.array(cutphitopo, mask=np.isnan(cutphitopo))

caxtopo = axtopo.imshow(masked_arraytopo, cmaptopo, interpolation='nearest',vmax=vmaxtopo,vmin=vmintopo)
#cax = ax.imshow(masked_array, cmap, interpolation='none',vmax=vmax,vmin=vmin)
divider = make_axes_locatable(axtopo)
ctopo = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(caxtopo, cax=ctopo)

axtopo.set_title('Topography (cropped)', fontsize=14)

############################

#fig.canvas.set_window_title('Vertical velocity and topo')

try:
    del ds
    del dstopo
except:
    pass

# Display the data
plt.tight_layout()
# ax.set_rasterized(True)
try:
    fig.savefig('{}_map.png'.format(outfile), format='PNG',dpi=180)
except:
    pass

#####################################################################################
# PLOT THE SCATTER PLOTS
#####################################################################################

############ Plot vertical against topo for each pixel

# Flatten the data into 1D arrays
xtopo = np.ravel(cutphitopo)
yvertvel = np.ravel(cutphi)

# Plot with ylims if specified
if yflag:
    fig2 = plt.figure(2,figsize=(12,8))
    axcompare = fig2.add_subplot(1,1,1)
    axcompare.scatter(xtopo,yvertvel,s=0.5,alpha=tp,color='black')
    if yflag:
        axcompare.set_ylim(ymin,ymax)
    axcompare.set_xlabel('z (m)', fontsize=14)
    axcompare.set_ylabel('LOS velocity (mm/yr)', fontsize=14)
    axcompare.set_title(outfile + ' with imposed y-axis limits', fontsize=14)

    #fig2.canvas.set_window_title('Cropped data + ylims')

    # Display the data
    plt.tight_layout()
    # ax.set_rasterized(True)
    try:
        fig2.savefig('{}_crop_ylims_plot.png'.format(outfile), format='PNG',dpi=180)
    except:
        pass

# Plot the cropped data on new plot without ylims
fig3 = plt.figure(3,figsize=(12,8))
axcomparebig = fig3.add_subplot(1,1,1)
axcomparebig.scatter(xtopo,yvertvel,s=0.5,alpha=tp,c=xtopo, edgecolors='none')

valid_mask = np.isfinite(xtopo) & np.isfinite(yvertvel)
xtopo_clean = xtopo[valid_mask]
yvertvel_clean = yvertvel[valid_mask]

params = np.polyfit(xtopo_clean, yvertvel_clean, 1)  # Le '1' indique un fit linéaire (degré 1)
x = np.linspace(np.min(xtopo), np.max(xtopo), 1000)
y = params[0] * x + params[1]

axcomparebig.plot(x, y, color='r', label=f'a.x + b (a = {params[0]:.4f}, b = {params[1]:.4f})')
print(params)

#plt.colorbar()
axcomparebig.set_xlabel('z (m)', fontsize=14)
axcomparebig.set_ylabel('LOS velocity (mm/yr)', fontsize=14)
axcomparebig.set_title(outfile, fontsize=14)

#fig3.canvas.set_window_title('Cropped data')

# Display the data
plt.tight_layout()
# ax.set_rasterized(True)
try:
    fig3.savefig('{}_crop_plot.png'.format(outfile), format='PNG',dpi=180)
except:
    pass

########### Show all plots

plt.show()



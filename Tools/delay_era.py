#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine (CRPG)
################################################################################


"""\
delay_era.py
-------------
Generates a plot of delay vs. time between two elevations.

Usage: 
    delay_era.py --filespath=<path> --alt_max=<value> --alt_min=<value> [--reg_type=<value>]

Options:
    -h, --help             Show this help message.
    --filespath=<path>     Path to the folder created by the NSBAS function nsb_extract_erai.pl (e.g., 58_10_6h).
    --alt_max=<value>      Maximum altitude.
    --alt_min=<value>      Minimum altitude.
    --reg_type=<value>     Type of regretion (lin_amp/amp), default : amp
"""

print()
print()
print('Author: Hugo Watine')
print()
print()

from datetime import datetime
import os, sys
import numpy as np
from osgeo import gdal, osr
import pandas as pd
from matplotlib import pyplot as plt

try:
    from nsbas import docopt
except:
    import docopt

import docopt
arguments = docopt.docopt(__doc__)

### Fonction ###

def calcul_delay(path, alt_max, alt_min):
    data = np.loadtxt(path, comments='#', unpack=True)
    alt = data[0]
    delay_m = data[3]

    index_max = np.where(alt == alt_max)[0][0]
    index_min = np.where(alt == alt_min)[0][0]

    delay = delay_m[index_min] - delay_m[index_max]

    return delay

def fit_linear_model(x, data, RMS, reg_type='amp'):

    if reg_type == 'amp':
        G = np.vstack([np.ones_like(x), np.cos(2*np.pi*x), np.sin(2*np.pi*x)]).T
    if reg_type == 'lin_amp':
        G = np.vstack([x, np.ones_like(x), np.cos(2*np.pi*x), np.sin(2*np.pi*x)]).T
    
    Cd = np.diag(RMS ** 2)  # Matrice de covariance basée sur RMS

    # Calcul du paramètre de régression a et b en minimisant les moindres carrés pondérés
    params = np.dot(np.linalg.inv(np.dot(np.dot(G.T, np.linalg.inv(Cd)), G)),
                    np.dot(np.dot(G.T, np.linalg.inv(Cd)), data))

    return params

path = arguments["--filespath"]
alt_max = arguments["--alt_max"]
alt_min = arguments["--alt_min"]
if arguments["--filespath"] ==  None or arguments["--alt_max"] ==  None or ["--alt_min"] ==  None:
    print('!! Check input parapeters !!')
    sys.exit()

if arguments["--reg_type"] ==  None:
    reg_type = 'amp'
else:
    reg_type= arguments["--reg_type"]

# Affichage de la position des données
lat, lon, alt = list(map(float, open(path+'/lat_lon_z_clean').readline().split(None, 3)[0:3]))
print('---------------')
print(f"Latitude : {lat}")
print(f"Longitude : {lon}")
print(f"Altitude : {alt}")
print('---------------')
print('')

#Extraction des fichiers .del
del_files = [f for f in os.listdir(path) if f.endswith('.del')]

# Lecture de chaques fichier .del pour en extraire la date et le delais entre les deux alitudes en entrée
date = []
DOY = []
delay = []

alt_max_round = round(float(alt_max) / 100) * 100
alt_min_round = round(float(alt_min) / 100) * 100

print(f'** Delay calculation between {alt_min_round} and {alt_max_round} m **')
print()

for file in del_files:
    formatted_dates = datetime.strptime(str(file.split('.')[0]), "%Y%m%d")
    year_fraction = formatted_dates.year + (formatted_dates.timetuple().tm_yday - 1) / 365.25
    doy = (formatted_dates.timetuple().tm_yday - 1) / 365.25
    
    date.append(year_fraction)
    DOY.append(doy)
    delay.append(calcul_delay(path+f'/{file}', alt_max_round, alt_min_round))

# Mettre la première date à 0
index = np.argmin(date)
delay = [d - delay[index]for d in delay]

# Affichage 

date = np.array(date)
DOY = np.array(DOY)
delay = np.array(delay) * 1000 # m to mm

sorted_indices = np.argsort(date)
date = date[sorted_indices]
DOY = DOY[sorted_indices]
delay = delay[sorted_indices]

fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True, gridspec_kw={'width_ratios': [2, 1]})

# ---- Graphique principal (delay vs date) ----
axes[0].plot(date, delay, 'o', c='black', alpha=0.6)

params = fit_linear_model(date, delay, np.ones_like(delay), reg_type)
x = np.linspace(np.min(date), np.max(date), 1000)

if reg_type == 'lin_amp':   
    y = params[0] * x + params[1] + np.sqrt(params[2]**2 + params[3]**2) * np.cos(2*np.pi*x + np.arctan2(-params[3], params[2]))
    axes[0].plot(x, y, color='r', label=f'a.x + b + c.cos(2.pi.x + phi) ; (a = {params[0]:.4f}, b = {params[1]:.4f}, c = {np.sqrt(params[2]**2 + params[3]**2):.4f}, phi = {np.arctan2(-params[3],params[2]):.4f})')
if reg_type == 'amp':
    y = np.sqrt(params[1]**2 + params[2]**2) * np.cos(2*np.pi*x + np.arctan2(-params[2], params[1])) + params[0]
    axes[0].plot(x, y, color='r', label=f'a.cos(2.pi.x + phi) + b ;  (a = {np.sqrt(params[1]**2 + params[2]**2):.4f}, phi = {np.arctan2(-params[2],params[1]):.4f}, b = {params[0]:.4f})')

axes[0].set_xlabel("Date")
axes[0].set_ylabel("Tropospheric delay (mm)")
axes[0].legend()
axes[0].set_title("Delay vs Date")

# ---- Second graphique (delay vs DOY) ----

axes[1].plot(DOY, delay, 'o', c='black', alpha=0.6)


print('delay std =', np.std(delay)) 
params = fit_linear_model(DOY, delay, np.ones_like(delay), reg_type)
x = np.linspace(np.min(DOY), np.max(DOY), 1000)

if reg_type == 'lin_amp':   
    y = params[0] * x + params[1] + np.sqrt(params[2]**2 + params[3]**2) * np.cos(2*np.pi*x + np.arctan2(-params[3], params[2]))
    axes[1].plot(x, y, color='r', label=f'a.x + b + c.cos(2.pi.x + phi) ; (a = {params[0]:.4f}, b = {params[1]:.4f}, c = {np.sqrt(params[2]**2 + params[3]**2):.4f}, phi = {np.arctan2(-params[3],params[2]):.4f})')
if reg_type == 'amp':
    y = np.sqrt(params[1]**2 + params[2]**2) * np.cos(2*np.pi*x + np.arctan2(-params[2], params[1])) + params[0]
    axes[1].plot(x, y, color='r', label=f'a.cos(2.pi.x + phi) + b ;  (a = {np.sqrt(params[1]**2 + params[2]**2):.4f}, phi = {np.arctan2(-params[2],params[1]):.4f}, b = {params[0]:.4f})')

axes[1].set_xlabel("DOY (Day of Year)")
axes[1].set_title("Delay vs DOY")
#axes[1].legend()

# Ajustement de la mise en page
plt.tight_layout()
plt.show()

# plt.figure(figsize=(8, 6))
# months = [int((d % 1) * 12) + 1 for d in DOY]
# df = pd.DataFrame({'Month': months, 'Delay': delay})

# box = df.boxplot(column='Delay', by='Month', grid=False)

# plt.xlabel("Month")
# plt.ylabel("Delay (mm)")
# plt.title("Distribution of Delay by Month")
# plt.suptitle("")  # Supprimer le titre automatique ajouté par pandas

# # Style global plus agréable
# plt.xticks(rotation=45)
# plt.tight_layout()
# plt.grid(axis='y', linestyle='--', alpha=0.7)
# plt.show()




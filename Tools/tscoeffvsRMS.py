#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tscoeffvsRMS.py
------------

Usage: tscoeffvsRMS.py <coeff>  <ts_coeff> <RMSinterfero>

Options:
-h | --help         Show this screen
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from cmcrameri import cm
from datetime import datetime, timedelta
import matplotlib.gridspec as gridspec
import sys

try:
    from nsbas import docopt
except:
    import docopt


def fit_linear_model(x, data, RMS=None, reg_type='lin'):
    if reg_type == 'lin':
        G = np.vstack([x, np.ones_like(x)]).T  
    if reg_type == 'amp':
        G = np.vstack([x, np.ones_like(x), np.cos(2*np.pi*x), np.sin(2*np.pi*x)]).T
    if reg_type == 'acc':
        G = np.vstack([x**2, x, np.ones_like(x), np.cos(2*np.pi*x), np.sin(2*np.pi*x)]).T
    if reg_type =='ampi':
        G = np.vstack([x, np.ones_like(x), np.cos(2*np.pi*x), np.sin(2*np.pi*x), x*np.cos(2*np.pi*x), x*np.sin(2*np.pi*x)]).T
            
    if RMS ==None:
        Cd = np.eye(len(x))
    else :
        Cd = np.diag(RMS ** 2)  # Matrice de covariance basée sur RMS

    params, *_ = np.linalg.lstsq(G, data, rcond=None)
    return params

 #   # Calcul du paramètre de régression a et b en minimisant les moindres carrés pondérés
 #   params = np.dot(np.linalg.inv(np.dot(np.dot(G.T, np.linalg.inv(Cd)), G)),
 #                   np.dot(np.dot(G.T, np.linalg.inv(Cd)), data))
#
#    return params

def dec2date(times):
    times = np.atleast_1d(times)
    dates = []

    for t in times:
        year = int(np.floor(t))
        dec = t - year
        day_of_year = int(np.round(dec * 365.1))

        date = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)
        dates.append(date.strftime('%Y%m%d'))
    return dates

def date2dec(dates):
    dates  = np.atleast_1d(dates)
    times = []
    for date in dates:
        x = datetime.strptime('{}'.format(date),'%Y%m%d')
        dec = float(x.strftime('%j'))/365.1
        year = float(x.strftime('%Y'))
        times.append(year + dec)
    return times

arguments = docopt.docopt(__doc__)
path_coeff = arguments["<coeff>"]
path_tscoeff = arguments["<ts_coeff>"]
path_RMS = arguments["<RMSinterfero>"]

date1, date2, coeff = np.loadtxt(path_coeff, comments="#", unpack=True, dtype='i,i,f', usecols=(0,1,2))
datedec, tscoeff = np.loadtxt(path_tscoeff, comments="#", unpack=True, dtype='f,f', usecols=(0,1))
date1RMS, date2RMS, RMSinterfero = np.loadtxt(path_RMS, comments="#", unpack=True, dtype='i,i,f', usecols=(1,2,3))

coeff_pairs = np.array([f"{d1}_{d2}" for d1, d2 in zip(date1, date2)])
RMS_pairs   = {f"{d1}_{d2}": rms for d1, d2, rms in zip(date1RMS, date2RMS, RMSinterfero)}

common_pairs = set(coeff_pairs) & set(RMS_pairs.keys())
only_in_coeff = set(coeff_pairs) - set(RMS_pairs.keys())
only_in_RMS   = set(RMS_pairs.keys()) - set(coeff_pairs)

# --- Affichage pour debug ---
print(f"✅ Nombre de paires communes : {len(common_pairs)}")
print(f"❌ Paires dans coeff mais pas dans RMS: {len(only_in_coeff)}")
print(f"❌ Paires dans RMS mais pas dans coeff: {len(only_in_RMS)}")

mask = np.array([p not in only_in_coeff for p in coeff_pairs])
date1 = date1[mask]
date2 = date2[mask]
coeff = coeff[mask]
RMSinterfero = np.array([RMS_pairs[p] for p in coeff_pairs[mask]])

datedec1 = np.array(date2dec(date1))
datedec2 = np.array(date2dec(date2))

params = fit_linear_model(datedec, tscoeff, RMS=None, reg_type='amp')
a, b, c, phi = params[0], params[1], np.sqrt(params[2]**2 + params[3]**2), np.arctan2(-params[3],params[2])
coeff1 = a * datedec1 + b + c*np.cos(2*np.pi*datedec1+phi)
coeff2 = a * datedec2 + b + c*np.cos(2*np.pi*datedec2+phi)

coeff_reconstruct = coeff2 - coeff1
residual = coeff - coeff_reconstruct

fig = plt.figure(figsize=(12,5))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
ax_left = fig.add_subplot(gs[0, 0])
ax_right = fig.add_subplot(gs[0, 1])

x = np.linspace(np.min(datedec), np.max(datedec), 500)
y =  params[0] * x + params[1] + np.sqrt(params[2]**2 + params[3]**2)*np.cos(2*np.pi*x+np.arctan2(-params[3],params[2]))

ax_left.plot(x, y, color='r', label=f'a.x + b + c.cos(2.pi.x + phi) (a = {params[0]:.4f}, b = {params[1]:.4f}, c = {np.sqrt(params[2]**2 + params[3]**2):.4f}, phi = {np.arctan2(-params[3],params[2]):.4f})')
ax_left.plot(datedec, tscoeff, 'ko')
ax_left.set_xlabel("Année")
ax_left.set_ylabel("Coefficient")
ax_left.set_title("Série temporelle")
ax_left.legend()
ax_left.grid(True)

ax_right.plot(RMSinterfero, abs(residual), 'ko')
ax_right.axhline(0, color='r', linestyle='--')
ax_right.set_xlabel("RMS")
ax_right.set_ylabel("Coeff - Coeff reconstruit")
ax_right.set_title("Résidus")
ax_right.grid(True)


print(len(date1), len(date2), len(coeff), len(coeff_reconstruct), len(residual), len(RMSinterfero))
sys.exit()
out = np.column_stack([date1, date2, coeff, coeff_reconstruct, residual, RMSinterfero])
np.savetxt(
    "coeff_residuals.txt",
    out,
    fmt="%8d %8d % .6f % .6f % .6f % .6f",
    header="date1 date2 coeff coeff_reconstruct residual RMS",
    comments=""
)

plt.show()





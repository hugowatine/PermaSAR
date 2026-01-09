#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################


"""\
invert_phi.py
-------------
Time series inversion: solve phase/noise/spectrum for each acquisition dates using phase/noise/spectrum of all interferograms
Solve problem : phi(AB) = phi(B) - phi(A) 

Usage: invert_phi.py [--datesfile=<path>] --input=<path> --output=<path> [--noise=<yes/no>]  [--cons=<yes/no>] [--plot] [--csv=<path>] [--rad2mm=<value>]
invert_phi.py  -h | --help

Options:
-h --help           Show this screen.
--datesfile PATH    list images file [default: baseline.rsc]
--input PATH        Infile containing the phase for each pair of int.
--output PATH       Outfile 
--noise PATH	    If yes, solve problem : sigma(AB)^2 = sigma(A)^2 + sigma(B)^2 [default: no]
--cons VALUE        Add postive constrain to the inversion
--plot    Plot results
"""

# docopt (command line parser)
import docopt

import scipy.linalg as lst
import scipy.optimize as opt
import numpy as np
import math
import pandas as pd
from datetime import datetime
import sys

def consInvert(A,b,sigmad=1,ineq=[None,None], cond=1.0e-10, iter=250,acc=1e-06):
    '''Solves the constrained inversion problem.

    Minimize:
    
    ||Ax-b||^2

    Subject to:
    Ex >= f
    '''
    
    Ain = A
    bin = b

    if Ain.shape[0] != len(bin):
        raise ValueError('Incompatible dimensions for A and b')

    Ein = ineq[0]
    fin = ineq[1]

    if Ein is not None:
        if Ein.shape[0] != len(fin):
            raise ValueError('Incompatible shape for E and f')
        if Ein.shape[1] != Ain.shape[1]:
            raise ValueError('Incompatible shape for A and E')

    ####Objective function and derivative
    _func = lambda x: np.sum(((np.dot(Ain,x)-bin)/sigmad)**2)
    _fprime = lambda x: 2*np.dot(Ain.T/sigmad, (np.dot(Ain,x)-bin)/sigmad)

    ######Inequality constraints and derivative
    if Ein is not None:
        _f_ieqcons = lambda x: np.dot(Ein,x)-fin
        _fprime_ieqcons = lambda x: Ein

    ######Actual solution of the problem
    temp = lst.lstsq(Ain,bin,cond=cond)   ####Initial guess.
    x0 = temp[0]

    if Ein is None:
        res = temp
    else:
        res = opt.fmin_slsqp(_func,x0,f_ieqcons=_f_ieqcons,fprime=_fprime, fprime_ieqcons=_fprime_ieqcons,iter=iter,full_output=True,acc=acc)
        if res[3] != 0:
            print('Exit mode %d: %s \n'%(res[3],res[4]))

    fsoln = res[0]
    return fsoln

# read arguments
arguments = docopt.docopt(__doc__)
liste_int = arguments["--input"]
outfile = arguments["--output"]
if arguments["--cons"] == None:
    cons = 'no'
else:
    cons = arguments["--cons"]
if arguments["--datesfile"] ==  None:
    basefile = "baseline.rsc"
    im,imd=np.loadtxt(basefile,comments="#",usecols=(0,4),unpack=True,dtype='i,f')
else:
    basefile=arguments["--datesfile"]
    im,imd=np.loadtxt(basefile,comments="#",usecols=(0,4),unpack=True,dtype='i,f')
if arguments["--noise"] == None:
    noise = 'no'
else:
    noise = arguments["--noise"]

#Data loading
print("int list=",liste_int)
date1,date2,spint=np.loadtxt(liste_int,comments="#",unpack=True,dtype='i,i,f', usecols=(0,1,2))

mask = ~np.isnan(spint) & ~np.isinf(spint)
date1, date2, spint = date1[mask], date2[mask], spint[mask]

kmax=len(date1)
print("number of interferogram: ",kmax)


print("image list=",basefile)
nmax=len(imd)
print("number of image: ",nmax)

#build G
G=np.zeros((kmax+1,nmax))
if noise=='yes':
  for k in range((kmax)):
    for n in range((nmax)):
        if (date1[k]==im[n]): 
          G[k,n]=1
        elif (date2[k]==im[n]):
          G[k,n]=1
else:
  for k in range((kmax)):
    for n in range((nmax)):
        if (date1[k]==im[n]): 
          G[k,n]=-1
        elif (date2[k]==im[n]):
          G[k,n]=1
# ini phi first image to 0 
G[-1,0]=1

#build d
d=np.zeros((kmax+1))
d[:kmax]=spint

print()
# Constrain
if cons=='yes':
    print("Add positive constrain to the inversion")
    f = np.zeros(nmax)
    E = np.diag(np.ones(nmax))
    
    print("Inversion....")
    sp = consInvert(G,d,ineq=[E,f])

else:
    print("Inversion....")
    sp = consInvert(G,d)

## check resolution of the inversion
#Res = np.diag(np.dot(np.dot(G,np.linalg.pinv(np.dot(G.T,G))),G.T))
#print('date , phase, res: ')
#for n in range(len(imd)):
#    print(imd[n], sp[n], Res[n])
#print()

#for n in range(len(imd)):
#    print(imd[n], sp[n])
#print()

# save in output file
np.savetxt(outfile, np.vstack([imd,sp]).T, fmt='%.6f')

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm
from datetime import datetime
from scipy.signal import medfilt

do_plot = bool(arguments["--plot"])

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

    # Calcul du paramètre de régression a et b en minimisant les moindres carrés pondérés
    params = np.dot(np.linalg.inv(np.dot(np.dot(G.T, np.linalg.inv(Cd)), G)),
                    np.dot(np.dot(G.T, np.linalg.inv(Cd)), data))

    return params

if do_plot:
    rad2mm = float(arguments.get("--rad2mm") or 1.)
    csvfile = arguments.get("--csv")
    if csvfile is not None:
        df = pd.read_csv(csvfile)

        # Détection automatique des colonnes
        # (par ex. "date", "system:time_start", "Day", "Time", etc.)
        col_date = df.columns[0]
        col_temp = df.columns[1]

        # Conversion date -> datetime
        df[col_date] = pd.to_datetime(df[col_date], errors='coerce')

        # Retire les lignes invalides
        df = df.dropna(subset=[col_date, col_temp])

        # Conversion datetime -> année décimale
        def to_decimal_year(dt):
            year = dt.year
            start = datetime(year, 1, 1)
            end = datetime(year + 1, 1, 1)
            return year + (dt - start).total_seconds() / (end - start).total_seconds()

        date_csv = np.array([to_decimal_year(d) for d in df[col_date]])
        temp_csv = df[col_temp].to_numpy() - 273.15
    else:
        date_csv, temp_csv = None, None


    if len(imd) == 0:
        print("⚠️  Aucun coefficient valide pour tracer le graphique.")
    else:
        import matplotlib.gridspec as gridspec

        if rad2mm != None:
            sp = sp * rad2mm

        # Figure avec 1 ligne, 2 colonnes égales
        fig = plt.figure(figsize=(12,5))
        gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
        ax_left = fig.add_subplot(gs[0, 0])
        ax_right = fig.add_subplot(gs[0, 1])

        # ----- Graphe gauche : série temporelle -----


        params = fit_linear_model(imd, sp, RMS=None, reg_type='amp')
        x = np.linspace(np.min(imd), np.max(imd), 500)

        y =  params[0] * x + params[1] + np.sqrt(params[2]**2 + params[3]**2)*np.cos(2*np.pi*x+np.arctan2(-params[3],params[2]))
        ax_left.plot(x, y, color='r', label=f'a.x + b + c.cos(2.pi.x + phi) (a = {params[0]:.4f}, b = {params[1]:.4f}, c = {np.sqrt(params[2]**2 + params[3]**2):.4f}, phi = {np.arctan2(-params[3],params[2]):.4f})')
        ax_left.plot(imd, sp, 'ko')
        ax_left.set_xlabel("Année")
        ax_left.set_ylabel("Coefficient")
        ax_left.set_title("Série temporelle")
        ax_left.legend()
        ax_left.grid(True)

        # ----- Graphe droit : série enroulée sur 1 an -----
        # convertir années décimales en mois
        months = imd % 1
        import datetime
        dates = [datetime.datetime(2019, 1, 1) + datetime.timedelta(days=float(v*365))
         for v in months]
        ax_right.plot(dates, sp, 'ko', markersize=4)
        ax_right.xaxis.set_major_formatter(mdates.DateFormatter("%m"))
        ax_right.xaxis.set_major_locator(mdates.MonthLocator())
        plt.gcf().autofmt_xdate()

        if date_csv is not None:
            months_csv = date_csv % 1
            ax_temp2 = ax_right.twinx()

            params = fit_linear_model(months_csv, temp_csv, RMS=None, reg_type='amp')
            #x = np.linspace(np.min(months_csv), np.max(months_csv), 1000)
            #y =  params[0] * months_csv + params[1] + np.sqrt(params[2]**2 + params[3]**2)*np.cos(2*np.pi*months_csv+np.arctan2(-params[3],params[2]))
            #ax_temp2.plot(months_csv, y, 'r--', alpha=0.6, label=f'a.x + b + c.cos(2.pi.x + phi) (a = {params[0]:.4f}, b = {params[1]:.4f}, c = {np.sqrt(params[2]**2 + params[3]**2):.4f}, phi = {np.arctan2(-params[3],params[2]):.4f})')
            
            ax_temp2.plot(months_csv, temp_csv, 'r.', alpha=0.6)
            ax_temp2.set_ylabel("Température (K)", color='r')
            ax_temp2.tick_params(axis='y', labelcolor='r')
        

        ax_right.set_xlabel("Mois")
        #ax_right.set_xlim([0., 1.])

    #     ax_right.set_xticks([
    # 0.00000,  # January
    # 0.08493,  # February
    # 0.16164,  # March
    # 0.24727,  # April
    # 0.32877,  # May
    # 0.41027,  # June
    # 0.49315,  # July
    # 0.57534,  # August
    # 0.65753,  # September
    # 0.73973,  # October
    # 0.82192,  # November
    # 0.90411   # December
    # ])
        #ax_right.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
        ax_right.set_title("Cycle annuel")
        ax_right.grid(True)

        plt.tight_layout()
        plt.savefig(outfile + ".pdf", dpi=300, bbox_inches='tight')
        plt.show()

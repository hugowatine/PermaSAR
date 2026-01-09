#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Hugo Watine (CRPG)
################################################################################

"""\
flatsim_process_card.py
-----------------------
Compare 2 TS analyses with the first order data quality:
RMSpixel, RMSdate, RMSint, model coeff, int weight, network.

Usage:
    flatsim_process_card.py --ts_path=<path> [--ts2_path=<path>]
    flatsim_process_card.py -h | --help

Options:
    --ts_path=<path>    Path to the first TS folder
    --ts2_path=<path>   Path to the second TS folder for comparison
    -h --help           Show this help message
"""

import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal
import os
import docopt
import matplotlib.dates as mdates
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import pydot
from matplotlib.dates import date2num
from datetime import datetime
import sys
import glob
import matplotlib as mpl
from datetime import datetime, timedelta

arguments = docopt.docopt(__doc__)
ts_path = arguments["--ts_path"]
ts2_path = arguments.get("--ts2_path")

def plot_network(ax, pair_file, baseline_file, use_weight=False, flatsim=False):
    """
    Dessine le réseau d'interférogrammes dans un axes Matplotlib
    avec colorimétrie et affichage similaire à nsb_plot_interferograms_network.
    """

    # Charger le graphe
    graph = pydot.Dot("interferogram_network", graph_type="digraph")

    # Ajouter les noeuds (baseline)
    for line in open(baseline_file, "r"):
        d, bp = line.split()[0], float(line.split()[1])
        graph.add_node(pydot.Node(d, label=d, bperp=bp))

    # Ajouter les edges (paires)
    pairs, weights = [], []
    for line in open(pair_file, "r"):
        p = line.split()
        if len(p) < 2: 
            continue
        if flatsim:
            d1, d2 = p[1], p[2]
            w = float(p[3]) if (use_weight and len(p) >= 4) else None
        else:
            d1, d2 = p[0], p[1]
            w = float(p[2]) if (use_weight and len(p) >= 3) else None
        graph.add_edge(pydot.Edge(d1, d2))
        pairs.append((d1,d2))
        weights.append(w)

    # Normalisation couleurs
    if use_weight and any(w is not None for w in weights):
        valid_weights = [w for w in weights if w is not None]
        vmin, vmax = np.min(valid_weights), np.max(valid_weights)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.cm.managua
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
    else:
        sm = None

    # Tracer les edges
    for i, edge in enumerate(graph.get_edges()):
        n1 = graph.get_node(edge.get_source())[0]
        n2 = graph.get_node(edge.get_destination())[0]
        x0 = date2num(datetime.strptime(n1.get_label(), "%Y%m%d"))
        y0 = float(n1.get_attributes()['bperp'])
        x1 = date2num(datetime.strptime(n2.get_label(), "%Y%m%d"))
        y1 = float(n2.get_attributes()['bperp'])
        dx, dy = x1 - x0, y1 - y0
        color = cmap(norm(weights[i])) if (use_weight and weights[i] is not None) else 'black'
        ax.arrow(x0, y0, dx, dy, linewidth=0.5, color=color, alpha=0.5, length_includes_head=True)

    # Tracer les noeuds
    xs, ys = [], []
    for node in graph.get_nodes():
        xs.append(date2num(datetime.strptime(node.get_label(), "%Y%m%d")))
        ys.append(float(node.get_attributes()['bperp']))
    ax.plot(xs, ys, "o", markersize=4, color="dodgerblue", mec="black", alpha=0.8, picker=5)

    # Axes et labels
    ax.set_xlabel("Acquisition Date")
    if flatsim == False:
        ax.set_ylabel("Perpendicular Baseline (m)")
    ax.set_title("Interferogram network")
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%Y/%m/%d"))
    ax.tick_params(axis="x", rotation=45)

    # Ajouter colorbar si weights
    if sm:
        plt.colorbar(sm, ax=ax, label="Weight value")

def open_gdal(file, band=1, supp_ndv=None, complex=False):
    """
    Use GDAL to open band as real value or complex interferogram.
    Returns (data, Xsize, Ysize)
    If complex=True: data = [amplitude, phase]
    """

    if not os.path.isfile(file):
        raise FileNotFoundError('File does not exists: {}'.format(file))
    ds = gdal.Open(file)

    band = ds.GetRasterBand(band)
    ndv = band.GetNoDataValue()
    data = band.ReadAsArray()
    Xsize = ds.RasterXSize
    Ysize = ds.RasterYSize
    if ndv is not None and ndv != np.nan:
        data[data==ndv] = np.nan
    if supp_ndv is not None and supp_ndv != np.nan:
        data[data==supp_ndv] = np.nan
    return data

def load_optional_file(primary, fallback):
    """Renvoie le fichier existant entre primary et fallback, ou None."""
    if os.path.exists(primary):
        return primary
    if os.path.exists(fallback):
        return fallback
    return None


def opents_data(folder_path):
    """Charge les fichiers RMS d'un dossier et renvoie un dictionnaire."""
    d = {}
    # --- RMSdate ---
    f = load_optional_file(
        os.path.join(folder_path, 'RMSdate'),
        os.path.join(folder_path, 'RMSdate.txt')
    )
    if f:
        date, rmsdate = np.loadtxt(f, comments='#', usecols=(1, 2),
                                   unpack=True, dtype='i,f')
        d['date'] = [datetime.strptime(str(x), "%Y%m%d") for x in date]
        d['rmsdate'] = rmsdate

    # --- RMSinterfero ---
    f = load_optional_file(
        os.path.join(folder_path, 'RMSinterfero'),
        os.path.join(folder_path, 'RMSinterfero.txt')
    )
    if f:
        nbint, rmsint, d1, d2 = np.loadtxt(f, comments='#',
                                           usecols=(0,3,1,2),
                                           unpack=True, dtype='f,f,i,i')
        d['nbint'], d['rmsint'], d['d1'], d['d2'] = nbint, rmsint, d1, d2

    # --- RMSmap (RMSpixel ou TIFF) ---
    rmspix = os.path.join(folder_path, 'RMSpixel')
    lec = os.path.join(folder_path, 'lect.in')
    tiff = os.path.join(folder_path, 'CNES_Net_geo_8rlks.tiff')

    if os.path.exists(rmspix) and os.path.exists(lec):
        ncols, nlines = map(int, open(lec).readline().split()[:2])
        rmsmap = np.fromfile(rmspix, dtype=np.float32)[:nlines*ncols]
        d['rmsmap'] = rmsmap.reshape((nlines, ncols)).astype(float)
        d['rmsmap'][d['rmsmap'] == 0] = np.nan
    elif os.path.exists(tiff):
        d['rmsmap'] = open_gdal(tiff)

    # --- Network ----
    faux = os.path.join(folder_path, '../2*/baseleline.rsc') # ../2*/baseleline.rsc

    return d

def align_rmsint(ts1_data, ts2_data):
    """
    Aligne les RMSinterfero de ts1 et ts2 sur les mêmes paires (d1,d2).
    Les valeurs manquantes sont remplacées par np.nan.
    """
    # Création d'un dictionnaire pour ts1 et ts2 : clé = (d1,d2)
    ts1_dict = {(int(d1), int(d2)): rms for d1, d2, rms in zip(ts1_data['d1'], ts1_data['d2'], ts1_data['rmsint'])}
    ts2_dict = {(int(d1), int(d2)): rms for d1, d2, rms in zip(ts2_data['d1'], ts2_data['d2'], ts2_data['rmsint'])}

    # Ensemble de toutes les paires uniques
    all_pairs = sorted(set(ts1_dict.keys()) | set(ts2_dict.keys()))

    # Créer des listes alignées
    rms1_aligned = []
    rms2_aligned = []

    for pair in all_pairs:
        rms1_aligned.append(ts1_dict.get(pair, np.nan))
        rms2_aligned.append(ts2_dict.get(pair, np.nan))

    # Retourner un dictionnaire avec les paires et les RMS alignés
    nbint_aligned = np.arange(len(all_pairs))
    return {
        'pairs': all_pairs,
        'rms1': np.array(rms1_aligned),
        'rms2': np.array(rms2_aligned),
        'nbint': nbint_aligned
    }


# Charger les données
ts1_data = opents_data(ts_path)
ts2_data = opents_data(ts2_path) if ts2_path else None

# Exemple d'utilisation
if ts2_data:
    aligned = align_rmsint(ts1_data, ts2_data)


fig = plt.figure(figsize=(18, 9))

# Colonnes égales
width_ratios  = [1, 1, 1, 1, 1, 1, 1, 1]
height_ratios = [1, 1, 1]

gs = gridspec.GridSpec(3, 8, width_ratios=width_ratios, height_ratios=height_ratios,
                       hspace=0.4, wspace=0.3)

# ZONES MODIFIÉES POUR QUE RMSmap AIT 2 COLONNES DE LARGE
zones = {
    1: (0, 2, 0, 2),  # RMSmap ts_path -> 2 colonnes (au lieu de 1)
    2: (0, 2, 2, 4),  # RMSmap ts2_path -> 2 colonnes
    3: (0, 1, 4, 6),  # RMSdate (inchangé, décalé)
    4: (1, 2, 4, 6),  # RMSinterfero
    5: (2, 3, 0, 2),  # Placeholder
    6: (2, 3, 2, 4),  # Placeholder
    7: (2, 3, 4, 6),  # Placeholder
    8: (0, 1, 6, 8),
    9: (1, 2, 6, 8)
}

# RMSmap ts_path
im_list = []
r0, r1, c0, c1 = zones[1]
ax1 = fig.add_subplot(gs[r0:r1, c0:c1])
im1 = ax1.imshow(ts1_data['rmsmap'], cmap=cm.viridis, origin='upper',
                 interpolation='nearest', vmin=0, vmax=1., alpha=0.9)
ax1.set_title('RMSmap Crpg')
im_list.append(im1)
ax1.set_xticks([]); ax1.set_yticks([])
for s in ax1.spines.values(): s.set_visible(False)

# RMSmap ts2
r0, r1, c0, c1 = zones[2]
ax2 = fig.add_subplot(gs[r0:r1, c0:c1])
if ts2_data:
    im2 = ax2.imshow(ts2_data['rmsmap'], cmap=cm.viridis, origin='upper',
                     interpolation='nearest', vmin=0, vmax=1., alpha=0.9)
    ax2.set_title('RMSmap Flatsim')
    im_list.append(im2)
    ax2.set_xticks([]); ax2.set_yticks([])
    for s in ax2.spines.values(): s.set_visible(False)

# COLORBAR UNIQUE
cbar = fig.colorbar(
    im_list[0],
    ax=[ax1, ax2] if ts2_data else [ax1],
    orientation='vertical',
    shrink=0.8,
    pad=0.1
)
cbar.set_label("RMS (0–1)")

# RMSdate
r0, r1, c0, c1 = zones[3]
ax3 = fig.add_subplot(gs[r0:r1, c0:c1])
if ts2_data:
    ax3.plot(ts2_data['date'], ts2_data['rmsdate'], 'ro', markersize=3, label='Flatsim', alpha=0.3)
ax3.plot(ts1_data['date'], ts1_data['rmsdate'], 'bo', markersize=3, label='Crpg', alpha=0.3)
ax3.set_xlabel('Date');  ax3.set_ylabel('RMSdate')
ax3.legend(); ax3.set_title('RMSdate')
ax3.set_ylim([0,2])

# RMSinterfero
r0, r1, c0, c1 = zones[4]
ax4 = fig.add_subplot(gs[r0:r1, c0:c1])
if ts2_data:
    aligned = align_rmsint(ts1_data, ts2_data)
    ax4.plot(aligned['nbint'], aligned['rms2'], 'rs-', markersize=3, label='Flatsim', alpha=0.3)
    ax4.plot(aligned['nbint'], aligned['rms1'], 'bs-', markersize=3, label='Crpg', alpha=0.6)
else:
    ax4.plot(ts1_data['nbint'], ts1_data['rmsint'], 'bs-', markersize=3, label='Crpg', alpha=0.6)

ax4.set_xlabel('Interfero')
ax4.set_ylabel('RMSint')
ax4.set_ylim([0,2])
ax4.legend()
ax4.set_title('RMSinterfero')

# Network plot

r0, r1, c0, c1 = zones[5]
ax_net = fig.add_subplot(gs[r0:r1, c0:c1])

pair_file = os.path.join(ts_path, "list_pair")
baseline_pattern = os.path.join(".", "2*/baseline.rsc")
baseline_matches = glob.glob(baseline_pattern)
baseline_file = baseline_matches[0] if baseline_matches else None

if os.path.exists(pair_file) and baseline_file and os.path.exists(baseline_file):
    plot_network(ax_net, pair_file, baseline_file, use_weight=True)
else:
    ax_net.text(0.5, 0.5, "network files missing", ha="center", va="center")
    ax_net.set_xticks([])
    ax_net.set_yticks([])

r0, r1, c0, c1 = zones[6]
ax_net = fig.add_subplot(gs[r0:r1, c0:c1])

if ts2_path is None:
    ax_net.text(0.5, 0.5, "ts2_path missing", ha="center", va="center")
    ax_net.set_xticks([]); ax_net.set_yticks([])
else:
    pair_file = os.path.join(ts2_path, 'RMSinterfero.txt')

    baseline_pattern = os.path.join(".", "2*/baseline.rsc")
    baseline_matches = glob.glob(baseline_pattern)
    baseline_file = baseline_matches[0] if baseline_matches else None

    if os.path.exists(pair_file) and baseline_file and os.path.exists(baseline_file):
        plot_network(ax_net, pair_file, baseline_file, use_weight=False, flatsim=True)
    else:
        ax_net.text(0.5, 0.5, "network files missing", ha="center", va="center")
        ax_net.set_xticks([]); ax_net.set_yticks([])

# temporal baseline
r0, r1, c0, c1 = zones[7]
ax7 = fig.add_subplot(gs[r0:r1, c0:c1])

dates_d1 = np.array([datetime.strptime(str(d), "%Y%m%d") for d in ts1_data['d1']])
dates_d2 = np.array([datetime.strptime(str(d), "%Y%m%d") for d in ts1_data['d2']])
time_baselines = np.array([(d2 - d1).days for d1, d2 in zip(dates_d1, dates_d2)])

# Histogramme
ax7.hist(time_baselines, bins=30, color='steelblue', alpha=0.7, edgecolor='black')
ax7.set_xlabel("Temporal baseline (days)")
#ax7.set_ylabel("Number of interferograms")
ax7.set_title("Temporal baseline histogram")
ax7.set_ylim([0,200])

baseline_range = np.linspace(0.1, time_baselines.max(), 200)  # éviter division par zéro
weight_curve = np.exp(-219 / baseline_range) + 0.1
weight_curve[weight_curve < 0.1] = 0.1
ax8 = ax7.twinx()
ax8.plot(baseline_range, weight_curve, 'r-', lw=2, label='weight exp(-0.6 / baseline) + 0.1')
ax8.set_ylabel('Weight', color='red')
ax8.tick_params(axis='y', colors='red')

def compute_model(params, x):
    """ Génère a.x + b + c*cos(2*pi*x + phi) """
    a, b, c2, s2 = params
    amp = np.sqrt(c2**2 + s2**2)
    phi = np.arctan2(-s2, c2)
    return a*x + b + amp * np.cos(2*np.pi*x + phi)

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

# ---- Lecture du fichier coef_date.txt ----
coef_file = os.path.join("INTERFERO", "coef_date.txt")
fallback_file = os.path.join("INTERFERO", "coeffs_ts.txt")
if os.path.exists(coef_file):
    data_coef = np.loadtxt(coef_file)
elif os.path.exists(fallback_file):
    data_coef = np.loadtxt(fallback_file)
imd = data_coef[:, 0]   # dates décimales
sp = data_coef[:, 1]    # valeurs associées

# --- Création de la figure (1 ligne, 2 colonnes) ---
r0, r1, c0, c1 = zones[8]   # emplacement du graphique 8
ax_left = fig.add_subplot(gs[r0:r1, c0:c1])
ax_left.set_ylim([-1.5,1.5])

r0, r1, c0, c1 = zones[9]   # emplacement du graphique 9
ax_right = fig.add_subplot(gs[r0:r1, c0:c1])
ax_right.set_ylim([-1.5,1.5])
# ---------- GRAPHE 8 : Série temporelle ----------
# Fit linéaire + cos (identique à ton code d’exemple)
params = fit_linear_model(imd, -sp, RMS=None, reg_type='amp')
x = np.linspace(np.min(imd), np.max(imd), 500)

y = (params[0] * x
     + params[1]
     + np.sqrt(params[2]**2 + params[3]**2)
       * np.cos(2*np.pi * x + np.arctan2(-params[3], params[2])))

ax_left.plot(x, y, 'r', label=f"a.x + b + c cos(2πx+φ)\n(a={params[0]:.3f}, c={np.sqrt(params[2]**2 + params[3]**2):.3f})")
ax_left.plot(imd, -sp, 'ko', markersize=4)

ax_left.set_xlabel("Année décimale")
#ax_left.set_ylabel("Coefficient")
ax_left.set_title("Série temporelle x -1")
ax_left.grid(True)
ax_left.legend()

# ---------- GRAPHE 9 : Cycle annuel ----------
months = imd % 1   # phase entre 0 et 1
dates = [datetime(2019,1,1) + timedelta(days=float(m*365)) for m in months]

ax_right.plot(dates, -sp, 'ko', markersize=4)
ax_right.set_title("Cycle annuel x -1")
ax_right.set_xlabel("Mois")
#ax_right.set_ylabel("Coefficient")

ax_right.xaxis.set_major_formatter(mdates.DateFormatter("%m"))
ax_right.xaxis.set_major_locator(mdates.MonthLocator())

ax_right.grid(True)

# Placeholders
for idx in []:
    r0, r1, c0, c1 = zones[idx]
    ax = fig.add_subplot(gs[r0:r1, c0:c1])
    ax.set_facecolor('lightgray')
    ax.text(0.5, 0.5, str(idx), ha='center', va='center', fontsize=16)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(True)

#plt.subplots_adjust(wspace=0.35, hspace=0.45)
#plt.tight_layout(pad=2.0)
# plt.subplots_adjust(
#     left=0.08,
#     right=0.97,
#     bottom=0.08,
#     top=0.93,
#     wspace=0.35,
#     hspace=0.45
# )
plt.tight_layout()
plt.show()


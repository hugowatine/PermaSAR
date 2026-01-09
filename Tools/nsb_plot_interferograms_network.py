#!/usr/bin/env python3
##########################################################################
#
#   This file is part of NSBAS.
#
#   NSBAS is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   NSBAS is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with NSBAS.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

"""\
nsb_plot_interferograms_network

Usage:
  nsb_plot_interferograms_network [-o <out_path>] [--weight=<yesno>] <pair_file> <baseline_file>
  nsb_plot_interferograms_network -h | --help

Options:
  -o OUT_PATH        Write figure to file.
  --weight=yes/no    Color edges using weight values from <pair_file> (default: no)
  -h --help          Show this screen.
"""

import string, os
from datetime import datetime

import numpy as np
import matplotlib as mpl
from nsbas import docopt
arguments = docopt.docopt(__doc__)

if arguments["-o"]:
    mpl.use('pdf')

import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.dates import date2num, num2date
import pydot

# Figure setup
fig = plt.figure(figsize=(13,6))
ax = fig.add_subplot(111)

# Load graph from <pair_file> and <baseline_file>
graph = pydot.Dot("interferogram_network", graph_type="digraph")

# Load baseline file (dates + perpendicular baselines)
for line in open(arguments["<baseline_file>"], "r"):
    lines = line.split()
    date, pbaseline = lines[0], lines[1]
    graph.add_node(pydot.Node(date.strip(), label=date.strip(), bperp=float(pbaseline)))

# Load pair file (master, slave [, weight])
pairs = []
weights = []

use_weight = arguments["--weight"] == "yes"

for line in open(arguments["<pair_file>"], "r"):
    parts = line.split()
    if len(parts) < 2:
        continue
    date1, date2 = parts[0].strip(), parts[1].strip()
    weight = float(parts[2]) if (use_weight and len(parts) >= 3) else None
    graph.add_edge(pydot.Edge(date1, date2))
    pairs.append((date1, date2))
    weights.append(weight)

# Draw nodes
x, y = [], []
for n in graph.get_nodes():
    x.append(date2num(datetime.strptime(n.get_label(), "%Y%m%d")))
    y.append(float(n.get_attributes()["bperp"]))
ax.plot(x, y, "o", color='dodgerblue', mec='black', markersize=4, picker=5)

# Prepare color mapping if weight is used
if use_weight and any(w is not None for w in weights):
    valid_weights = [w for w in weights if w is not None]
    wmin, wmax = np.min(valid_weights), np.max(valid_weights)
    cmap = plt.cm.Blues
    norm = mpl.colors.Normalize(vmin=wmin, vmax=wmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
else:
    sm = None

# Draw edges
for i, edge in enumerate(graph.get_edges()):
    master = graph.get_node(edge.get_source())[0]
    slave = graph.get_node(edge.get_destination())[0]
    x0 = date2num(datetime.strptime(master.get_label(), "%Y%m%d"))
    y0 = float(master.get_attributes()["bperp"])
    x1 = date2num(datetime.strptime(slave.get_label(), "%Y%m%d"))
    y1 = float(slave.get_attributes()["bperp"])
    dx, dy = x1 - x0, y1 - y0

    if use_weight and weights[i] is not None:
        color = sm.to_rgba(weights[i])
    else:
        color = 'black'

    ax.arrow(x0, y0, dx, dy, linewidth=0.5, color=color, alpha=0.5, length_includes_head=True)

# Register click usage to display date of nearest point
def onpick(event):
    global ax, an
    pt = event.artist
    ind = event.ind
    x = np.take(pt.get_xdata(), ind)[0]
    y = np.take(pt.get_ydata(), ind)[0]
    an.xy = (x, y)
    if hasattr(an, "xyann"): # matplotlib >= 1.4
        an.xyann = (x + abs(ax.get_xlim()[1]-ax.get_xlim()[0])/30,
                    y + abs(ax.get_ylim()[1]-ax.get_ylim()[0])/15)
    else:
        an.xytext = (x + abs(ax.get_xlim()[1]-ax.get_xlim()[0])/30,
                     y + abs(ax.get_ylim()[1]-ax.get_ylim()[0])/15)
    an.set_text(num2date(x).strftime("%Y/%m/%d"))
    an.set_visible(True)
    fig.canvas.draw()

an = ax.annotate("",
                 xy=(0, 0), xycoords="data",
                 xytext=(0, 0), textcoords="data",
                 arrowprops=dict(facecolor="black", width=1, shrink=0.3),
                 bbox=dict(boxstyle="round", fc="w"))
an.set_visible(False)
fig.canvas.mpl_connect("pick_event", onpick)

# Legend and labels
fig.suptitle("Interferogram network")
fig.autofmt_xdate()
ax.set_xlabel("Acquisition Date")
ax.set_ylabel("Perpendicular Baseline (m)")
ax.xaxis.set_major_formatter(dates.DateFormatter("%Y/%m/%d"))

if sm:
    cbar = plt.colorbar(sm, ax=ax, label="Weight value")

# Save or show
if arguments["-o"]:
    plt.savefig(arguments["-o"], bbox_inches="tight")
else:
    plt.show()


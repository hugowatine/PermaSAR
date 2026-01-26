#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Author        : Simon DAOUT (CRPG)
################################################################################

"""compute_shift_geomap.py
-------------------
Extracts spatial metadata from a files and calculates the offsets required to align a geocoded product with a reference geomap.

Usage: compute_shift_geomap.py [--geo=<file>] [--lut=<file>] [-h]

Options:
  --geo=<file>              Path to the geocoded metadata [default: CNES_CosENU_geo_8rlks.unw.rsc]
  --lut=<file>              Path to the Lookup Tables metadata [default: CNES_Lut_geo_b2_8rlks.unw.rsc]
  -h --help                 Show this screen and exit.
"""

import numpy as np
import docopt

def extraire_valeurs_rsc(fichier):
    valeurs = {}
    with open(fichier, 'r') as f:
        for ligne in f:
            ligne = ligne.strip()
            if ligne and not ligne.startswith('#'):
                parties = ligne.split()
                if len(parties) >= 2:
                    cle = parties[0]
                    val = ' '.join(parties[1:])
                    try:
                     valeurs[cle] = float(val)
                    except ValueError:
                        valeurs[cle] = val
    return valeurs

# Extraction .rsc
arguments = docopt.docopt(__doc__)

if arguments["--geo"] ==  None:
    geocoded_path = "CNES_CosENU_geo_8rlks.unw.rsc"
else:
    geocoded_path = arguments["--geo"]

if arguments["--lut"] == None:
    geomap_path = "CNES_Lut_geo_b2_8rlks.unw.rsc"
else:
    geomap_path = arguments["--lut"]

geocoded = extraire_valeurs_rsc(geocoded_path)
geomap = extraire_valeurs_rsc(geomap_path)

# Affichage des .rsc
print("Valeurs extraites du fichier geocoded :")
for cle, val in geocoded.items():
    print(f"{cle}: {val}")

print("\nValeurs extraites du fichier geomap :")
for cle, val in geomap.items():
    print(f"{cle}: {val}")

# Extraction des coordonnées
X_FIRST_geocoded = geocoded["X_FIRST"]
X_STEP_geocoded = geocoded["X_STEP"]
Y_FIRST_geocoded = geocoded["Y_FIRST"]
Y_STEP_geocoded = geocoded["Y_STEP"]
WIDTH_geocoded = geocoded["WIDTH"]
LENGTH_geocoded = geocoded["FILE_LENGTH"]

X_FIRST_geomap = geomap["X_FIRST"]
X_STEP_geomap = geomap["X_STEP"]
Y_FIRST_geomap = geomap["Y_FIRST"]
Y_STEP_geomap = geomap["Y_STEP"]
WIDTH_geomap = geomap["WIDTH"]
LENGTH_geomap = geomap["FILE_LENGTH"]

# Calcule des décalages
XMIN = (X_FIRST_geocoded - X_FIRST_geomap) / X_STEP_geomap
YMIN = (Y_FIRST_geocoded - Y_FIRST_geomap) / Y_STEP_geomap
XMAX = ((X_FIRST_geocoded + WIDTH_geocoded * X_STEP_geocoded) - X_FIRST_geomap) / X_STEP_geomap
YMAX = ((Y_FIRST_geocoded + LENGTH_geocoded * Y_STEP_geocoded) - Y_FIRST_geomap) / Y_STEP_geomap

#print(int(XMAX), int(YMAX))
#XMAX = WIDTH_geomap - ((X_FIRST_geomap + WIDTH_geomap * X_STEP_geomap) - (X_FIRST_geocoded + WIDTH_geocoded * X_STEP_geocoded))/X_STEP_geomap 
#YMAX = LENGTH_geomap - ((Y_FIRST_geomap + LENGTH_geomap * Y_STEP_geomap) - (Y_FIRST_geocoded + LENGTH_geocoded * Y_STEP_geocoded))/Y_STEP_geomap 
#print(int(XMAX), int(YMAX))
#print(Y_FIRST_geocoded)
#print(Y_FIRST_geomap)
#print(Y_FIRST_geocoded + WIDTH_geocoded * Y_STEP_geocoded)
#print(Y_FIRST_geomap + WIDTH_geomap * Y_STEP_geomap)

print(f"\nDécalages calculés : XMIN = {int(XMIN)}, YMIN = {int(YMIN)}, XMAX = {int(XMAX)}, YMAX = {int(YMAX)}")



#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

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
geocoded = extraire_valeurs_rsc("CNES_CosENU_geo_8rlks.unw.rsc")
geomap = extraire_valeurs_rsc("CNES_Lut_geo_b2_8rlks.unw.rsc")

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

X_FIRST_geomap = geomap["X_FIRST"]
X_STEP_geomap = geomap["X_STEP"]
Y_FIRST_geomap = geomap["Y_FIRST"]
Y_STEP_geomap = geomap["Y_STEP"]
WIDTH_geomap = geomap["WIDTH"]

# Calcule des décalages
XMIN = (X_FIRST_geocoded - X_FIRST_geomap) / X_STEP_geomap
YMIN = (Y_FIRST_geocoded - Y_FIRST_geomap) / Y_STEP_geomap

XMAX = ((X_FIRST_geocoded + WIDTH_geocoded * X_STEP_geocoded) - X_FIRST_geomap) / X_STEP_geomap
YMAX = ((Y_FIRST_geocoded + WIDTH_geocoded * Y_STEP_geocoded) - Y_FIRST_geomap) / Y_STEP_geomap
#print(Y_FIRST_geocoded)
#print(Y_FIRST_geomap)
#print(Y_FIRST_geocoded + WIDTH_geocoded * Y_STEP_geocoded)
#print(Y_FIRST_geomap + WIDTH_geomap * Y_STEP_geomap)

print(f"\nDécalages calculés : XMIN = {int(XMIN)}, YMIN = {int(YMIN)}, XMAX = {int(XMAX)}, YMAX = {int(YMAX)}")



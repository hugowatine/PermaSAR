#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""
batch_checkcoeff.py
--------------------
Script pour recuperer tout les coeffient de proportionalité avec le model dans tous les dossiers int_YYYYMMDD_YYYYMMDD.

Usage:
    batch_checkcoeff.py --path=<path> [--output=<path>] [--file_suffix=<value>] [--plot]
    batch_checkcoeff.py -h | --help

Options:
    --path=<path>           Répertoire contenant les dossiers à traiter (obligatoire)
    --output=<path>         Fichier résumé en sortie (défaut : path/coeffs_summary.txt)
    --file_suffix=<value>   Suffixe éventuel pour le fichier .stack D1-D2suffix.stack (défaut : D1-D2.stack)
    --plot                  Affiche un graphique du coeff en fonction de la baseline
    -h --help               Affiche ce message d'aide.
"""

import os
import re
import subprocess
import sys
import matplotlib.pyplot as plt
from docopt import docopt
from datetime import datetime
import matplotlib.dates as mdates
import numpy as np
from matplotlib.colors import LogNorm

def main():
    arguments = docopt(__doc__)
    base_dir = arguments["--path"]
    do_plot = arguments.get("--plot", False)

    if arguments["--file_suffix"] == None:
        suffix = None
    else : 
        suffix = arguments["--file_suffix"]

    if arguments["--output"] is None:
        output_file = os.path.join(base_dir, "coeffs_summary.txt")
    else:
        output_file = arguments["--output"]

    if base_dir is None:
        print("❌ Erreur : l'option --path est obligatoire.")
        sys.exit(1)

    # base_dir = os.path.normpath(os.path.join(os.getcwd(), base_dir))

    # if not os.path.isdir(base_dir):
    #     print(f"❌ Erreur : le chemin '{base_dir}' n'existe pas ou n'est pas un dossier.")
    #     sys.exit(1)

    pattern_dir = re.compile(r"^int_\d{8}_\d{8}$")

    if suffix == None:
        pattern_stack = re.compile(r"^(\d{8})-(\d{8})\.stack$")
    else:
        pattern_stack = re.compile(r"^(\d{8})-(\d{8})_testlin\.stack$")

    summary_lines = []
    baseline_days = []
    coeff_values = []
    date1_list = []

    print(f"\n🔍 Recherche des dossiers dans : {base_dir}\n")

    # Filtrer les dossiers valides
    all_subdirs = sorted([
        d for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d)) and pattern_dir.match(d)
    ])

    total = len(all_subdirs)

    if total == 0:
        print("❌ Aucun dossier 'int_YYYYMMDD_YYYYMMDD' trouvé dans ce répertoire.")
        sys.exit(1)

    for idx, subdir in enumerate(all_subdirs, start=1):
        full_path = os.path.join(base_dir, subdir)
        print(f"\n📁 Traitement du dossier {idx}/{total} : {subdir}")
        date1, date2 = subdir.split('_')[1], subdir.split('_')[2] 

        found_stack = False
        for f in os.listdir(full_path):
            #print("DEBUG fichier trouvé :", f)
            stack_match = pattern_stack.match(f)
            if stack_match:
                date1, date2 = stack_match.groups()
                stack_file = os.path.join(full_path, f)
                try:
                    with open(stack_file, "r") as file:
                        value = float(file.read().strip())
                    summary_lines.append(f"{date1}\t{date2}\t{value}\n")
                    found_stack = True
                    if do_plot:
                        d1 = datetime.strptime(date1, "%Y%m%d")
                        d2 = datetime.strptime(date2, "%Y%m%d")
                        baseline_days.append((d2 - d1).days)
                        coeff_values.append(value)
                        date1_list.append(date1)
                except Exception as e:
                    print(f"  ⚠️  Erreur lecture fichier {f} : {e}")

        if not found_stack:
            summary_lines.append(f"{date1}\t{date2}\tNaN\n")
            print(f"  ⚠️  Aucun fichier .stack trouvé dans {subdir}")

    # Écriture du fichier résumé
    with open(output_file, "w") as out_f:
        #out_f.write("Date1\tDate2\tCoeff\n")
        out_f.writelines(summary_lines)

    print(f"\n✅ Extraction terminée pour {total} dossier(s). Résultats sauvegardés dans '{output_file}'.\n")
    
    if do_plot:
        from scipy.signal import medfilt
        if len(baseline_days) == 0:
            print("⚠️  Aucun coefficient valide pour tracer le graphique.")
        else:
            dates1 = []
            for d in date1_list:
                try:
                    parsed_date = datetime.strptime(d, "%Y%m%d")
                    same_year_date = parsed_date.replace(year=2000)  # année factice
                    dates1.append(same_year_date)
                except Exception:
                    continue
            #dates1 = [datetime.strptime(d, "%Y%m%d") for d in date1_list]  # date1_list à créer
            
            coeffs = [abs(c) for c in coeff_values]
            #coeffs = [c for c in coeff_values]

            sorted_data = sorted(zip(dates1, coeffs, baseline_days), key=lambda x: x[0])
            dates1, coeffs, baselines = zip(*sorted_data)

            if len(coeffs) >= 3:
                coeffs_filtered = medfilt(coeffs, kernel_size=31)
            else:
                coeffs_filtered = coeffs
            vmin = np.percentile(baseline_days, 2)
            vmax = np.percentile(baseline_days, 98)
            norm = LogNorm(vmin=max(vmin, 1e-6), vmax=vmax)
            plt.figure(figsize=(10,5))
            sc = plt.scatter(dates1, coeffs, c=baseline_days, cmap='gist_ncar', edgecolor='k')
            plt.plot(dates1, coeffs_filtered, 'r-')
            plt.xlabel("Mois")
            plt.ylabel("Coefficient")
            plt.title("Coefficient en fonction du mois")
            plt.grid(True)
            plt.colorbar(sc, label="Baseline temporelle (jours)")
        
            # Formater l'axe des X pour des dates lisibles
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            plt.gcf().autofmt_xdate()  # rotation des dates
            plt.show()
if __name__ == "__main__":
    main()

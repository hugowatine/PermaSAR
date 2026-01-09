#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""
batch_unw_nodata0.py
--------------------
Script pour lancer automatiquement unw_nodata0.py sur tous les dossiers int_YYYYMMDD_YYYYMMDD.

Usage:
    batch_unw_nodata0.py --path=<path> [--nproc=<n>] [--list=<file>] [--prefix=<prefix>] [--suffix=<suffix>]
    batch_unw_nodata0.py -h | --help

Options:
    --path=<path>    Répertoire contenant les dossiers à traiter (obligatoire)
    --list=<file>    Fichier texte contenant les paires à traiter (col1=YYYYMMDD col2=YYYYMMDD)
    --prefix=<prefix> Préfixe des fichiers IFG [default: CNES_InU_geo]
    --suffix=<suffix>   Suffixe des fichiers IFG [default: _era_8rlks.tiff]
    --nproc=<n>        Nombre de processus parallèles [default: 4]
    -h --help        Affiche ce message d'aide.
"""

import os
import re
import subprocess
import sys
from docopt import docopt
from concurrent.futures import ThreadPoolExecutor, as_completed

def read_list_file(list_file):
    """Lit le fichier list et retourne une liste des dossiers à traiter."""
    pairs = []
    with open(list_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                date1, date2 = parts[0], parts[1]
                if re.match(r"^\d{8}$", date1) and re.match(r"^\d{8}$", date2):
                    pairs.append(f"int_{date1}_{date2}")
    return pairs

def process_subdir(unw_script, base_dir, subdir, ifg_start, ifg_end):
    """Traite un sous-dossier donné."""
    full_path = os.path.join(base_dir, subdir)
    print(f"\n📁 Traitement du dossier : {subdir}")

    ifg = None
    for f in os.listdir(full_path):
        if f.startswith(ifg_start) and f.endswith(ifg_end):
            ifg = os.path.join(full_path, f)

    if ifg:
        cmd = [
            "python3",
            unw_script,
            f"--unw={ifg}",
        ]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Erreur dans {subdir} : {e}")
    else:
        print(f"  ⚠️  Fichiers ifg manquants dans {subdir}, traitement ignoré.")

def main():
    arguments = docopt(__doc__)

    base_dir = arguments["--path"]
    list_file = arguments["--list"]
    ifg_start = arguments["--prefix"] or "CNES_InU_geo"
    ifg_end = arguments["--suffix"] or "_era_8rlks.tiff"
    nproc = int(arguments.get("--nproc", 4))

    if base_dir is None:
        print("❌ Erreur : l'option --path est obligatoire.")
        sys.exit(1)

    base_dir = os.path.abspath(base_dir)

    if not os.path.isdir(base_dir):
        print(f"❌ Erreur : le chemin '{base_dir}' n'existe pas ou n'est pas un dossier.")
        sys.exit(1)

    unw_script = os.path.expanduser("~/PermaSAR/Tools/unw_nodata0.py")
    pattern = re.compile(r"^int_\d{8}_\d{8}$")

    print(f"\n🔍 Recherche des dossiers dans : {base_dir}\n")

    if list_file:
        if not os.path.isfile(list_file):
            print(f"❌ Erreur : le fichier '{list_file}' n'existe pas.")
            sys.exit(1)
        all_subdirs = read_list_file(list_file)
    else:
        all_subdirs = sorted([
            d for d in os.listdir(base_dir)
            if os.path.isdir(os.path.join(base_dir, d)) and pattern.match(d)
        ])

    total = len(all_subdirs)

    if total == 0:
        print("❌ Aucun dossier 'int_YYYYMMDD_YYYYMMDD' trouvé dans ce répertoire.")
        sys.exit(1)

    # Exécution parallèle
    with ThreadPoolExecutor(max_workers=nproc) as executor:
        futures = [
            executor.submit(process_subdir, unw_script, base_dir, subdir, ifg_start, ifg_end)
            for subdir in all_subdirs
        ]

        for future in as_completed(futures):
            future.result()

    print(f"\n✅ Traitement terminé pour {total} dossier(s).\n")

if __name__ == "__main__":
    main()

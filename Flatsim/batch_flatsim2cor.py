#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""
batch_flatsim2cor.py
--------------------
Script pour lancer automatiquement flatsim2cor.py sur tous les dossiers int_YYYYMMDD_YYYYMMDD.

Usage:
    batch_flatsim2cor.py --path=<path> [--suffix=<value>] [--prefix=<value>][--nproc=<n>] [--list_int=<file>]
    batch_flatsim2cor.py -h | --help

Options:
    --path=<path>      Répertoire contenant les dossiers à traiter (obligatoire)
    --list_int=<file>  Fichier texte contenant "date1 date2" par ligne (optionnel)
    --prefix=<value>   Prefix of the data at to be process $prefix$date1-$date2$suffix$ 
    --suffix=<value>   Suffix of the data at the starting of the processes $prefix$date1-$date2$suffix$ 
    --nproc=<n>        Nombre de processus parallèles [default: 4]
    -h --help          Affiche ce message d'aide.
"""

import os
import re
import subprocess
import sys
import shutil
from docopt import docopt
from concurrent.futures import ThreadPoolExecutor, as_completed


def process_subdir(flatsim_script, base_dir, subdir, prefix, suffix):
    full_path = os.path.join(base_dir, subdir)
    print(f"\n📁 Traitement du dossier : {subdir}")

    coh = None
    for f in os.listdir(full_path):
        if f.startswith(prefix) and f.endswith(suffix):
            coh = os.path.join(full_path, f)

    if coh:
        cmd = [
            "python3",
            flatsim_script,
            f"--coh={coh}"
        ]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            #subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Erreur dans {subdir} : {e}")
    else:
        print(f"  ⚠️  Fichier coh manquant dans {subdir}, traitement ignoré.")


def read_list_int(file_path):
    """Lit le fichier list_int et retourne un set des dossiers à traiter."""
    wanted = set()
    with open(file_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            d1, d2 = parts[0], parts[1]
            folder = f"int_{d1}_{d2}"
            wanted.add(folder)
    return wanted


def main():
    arguments = docopt(__doc__)

    base_dir = arguments["--path"]
    list_int_file = arguments.get("--list_int", None)
    nproc = int(arguments.get("--nproc", 4))
    if arguments['--prefix'] == None:
        prefix = 'CNES_Coh_geo'
    else:
        prefix = arguments['--prfix']
    if arguments['--suffix'] == None:
        suffix = '_8rlks.tiff'
    else:
        suffix = arguments['--suffix']

    if base_dir is None:
        print("❌ Erreur : l'option --path est obligatoire.")
        sys.exit(1)

    base_dir = os.path.abspath(base_dir)

    if not os.path.isdir(base_dir):
        print(f"❌ Erreur : le chemin '{base_dir}' n'existe pas ou n'est pas un dossier.")
        sys.exit(1)

    path_dirs = os.environ["PATH"].split(":")
    flatsim_script = shutil.which("flatsim2cor.py")
    pattern = re.compile(r"^int_\d{8}_\d{8}$")

    print(f"\n🔍 Recherche des dossiers dans : {base_dir}\n")

    # Dossiers existants
    all_subdirs = sorted([
        d for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d)) and pattern.match(d)
    ])

    if list_int_file:
        print(f"📄 Filtrage via la liste : {list_int_file}")
        wanted = read_list_int(list_int_file)

        # On garde uniquement les dossiers présents dans wanted
        filtered_subdirs = [d for d in all_subdirs if d in wanted]

        missing = wanted - set(all_subdirs)
        if missing:
            print("⚠️  Les dossiers suivants sont demandés mais inexistants :")
            for m in sorted(missing):
                print(f"   - {m}")

        subdirs_to_process = filtered_subdirs
    else:
        subdirs_to_process = all_subdirs

    total = len(subdirs_to_process)
    if total == 0:
        print("❌ Aucun dossier à traiter.")
        sys.exit(1)

    # Exécution parallèle
    with ThreadPoolExecutor(max_workers=nproc) as executor:
        futures = [
            executor.submit(process_subdir, flatsim_script, base_dir, subdir, prefix, suffix)
            for subdir in subdirs_to_process
        ]

        for future in as_completed(futures):
            future.result()

    print(f"\n✅ Traitement terminé pour {total} dossier(s).\n")


if __name__ == "__main__":
    main()


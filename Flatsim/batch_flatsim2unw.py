#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""
batch_flatsim2unw.py
--------------------
Script pour lancer automatiquement flatsim2unw.py sur tous les dossiers int_YYYYMMDD_YYYYMMDD.

Usage:
    batch_flatsim2unw.py --path=<path> [--prefix=<prefix>] [--suffix=<suffix>] [--nproc=<n>] [--list_int=<file>]
    batch_flatsim2unw.py -h | --help

Options:
    --path=<path>      Répertoire contenant les dossiers à traiter (obligatoire)
    --prefix=<prefix>  Préfixe des fichiers IFG [default: CNES_InU_geo]
    --suffix=<suffix>  Suffixe des fichiers IFG [default: _era_8rlks.tiff]
    --list_int=<file>  Fichier texte contenant "date1 date2" par ligne (optionnel)
    --nproc=<n>        Nombre de processus parallèles [default: 4]
    -h --help          Affiche ce message d'aide.
"""

import os
import re
import subprocess
import sys
from docopt import docopt
from concurrent.futures import ThreadPoolExecutor, as_completed


def process_subdir(flatsim_script, base_dir, subdir, ifg_start, ifg_end):
    full_path = os.path.join(base_dir, subdir)
    print(f"\n📁 Traitement du dossier : {subdir}")

    ifg = None
    coh = None
    for f in os.listdir(full_path):
        if f.startswith(ifg_start) and f.endswith(ifg_end):
            ifg = os.path.join(full_path, f)
        elif f.startswith("CNES_Coh_geo") and f.endswith("mask_8rlks.tiff"):
            coh = os.path.join(full_path, f)

    if ifg and coh:
        cmd = [
            "python3",
            flatsim_script,
            f"--ifg={ifg}",
            f"--coh={coh}"
        ]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Erreur dans {subdir} : {e}")
    else:
        print(f"  ⚠️  Fichiers ifg/coh manquants dans {subdir}, traitement ignoré.")


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
    ifg_start = arguments.get("--prefix") or "CNES_InU_geo"
    ifg_end = arguments.get("--suffix") or "_era_8rlks.tiff"
    nproc = int(arguments.get("--nproc", 4))

    if base_dir is None:
        print("❌ Erreur : l'option --path est obligatoire.")
        sys.exit(1)

    base_dir = os.path.abspath(base_dir)
    if not os.path.isdir(base_dir):
        print(f"❌ Erreur : le chemin '{base_dir}' n'existe pas ou n'est pas un dossier.")
        sys.exit(1)

    flatsim_script = os.path.expanduser("~/PermaSAR/Tools/flatsim2unw.py")
    pattern = re.compile(r"^int_\d{8}_\d{8}$")

    print(f"\n🔍 Recherche des dossiers dans : {base_dir}\n")

    # Dossiers existants
    all_subdirs = sorted([
        d for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d)) and pattern.match(d)
    ])

    if list_int_file:
        if not os.path.isfile(list_int_file):
            print(f"❌ Erreur : le fichier '{list_int_file}' n'existe pas.")
            sys.exit(1)

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
            executor.submit(process_subdir, flatsim_script, base_dir, subdir, ifg_start, ifg_end)
            for subdir in subdirs_to_process
        ]
        for future in as_completed(futures):
            future.result()

    print(f"\n✅ Traitement terminé pour {total} dossier(s).\n")


if __name__ == "__main__":
    main()




# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# ############################################
# # Author        : Hugo WATINE (CRPG)
# ############################################

# """
# batch_flatsim2unw.py
# --------------------
# Script pour lancer automatiquement flatsim2unw.py sur tous les dossiers int_YYYYMMDD_YYYYMMDD.

# Usage:
#     batch_flatsim2unw.py --path=<path> [--prefix=<prefix>] [--suffix=<suffix>] [--nproc=<n>]
#     batch_flatsim2unw.py -h | --help

# Options:
#     --path=<path>    Répertoire contenant les dossiers à traiter (obligatoire)
#     --prefix=<prefix> Préfixe des fichiers IFG [default: CNES_InU_geo]
#     --suffix=<suffix> Suffixe des fichiers IFG [default: _era_8rlks.tiff]
#     --nproc=<n>      Nombre de processus parallèles [default: 4]
#     -h --help        Affiche ce message d'aide.
# """

# import os
# import re
# import subprocess
# import sys
# from docopt import docopt
# from concurrent.futures import ThreadPoolExecutor, as_completed

# def read_list_file(list_file):
#     pairs = []
#     with open(list_file, "r") as f:
#         for line in f:
#             parts = line.strip().split()
#             if len(parts) >= 2:
#                 date1, date2 = parts[0], parts[1]
#                 if re.match(r"^\d{8}$", date1) and re.match(r"^\d{8}$", date2):
#                     pairs.append(f"int_{date1}_{date2}")
#     return pairs

# def process_subdir(flatsim_script, base_dir, subdir, ifg_start, ifg_end):
#     full_path = os.path.join(base_dir, subdir)
#     print(f"\n📁 Traitement du dossier : {subdir}")

#     ifg = None
#     coh = None
#     for f in os.listdir(full_path):
#         if f.startswith(ifg_start) and f.endswith(ifg_end):
#             ifg = os.path.join(full_path, f)
#         elif f.startswith("CNES_Coh_geo") and f.endswith("mask_8rlks.tiff"):
#             coh = os.path.join(full_path, f)

#     if ifg and coh:
#         cmd = [
#             "python3",
#             flatsim_script,
#             f"--ifg={ifg}",
#             f"--coh={coh}"
#         ]
#         # Exécution silencieuse
#         try:
#             subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
#         except subprocess.CalledProcessError as e:
#             print(f"  ❌ Erreur dans {subdir} : {e}")
#     else:
#         print(f"  ⚠️  Fichiers ifg/coh manquants dans {subdir}, traitement ignoré.")


# def main():
#     arguments = docopt(__doc__)

#     base_dir = arguments["--path"]
#     list_file = arguments.get("--list")
#     ifg_start = arguments.get("--prefix") or "CNES_InU_geo"
#     ifg_end   = arguments.get("--suffix") or "_era_8rlks.tiff"
#     nproc = int(arguments.get("--nproc", 4))

#     if base_dir is None:
#         print("❌ Erreur : l'option --path est obligatoire.")
#         sys.exit(1)

#     base_dir = os.path.abspath(base_dir)
#     if not os.path.isdir(base_dir):
#         print(f"❌ Erreur : le chemin '{base_dir}' n'existe pas ou n'est pas un dossier.")
#         sys.exit(1)

#     flatsim_script = os.path.expanduser("~/PermaSAR/Tools/flatsim2unw.py")
#     pattern = re.compile(r"^int_\d{8}_\d{8}$")

#     print(f"\n🔍 Recherche des dossiers dans : {base_dir}\n")

#     if list_file:
#         if not os.path.isfile(list_file):
#             print(f"❌ Erreur : le fichier '{list_file}' n'existe pas.")
#             sys.exit(1)
#         all_subdirs = read_list_file(list_file)
#     else:
#         all_subdirs = sorted([
#             d for d in os.listdir(base_dir)
#             if os.path.isdir(os.path.join(base_dir, d)) and pattern.match(d)
#         ])

#     total = len(all_subdirs)
#     if total == 0:
#         print("❌ Aucun dossier 'int_YYYYMMDD_YYYYMMDD' trouvé dans ce répertoire.")
#         sys.exit(1)

#     # Parallélisation avec ThreadPoolExecutor
#     with ThreadPoolExecutor(max_workers=nproc) as executor:
#         futures = [
#             executor.submit(process_subdir, flatsim_script, base_dir, subdir, ifg_start, ifg_end)
#             for subdir in all_subdirs
#         ]
#         for future in as_completed(futures):
#             future.result()  # récupère les exceptions

#     print(f"\n✅ Traitement terminé pour {total} dossier(s).\n")


# if __name__ == "__main__":
#     main()



# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# ############################################
# # Author        : Hugo WATINE (CRPG)
# ############################################

# """
# batch_flatsim2unw.py
# --------------------
# Script pour lancer automatiquement flatsim2unw.py sur tous les dossiers int_YYYYMMDD_YYYYMMDD.

# Usage:
#     batch_flatsim2int.py --path=<path> [--list=<file>] [--prefix=<prefix>] [--suffix=<suffix>]
#     batch_flatsim2int.py -h | --help

# Options:
#     --path=<path>    Répertoire contenant les dossiers à traiter (obligatoire)
#     --list=<file>    Fichier texte contenant les paires à traiter (col1=YYYYMMDD col2=YYYYMMDD)
#     --prefix=<prefix> Préfixe des fichiers IFG [default: CNES_InU_geo]
#     --suffix=<suffix>   Suffixe des fichiers IFG [default: _era_8rlks.tiff]
#     -h --help        Affiche ce message d'aide.
# """

# import os
# import re
# import subprocess
# import sys
# from docopt import docopt

# def read_list_file(list_file):
#     pairs = []
#     with open(list_file, "r") as f:
#         for line in f:
#             parts = line.strip().split()
#             if len(parts) >= 2:
#                 date1, date2 = parts[0], parts[1]
#                 if re.match(r"^\d{8}$", date1) and re.match(r"^\d{8}$", date2):
#                     pairs.append(f"int_{date1}_{date2}")
#     return pairs

# def main():
#     arguments = docopt(__doc__)

#     base_dir = arguments["--path"]
#     list_file = arguments["--list"]
#     ifg_start = arguments["--prefix"] or "CNES_InU_geo"
#     ifg_end   = arguments["--suffix"] or "_era_8rlks.tiff"
    
#     if base_dir is None:
#         print("❌ Erreur : l'option --path est obligatoire.")
#         sys.exit(1)

#     base_dir = os.path.abspath(base_dir)

#     if not os.path.isdir(base_dir):
#         print(f"❌ Erreur : le chemin '{base_dir}' n'existe pas ou n'est pas un dossier.")
#         sys.exit(1)

#     flatsim_script = os.path.expanduser("~/PermaSAR/Tools/flatsim2unw.py")
#     pattern = re.compile(r"^int_\d{8}_\d{8}$")

#     print(f"\n🔍 Recherche des dossiers dans : {base_dir}\n")

#     if list_file:
#         # Lire les paires spécifiées
#         if not os.path.isfile(list_file):
#             print(f"❌ Erreur : le fichier '{list_file}' n'existe pas.")
#             sys.exit(1)
#         all_subdirs = read_list_file(list_file)
#     else:
#         # Prendre tous les dossiers disponibles
#         all_subdirs = sorted([
#             d for d in os.listdir(base_dir)
#             if os.path.isdir(os.path.join(base_dir, d)) and pattern.match(d)
#         ])

#     total = len(all_subdirs)

#     if total == 0:
#         print("❌ Aucun dossier 'int_YYYYMMDD_YYYYMMDD' trouvé dans ce répertoire.")
#         sys.exit(1)

#     for idx, subdir in enumerate(all_subdirs, start=1):
#         full_path = os.path.join(base_dir, subdir)
#         print(f"\n📁 Traitement du dossier {idx}/{total} : {subdir}")

#         ifg = None
#         coh = None
#         for f in os.listdir(full_path):
#             if f.startswith(ifg_start) and f.endswith(ifg_end):
#                 ifg = os.path.join(full_path, f)
#             elif f.startswith("CNES_Coh_geo") and f.endswith(".tiff"):
#                 coh = os.path.join(full_path, f)
#         if ifg and coh:
#             cmd = [
#                 "python3",
#                 flatsim_script,
#                 f"--ifg={ifg}",
#                 f"--coh={coh}"
#             ]
#             print("  → Exécution :", " ".join(cmd))
#             try:
#                 subprocess.run(cmd, check=True)
#             except subprocess.CalledProcessError as e:
#                 print(f"  ❌ Erreur dans {subdir} : {e}")
#         else:
#             print(f"  ⚠️  Fichiers ifg/coh manquants dans {subdir}, traitement ignoré.")

#     print(f"\n✅ Traitement terminé pour {total} dossier(s).\n")

# if __name__ == "__main__":
#     main()

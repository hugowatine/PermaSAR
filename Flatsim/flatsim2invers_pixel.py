#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""
flatsim2invers_pixel.py
-------------
Runs a sequence of processing steps to convert FLATSIM outputs into a first time series analysis without correction:

1- Creates a readme.txt
2- Renames the DAUX and INT file to match the master name and INTERFERO
3- Generates a list of IFGs used in FLATSIM
4- Unzips the IFGs
5- Renames the IFG files to the format int_date1_date2
6- Applies a mask to the CNES_InU IFGs using a shapefile
7- Removes linear ramps (r + az) and referencing the IFGs
8- Launches `prep_invers_pixel.py`
9- Modifies `input_inv_send.py` (without correction)
10- Runs `invers_pixel`

Usage: flatsim2invers_pixel.py --track_path=<path> [--shapefile_path=<path>] [--step=<value>]

Options:
    -h --help               Show this screen.
    --track_path=<path>     Path to the track folder (e.g., A056_centre)
    --shapefile_path=<path> Path to the shapefile used for cropping [default: None]
    --step=<value>          Processing step at which to start [default: 1]
"""

print("\n\nAuthor: Hugo WATINE\n")

import os
import sys
import glob
import re
import subprocess
from osgeo import gdal, osr, gdalconst

try:
    from nsbas import docopt
except ImportError:
    import docopt

# Read arguments
arguments = docopt.docopt(__doc__)
track_path = arguments["--track_path"]

# Handle shapefile path
shapefile = arguments['--shapefile_path']
if shapefile in [None, "None"]:
    print('No shapefile --> no crop during the process')
    shapefile = None

# Handle step value
step = int(arguments['--step']) if arguments['--step'] else 1

print('-------')
print(f'Starting processing the track {track_path} at step {step}')
print('-------\n')

# Step 1 - Create readme.txt
if step == 1:
    readme_path = os.path.join(track_path, 'readme.txt')
    with open(readme_path, 'w') as f:
        f.write("This file contains information about the FLATSIM processing workflow.\n")
        f.write(f"Track path: {track_path}\n")
        f.write("Author: Hugo Watine (CRPG)\n\n")
        f.write("Step 1: Creation of readme.txt\n\n")
    print(f'1. readme.txt file created at {readme_path}')
    step = 2  # Automatically continue to next step
readme_path = os.path.join(track_path, 'readme.txt')  # Re-define in case step 1 was skipped


# Step 2 - Rename DAUX and INT folder to match master date

if step == 2:
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSBAS_DAUX*')) if os.path.isdir(d)]
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INT')) if os.path.isdir(d)]
    if not dint_dirs:
        print("2. No INT/ directory found")
    else :
        subprocess.run(['mv', dint_dirs[0], track_path + '/INTERFERO'])
        print(f'2. Renamed INT to INTERFERO')

    if not daux_dirs:
        print("2. No NSBAS_DAUX directory found (maybe only a .zip exists). Skipping step 2.")
    else:
        daux = daux_dirs[0]
        proc = os.path.join(daux, 'nsbas.proc')

        if not os.path.exists(proc):
            print(f"2. 'nsbas.proc' not found in {daux}")
            sys.exit(1)

        with open(proc) as f:
            match = re.search(r'SarMasterDir\s*=\s*(\d+)', f.read())
            if not match:
                print("2. No SarMasterDir found in nsbas.proc")
                sys.exit(1)
            date = match.group(1)

        new_path = os.path.join(track_path, date)
        if not os.path.exists(new_path):
            subprocess.run(['mv', daux, new_path])
            print(f'2. Renamed {os.path.basename(daux)} to {date}')
            with open(readme_path, 'a') as f:
                f.write(f"Step 2: Renamed {os.path.basename(daux)} to {date}\n")
                f.write(f"mv {daux} {new_path}\n")
        else:
            print(f'2. Target folder {new_path} already exists, skipping rename.')
    step = 3

# Step 3 - Generates a list of IFGs used in FLATSIM
if step == 3:
    dts_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSBAS_TS*')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    if not dts_dirs or not daux_dirs:
        print("Error no NSBAS_TS or SarMasterDir (DAUX) directory found")
        sys.exit(1)
    else:
        daux = daux_dirs[0] + '/interf_pair.rsc'
        dts = dts_dirs[0] + '/RMSinterfero.txt'

        script_path = os.path.expanduser("~/PermaSAR/Flatsim/diffIGF_flatsimVSNsbas.py")
        try:
            subprocess.run([
                "python3", script_path,
                f"--RMSinterfero={dts}",
                f"--inter_pairs={daux}"
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running comparison script: {e}")
            sys.exit(1)
    
        comparison_file = os.path.join(track_path, "comparaisonIFGs_Flatsim_NSBAS.txt")
        if not os.path.exists(comparison_file):
            print(f"Fichier de comparaison non trouvé : {comparison_file}")
            sys.exit(1)

        # Lire les paires à supprimer
        to_remove = []
        with open(comparison_file, "r") as f:
            lines = f.readlines()
            in_section = False
            for line in lines:
                if "Pairs only in inter_pair.rsc (NSBAS):" in line:
                    in_section = True
                    continue
                if in_section:
                    line = line.strip()
                    if not line or not re.match(r"\d{8}_\d{8}", line):
                        break
                    to_remove.append(tuple(line.split("_")))

        # Construire une copie filtrée du fichier interf_pair.rsc
        daux_copy = os.path.join(track_path, "filtered_inter_pair.rsc")
        with open(daux, "r") as f_in:
            all_lines = f_in.readlines()

        # Format des paires à supprimer : "YYYYMMDD YYYYMMDD"
        remove_set = {f"{start} {end}" for start, end in to_remove}

        filtered_lines = [line for line in all_lines if line.strip() not in remove_set]

        with open(daux_copy, "w") as f_out:
            f_out.writelines(filtered_lines)

        print(f"Fichier filtré écrit : {daux_copy} ({len(all_lines) - len(filtered_lines)} paires supprimées)")

        # Mise à jour du readme
        with open(readme_path, 'a') as f:
            f.write("\nStep 3: Filtrage des paires NSBAS absentes de FLATSIM\n")
            f.write(f"Fichier généré : {daux_copy} ({len(all_lines) - len(filtered_lines)} paires supprimées)\n")

    step =4

# Step 4 - Unzips the IFGs
if step == 4 :
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    if not dint_dirs:
        print("2. No INTERFERO/ directory found")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        print(f"4. Unzipping all files in {dint_dirs[0]}")
        subprocess.run("unzip '*'", shell=True, cwd=dint_dirs[0])
        
        with open(readme_path, 'a') as f:
            f.write("\nStep 4: Unzipped all files in INTERFERO/\n")
            f.write(f"cd {dint_dirs[0]}\nunzip '*'\n")


    








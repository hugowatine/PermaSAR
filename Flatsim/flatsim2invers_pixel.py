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
6- Applies a mask to the CNES_InU IFGs using a shapefile (e.g. WestTibet/TP_shapefile/DBATP_Polygon.shp)
7- Create .int, .cor, .unw
8- Extract the seasonal amplitude from the row flatsim cube and flat it
9- Create a list of ifgs to process from rox flatsim process
10- Remove a model to the .int
11- Change coherence to colinearity 
12- Unwarp .int
13- Add back the model
14- Removes linear ramps (r + az) and referencing the IFGs
15- Launches `prep_invers_pixel.py`
16- Modifies `input_inv_send.py` (without correction)
17- Runs `invers_pixel`
18- plot and save the output of invers_pixel

Usage: flatsim2invers_pixel.py --track_path=<path> --step=<value> [--shapefile_path=<path>] [--list_int=<path>] [--cutfile=<path>]

Options:
    -h --help               Show this screen.
    --track_path=<path>     Path to the track folder (e.g., A056_centre)
    --shapefile_path=<path> Path to the shapefile used for cropping [default: None]
    --list_int=<path>       Path to the list of int to be process [default: filtered_inter_pair.rsc]
    --cutfile=<path>        Path to the cutfile used for unwrapping (default: 'CNES_DEM_geo_8rlks_slope_clean_shifted.hgt')
    --step=<value>          Processing step (e.g. 1,2,3)
"""

print("\n\nAuthor: Hugo WATINE\n")

import os
import sys
import glob
import re
import subprocess
from osgeo import gdal, osr, gdalconst
gdal.UseExceptions()

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

list_int = arguments['--list_int']
if list_int == None:
    list_int = 'filtered_inter_pair.rsc'

# Handle step value
step = [int(s) for s in re.split(r'[,\s]+', arguments['--step'].strip()) if s]

print('-------')
print(f'Processing the track {track_path} with step {step}')
print('-------\n')

# Step 1 - Create readme.txt
if 1 in step:
    readme_path = os.path.join(track_path, 'readme.txt')
    with open(readme_path, 'w') as f:
        f.write("This file contains information about the FLATSIM processing workflow.\n")
        f.write(f"Track path: {track_path}\n")
        f.write("Author: Hugo Watine (CRPG)\n\n")
        f.write("Step 1: Creation of readme.txt\n\n")
    print(f'1. readme.txt file created at {readme_path}')

readme_path = os.path.join(track_path, 'readme.txt')  # Re-define in case step 1 was skipped


# Step 2 - Rename DAUX and INT folder to match master date

if 2 in step:
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

# Step 3 - Generates a list of IFGs used in FLATSIM
if 3 in step:
    dts_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSBAS_TS*')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    if not dts_dirs or not daux_dirs:
        print("Error no NSBAS_TS or SarMasterDir (DAUX) directory found")
        sys.exit(1)
    else:
        daux = daux_dirs[0] + '/interf_pair.rsc'
        dts = dts_dirs[0] + '/RMSinterfero.txt'

        comparaison_file = os.path.join(track_path, "comparaisonIFGs_Flatsim_NSBAS.txt")

        script_path = os.path.expanduser("~/PermaSAR/Flatsim/diffIGF_flatsimVSNsbas.py")
        try:
            subprocess.run([
                "python3", script_path,
                f"--RMSinterfero={dts}",
                f"--inter_pairs={daux}",
                f"--output={comparaison_file}"
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running comparison script: {e}")
            sys.exit(1)
    

        
        if not os.path.exists(comparaison_file):
            print(f"Fichier de comparaison non trouvé : {comparaison_file}")
            sys.exit(1)

        # Lire les paires à supprimer
        to_remove = []
        with open(comparaison_file, "r") as f:
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

# Step 4 - Unzips the IFGs
if 4 in step :
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

# 5- Renames the IFG files to the format int_date1_date2
if 5 in step :
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]

    if not dint_dirs:
        print("2. No INTERFERO/ directory found")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        print(f"5. Renames the IFG folders inside {dint} to the format int_date1_date2")
        script_path = os.path.expanduser("/home/hugo/PermaSAR/Flatsim/rename_int_folders.sh")
        subprocess.run([script_path, dint], check=True)

        with open(readme_path, 'a') as f:
            f.write("\nStep 5: Renamed all IFG folders in INTERFERO/\n")
            f.write(f"{script_path} {dint}\n")

# compter le nombre de fihcier .zip ls -1 *.zip | wc -l
# compter les dossier int_date1_date2 ls -d int_* 2>/dev/null | wc -l
# MEttre les zip problématique dans un fichier txt 
# for z in *.zip; do d="int_$(echo "$z" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{4}-[0-9]{2}-[0-9]{2}' | sed 's/-//g')";     [ ! -d "$d" ] && echo "$z"; done > zip_sans_dossier.txt
# -> différences = besoin de dezzpiper 

# 6- Applies a mask to the CNES_InU IFGs using a shapefile #--shapefile=./TP_shapefile/DBATP_Polygon.shp
if 6 in step :
    if shapefile == None:
        print('6. No input shapefile, exit')
        sys.exit(1)

    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]

    if not dint_dirs:
        print("6. No INTERFERO/ directory found")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        script_path = os.path.expanduser("/home/hugo/PermaSAR/Flatsim/folder_maskshp.sh")
        nproc = 15

        print(f"6. Applies a mask to the CNES_InW IFGs using a shapefile")
        file2mask = "CNES_InW_geo_*_sd_era_8rlks.tiff"
        subprocess.run([script_path, dint, shapefile, file2mask, str(nproc)], check=True)
        with open(readme_path, 'a') as f:
            f.write("\nStep 6: Applies a mask to the CNES_InW IFGs and CNES_Coh using a shapefile\n")
            f.write(f"{script_path} {dint} {shapefile} {file2mask} {nproc}\n")

        print(f"6. Applies a mask to the CNES_Coh IFGs using a shapefile")
        file2mask = "CNES_Coh_geo_*_8rlks.tiff"
        subprocess.run([script_path, dint, shapefile, file2mask, str(nproc)], check=True)
        with open(readme_path, 'a') as f:
            f.write(f"{script_path} {dint} {shapefile} {file2mask} {nproc}\n")

        # Pour trouver les sous dossier ou ca ne marche pas : "find . -type d -print0 | while IFS= read -r -d '' dir; do     if ! find "$dir" -maxdepth 1 -type f -name "CNES_InW_geo_*_*_sd_era_mask_8rlks.tiff" | grep -q .; then         echo "$dir";     fi; done"
        # find . -maxdepth 1 -type d -name "int_*" | while read -r dir; do     if ! find "$dir" -maxdepth 1 -type f -name "CNES_InW_geo_*_*_sd_era_mask_8rlks.tiff" | grep -q .; then         echo "Dossier sans fichier CNES_InW_geo_*_*_sd_era_mask_8rlks.tiff : $dir";     fi; done

# 7- Create .int, .cor, .unw
if 7 in step :
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    if not dint_dirs:
        print("7. No INTERFERO/ directory found")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        nproc = 20

        try:
             subprocess.run(["python3", '/home/hugo/PermaSAR/Flatsim/batch_flatsim2int.py', f"--path={dint}", f"--nproc={nproc}"], check=True)
        except subprocess.CalledProcessError as e:
             print(f"Error running comparison script: {e}")
             sys.exit(1)
        try:
            subprocess.run(["python3", '/home/hugo/PermaSAR/Flatsim/batch_flatsim2cor.py', f"--path={dint}", f"--nproc={nproc}"], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running comparison script: {e}")
            sys.exit(1)
        try:
             subprocess.run(["python3", '/home/hugo/PermaSAR/Flatsim/batch_flatsim2unw.py', f"--path={dint}", f"--nproc={nproc}"], check=True)
        except subprocess.CalledProcessError as e:
             print(f"Error running comparison script: {e}")
             sys.exit(1)

        for d in glob.glob(os.path.join(dint, 'int_*_*')):
            if not os.path.isdir(d):
                continue
            for pattern in ['*_mask_*.cor', '*_mask_*.cor.rsc']:
                for f in glob.glob(os.path.join(d, pattern)):
                    if not os.path.exists(f):
                        continue
                    new_name = f.replace('_mask', '')
                    os.rename(f, new_name)
        
        with open(readme_path, 'a') as f:
            f.write("\nStep 7: Create .int, .cor, .unw/\n")
            f.write(f"~/PermaSAR/Flatsim/batch_flatsim2int.py --path={dint}\n")
            f.write(f"~/PermaSAR/Flatsim/batch_flatsim2cor.py --path={dint}\n")
            f.write(f"~/PermaSAR/Flatsim/batch_flatsim2unw.py --path={dint}\n")
# Tej les _mask du .cor et .cor.rsc
# for d in int_*_*; do [ -d "$d" ] || continue; for f in "$d"/*_mask_*.cor; do [ -e "$f" ] || continue; mv "$f" "${f/_mask/}"; done; done
# for d in int_*_*; do [ -d "$d" ] || continue; for f in "$d"/*_mask_*.cor.rsc; do [ -e "$f" ] || continue; mv "$f" "${f/_mask/}"; done; done

# for d in int_*_*; do [ -d "$d" ] || continue; for f in "$d"/*.cor; do [ -e "$f" ] || continue; mv "$f" "$d/coh_$(basename "$f")"; done; done
# for d in int_*_*; do [ -d "$d" ] || continue; for f in "$d"/*.cor.rsc; do [ -e "$f" ] || continue; mv "$f" "$d/coh_$(basename "$f")"; done; done
        
#for f in filtSWc_*_mask_*rlks.int; do [ -e "$f" ] || continue cp "$f" "${f/_mask_/_mask_testwmedian_}" done;

# 8- Extract the seasonal amplitude from the row flatsim cube and flat it

if 8 in step :
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    drow_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSB*')) if os.path.isdir(d)]

    if not daux_dirs:
        print("8. No auxiliary directory found")
        sys.exit(1)
    elif not drow_dirs:
        print("8. No ROW NSBAS directory found")
        sys.exit(1)
    else : 
        daux = daux_dirs[0]
        drow = drow_dirs[0]

        print(f"8. Extract the seasonal amplitude from {drow} and flat it")

        baseline_rsc = os.path.join(daux, "baseline.rsc")
        list_images_txt = os.path.join(daux, "list_images.txt")
        rmsdate_txt = os.path.join(drow, "RMSdate.txt")
        inaps_txt = os.path.join(drow, "inaps.txt")

        try:
            subprocess.run(f"awk '{{print 0, $1, 0, $5, 0, $2}}' {baseline_rsc} > {list_images_txt}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running comparison script: {e}")
            sys.exit(1)
        try:
            subprocess.run(f"cp {list_images_txt} {drow}/. ", shell=True, check=True)
            list_images_txt = os.path.join(drow, "list_images.txt")
        except subprocess.CalledProcessError as e:
            print(f"Error cpoy {list_images_txt}: {e}")
            sys.exit(1)
        try:
            subprocess.run(f"awk 'NR==FNR {{a[$2]; next}} ($2 in a)' {rmsdate_txt} {list_images_txt} > tmp_list_images.txt && mv tmp_list_images.txt {list_images_txt}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error filtering list_images.txt: {e}")
            sys.exit(1)
        try:
            subprocess.run(f"awk '{{print $3}}' {rmsdate_txt} > {inaps_txt}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error creating inaps.txt: {e}")
            sys.exit(1)

        with open(readme_path, 'a') as f:
            f.write("\nStep 8: Extract the seasonal amplitude from {drow}\n")
            f.write("Aux folder: awk '{print 0, $1, 0, $5, 0, $2}' baseline.rsc > list_images.txt\n")
            f.write("NSBAS folder: cp ../aux/list_images.txt .\n")
            f.write("NSBAS folder: awk 'NR==FNR {a[$2]; next} !($2 in a) {print $2}' RMSdate.txt list_images.txt\n")
            f.write("--> supprimer certaines lignes dans list_images.txt\n")
            f.write("awk '{ print $3 }' RMSdate.txt > inaps.txt\n")
            f.write("invers_disp2coef.py --cube=CNES_DTs_geo_8rlks.tiff --linear=yes --seasonal=yes --niter=2 --aps=inaps.txt --plot=no --list_images=list_images.txt --nproc=15\n")

        cube = "CNES_DTs_geo_8rlks.tiff"
        list_images_txt = "list_images.txt"
        inaps_txt = "inaps.txt"

        try:
            subprocess.run(f"invers_disp2coef.py --cube={cube} --linear=yes --seasonal=yes --niter=2 --aps={inaps_txt} --plot=no --list_images={list_images_txt} --nproc=15", shell=True, check=True, cwd=drow)
        except subprocess.CalledProcessError as e:
            print(f"Error runing invers_disp2coef.py : {e}")
            sys.exit(1)
        
        #find . -maxdepth 1 -type f -newermt "2025-12-04 16:58:00" -exec cp {} ./A012_sud/NSBAS_TS-PKG_S1_TIBET-HIM-A012SUD-VV-2014-2022_IW123_2015-06-03_2022-05-21/ \;

# 9 - Modify nsbas.proc to nsbas_col.proc for interferogram processing
if 9 in step:
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]

    if not daux_dirs:
        print("9. No auxiliary directory found.")
        sys.exit(1)

    else:
        daux = daux_dirs[0]
        print(f"9. Modifying nsbas.proc to nsbas_col.proc for interferogram processing...")
        nsbas_proc = os.path.join(daux, "nsbas.proc")
        nsbas_col_proc = os.path.join(daux, "nsbas_col.proc")

        try:
            subprocess.run(f"cp {nsbas_proc} {nsbas_col_proc}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error copying {nsbas_proc}: {e}")
            sys.exit(1)

        print(f"9. Please open CNES_Net_geo_8rlks.tiff (band 4) and CNES_TCoh_radar_2rlks to found a new referecence site (200x200) and seed for unwarping (1x1) ")
        # scp hugo@193.54.29.81:/data2/FLATSIM/WestTibet/D063_sud/NSBAS_TS-PKG_S1_TIBET-HIM-D063SUD-VV-2014-2022_IW123_2014-10-22_2022-05-25/CNES_Net_geo_8rlks.tiff /Users/hugowatine/Desktop/PHD/These/TIBET/WestTibet/D063_sud/

        with open(nsbas_col_proc, "r+") as f:
            lines = f.readlines()
            f.seek(0)
            for line in lines:
                stripped = line.strip()
                if stripped.startswith("Rlooks_int"):
                    f.write("Rlooks_int = 8\n")
                elif stripped.startswith("filterstyle"):
                    f.write("filterstyle = SWc\n")
                elif stripped.startswith("SWamplim"):
                    f.write("SWamplim = 0.05\n")
                elif stripped.startswith("unw_seedx"):
                    seedx = input(f"Enter the new value of seedx for geo unwarping (current line: {stripped}): ")
                    f.write(f"seedx = {seedx}\n")
                elif stripped.startswith("unw_seedy"):
                    seedy = input(f"Enter the new value of seedy for geo unwarping (current line: {stripped}): ")
                    f.write(f"seedy = {seedy}\n")
                elif stripped.startswith("RefLeft"):
                    ref_left = input(f"Enter the new value of RefLeft for geo referencing e.g. ~like seedx (current line: {stripped}): ")
                    f.write(f"RefLeft = {ref_left}\n")
                elif stripped.startswith("RefTop"):
                    ref_top = input(f"Enter the new value of RefTop for geo referencing e.g. ~like seedy (current line: {stripped}): ")
                    f.write(f"RefTop = {ref_top}\n")
                else:
                    f.write(line)
            f.truncate()

        with open(readme_path, 'a') as f:
            f.write("\nStep 9: Modify nsbas.proc to nsbas_col.proc for IFG processing\n")
            f.write("Copied nsbas.proc → nsbas_col.proc\n")
            f.write("Updated parameters:\n")
            f.write(" - Rlooks_int = 8\n")
            f.write(" - filterstyle = SWc\n")
            f.write(" - SWamplim = 0.05\n")
            f.write(f" - unw_seedx = {seedx}\n")
            f.write(f" - unw_seedy = {seedy}\n")
            f.write(f" - RefLeft = {ref_left}\n")
            f.write(f" - RefTop = {ref_top}\n")
        
# 10- filter the .int with coherence
if 10 in step:
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    
    if not dint_dirs:
        print("10. No INTERFERO/ directory found")
        sys.exit(1)
    elif not daux_dirs:
        print("10. No auxiliary directory found.")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        daux = daux_dirs[0]
        nsbas_col_proc = os.path.relpath(os.path.join(daux, 'nsbas_col.proc'), start=track_path)
        list_int = list_int
        job = 'filterSW'
        prefix = 'CNES_InW_geo_'
        suffix = '_sd_era_mask'
        #suffix = '_sd_era'
        nproc=20

        try:
            subprocess.run(f"/home/hugo/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} -f", shell=True, check=True, cwd=track_path)
        except subprocess.CalledProcessError as e:
            print(f"Error runing nsb_filtflatunw.py : {e}")
            sys.exit(1)

        with open(readme_path, 'a') as f:
            f.write("\nStep 10: filter the .int\n")
            f.write(f"~/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc}\n")

# 11- remove model 
if 11 in step:
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    drow_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSB*')) if os.path.isdir(d)]

    if not dint_dirs:
        print("11. No INTERFERO/ directory found")
        sys.exit(1)
    elif not daux_dirs:
        print("11. No auxiliary directory found.")
        sys.exit(1)
    elif not drow_dirs:
        print("11. No ROW NSBAS directory found")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        daux = daux_dirs[0]
        drow = drow_dirs[0]
        nsbas_col_proc = os.path.relpath(os.path.join(daux, 'nsbas_col.proc'), start=track_path)
        list_int = list_int
        job = 'flat_model'
        prefix = 'CNES_InW_geo_'
        suffix = '_sd_era_mask'
        nproc=10
        #model = os.path.relpath(os.path.join(drow, 'stack_1001_0401_lowRMS_mask_HP158_interpolate_LP2_inverted.tiff'), start=track_path)
        model = os.path.relpath(os.path.join(drow, 'ampwt_coeff_mask_filter400_interpolate_filter2_quad.tiff'), start=track_path)
        #model = os.path.relpath(os.path.join(drow, 'lin_model_filter400_interpolate_filter2_quad.tif'), start=track_path)

        try:
            subprocess.run(f"/home/hugo/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} --model={model} -f", shell=True, check=True, cwd=track_path, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing nsb_filtflatunw.py : {e}")
            sys.exit(1)

        with open(readme_path, 'a') as f:
            f.write("\nStep 11: remove model \n")
            f.write(f"~/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} --model={model}\n")
# /home/hugo/PermaSAR/Tools/flatten_stack.py CNES_InW_geo_D1-D2_sd_era_mask_8rlks.int filtSWc_CNES_InW_geo_D1-D2_sd_era_mask_8rlks.int /data2/FLATSIM/WestTibet/A085_sud/NSBAS_TS-PKG_S1_TIBET-HIM-A085SUD-VV-2014-2022_IW123_2015-05-03_2022-05-26/ampwt_coeff_mask_filter400_interpolate_filter2_quad.tiff --nreg=60 --thresh_amp=0.1 --thresh_cohreg=0.6 --thresh_model=0.4
            
# 12- Generate time series of the coefficients between the model and the interferograms.
if 12 in step:
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]

    if not dint_dirs:
        print("12. No INTERFERO/ directory found")
        sys.exit(1)
    elif not daux_dirs:
        print("12. No auxiliary directory found.")
        sys.exit(1)
    else : 
        dint = os.path.abspath(dint_dirs[0])
        daux = daux_dirs[0]
        baseline = daux + '/baseline.rsc'
        input = dint + '/coeffs_summary.txt'
        output = dint + '/coeffs_ts.txt'

        try:
            subprocess.run(f"/home/hugo/PermaSAR/Flatsim/batch_flatsimcheckcoeff.py --path={dint}", shell=True, check=True, cwd=track_path)
        except subprocess.CalledProcessError as e:
            print(f"Error runing batch_flatsimcheckcoeff.py : {e}")
            sys.exit(1)
        
        try:
            subprocess.run(f"/home/hugo/PermaSAR/TimeSeries/invert_phi.py --datesfile={baseline} --input={input} --output={output} --rad2mm=-1 --plot", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error runing batch_flatsimcheckcoeff.py : {e}")
            sys.exit(1)


        with open(readme_path, 'a') as f:
            f.write("\nStep 12: Generate time series of the coefficients between the model and the interferograms \n")
            f.write(f"~/PermaSAR/Flatsim/batch_flatsimcheckcoeff.py --path={dint}\n")
            f.write(f"~/PermaSAR/TimeSeries/invert_phi.py --datesfile={baseline} --input={input} --output={output} --plot=yes\n")

# 13- Change coherence to colinearity
if 13 in step:
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    drow_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSB*')) if os.path.isdir(d)]

    if not dint_dirs:
        print("13. No INTERFERO/ directory found")
        sys.exit(1)
    elif not daux_dirs:
        print("13. No auxiliary directory found.")
        sys.exit(1)
    elif not drow_dirs:
        print("13. No ROW NSBAS directory found")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        daux = daux_dirs[0]
        drow = drow_dirs[0]
        nsbas_col_proc = os.path.relpath(os.path.join(daux, 'nsbas_col.proc'), start=track_path)
        list_int = list_int
        job = 'colin'
        prefix = 'CNES_InW_geo_'
        suffix = '_sd_era_mask_nomodel'
        nproc=10

        try:
            subprocess.run(f"/home/hugo/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} -f", shell=True, check=True, cwd=track_path, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing nsb_filtflatunw.py : {e}")
            sys.exit(1)

        with open(readme_path, 'a') as f:
            f.write("\nStep 13: Change coherence to colinearity \n")
            f.write(f"~PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} -f\n")

            
# 14- filter the .int with colinearity
if 14 in step:
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    
    if not dint_dirs:
        print("14. No INTERFERO/ directory found")
        sys.exit(1)
    elif not daux_dirs:
        print("14. No auxiliary directory found.")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        daux = daux_dirs[0]
        nsbas_col_proc = os.path.relpath(os.path.join(daux, 'nsbas_col.proc'), start=track_path)
        list_int = list_int
        job = 'filterSW'
        prefix = 'col_'
        suffix = '_sd_era_mask_nomodel'
        #suffix = '_sd_era_mask'
        nproc=20

        # try:
        #     subprocess.run(["python3", '/home/hugo/PermaSAR/Flatsim/batch_col2cor.py', f"--path={dint}", f"--nproc={nproc}"], check=True)
        # except subprocess.CalledProcessError as e:
        #     print(f"Error running comparison script: {e}")
        #     sys.exit(1)

        try:
            subprocess.run(f"/home/hugo/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} -f", shell=True, check=True, cwd=track_path)
        except subprocess.CalledProcessError as e:
            print(f"Error runing nsb_filtflatunw.py : {e}")
            sys.exit(1)

        with open(readme_path, 'a') as f:
            f.write("\nStep 14: filter the .int with colinearity\n")
            f.write(f"~PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} -f\n")

# 15- Generate cutfile
            


# 15- Unwarp .int
if 15 in step:
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    drow_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSB*')) if os.path.isdir(d)]
    
    if not dint_dirs:
        print("15. No INTERFERO/ directory found")
        sys.exit(1)
    elif not daux_dirs:
        print("15. No auxiliary directory found.")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        daux = daux_dirs[0]
        drow = drow_dirs[0]
        nsbas_col_proc = os.path.relpath(os.path.join(daux, 'nsbas_col.proc'), start=track_path)
        if arguments['--cutfile'] == None:
            #cutfile = os.path.relpath(os.path.join(daux, 'CNES_DEM_geo_8rlks_slope_clean_shifted.hgt'), start=track_path)
            cutfile = os.path.relpath(os.path.join(drow, 'proxy.hgt'), start=track_path)
        else:
            cutfile = arguments['--cutfile']
        list_int = list_int
        job = 'unwrapping'
        prefix = 'col_'
        suffix = '_sd_era_mask_nomodel'
        #suffix = '_sd_era_mask'
        #suffix = '_sd_era_nomodel'
        nproc=20

        try:
            subprocess.run(f"/home/hugo/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} --cutfile={cutfile} -f", shell=True, check=True, cwd=track_path)
        except subprocess.CalledProcessError as e:
            print(f"Error runing nsb_filtflatunw.py : {e}")
            sys.exit(1)

        with open(readme_path, 'a') as f:
            f.write("\nStep 15: Unwarp .int\n")
            f.write(f"~/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} --cutfile={cutfile} -f\n")

# 16- Add back the model
if 16 in step:
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    drow_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSB*')) if os.path.isdir(d)]

    if not dint_dirs:
        print("16. No INTERFERO/ directory found")
        sys.exit(1)
    elif not daux_dirs:
        print("16. No auxiliary directory found.")
        sys.exit(1)
    elif not drow_dirs:
        print("16. No ROW NSBAS directory found")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        daux = daux_dirs[0]
        drow = drow_dirs[0]
        nsbas_col_proc = os.path.relpath(os.path.join(daux, 'nsbas_col.proc'), start=track_path)
        list_int = list_int
        job = 'add_model_back'
        prefix = 'col_'
        suffix = '_sd_era_mask_nomodel'
        #suffix = '_sd_era_nomodel'
        nproc=20
        model = os.path.relpath(os.path.join(drow, 'ampwt_coeff_mask_filter400_interpolate_filter2_quad.tiff'), start=track_path)
        #model = os.path.relpath(os.path.join(drow, 'ampwt_coeff_mask_filter400_interpolate_filter2.tiff'), start=track_path)

        try:
            subprocess.run(f"/home/hugo/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} --model={model} -f", shell=True, check=True, cwd=track_path, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing nsb_filtflatunw.py : {e}")
            sys.exit(1)

        with open(readme_path, 'a') as f:
            f.write("\nStep 16: Add back the model \n")
            f.write(f"~/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix={prefix} --suffix={suffix} --jobs={job} --list_int={list_int} {nsbas_col_proc} --nproc={nproc} --model={model} -f\n")

# 17- Removes linear ramps (r + az) and referencing the IFGs
if 17 in step:
    dint_dirs = [d for d in glob.glob(os.path.join(track_path, 'INTERFERO')) if os.path.isdir(d)]
    daux_dirs = [d for d in glob.glob(os.path.join(track_path, '2*')) if os.path.isdir(d)]
    drow_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSB*')) if os.path.isdir(d)]

    if not dint_dirs:
        print("17. No INTERFERO/ directory found")
        sys.exit(1)
    elif not daux_dirs:
        print("17. No auxiliary directory found.")
        sys.exit(1)
    elif not drow_dirs:
        print("17. No ROW NSBAS directory found")
        sys.exit(1)
    else : 
        dint = dint_dirs[0]
        daux = daux_dirs[0]
        drow = drow_dirs[0]
        nsbas_col_proc = os.path.join(daux, 'nsbas_col.proc')
        list_int = 'filtered_inter_pair.rsc'
        prefix = 'filt_col_'
        #suffix = '_sd_era_mask_cohcol'
        suffix = '_sd_era_mask'
        flat = 3
        nproc=10
        estim = 'yes'
        suffix_output = '_corr'
        cohpixel = 'yes'
        threshold_coh = 0.3
        tsinv = 'yes'
        ref = os.path.relpath(os.path.join(daux, 'CNES_DEM_geo_8rlks_slope_clean.unw'), start=track_path)
        perc = 95
        rlook = 8

        with open(nsbas_col_proc) as f:
            vals = {k: int(v) for k, v in re.findall(r'(RefLeft|RefTop|RefWidth|RefLength)\s*=\s*(\d+)', f.read())}

        ref_zone = (vals['RefTop'], vals['RefTop'] + vals['RefWidth'], vals['RefLeft'], vals['RefLeft'] + vals['RefLength'])
        ref_zone = f"{ref_zone[0]},{ref_zone[1]},{ref_zone[2]},{ref_zone[3]}"

        print('ref_zone : ', ref_zone)

        try:
            subprocess.run(f"/home/hugo/PermaSAR/Flatsim/flatsim_invert_ramp_topo_unw.py --ref_zone={ref_zone} --int_list={list_int} --int_path=./INTERFERO/ --prefix={prefix} --suffix={suffix} --flat={flat} --tsinv={tsinv} --estim={estim} --perc={perc} --suffix_output={suffix_output} --ref={ref} --cohpixel={cohpixel} --threshold_coh={threshold_coh} --rlook={rlook} --nproc={nproc}", shell=True, check=True, cwd=track_path)
        except subprocess.CalledProcessError as e:
            print(f"Error runing nsb_filtflatunw.py : {e}")
            sys.exit(1)

        with open(readme_path, 'a') as f:
            f.write("\nStep 17: Removes linear ramps (r + az) and referencing the IFGs \n")
            f.write(f"~/PermaSAR/Flatsim/flatsim_invert_ramp_topo_unw.py --ref_zone={ref_zone} --int_list={list_int} --int_path={dint} --prefix={prefix} --suffix={suffix} --flat={flat} --tsinv={tsinv} --estim={estim} --perc={perc} --suffix_output={suffix_output} --ref={ref} --cohpixel={cohpixel} --threshold_coh={threshold_coh} --rlook={rlook} --nproc={nproc}\n")

# 18- Generate the mask
if 18 in step:
    drow_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSB*')) if os.path.isdir(d)]
    dts_dirs = [d for d in glob.glob(os.path.join(track_path, 'ts_corr_bt_bpAPS')) if os.path.isdir(d)]

    if not dts_dirs:
        print("17. No ts_corr_bt_bpAPS/ directory found")
        sys.exit(1)
    elif not drow_dirs:
        print("17. No ROW NSBAS directory found")
        sys.exit(1)
    else : 
        drow = drow_dirs[0]
        dts = dts_dirs[0]
        
        track_name = dts.split('/')[-2]

        # RMS
        print("generate RMSpixel.tiff")
        try:
            subprocess.run(f"~/PermaSAR/Tools/r4totif.py --infile=RMSpixel --outfile={track_name}_RMSpixel.tiff --ref_file=../{drow.split('/')[-1]}/CNES_MV-LOS_geo_8rlks.tiff", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
            subprocess.run(f"gdal_edit.py -a_nodata 0.0 {track_name}_RMSpixel.tiff", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing r4totif.py : {e}")
            sys.exit(1)
        
        # biais
        print("generate CNES_Net_geo_8rlks_bias_mask.tiff")
        try:
            subprocess.run(f"gdal_translate -b 5 CNES_Net_geo_8rlks.tiff {track_name}_CNES_Net_geo_8rlks_bias.tiff", shell=True, check=True, cwd=drow, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing gdal_translate : {e}")
            sys.exit(1)
        try:
            subprocess.run(f"~/PermaSAR/Tools/raster_maskshp.py --raster={track_name}_CNES_Net_geo_8rlks_bias.tiff --shapefile=../../TP_shapefile/DBATP_Polygon.shp --output={track_name}_CNES_Net_geo_8rlks_bias_mask.tiff", shell=True, check=True, cwd=drow, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing raster_maskshp : {e}")
            sys.exit(1)
        try:
            subprocess.run(f"cp {track_name}_CNES_Net_geo_8rlks_bias_mask.tiff ../ts_corr_bt_bpAPS/", shell=True, check=True, cwd=drow)
        except subprocess.CalledProcessError as e:
            print(f"Error runing cp : {e}")
            sys.exit(1)
        # gdal_calc.py -A CNES_Net_geo_8rlks_bias_mask.tiff --outfile=CNES_Net_geo_8rlks_bias_mask_mmyr.tiff --calc="(A*365*(-4.413824938174363)/24).astype('float32')" --overwrite
        # gdal_calc.py -A CNES_Net_geo_8rlks_bias.tiff --outfile=CNES_Net_geo_8rlks_bias_mask.tiff --calc="((A <= 0.015) & (A >= -0.015)).astype('uint8')" --overwrite

        # Mean miscolsure
        print("generate mean_misclosure.tiff")
        try:
            subprocess.run(f"~/PermaSAR/Flatsim/flatsim_compute_mean_misclosure.py --output={track_name}_mean_misclosure --weight", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing gdal_translate : {e}")
            sys.exit(1)

        # bt moyenne
        print("generate average_bt.tiff")
        try:
            subprocess.run("~/PermaSAR/Tools/count_unw_pixel.py", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing count_unw_pixel.py : {e}")
            sys.exit(1)
        try:
            subprocess.run(f"~/PermaSAR/Tools/r4totif.py --infile=average_bt.r4 --outfile={track_name}_average_bt.tiff --ref_file=../{drow.split('/')[-1]}/CNES_MV-LOS_geo_8rlks.tiff", shell=True, check=True, cwd=dts,  stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error runing r4totif.py : {e}")
            sys.exit(1)

        ref = gdal.Open(f"{dts}/{track_name}_RMSpixel.tiff")
        target = gdal.Open(f"{dts}/{track_name}_mean_misclosure.tiff", gdal.GA_Update) 
        target.SetProjection(ref.GetProjection())
        target.SetGeoTransform(ref.GetGeoTransform())
        target = None
        ref = None
        # gdal_calc.py -A average_bt.tiff --outfile=average_bt_mask.tiff --calc="(A > 0.4).astype('uint8')" --overwrite
        
        # Telechargement 
        # scp hugo@193.54.29.81:/data2/FLATSIM/WestTibet/D121_sud/ts_corr_bt_bpAPS/D121_sud_*.tiff /Users/hugowatine/Desktop/PHD/These/TIBET/WestTibet/D121_sud/

        # with open(readme_path, 'a') as f:
        #     f.write("\nStep 17: Removes linear ramps (r + az) and referencing the IFGs \n")
        #     f.write(f"~PermaSAR/Flatsim/flatsim_invert_ramp_topo_unw.py --ref_zone={ref_zone} --int_list={list_int} --int_path={dint} --prefix={prefix} --suffix={suffix} --flat={flat} --tsinv={tsinv} --estim={estim} --perc={perc} --suffix_output={suffix_output} --ref={ref} --cohpixel={cohpixel} --threshold_coh={threshold_coh} --rlook={rlook} --nproc={nproc}\n")

if 19 in step:
    drow_dirs = [d for d in glob.glob(os.path.join(track_path, 'NSB*')) if os.path.isdir(d)]
    dts_dirs = [d for d in glob.glob(os.path.join(track_path, 'ts_corr_bt_bpAPS')) if os.path.isdir(d)]

    if not dts_dirs:
        print("17. No ts_corr_bt_bpAPS/ directory found")
        sys.exit(1)
    elif not drow_dirs:
        print("17. No ROW NSBAS directory found")
        sys.exit(1)
    else : 
        drow = drow_dirs[0]
        dts = dts_dirs[0]
        
        track_name = dts.split('/')[-2]

        # RMS
        print("generate amp_model.tiff")
        try:
            subprocess.run(f"~/PermaSAR/Tools/r4totif.py --infile=amp_modele --outfile={track_name}_amp_modele.tiff --ref_file=../{drow.split('/')[-1]}/CNES_MV-LOS_geo_8rlks.tiff", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
            subprocess.run(f"gdal_edit.py -a_nodata 0.0 {track_name}_amp_modele.tiff", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
            subprocess.run(f"gdal_calc.py -A {track_name}_amp_modele.tiff --calc='A*4.413824938174363' --type=Float32 --outfile={track_name}_amp_modele_mm.tiff --overwrite", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error generating amp_model.tiff : {e}")
            sys.exit(1)

        print("generate slope_corr_var.tiff")
        try:
            subprocess.run(f"~/PermaSAR/Tools/r4totif.py --infile=slope_corr_var --outfile={track_name}_slope_corr_var.tiff --ref_file=../{drow.split('/')[-1]}/CNES_MV-LOS_geo_8rlks.tiff", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
            subprocess.run(f"gdal_edit.py -a_nodata 0.0 {track_name}_slope_corr_var.tiff", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
            subprocess.run(f"gdal_calc.py -A {track_name}_slope_corr_var.tiff --calc='A*-4.413824938174363' --type=Float32 --outfile={track_name}_slope_corr_var_mmyr.tiff --overwrite", shell=True, check=True, cwd=dts, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Error generating slope_corr_var.tiff : {e}")
            sys.exit(1)
        # scp hugo@193.54.29.81:"/data2/FLATSIM/WestTibet/*/ts_corr_bt_bpAPS/*mm*.tiff"  //Users/hugowatine/Desktop/PHD/These/TIBET/WestTibet/

#../batch_unw_nodata0.py --path=./INTERFERO/ --list=filtered_inter_pair.rsc --prefix=filt_col_ --suffix=_sd_era_mask_corr_8rlks.unw
#prep_invers_pixel.py --int_path=INTERFERO/ --outputdir=./ts_nocorr --int_list=filtered_inter_pair.rsc --dates_list=./20190516/baseline.rsc --prefix=filt_col_ --suffix=_sd_era_mask_corr_8rlks
#vi input_inv_send
#0     #   mask pixels with large RMS misclosure  (y=0;n=1)
#4    #  threshold for the mask on RMS misclosure (in same unit as input files)
#4      #  range and azimuth downsampling (every n pixel)
#0      #  iterations to correct unwrapping errors (y:nb_of_iterations,n:0)
#3      #  iterations to weight pixels of interferograms with large residual? (y:nb_of_iterations,n:0)
# invers_pixel < input_inv_send
#../../quality_card_int.py
            
# RELANCER TRAITEMENT RAPIDEMENT
# ~/PermaSAR/Flatsim/flatsim_nsb_filtflatunw.py --prefix=CNES_InW_geo_ --suffix=_sd_era_mask_testwmedian_nomodel --jobs=colin/filterSW/unwrapping/add_model_back --model=NSBAS_TS-PKG_S1_TIBET-HIM-A158CENTRE-VV-2014-2022_IW123_2014-10-16_2022-05-31/ampwt_coeff_mask_filter400_interpolate_filter2_quad.tiff --list_int=filtered_inter_pair.rsc 20190417/nsbas_col.proc --nproc=20


#nsb_filtflatunw.py --prefix=CNES_InW_geo_ --suffix=_sd_era --jobs=filterSW --list_int=filtered_inter_pair.rsc ./20190407/nsbas_col.p
#roc --nproc=10
# 10- Remove a model to the .int
#if 9 in step:


# 11- Change coherence to colinearity 
# 12- Unwarp .int
# 13- Add back the model
# 14- Removes linear ramps (r + az) and referencing the IFGs
# 15- Launches `prep_invers_pixel.py`
# 16- Modifies `input_inv_send.py` (without correction)
# 17- Runs `invers_pixel`
# 18- plot and save the output of invers_pixel


#awk '$4<0.5 {print $2, $3}' RMSinterfero > ../interf_pair_lowRMS.rsc
# awk '{
#     cmd1 = "date -d " $2 " +%s"
#     cmd1 | getline t1
#     close(cmd1)
#     cmd2 = "date -d " $3 " +%s"
#     cmd2 | getline t2
#     close(cmd2)
#     bt = (t2 - t1) / 86400
#     printf "%s\t%s\t%.6f\t%d\n", $2, $3, $4, bt
# }' RMSinterfero | LC_NUMERIC=C sort -k4,4n -k3,3g > RMSinterfero_bt

# awk 'NR==FNR && $3>0.9 && $4<80 {bad[$1,$2]; next}
#      NR!=FNR && !(($1,$2) in bad)' RMSinterfero_bt ../filtered_inter_pair_lowRMS.rsc > ../filtered_inter_pair_lowRMS_clean.rsc


### Lancement avec correction en poids

# Poid par date 
# awk 'NR==FNR { a[$1] = $2; next }
#     $2 in a { inv = (a[$2] != 0 ? 1/a[$2] : 0); print $2, $4, inv }' ../rms_corr_modif2date.txt list_dates > quality_bp_APS.txt
#



# Final ts
# ~/PermaSAR/Flatsim/flatsim_prep_invers_pixel.py --int_path=INTERFERO/ --outputdir=./ts_corr_bt_bpAPS --int_list=filtered_inter_pair_lowRMS_lowRMS.rsc --dates_list=./20190516/baseline.rsc --prefix=filt_col_ --suffix=_sd_era_mask_corr_8rlks --Bc=0.6,_
# 0      #  Weigthing by image quality (y:0,n:1) ? (then read quality in the list of input images)
# 8      #  iterations to correct unwrapping errors (y:nb_of_iterations,n:0)
# cp rms_corr.txt rms_corr_modif.txt
# vi rms_corr_modif.txt --> Supprimer première ligne 
# ~/PermaSAR/TimeSeries/invert_phi.py --input=rms_corr_modif.txt --output=rms_corr_modif2date.txt --plot --rad2mm=1 --noise=yes --datesfile=./20190720/baseline.rsc
# ~/PermaSAR/Flatsim/flatsim_quality_listdate.py --list_date=list_dates --aps_date=../rms_corr_modif2date.txt --output=list_dates_quality.txt
# cp list_dates_quality.txt list_dates
# invers_pixel < input_inv_send


# preview int or unw  
# gdal_translate -b 2 CNES_DEM_radar_8rlks.tiff CNES_DEM_radar_8rlks_b2.tiff
# ~/PermaSAR/Tools/flatsim2unw.py --ifg=CNES_DEM_radar_8rlks_b2.tiff --coh=CNES_DEM_radar_8rlks_b2.tiff
# nsb_geocode.pl CNES_Lut_geo_8rlks.trans CNES_DEM_radar_8rlks_b2.unw CNES_DEM_geo_8rlks.unw
# gdal_translate -b 1 -of GTiff CNES_DEM_geo_8rlks.unw CNES_DEM_geo_8rlks.tif
#~/PermaSAR/Tools/shift_raster.py CNES_DEM_geo_8rlks.tif --shift_y=1
#~/PermaSAR/Tools/tif2hgt.py --tif=CNES_DEM_geo_8rlks_shifted.tif
# preview_int.py --outputdir=./INT_JPG --radar=./20190407/CNES_DEM_geo_8rlks_shifted.hgt --int_list=newprocessing_inter_pair.rsc --prefix=filtSWc_CNES_InW_geo_ --suffix=_sd_era --rlook=8 --int_path=INTERFERO --nproc=10
# preview_int.py --outputdir=./INT_NOMODEL_JPG --radar=./20190407/CNES_DEM_geo_8rlks_shifted.hgt --int_list=newprocessing_inter_pair.rsc --prefix=filt_col_ --suffix=_sd_era_nomodel --rlook=8 --int_path=INTERFERO --nproc=10
# preview_unw.py --outputdir=./UNW_JPG --radar=./20190511/CNES_DEM_geo_8rlks_shifted.hgt --int_list=filtered_inter_pair.rsc --prefix=filt_col_ --suffix=_sd_era_mask_corr --rlook=8 --int_path=INTERFERO --nproc=10          


# Mask
# ~/PermaSAR/Tools/count_unw_pixel.py
# ~/PermaSAR/Tools/r4totif.py --infile=average_bt.r4 --outfile=average_bt.tiff --ref_file=../NSBAS_TS-PKG_S1_TIBET-HIM-A012SUD-VV-2014-2022_IW123_2015-06-03_2022-05-21/CNES_MV-LOS_geo_8rlks.tiff           
# gdal_calc.py -A average_bt.tiff --outfile=average_bt_mask.tiff --calc="(A > 0.4).astype('uint8')" --overwrite
            
# gdal_translate -b 5 CNES_Net_geo_8rlks.tiff CNES_Net_geo_8rlks_bias.tiff
# raster_maskshp.py --raster=CNES_Net_geo_8rlks_bias.tiff --shapefile=../../TP_shapefile/DBATP_Polygon.shp --output=CNES_Net_geo_8rlks_bias_mask.tiff
# gdal_calc.py -A CNES_Net_geo_8rlks_bias_mask.tiff --outfile=CNES_Net_geo_8rlks_bias_mask_mmyr.tiff --calc="(A*365*(-4.413824938174363)/24).astype('float32')" --overwrite
# cp CNES_Net_geo_8rlks_bias_mask_mmyr.tiff ../ts_corr_bt_bpAPS/
# gdal_calc.py -A CNES_Net_geo_8rlks_bias.tiff --outfile=CNES_Net_geo_8rlks_bias_mask.tiff --calc="((A <= 0.015) & (A >= -0.015)).astype('uint8')" --overwrite


# retirer la tecto

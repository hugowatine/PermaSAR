#!/bin/bash

BASE_DIR=$1
SHAPEFILE=$2
NAME=$3
NPROC=${4:-4}
LIST_FILE=""

# vérifier si un argument --list_int est fourni
for arg in "$@"; do
    [[ $arg == --list_int=* ]] && LIST_FILE="${arg#*=}"
done

[ -z "$BASE_DIR" ] || [ -z "$SHAPEFILE" ] || [ -z "$NAME" ] && { echo "Usage: $0 <base_dir> <shapefile> <name> [nproc] [--list_int=<file>]"; exit 1; }

MASK_SCRIPT="/home/hugo/PermaSAR/Tools/raster_maskshp.py"
declare -A VALID_DIRS
[ -n "$LIST_FILE" ] && while read -r d1 d2 _; do [[ "$d1" =~ ^[0-9]{8}$ && "$d2" =~ ^[0-9]{8}$ ]] && VALID_DIRS["int_${d1}_${d2}"]=1; done < "$LIST_FILE"

find "$BASE_DIR" -type d -name "int_*" | while read -r subdir; do
    dirn=$(basename "$subdir")
    [ -n "$LIST_FILE" ] && [ -z "${VALID_DIRS[$dirn]}" ] && continue
    find "$subdir" -maxdepth 1 -type f -name "$NAME"
done | xargs -I {} -P "$NPROC" bash -c '
raster={}; dirn=$(dirname "$raster"); fname=$(basename "$raster")
out="${fname/_8rlks/_mask_8rlks}"; outpath="$dirn/$out"
echo "→ Masking $raster"
python3 '"$MASK_SCRIPT"' --raster="$raster" --shapefile="'"$SHAPEFILE"'" --output="$outpath"
'

echo "✅ All rasters processed!"





# #!/bin/bash

# # Arguments
# BASE_DIR=$1
# SHAPEFILE=$2
# NAME=$3
# NPROC=${4:-4}   # Nombre de processus parallèles (par défaut 4)
# MASK_SCRIPT_PATH="/home/hugo/PermaSAR/Tools/raster_maskshp.py"

# if [ -z "$BASE_DIR" ] || [ -z "$SHAPEFILE" ] || [ -z "$NAME" ]; then
#   echo "Usage: $0 <base_directory> <shapefile_path> <name_pattern> [nproc]"
#   exit 1
# fi

# echo "Starting batch masking in: $BASE_DIR"
# echo "Using shapefile: $SHAPEFILE"
# echo "To mask: $NAME"
# echo "Using $NPROC parallel processes"

# find "$BASE_DIR" -type d -name "int_*" | while read -r subdir; do
#     find "$subdir" -maxdepth 1 -type f -name "$NAME"
# done | xargs -I {} -P "$NPROC" bash -c '
# raster_path={}
# dirname=$(dirname "$raster_path")
# filename=$(basename "$raster_path")
# output_name="${filename/_8rlks/_mask_8rlks}"
# output_path="$dirname/$output_name"

# python3 '"$MASK_SCRIPT_PATH"' --raster="$raster_path" --shapefile="'"$SHAPEFILE"'" --output="$output_path"
# '

# echo "All rasters processed!"
# find "$BASE_DIR" -type d -name "int_*" | while read -r subdir
#     find "$subdir" -maxdepth 1 -type f -name "$NAME"
# done | xargs -I {} -P "$NPROC" bash -c '

# raster_path="{}"
# dirname=$(dirname "$raster_path")
# filename=$(basename "$raster_path")
# output_name="${filename/_sd_era_8rlks/_sd_era_mask_8rlks}"
# output_path="$dirname/$output_name"

# python3 '"$MASK_SCRIPT_PATH"' --raster="$raster_path" --shapefile="'"$SHAPEFILE"'" --output="$output_path"
# '

# echo "All rasters processed!"



# #!/bin/bash

# # Arguments
# BASE_DIR=$1
# SHAPEFILE=$2
# NAME=$3
# MASK_SCRIPT_PATH="/home/hugo/PermaSAR/Tools/raster_maskshp.py"


# if [ -z "$BASE_DIR" ] || [ -z "$SHAPEFILE" ]; then
#   echo "Usage: $0 <base_directory> <shapefile_path>"
#   exit 1
# fi

# echo "Starting batch masking in: $BASE_DIR"
# echo "Using shapefile: $SHAPEFILE"
# echo "To mask: $NAME"

# # Trouve tous les sous-dossiers int_YYYYMMDD_YYYYMMDD
# find "$BASE_DIR" -type d -name "int_*" | while read -r subdir; do
#   echo "Searching in: $subdir"
  
#   # Trouve tous les fichiers matching NAME
#   find "$subdir" -maxdepth 1 -type f -name $NAME | while read -r raster_path; do
#     #echo "→ Found raster: $raster_path"
    
#     filename=$(basename "$raster_path")
#     dirname=$(dirname "$raster_path")

#     # Nouveau nom de fichier de sortie avec _mask_
#     output_name="${filename/_sd_era_8rlks/_sd_era_mask_8rlks}"
#     #output_path="$dirname/$output_name"
#     output_path="$subdir/$output_name"
    
#     # Lancer le masque
#     #echo "Masking and saving to: $output_path"
#     python3 "$MASK_SCRIPT_PATH" --raster="$raster_path" --shapefile="$SHAPEFILE" --output="$output_path"

#     #echo "Done with: $filename"
#     echo
#   done
# done

# echo "All rasters processed!"


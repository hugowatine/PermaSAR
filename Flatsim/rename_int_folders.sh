#!/bin/bash

# Chemin vers le dossier passé en argument
input_dir="$1"

# Vérifie que le dossier existe
if [ ! -d "$input_dir" ]; then
  echo "Le dossier spécifié n'existe pas."
  exit 1
fi

# Parcours des sous-dossiers
for folder in "$input_dir"/NSBAS_INT-PKG_S1*; do
  # Vérifie que c'est un dossier
  if [ -d "$folder" ]; then
    # Extraction des deux dates
    basename=$(basename "$folder")
    if [[ "$basename" =~ ([0-9]{4})-([0-9]{2})-([0-9]{2})_([0-9]{4})-([0-9]{2})-([0-9]{2}) ]]; then
      date1="${BASH_REMATCH[1]}${BASH_REMATCH[2]}${BASH_REMATCH[3]}"
      date2="${BASH_REMATCH[4]}${BASH_REMATCH[5]}${BASH_REMATCH[6]}"
      newname="int_${date1}_${date2}"
      mv "$folder" "$input_dir/$newname"
      echo "Renommé : $basename -> $newname"
    fi
  fi
done


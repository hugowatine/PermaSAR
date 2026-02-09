#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# Extraction des coordonnées de l'étendue spatiale
xmin=$(gdalinfo TCoh_geo_8rlks.tiff | grep "Lower Left" | awk '{print $4}' | tr -d '(),')
ymin=$(gdalinfo TCoh_geo_8rlks.tiff | grep "Lower Left" | awk '{print $5}' | tr -d '(),')
xmax=$(gdalinfo TCoh_geo_8rlks.tiff | grep "Upper Right" | awk '{print $4}' | tr -d '(),')
ymax=$(gdalinfo TCoh_geo_8rlks.tiff | grep "Upper Right" | awk '{print $5}' | tr -d '(),')

# Extraction de la résolution
xres=$(gdalinfo TCoh_geo_8rlks.tiff | grep "Pixel Size" | awk '{print $4}' | tr -d '(),' | awk -F'-' '{print $2}') 

echo $xmin
echo $ymin 
echo $xmax 
echo $ymax 
echo $xres 
echo $yres
gdalwarp -te $xmin $ymin $xmax $ymax -tr $xres -$xres -r cubicspline -of GTiff NDSI.tif NDSI_reprojected.tif
gdalwarp -te $xmin $ymin $xmax $ymax -tr $xres -$xres -r cubicspline -of GTiff NDWI.tif NDWI_reprojected.tif

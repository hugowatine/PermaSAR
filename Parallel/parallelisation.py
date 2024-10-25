import numpy as np
from osgeo import gdal, osr
import matplotlib.pyplot as plt
import time
import multiprocessing as mp
import sys

#import multiprocessing.shared_memory as shm

# Fonction pour ajuster une régression linéaire
def linear_regression(y):
    global formatted_dates
    """Renvoie la pente (slope) d'une régression linéaire de y par rapport à x."""
    A = np.vstack([formatted_dates, np.ones(len(formatted_dates))]).T
    m, _ = np.linalg.lstsq(A, y, rcond=None)[0]  # m est la pente
    return m

def read_band(band_idx, file_path, line, ncol, num_lines):
    ds = gdal.Open(file_path)  # Chaque processus ouvre son propre dataset
    band = ds.GetRasterBand(band_idx)
    return band.ReadAsArray(0, line, ncol, num_lines)

def init(arr):
    global formatted_dates
    formatted_dates=arr

if __name__ == '__main__':

    # Fichier :
    #fcube = "/Users/hugowatine/Desktop/Hugo/These/TIBET/Script_Hugo/Parallelisation/CNES_DTs_geo_8rlks_cropped.tif"
    #fcube = "/Users/hugowatine/Desktop/Hugo/These/TIBET/Script_Hugo/Parallelisation/CNES_DTs_geo_8rlks_cropped1000x2000.tif"
    fcube = "/Users/hugowatine/Desktop/Hugo/These/TIBET/Script_Hugo/Parallelisation/CNES_DTs_geo_8rlks.tiff"
    fliste_image = "/Users/hugowatine/Desktop/Hugo/These/TIBET/Script_Hugo/Parallelisation/list_images.txt"

    outputf = "/Users/hugowatine/Desktop/Hugo/These/TIBET/Script_Hugo/Parallelisation/slope_raster.tiff"
    num_processes = int(sys.argv[1])
    block_size = int(sys.argv[2])

    # Lecture de liste images:
    dates, idates,rms=np.loadtxt(fliste_image, comments='#', usecols=(1,3,4), unpack=True, dtype='f,f,f')
    formatted_dates = idates


    #Lecture du cube:
    ds = gdal.Open(fcube)
    ncol, nlines = ds.RasterXSize, ds.RasterYSize
    N = ds.RasterCount

    # Initialiser un tableau pour stocker les coefficients linéaires
    slope_array = np.zeros((nlines, ncol), dtype=np.float32)


    # Parcourir chaque pixel
    t1 = time.time()

    with mp.Pool(processes=num_processes, initializer=init, initargs=(formatted_dates,)) as pool:
        for line in range(0, nlines, block_size):
            end_line = min(line + block_size, nlines)
            block_size_local = end_line - line
            if line%20 == 0:
                print(line)
            # Extraire la série temporelle pour chaque pixel
            t2 = time.time()
            #line_time_series = np.array([ds.GetRasterBand(b+1).ReadAsArray(0, line, ncol, 1) for b in range(N)])
            results = [pool.apply_async(read_band, args=(i+1, fcube, line, ncol, end_line - line)) for i in range(N)]
            bands_data = [result.get() for result in results]
            temp_array = np.concatenate([np.expand_dims(band, axis=0) for band in bands_data], axis=0)

            t21 = time.time()
            line_time_series = np.array(temp_array)
            t22 = time.time()
            line_time_series = np.reshape(line_time_series, (N,ncol*block_size_local))

            t3 = time.time()
            line_time_series = np.squeeze(line_time_series)
            t4 = time.time()

            # Diviser le tableau en sous-tableaux selon la troisième dimension
            chunks = np.array_split(line_time_series, num_processes, axis=1)
            t5 = time.time()

            # Créer un pool de processus


                # Appliquer la fonction sur chaque portion en parallèle
            results = pool.map(linear_regression, chunks)
            t6 = time.time()
            
            if line%20 == 0:
                print(t21-t2)
                print(t22-t21)
                print(t3-t22)
                print(t3-t2)
                print(t4-t3)
                print(t5-t4)
                print(t6-t5)

            slope = np.concatenate(results, axis=0)
            slope = np.reshape(slope, (block_size_local, ncol))
            # Stocker la pente dans l'array
            slope_array[line:end_line] = slope

            # Effectuer la régression linéaire sur les dates
            #slope = linear_regression(formatted_dates, pixel_time_series)

    print(f"Elapsed time : {time.time()-t1} s")

    # Écriture d'un nouveau raster avec les coefficients linéaires (slopes)
    driver = gdal.GetDriverByName('GTiff')
    out_raster = driver.Create(outputf, ncol, nlines, 1, gdal.GDT_Float32)

    # Copier les géoréférencements et les projections du raster d'origine
    out_raster.SetGeoTransform(ds.GetGeoTransform())
    out_raster.SetProjection(ds.GetProjection())

    # Écrire les valeurs des pentes dans la première bande du nouveau raster
    out_band = out_raster.GetRasterBand(1)
    out_band.WriteArray(slope_array)

    # Nettoyage
    out_band.FlushCache()
    out_raster = None
    ds = None

    # Affichage du raster de pente
    def plot_raster(file_path):
        # Ouvrir le raster
        ds = gdal.Open(file_path)
        
        # Lire les données de la première bande
        slope_data = ds.GetRasterBand(1).ReadAsArray()
        
        # Créer une figure et afficher le raster
        plt.figure(figsize=(10, 6))
        plt.imshow(slope_data, cmap='RdYlBu', interpolation='none', vmax=np.nanpercentile(slope_data, 98), vmin=np.nanpercentile(slope_data,2))
        plt.colorbar(label="Slope Coefficient")
        plt.title('Slope Raster')
        plt.xlabel('Column Index')
        plt.ylabel('Row Index')
        plt.show()

    # Afficher le raster
    plot_raster(outputf)
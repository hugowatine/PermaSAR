#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
# Author        : Hugo WATINE (CRPG)
############################################

"""\
invers_disp2coef_par_Flatsim.py
-------------
Test décomposition pour produit flatsim dans le cadre du projet de thèse.

Usage: invers_disp2coef_par_Flatsim.py [--cube=<path>] [--lectfile=<path>] [--list_images=<path>] [--aps=<path>] [--rms=<path>]\
[--rmspixel=<path>] [--threshold_rms=<value>] [--linear=<yes/no>] [--seasonal=<yes/no>] [--seasonal_increase=<yes/no>] [--acceleration=<yes/no>] [--imref=<value>]\
[--nproc=<nb_cores>] [--niter=<value>] [--plot=<yes/no>] [--ineq=<yes/no>] [--cond=<value>] [--block_size=<value>]

Options:
-h --help           Show this screen.
--cube=<path>           Path to time series displacements cube file [default: no]
--lectfile=<path>       Path to the lect.in file. Simple text file containing width and length and number of images of the time series cube (output of invers_pixel). By default the program will try to find an .hdr file. [default: lect.in].
--list_images=<path>    Path to list images file. text file made of 5 columns containing for each images 1) number 2) Doppler freq (not read) 3) date in YYYYMMDD format 4) numerical date 5) perpendicular baseline [default: images_retenues].
--aps=<path>            Path to the APS file. Text file with two columns (dates, APS) giving an input APS to each dates. By default, no weigthing if no empirical estimation or misfit from the spatial estimation used as input uncertianties [default: None].
--rms=<path>            Path to the RMS file. Text file with two columns (dates, RMS) giving an input RMS to each dates. By default, no weigthing if no empirical estimation or misfit from the spatial estimation used as input uncertianties [default: None].
--rmspixel=<path>       Path to the RMS map. Map in r4 or tiff format that gives an error for each pixel (e.g RMSpixel, output of invers_pixel) [default: None].
--threshold_rms=<value> Threshold on rmspixel argument for empitical spatial estimations [default: 1.]
--niter=<value>         Number of iterations. At the first iteration, image uncertainties is given by aps file or misfit spatial iteration, while for the next itarations, uncertainties are equals to the global RMS previous temporal decomposition [default: 0].
--linear=<yes/no>       Add a linear function in the inversion [default:yes]
--acceleration=<yes/no> Add an accelration fonction inn the inversion [default: no]
--seasonal=<yes/no>       If yes, add seasonal terms in the decomposition [default: no]
--seasonal_increase=<yes/no>   If yes, add seasonal terms function of time in the inversion
--imref=<value>         Reference image number [default: 1]
--plot=<yes/no>         Display plots [default: no]
--nproc=<nb_cores>        Use <nb_cores> local cores to create delay maps [Default: 4]
--ineq=<yes/no>           If yes, sequential least-square optimisation. If no, SVD inversion with mask on eigenvalues smaller than --cond value. If postseimsic functions, add inequality constraints in the inversion. Use least square results without post-seismic functions as a first guess to iterate the inversion. Then, force postseismic to be the same sign and inferior than steps steps of the first guess [default: yes].
--cond=<value>            Condition value for optimization: Singular value smaller than cond are considered zero [default: 1e-6]
--block_size=<value>       Number of line opend at each steaps for the temporal decomposition  
"""

print()
print('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
print('#'                                                                 '#')
print('#         Linear Inversion of InSAR time series displacements       #')
print('#                  with a decomposition in time                     #')
print('#'                                                                 '#')
print('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
print()

import numpy as np
from osgeo import gdal, osr, gdalconst
import matplotlib.pyplot as plt
import time
import multiprocessing as mp
import sys
from os import path, environ, getcwd
import os
from numpy.lib.stride_tricks import as_strided
import math

try:
    from nsbas import docopt
except:
    import docopt

gdal.UseExceptions()

import logging

logging.basicConfig(level=logging.INFO,\
    format='line %(lineno)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('invers_disp2coef.log')
start_time_1 = time.time()
################################
# Fonction
################################
def checkinfile(file):
    if path.exists(file) is False:
        logger.info("File: {0} not found in {1}, Exit !".format(file,getcwd()))
        sys.exit()

def write_envi_hdr(filename, shape, dtype='float32', interleave='bip'):
    """
    Crée un fichier ENVI .hdr à partir d'un fichier.
    
    Parameters:
    - filename: nom du fichier sans extension (.hdr sera ajouté)
    - shape: tuple (lines, samples, bands)
    - dtype: type numpy ('float32', 'int16', ...)
    - interleave: 'bsq', 'bil' ou 'bip'
    """

    if len(shape) == 3:
        lines, samples, bands = shape
    else:
        lines, samples = shape
        bands = 1
    dtype_map = {
        'uint8': 1,
        'int16': 2,
        'int32': 3,
        'float32': 4,
        'float64': 5,
        'complex64': 6,
        'complex128': 9,
        'uint16': 12,
        'uint32': 13,
        'int64': 14,
        'uint64': 15,
    }

    if dtype not in dtype_map:
        raise ValueError(f"Unsupported data type for ENVI: {dtype}")
    
    hdr_content = f"""ENVI
samples = {samples}
lines   = {lines}
bands   = {bands}
data type = {dtype_map[dtype]}
interleave = {interleave}
byte order = 0
"""

    with open(filename + '.hdr', 'w') as f:
        f.write(hdr_content)

def consInvert(A, b, sigmad, ineq='yes', cond=1e-6, iter=60, acc=5e-4, equality=False):
    """
    Résout Ax ≈ b sous contraintes d'inégalité (et éventuellement d’égalité).

    Retourne :
        fsoln : vecteur de solution
        sigmam : incertitudes (diag de la covariance)
    """
    global indexpo, indexco, indexpofull, pos, indexseas, indexseast
    if A.shape[0] != len(b):
        raise ValueError('Dimensions incompatibles pour A et b')

    if ineq == 'no':
        W = np.diag(1.0 / sigmad)
        fsoln = np.linalg.lstsq(W @ A, W @ b, rcond=cond)[0]
    else:
        # Initialisation
        if indexpo is not None and len(indexpo) > 0:
            Ain = np.delete(A, indexpo, axis=1)
            Win = np.diag(1.0 / np.delete(sigmad, indexpo))
            mtemp = np.linalg.lstsq(Win @ Ain, Win @ b, rcond=cond)[0]

            # Réinsérer les post-sismiques
            for z in range(len(indexpo)):
                mtemp = np.insert(mtemp, indexpo[z], 0.0)
            minit = np.copy(mtemp)
        else:
            W = np.diag(1.0 / sigmad)
            minit = np.linalg.lstsq(W @ A, W @ b, rcond=cond)[0]

        # Définir les bornes
        n = len(minit)
        mmin = -np.ones(n) * np.inf
        mmax = np.ones(n) * np.inf

        if indexpo is not None and indexco is not None:
            for i in range(len(indexco)):
                ico = int(indexco[i])
                ipo = int(indexpofull[i])
                if pos[i] > 0. and minit[ico] > 0.:
                    mmin[ipo], mmax[ipo] = 0, np.inf
                    mmin[ico], mmax[ico] = 0, minit[ico]
                elif pos[i] > 0. and minit[ico] < 0.:
                    mmin[ipo], mmax[ipo] = -np.inf, 0
                    mmin[ico], mmax[ico] = minit[ico], 0

        bounds = list(zip(mmin, mmax))

        # Fonction à minimiser
        def _func(x):
            return np.sum(((A @ x - b) / sigmad) ** 2)

        def _fprime(x):
            return 2 * A.T @ ((A @ x - b) / sigmad**2)

        # Contraintes d'égalité
        if equality:
            def eq_cond(x):
                return (x[indexseast + 1] / x[indexseast]) - (x[indexseas + 1] / x[indexseas])
            res = opt.fmin_slsqp(_func, minit, bounds=bounds, fprime=_fprime,
                                 eqcons=[eq_cond], iter=iter, acc=acc,
                                 full_output=True, iprint=0)
        else:
            res = opt.fmin_slsqp(_func, minit, bounds=bounds, fprime=_fprime,
                                 iter=iter, acc=acc, full_output=True, iprint=0)

        fsoln, fx, its, imode, msg = res
        if imode != 0:
            #logger.warning("SLSQP did not converge: {msg}")
            fsoln = minit  # ou np.full_like(minit, np.nan)

    # Calcul de l'incertitude
    try:
        varx = np.linalg.pinv(A.T @ A)
        res2 = np.sum((b - A @ fsoln) ** 2)
        scale = 1. / (A.shape[0] - A.shape[1])
        sigmam = np.sqrt(scale * res2 * np.diag(varx))
    except np.linalg.LinAlgError:
        sigmam = np.full(A.shape[1], np.nan)

    return fsoln, sigmam

def compute_auto_block_size(new_lines, new_cols, N, dtype='float32', target_memory_MB=200):
    """
    Détermine une taille de bloc (nb de lignes) pour rester dans target_memory_MB
    """
    bytes_per_element = np.dtype(dtype).itemsize
    bytes_per_pixel = N * bytes_per_element
    bytes_per_line = new_cols * bytes_per_pixel

    target_bytes = target_memory_MB * 1024 * 1024
    block_size = target_bytes // bytes_per_line

    # Toujours >= 1 et <= new_lines
    return max(1, min(int(block_size), new_lines))

################################
# Read Arguments
################################

arguments = docopt.docopt(__doc__)

if arguments["--plot"] == 'yes' :
    plot = 'yes'
else:
    plot = 'no'

# Data
if arguments["--cube"] ==  None:
    arguments["--cube"] = "depl_cumule"
if arguments["--lectfile"] ==  None:
    arguments["--lectfile"] = "lect.in"
if arguments["--list_images"] ==  None:
    try:
        checkinfile("images_retenues")
        arguments["--list_images"] = "images_retenues"
    except:
        checkinfile("list_images.txt")
        arguments["--list_images"] = "list_images.txt"

if arguments["--rmspixel"] ==  None:
    arguments["--rmspixel"] = None

# Decomp fct
if arguments["--linear"] ==  None:
    arguments["--linear"] = 'yes'
if arguments["--seasonal"] ==  None:
    arguments["--seasonal"] = 'no'
if arguments["--seasonal_increase"] ==  None:
    arguments["--seasonal_increase"] = 'no'
if arguments["--acceleration"] == None:
    arguments["--acceleration"] = 'no'

if arguments["--niter"] ==  None:
    arguments["--niter"] = 1

if arguments["--ineq"] ==  None:
    arguments["--ineq"] = 'no'

if arguments["--nproc"] ==  None:
    nproc = 5
else:
    nproc = int(arguments["--nproc"])

if arguments["--imref"] ==  None:
    imref = 0
elif int(arguments["--imref"]) < 1:
    logger.warning('--imref must be between 1 and Nimages')
else:
    imref = int(arguments["--imref"]) - 1

if arguments["--cond"] ==  None:
    arguments["--cond"] = 1e-6

################################
# INITIALISATION
################################
print()
print('## Initialisation ##')
print()

# Lecture des dates
checkinfile(arguments["--list_images"])
logger.debug('Load list of dates file: {}'.format(arguments["--list_images"]))

_ , idates, _ =np.loadtxt(arguments["--list_images"], comments='#', usecols=(1,3,4), unpack=True, dtype='f,f,f')

datemin, datemax = int(np.min(idates)), int(np.max(idates)) + 1
dmin = str(int(np.min(idates))) + '0101'
dmax = str(int(np.max(idates))+1) + '0101'

# Lecture des données

checkinfile(arguments["--cube"])  
try:
    ds = gdal.Open(arguments["--cube"])
    if ds is None:
        raise FileNotFoundError
except:
    print(f'.hdr file time series cube {arguments["--cube"]} not found.')
    print(f'Trying to create it using {arguments["--lectfile"]}...')
    # Ouverture de lectfile
    with open(arguments["--lectfile"]) as f:
        lines = f.readlines()
        ncol, nlines = map(int, lines[0].split()[:2])
        N = int(lines[3].strip())
    # Création du fichier .hdr
    write_envi_hdr(arguments["--cube"], shape=(nlines, ncol, N), dtype='float32', interleave='bip')
    print(arguments["--cube"] + '.hdr created')

    ds = gdal.Open(arguments["--cube"])
    if not ds:
        sys.exit(f"ERROR: Unable to open the file {arguments['--cube']}. Please check its format or contents.")

# Lecture des dimensions depuis le fichier ouvert
ncol, nlines, N = ds.RasterXSize, ds.RasterYSize, ds.RasterCount
band_number = list(range(N))

# Affichage des infos
print('\nTime series cube successfully opened.')
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", nlines, 'x', ncol, 'x', N)


# ETAPE DE CROP POSSIBLE ICI 

#######################################################
# Save new lect.in file
#######################################################

fid = open('lect_ts.in','w')
np.savetxt(fid, (ncol,nlines,N),fmt='%6i',newline='\t')
fid.close()


#######################################################
# Create functions of decomposition
######################################################

print()
print('## Create functions of decomposition ##')
print()

class pattern:
    def __init__(self,name,reduction,date):
        self.name=name
        self.reduction=reduction
        self.date=date

    def info(self):
        print(self.name, self.date)

class reference(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
    def g(self,t):
        return np.ones((t.size))

class linear(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=(t-self.to)
        return func

class acceleration(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func= (t-self.to)**2
        return func

class cosvar(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=math.cos(2*math.pi*(t[i]-self.to))
        return func

class sinvar(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=math.sin(2*math.pi*(t[i]-self.to))
        return func

class sint(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=(t[i]-self.to)*math.sin(2*math.pi*(t[i]-self.to))
        return func

class cost(pattern):
    def __init__(self,name,reduction,date):
        pattern.__init__(self,name,reduction,date)
        self.to=date

    def g(self,t):
        func=np.zeros(t.size)
        for i in range(t.size):
            func[i]=(t[i]-self.to)*math.cos(2*math.pi*(t[i]-self.to))
        return func

basis = [
    reference(name='reference', date=datemin, reduction='ref'),
]

index = len(basis)
iteration = False

if arguments["--linear"] == 'yes':
    indexinter=index
    basis.append(linear(name='linear', reduction='lin', date=datemin))
    index += 1

if arguments["--acceleration"] == 'yes':
    indexacc=index
    basis.append(acceleration(name='acceleration', reduction='acc', date=datemin))
    index += 1

if arguments["--seasonal"] == 'yes':
    indexseas = index
    basis.append(cosvar(name='seas. var (cos)', reduction='cos', date=datemin))
    basis.append(sinvar(name='seas. var (sin)', reduction='sin', date=datemin))
    index += 2

if arguments["--seasonal_increase"] == 'yes':
    indexseast = index
    basis.append(cost(name='increased seas. var (cos)', reduction='cost', date=datemin))
    basis.append(sint(name='increased seas. var (sin)', reduction='sint', date=datemin))
    index += 2

equality = False
if arguments["--seasonal_increase"] == 'yes' and arguments["--seasonal"] == 'yes':
    equality = True
    arguments["--ineq"] = 'yes'

Mbasis=len(basis)
print('Number of basis functions: {}'.format(Mbasis))

# Initialisation des matrices à NaN pour chaque fonction de base
for l in range(len(basis)):
    basis[l].m = np.ones((nlines, ncol)) * float('NaN')
    basis[l].sigmam = np.ones((nlines, ncol)) * float('NaN')

#######################################################
# Initialize aps and rms
######################################################
print()
print('## Initialize aps and rms ##')
print()

# aps
if  arguments["--aps"] is None:
    in_aps = np.ones((N), dtype='f') # no weigthing for the first itertion
    fimages =  arguments["--aps"]
    logger.info('APS are set to 1.')
else:
    fimages =  arguments["--aps"]
    try:
        in_aps = np.loadtxt(fimages, unpack=True, comments='#', usecols=(2), dtype='f')
    except:
        logger.warning('APS file is in decrepicated format, requiered two columns text file`')
        in_aps = np.loadtxt(fimages, comments='#', dtype='f')
    logger.info('Input APS: {}'.format(in_aps))
    logger.info('Set very low values to the 2 percentile to avoid overweighting...')
    min_aps= np.nanpercentile(in_aps,2)
    index = np.flatnonzero(in_aps<min_aps)
    in_aps[index] = min_aps
    #in_aps = in_aps[indexd] # A modifier si on supprime des dates

# rms
if  arguments["--rms"] is None:
    in_rms = np.ones((N), dtype='f') # no weigthing for the first itertion
    logger.info('RMS are set to 1.')
else:
    fimages =  arguments["--rms"]
    try:
        in_rms = np.loadtxt(fimages, unpack=True, comments='#', usecols=(2), dtype='f')
    except:
        logger.warning('RMS file is in decrepicated format, requiered two columns text file`')
        in_rms = np.loadtxt(fimages, comments='#', dtype='f')
    logger.info('Input RMS: {}'.format(in_rms))
    logger.info('Set very low values to the 2 percentile to avoid overweighting...')
    min_rms= np.nanpercentile(in_rms,2)
    index = np.flatnonzero(in_rms < min_rms)
    in_rms[index] = min_rms
    #in_rms = in_rms[indexd] # A modifier si on supprime des dates

## initialize input uncertainties
in_sigma = in_rms * in_aps
logger.info('Input uncertainties (aps * rms): {}'.format(in_sigma))

#######################################################
# Time Decomposition
######################################################
print()
print('---------------')
print('Time Decomposition')
print('---------------')
print()
from scipy import optimize as opt

def consInvert(A, b, sigmad, ineq='yes', cond=1e-6, iter=60, acc=5e-4, equality=False):
    """
    Résout Ax ≈ b avec ou sans contraintes d'inégalité, et éventuellement d’égalité.
    
    Entrées :
        - A : matrice design (N x M)
        - b : observations (N,)
        - sigmad : incertitudes sur b (N,)
        - ineq : 'yes' pour activer contraintes (positives par défaut)
        - cond : valeur seuil pour pseudo-inversion
        - equality : booléen, active contrainte d’égalité sur seasonal
    
    Sorties :
        - fsoln : solution estimée
        - sigmam : incertitudes sur les paramètres
    """
    global indexseas, indexseast
    
    if A.shape[0] != len(b):
        raise ValueError('Dimensions incompatibles pour A et b')

    if ineq == 'no':
        W = np.diag(1.0 / sigmad)
        fsoln = np.linalg.lstsq(W @ A, W @ b, rcond=cond)[0]
    else:
        W = np.diag(1.0 / sigmad)
        minit = np.linalg.lstsq(W @ A, W @ b, rcond=cond)[0]
        
        # Définir les bornes
        n = len(minit)

        mmin = -np.ones(n) * np.inf
        mmax = np.ones(n) * np.inf

        bounds = list(zip(mmin, mmax))
        # Fonction à minimiser
        def _func(x):
            return np.sum(((A @ x - b) / sigmad) ** 2)

        def _fprime(x):
            return 2 * A.T @ ((A @ x - b) / sigmad**2)
        
        # Contraintes d'égalité
        if equality:
            def eq_cond(x):
                return (x[indexseast + 1] / x[indexseast]) - (x[indexseas + 1] / x[indexseas])
            res = opt.fmin_slsqp(_func, minit, bounds=bounds, fprime=_fprime,
                                    eqcons=[eq_cond], iter=iter, acc=acc,
                                    full_output=True, iprint=0)
        else:
            res = opt.fmin_slsqp(_func, minit, bounds=bounds, fprime=_fprime,
                                    iter=iter, acc=acc, full_output=True, iprint=0)
            
        fsoln, fx, its, imode, msg = res
        if imode != 0:
            #logger.warning("SLSQP did not converge: {msg}")
            fsoln = minit  # ou np.full_like(minit, np.nan)

    # Calcul de l'incertitude
    try:
        varx = np.linalg.pinv(A.T @ A)
        res2 = np.sum((b - A @ fsoln) ** 2)
        scale = 1. / (A.shape[0] - A.shape[1])
        sigmam = np.sqrt(scale * res2 * np.diag(varx))
    except np.linalg.LinAlgError:
        sigmam = np.full(A.shape[1], np.nan)

    return fsoln, sigmam


def temporal_decomp(disp, dates, sigma, cond, ineq, equality):
    """
    Décomposition temporelle pour une série temporelle donnée.

    Entrées :
        - disp : série temporelle du pixel (longueur N)
        - dates : tableau des dates numériques (longueur N)
        - sigma : incertitudes associées à chaque date (longueur N)
        - cond, ineq, equality : paramètres d’inversion

    Sorties :
        - m : coefficients des fonctions de base
        - sigmam : incertitudes sur les coefficients
        - mdisp : série temporelle reconstruite
    """

    N = len(disp)
    M = len(basis)  # nombre de fonctions de base

    mdisp = np.full(N, np.nan, dtype=np.float32)
    m = np.full(M, np.nan, dtype=np.float32)
    sigmam = np.full(M, np.nan, dtype=np.float32)

    # Indices valides (non NaN)
    k = np.flatnonzero(~np.isnan(disp))
    kk = len(k)

    if kk > N / 6:
        tabx = dates[k]
        taby = disp[k]

        G = np.zeros((kk, M), dtype=np.float32)

        for l in range(M):
            G[:, l] = basis[l].g(tabx)

        # Si tu as des kernels supplémentaires :
        # for l in range(Mker):
        #     G[:, Mbasis + l] = kernels[l].g(tabx)

        m, sigmam = consInvert(G, taby, sigma[k], cond=cond, ineq=ineq, equality=equality)

        # reconstruction du modèle (forward model)
        mdisp[k] = G @ m

    return m, sigmam, mdisp

def paralelisation(args):
    """
    Applique une décomposition temporelle pixel par pixel sur un bloc 3D.

    Entrée:
        args = (data_block, in_sigma, cond, ineq, equality)

    Sortie:
        output_array: 2D numpy array (nlines_chunk, ncols_chunk)
                      contenant un résultat scalaire par pixel (ex: m[0])
    """
    data_block, idates, in_sigma, cond, ineq, equality, M = args
    nband, block_lines, block_cols = data_block.shape

    models = np.zeros((block_lines, block_cols, nband),dtype=np.float32)
    m_array = np.full((block_lines, block_cols, M), np.nan, dtype=np.float32)
    sigmam_array = np.full((block_lines, block_cols, M), np.nan, dtype=np.float32)

    for i in range(block_lines):
        for j in range(block_cols):
            ts = data_block[:, i, j]  # série temporelle du pixel
            if np.sum(~np.isnan(ts)) <= 1:  # facultatif : éviter de traiter ts vides
                continue
            try:
                m, sigmam, models[i,j,:] = temporal_decomp(ts, idates, in_sigma, cond, ineq, equality)
                m_array[i, j] = m
                sigmam_array[i, j] = sigmam
            except Exception as e:
                logger.warning(f"Erreur lors du lancement de la decomposition temporelle pour le pixel ({i}, {j}) : {str(e)}")

                continue

    result = m_array, sigmam_array, models
    return result

M = len(basis)
output_param = np.full((nlines, ncol, M), np.nan, dtype=np.float32)
output_sigparam = np.full((nlines, ncol, M), np.nan, dtype=np.float32)

if arguments["--block_size"] ==  None:
    block_size = compute_auto_block_size(nlines, ncol, N, dtype='float32', target_memory_MB=200)
else :
    block_size = arguments["--block_size"]

logger.info('Block size for parallelisation: {} '.format(block_size))
logger.info('Number of cores for parallelisation: {} '.format(nproc))

band_number=list(range(1, N+1))

# A tej après
end_time_1 = time.time()
print(f"Total time : {end_time_1 - start_time_1:.3f} secondes")

for ii in range(int(arguments["--niter"])):
    print()
    print('---------------')
    print('iteration: {}'.format(ii+1))
    print('---------------')

    output_models = np.zeros((nlines,ncol,N),dtype=np.float32)

    start_time_1 = time.time()
    for start_line in range(0, nlines, block_size):
        start_time_2 = time.time()
        # Ouverture de block_size ligne
        current_block_size = min(block_size, nlines - start_line)
        block_cube = ds.ReadAsArray(0, start_line, ncol, current_block_size, band_list=band_number)
        # remplacement des nodata en np.nan
        mask = (block_cube == 9990) | (block_cube == 9999)
        block_cube[mask] = np.nan
        # Referencement a imref
        cst = np.copy(block_cube[imref, :, :])
        cst[np.isnan(cst)] = 0.0
        for l in range(N):
            block_cube[l, :, :] = block_cube[l, :, :] - cst
            if l != imref:
                mask = (block_cube[l, :, :] == 0.0)
                block_cube[l, :, :][mask] = np.nan

        chunks = np.array_split(block_cube, nproc, axis=1)
        del block_cube
        args_list = [(chunk, idates, in_sigma, arguments["--cond"], arguments["--ineq"], equality, M) for chunk in chunks]
        with mp.Pool(processes=nproc) as pool:
            results = pool.map(paralelisation, args_list)

        line_start = start_line  # ligne de départ dans output_param pour ce bloc global

        for i, (m_array, sigmam_array, model_array) in enumerate(results):
            n_lines_chunk = m_array.shape[0]
            n_cols_chunk = m_array.shape[1]
            
            output_param[line_start:line_start+n_lines_chunk, 0:n_cols_chunk, :] = m_array
            output_sigparam[line_start:line_start+n_lines_chunk, 0:n_cols_chunk, :] = sigmam_array
            output_models[line_start:line_start+n_lines_chunk, 0:n_cols_chunk, :] = model_array
            
            line_start += n_lines_chunk
        end_time_2 = time.time()
        logger.info(f"Lignes {start_line + 1} to {start_line + block_size} done - time : {end_time_2 - start_time_2:.3f} secondes")        
    # data = output_models[:, :, 70]
    # data_for_percentile = np.where(data == 0.00, np.nan, data)
    # vmin = np.nanpercentile(data_for_percentile, 2)
    # vmax = np.nanpercentile(data_for_percentile, 98)

    # plt.imshow(data, cmap='viridis', vmin=vmin, vmax=vmax)
    # plt.colorbar()
    # plt.title('Carte bande 70 avec échelle percentile 2-98')
    # plt.show()
        
    # compute RMSE
    # remove outiliers

    output_models = np.transpose(output_models, (2, 0, 1))
    index = np.logical_or(output_models>9999., output_models<-9999)
    output_models[index] = 0.

    # Open the full initial cube 
    #start_time_3 = time.time()
    cube = ds.ReadAsArray(0, 0, ncol, nlines, band_list=band_number)
    #end_time_3 = time.time()
    #print(f"Temps ouverture complet raster : {end_time_2 - start_time_2:.3f} secondes")

    # Copute the residual
    squared_diff = (np.nan_to_num(cube,nan=0) - np.nan_to_num(output_models, nan=0))**2
    del cube, output_models
    res = np.sqrt(np.nanmean(squared_diff, axis=(1, 2))**2)  

    # remove low res to avoid over-fitting in next iter
    min_res= np.nanpercentile(res,2)
    index = np.flatnonzero(res < min_res)
    res[index] = min_res

    print('Residuals:\n', res)
    #print('Dates      Residuals  ')
    #for l in range(N):
    #    print (idates[l], res[l])

    np.savetxt('aps_{}.txt'.format(ii), res.T, fmt=('%.6f'))
    # set apsf is yes for next iteration
    arguments["--aps"] == 'yes'
    # update aps for next iterations taking into account in_aps and the residues of the last iteration
    in_sigma = res * in_aps * in_rms

end_time_1 = time.time()
print(f"Total time : {end_time_1 - start_time_1:.3f} secondes")

#######################################################
# Save functions in binary file
#######################################################

def save_geotiff(array, outname, ncols, nlines, gt, proj, nodata=-9999.0):
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outname, ncols, nlines, 1, gdal.GDT_Float32)
    band = ds.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(nodata)
    ds.SetGeoTransform(gt)
    ds.SetProjection(proj)
    band.FlushCache()
    del ds

def save_r4(array, outname):
    with open(outname, 'wb') as f:
        array.astype('float32').tofile(f)

# if arguments["--geotiff"] is not None:
#     for l in range((M)):
#         outname = 'test_{}_coeff.tif'.format(basis[l].reduction)
#         logger.info('Save: {}'.format(outname))
#         save_geotiff(basis[l].m, outname, ncol, nlines, gt, proj)

#         outname = 'test_{}_sigcoeff.tif'.format(basis[l].reduction)
#         logger.info('Save: {}'.format(outname))
#         save_geotiff(basis[l].sigmam, outname, ncol, nlines, gt, proj)


for l in range((M)):
    outname = 'test_{}_coeff.r4'.format(basis[l].reduction)
    logger.info('Save: {}'.format(outname))
    save_r4(output_param[:, :, l], outname)

    outname = 'test_{}_sigcoeff.r4'.format(basis[l].reduction)
    logger.info('Save: {}'.format(outname))
    save_r4(output_sigparam[:, :, l], outname)

#######################################################
# Compute Amplitude and phase seasonal
#######################################################
    
if arguments["--seasonal"]  == 'yes':
    cosine = output_param[:, :, indexseas]
    sine = output_param[:, :, indexseas+1]
    amp = np.sqrt(cosine**2+sine**2)
    phi = np.arctan2(sine,cosine)

    sigcosine = output_sigparam[:, :, indexseas]
    sigsine = output_sigparam[:, :, indexseas+1]
    sigamp = np.sqrt(sigcosine**2+sigsine**2)
    sigphi = (sigcosine*abs(sine)+sigsine*abs(cosine))/(sigcosine**2+sigsine**2)

    logger.info('Save: {}'.format('ampwt_coeff.r4'))
    save_r4(amp, 'ampwt_coeff.r4')
    logger.info('Save: {}'.format('ampwt_sigcoeff.r4'))
    save_r4(sigamp, 'ampwt_sigcoeff.r4')
    logger.info('Save: {}'.format('phiwt_coeff.r4'))
    save_r4(phi, 'phipwt_coeff.r4')
    logger.info('Save: {}'.format('phiwt_sigcoeff.r4'))
    save_r4(phi, 'phipwt_sigcoeff.r4')

if arguments["--seasonal_increase"]  == 'yes':
    cosine = output_param[:, :, indexseast+2]
    sine = output_param[:, :, indexseast+3]
    amp = np.sqrt(cosine**2+sine**2)
    phi = np.arctan2(sine,cosine)

    sigcosine = output_sigparam[:, :, indexseast+2]
    sigsine = output_sigparam[:, :, indexseast+3]
    sigamp = np.sqrt(sigcosine**2+sigsine**2)
    sigphi = (sigcosine*abs(sine)+sigsine*abs(cosine))/(sigcosine**2+sigsine**2)

    logger.info('Save: {}'.format('ampwt_increase_coeff.r4'))
    save_r4(amp, 'ampwt_coeff.r4')
    logger.info('Save: {}'.format('ampwt_increase_sigcoeff.r4'))
    save_r4(sigamp, 'ampwt_sigcoeff.r4')
    logger.info('Save: {}'.format('phiwt_increase_coeff.r4'))
    save_r4(phi, 'phipwt_coeff.r4')
    logger.info('Save: {}'.format('phiwt_increase_sigcoeff.r4'))
    save_r4(phi, 'phipwt_sigcoeff.r4')

if arguments["--acceleration"] == 'yes':
    # ax2 + bx + c
    t = np.linspace(idates[0], idates[-1], 100)
    t10perc = np.nanpercentile(t, 10)
    t50perc = np.nanpercentile(t, 50)
    t90perc = np.nanpercentile(t, 90)
    print(t10perc, t50perc, t90perc)

    b = output_param[:, :, indexinter]
    a = output_param[:, :, indexacc]

    sig_b =  output_sigparam[:, :, indexinter]
    sig_a =  output_sigparam[:, :, indexacc]

    lin_moy = a*(0 + idates[-1] - idates[0])/2 + b
    lin_moy3ans = a*3 + b
    lin_10perc = a*(t10perc - idates[0]) + b
    lin_50perc = a*(t50perc - idates[0]) + b
    lin_90perc = a*(t90perc - idates[0]) + b

    #siglin_moy = np.sqrt(sig_a**2 * ((0 + idates[-1] - idates[0])/2)**2 + sig_b**2)

    logger.info('Save: {}'.format('linmoy_coeff.r4'))
    save_r4(lin_moy, 'linmoy_coeff.r4')
    logger.info('Save: {}'.format('linmoy3ans_coeff.r4'))
    save_r4(lin_moy3ans, 'linmoy3ans_coeff.r4')
    logger.info('Save: {}'.format('linmoy10perc_coeff.r4'))
    save_r4(lin_10perc, 'linmoy10perc_coeff.r4')
    logger.info('Save: {}'.format('linmoy50perc_coeff.r4'))
    save_r4(lin_50perc, 'linmoy50perc_coeff.r4')
    logger.info('Save: {}'.format('linmoy90perc_coeff.r4'))
    save_r4(lin_90perc, 'linmoy90perc_coeff.r4')

    #logger.info('Save: {}'.format('linmoy_sigcoeff.r4'))
    #save_r4(siglin_moy, 'linmoy_coeff_sigcoeff.r4')


sys.exit()
 

image_to_show = output_param[:, :, 1]
vmin, vmax = np.nanpercentile(image_to_show, [2, 98])

plt.figure(figsize=(10, 6))
plt.imshow(image_to_show, cmap='rainbow', interpolation='none', vmin=vmin, vmax=vmax)
plt.colorbar(label='Valeur')
plt.title('Affichage de la vitesse de output_param')
plt.xlabel('Colonnes')
plt.ylabel('Lignes')

image_to_show = np.sqrt(output_param[:, :, 2]**2 + output_param[:, :, 3]**2)
vmin, vmax = np.nanpercentile(image_to_show, [2, 98])

plt.figure(figsize=(10, 6))
plt.imshow(image_to_show, cmap='rainbow', interpolation='none', vmin=vmin, vmax=vmax)
plt.colorbar(label='Valeur')
plt.title("Affichage de l'amplitude de output_param")
plt.xlabel('Colonnes')
plt.ylabel('Lignes')

plt.show()

print(f"Temps total : {end_time_1 - start_time_1:.3f} secondes")






sys.exit()





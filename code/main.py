import effk
import effkio
import numpy as np
import scipy.io as scio

# inputs
windows = [((10,5),(50,20)),((40,10),(55,32))] # ((x0,y0),(x1,y1)) corner coordinates of windows to evaluate over
dxy = 1.0 # grid cell size in [m]
fname_condmatrix, vname_condmatrix = 'inverted_transmissivity_field', 'Ks_plot'
#fname_condmatrix, vname_condmatrix = 'periodic_pattern', 'Ks_plot'

# conversion from conductivity to permeability
rho, mu, g = 1000., 0.001, 9.80665
condmatrix = scio.loadmat(fname_condmatrix)[vname_condmatrix].astype(float) # log10([m s-1])
kmatrix = np.multiply(np.power(10, condmatrix), mu/(rho*g)) # [m2]

# run
Lx, Ly, results = effk.effk(kmatrix, dxy, windows)

# screen and file output
effkio.screen_out(results)
effkio.file_out(results)
effkio.ellipse_plot(results, Lx, Ly, fname='keff_windows.png')
effkio.ellipse_plot(results, Lx, Ly, fname='keff_windows.pdf')

# prompt for exit
effkio.prompt_for_exit()
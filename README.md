## Permeability upscaling in two dimensions

Two-dimensional permeabilty upscaling tool for block-grids based on python and fipy
using Matlab data structures.

### Usage 
Our basic imports are
```
import effk
import effkio
import numpy as np
import scipy.io as scio
```
We define windows to evaluate upscaled permeability over using corner coordianates, the lattice size and
the name of the Matlab matrix file and the variable name the data are stored under.
```
# ((x0,y0),(x1,y1)) corner coordinates of windows to evaluate over
windows = [((10,5),(50,20)),((40,10),(55,32))] 
dxy = 1.0 # grid cell size in [m]
fname_condmatrix, vname_condmatrix = 'inverted_transmissivity_field', 'Ks_plot'
```
Data loading vie scipy functions, eventual unit conversion
```
# conversion from conductivity to permeability
rho, mu, g = 1000., 0.001, 9.80665
condmatrix = scio.loadmat(fname_condmatrix)[vname_condmatrix].astype(float) # log10([m s-1])
kmatrix = np.multiply(np.power(10, condmatrix), mu/(rho*g)) # [m2]
```
Call the solver
```
Lx, Ly, results = effk.effk(kmatrix, dxy, windows)
```
Results IO
```
effkio.screen_out(results)
effkio.file_out(results)
effkio.ellipse_plot(results, Lx, Ly, fname='keff_windows.png')
```
See `main.py` for a complete analysis example using data in `doc`.

## Examples
For a conductivity field like 
<p align="left">
  <img src="https://raw.githubusercontent.com/plang85/PyUp/master/doc/permeability.png" height="400">
  <br/>
</p>
And the upscaled permeability for the windows specified illustrated as ellipses
<p align="left">
  <img src="https://raw.githubusercontent.com/plang85/PyUp/master/doc/keff_windows.png" height="400">
  <br/>
</p>

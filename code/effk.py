from fipy import *
import numpy as np
import effkio
from scipy.stats import hmean


def effk(kmatrix, dxy, windows):
    """
    Computes the upscaled permeability tensor over a list of windows within the provided permeability map.

    Parameters
    ----------
    kmatrix: ndarray
        Two-dimensional array containing permeability [m2] values for each grid cell.
    dxy: float
        Uniform grid cell length in x- and y-direction [m].
    windows: list
        A list of tuples of the form ((x0,y0),(x1,y1)) in [m] that specify the windows permeability is upscaled over.

    Returns
    -------
    Lx, Ly: float
        Model length [m] in x- and y-direction.
    result_data:
        A list of [x0, y0, x1, y1, kface_w_mean, kface_w_hmean, kxx, kxy, kyy, kmin, kmax, eigw, eigv] per window.

    Todo
    ----
		Outsource.
    Results class.
    """
    # generating grid
    dim = 2
    print 'generating grid, discretizing variables, setting up BC\'s...'
    ny, nx = kmatrix.shape
    Lx, Ly = nx*dxy, ny*dxy
    mesh = Grid2D(dx=dxy, dy=dxy, nx=nx, ny=ny)

    # dicretized varibles
    p = [ CellVariable(name="fluid pressure", mesh=mesh, value=0.) for i in range(dim) ]
    k = CellVariable(name="permeability", mesh=mesh, value=0.)
    k.setValue(np.ravel(kmatrix))
    if k.min().value <= 0.:
        raise ValueError('erroneous data:permeability between %.2e and %.2e [m2]' % (k.min().value, k.max().value))
    kface = k.harmonicFaceValue # harmonic average at faces for FV diffusion solution
    print '\tplotting...'
    effkio.fipy_scalar_field_plot( k, nx, ny, Lx, Ly, fname = 'permeability.png', log = 1, label = r'Permeability [log$_{10}$(m$^2$)]' )

    # BC's in order x-I[0],y-II[1]
    fin, fout = [mesh.facesLeft,mesh.facesTop], [mesh.facesRight,mesh.facesBottom]
    pin, pout = 1.0, -1.0

    # pressure solution & postprocessing 
    print 'solving steady-state pressure, post-processing velocities and pressure-gradients...'
    pgrad = []
    vel = []
    for i in range(dim): # x-I[0],y-II[1]
        p[i].constrain(pin, fin[i])
        p[i].constrain(pout, fout[i])
        DiffusionTerm(kface).solve(var=p[i])
        pgrad.append(p[i].faceGrad)
        vel.append(pgrad[i]*kface)
        print '\tplotting...'
        fname = str('pressure_' + str(i) + '.png')
        effkio.fipy_scalar_field_plot( p[i], nx, ny, Lx, Ly, fname = fname, label = 'Pressure [Pa]' )


    # evaluating keff over windows
    print 'evaluating keff over windows...'
    X,Y = mesh.faceCenters
    result_data = []
    for i in range(len(windows)):
        x0, y0 = windows[i][0] # bounding box lower left corner
        x1, y1 = windows[i][1] # bounding box upper right corner
        where_w = ( mesh.interiorFaces & (X <= x1) & (X >= x0) & (Y <= y1) & (Y >= y0) )
        kface_w = kface[where_w.value]
        kface_w_mean = numerix.mean(kface_w)
        kface_w_hmean = hmean(kface_w.value)
        vavg = [] # will contain (<vx>, <vy>) for I[0] and II[1]
        pgradavg = []# will contain (<pgradx>, <pgrady>) for I[0] and II[1]
        for j in range(dim): # x(I) and y(II)
            # windowing applied herein
            vavg.append( (np.mean(vel[j].value[0][where_w.value]), np.mean(vel[j].value[1][where_w.value])) )
            pgradavg.append( (np.mean(pgrad[j].value[0][where_w.value]),np.mean(pgrad[j].value[1][where_w.value])) )
        # solving overdetermined problem
        A = np.array( ( (pgradavg[0][0],pgradavg[0][1],0.,0.), (0.,0.,pgradavg[0][0],pgradavg[0][1]), 
                        (pgradavg[1][0],pgradavg[1][1],0.,0.), (0.,0.,pgradavg[1][0],pgradavg[1][1]), 
                        (0.,1.,-1.,0.) ) )
        b = np.array([vavg[0][0],vavg[0][1],vavg[1][0],vavg[1][1],0.])
        a = np.dot( A.T, A )
        b = np.dot( A.T, b ).T
        x = np.linalg.solve( a, b )
        kxx, kxy, kyx, kyy = x
        eigw, eigv = np.linalg.eig( np.array( [[kxx,kxy],[kxy,kyy]] ) )
        kmin, kmax = np.min(eigw), np.max(eigw)
        # book keeping
        result_data += [ (x0, y0, x1, y1, kface_w_mean, kface_w_hmean, kxx, kxy, kyy, kmin, kmax, eigw, eigv) ]    
    return Lx, Ly, result_data




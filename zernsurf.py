import numpy as np

def tracezern(rays, coeff, rorder,aorder,arrsize,rad,eliminate='nan',maxiter=10):
    '''
    Function to trace set of rays to Zernike surface
    Inputs are position and cosines of rays
    Array of Zernike coefficients
    Output is position of rays at surface
    and new cosines after reflection
    '''
    
    # Temp list so inputs are not modified
    raycopy = []
    for lst in rays:
        raycopy.append(lst.copy())
    opd,x,y,z,l,m,n,ux,uy,uz = raycopy
    
    delta = np.ones(len(x)) * 100
    i = 0
    
    while (i < maxiter):
        
        rho = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y,x)
        
        Frho = np.zeros(len(x))
        Ftheta = np.zeros(len(x))
import numpy as np
import matplotlib.pyplot as plt
import prtp
import prtp.conicsolve
from prtp.Rays import Rays
from prtp.WolterOptic import WolterTypeOne
from prtp.Grating import Grating
from prtp.Instrument import Instrument
from prtp.Modification import Modification
from prtp.Sources import Subannulus
import astropy.units as u
import prtp.analyses as analyses

def islope(pos,focus,y=0):
    '''
    pos - initial position of the groove (could eventually be calculated using the index, pos = init_sep * index
    focus - distance to the focus from bottom of wafer
    y - The intial y-position of the point, will almost always be 0
    
    returns: inverse slope (i.e; run/rise). Multiplying this by the change in height will give you the change in x
    '''
    return pos / (focus-y)

def propagate(posns,dist,islopes=[]):
    '''
    posns - initial position of the groove positions in x
    dist - the distance by which you want to move the grooves
    
    returns - a list of the new groove positions
    '''
    return posns - islopes*dist

def fracture(posns,nearest=0.3):
    '''
    posns - The list of current groove positions
    nearest - What you are rounding to. By default, rounds to the nearest 3 Angstroms (0.3 nm)
    
    returns - a new array, with each value rounded to the nearest specified value
    '''
    return np.round(posns / nearest) * nearest


# Define Wolter-I parameters.
r0 = (165. * u.mm)  # Radius at the intersection plane of primary/secondary.
z0 = 3500. * u.mm 
mirror_len = 100. * u.mm
mirror_sep = 5. * u.mm

d = 160. * u.nm  # Groove period @ 3300 mm [nm]
L = 3250. * u.mm # Center of grating [mm]
d *= L  / (3300 * u.mm)  # Find groove period at center of grating.

# Define inner & outer radii to create rays.
rp_front = prtp.conicsolve.primrad(z0 + mirror_sep/2 + mirror_len, r0, z0)
rp_back = prtp.conicsolve.primrad(z0 + mirror_sep/2, r0, z0)

# Define initial rays in subannulus.
source = Subannulus(50000,rp_back, rp_front, np.radians(30.)*u.rad,wave=0.98903*u.nm,order=-5)

# Define Wolter Optic
wolter = WolterTypeOne(r0=r0,z0=z0,beckmann_scatter=True,ripple=1.48e-5)

# Define Grating (values from old code)
grat = Grating(0.*u.mm,151.86466758*u.mm,3247.56521956*u.mm,
            0.,0.99978888,-0.0205472,
            -0.01518378,-0.02054483,-0.99967363,
            l=100*u.mm,w=100*u.mm,d=d,radial=True,fdist=3250*u.mm)


def pfunc(rays):
    
    ystep = (4.5 * u.um).to(u.mm).value
    xstep = (.5 * u.nm).to(u.mm).value
    
    # Get position of each photon on the Grating's surface
    x,y = grat.getPosns(rays)
    # Translate origin so it is now in the bottom center of the grating
    y += (grat.l / 2).value
    # Only consider the right-hand side of the grating (the left side is just mirrored)
    x = np.abs(x)
    
    # Periods observed by each photon will eventually go into this list
    ds = []
    
    for i in range(len(rays)):
        
        # Get the position of the current photon
        posns = x[i]
        # Get the inverse slope of the groove this photon experiences
        isl = islope(posns,(grat.fdist.value+50),y=y[i])
        # Propagate this photon back down to the baseline
        posns = propagate(posns,-y[i],isl)
        # Make a second groove exactly 160 nm away
        posns2 = posns + (160 * u.nm).to(u.mm).value
        # Find its inverse slope
        isl2 = islope(posns2,grat.fdist.value+50)
        
        # Find how many 1.63 micron subfields we will need to traverses
        numsteps = (y[i]//ystep)
        
        # Propagate 
        currentxs = propagate(posns,numsteps*ystep,isl)
        currentxs2 = propagate(posns2,numsteps*ystep,isl2)
        
        ## Comment these two fracture lines to approximate a true radial grating
        currentxs = fracture(currentxs,xstep)
        currentxs2 = fracture(currentxs2,xstep)
        
        ## Note: the commented lines here are not needed currently, may be needed in the future
        nextxs = propagate(posns,(numsteps+1)*ystep,isl)
        nextxs2 = propagate(posns2,(numsteps+1)*ystep,isl2)
        ## Comment these two fracture lines to approximate a true radial grating
        nextxs = fracture(nextxs,xstep)
        nextxs2 = fracture(nextxs2,xstep)
        
        deltax = nextxs - currentxs
        deltax2 = nextxs2 - currentxs2
        
        ydisp = (y[i] % ystep) / ystep
        
        newposns = currentxs + ydisp*deltax
        newposns2 = currentxs2 + ydisp*deltax2
        
        d = newposns2 - newposns
        
        d += np.random.normal(0,.2*u.nm.to(u.mm))
        
        ds.append(d)
    
    # Add units to the groove periods
    ds = np.array(ds) * u.mm
    
    # These next four lines find the period at the center given the period at the photon's position.
    x,y = grat.getPosns(rays)
    dist = np.sqrt(x**2 + (grat.fdist.value - y)**2)
    ds /= dist
    ds *= grat.fdist.value
    
    ## To get a perfect radial grating, uncomment this line:
    ds = np.ones(len(rays)) * 160 * u.nm + (np.random.normal(0,.0,len(rays)) * u.nm)
    
    ## Uncomment for scatter plot of grating periods
    # plt.figure()
    # plt.scatter(x,y,c=ds,s=1.)
    # plt.colorbar()
    # plt.show()
    ##
    
    
    return ds.to(u.nm).value


grat.periodfunction = pfunc

# Initialize the instrument and add components
i = Instrument(source)
i.addComponent(wolter)
i.addComponent(grat)

# Simulate the Rays through the instrument
i.simulate()

# Access the final rays
rays = i.getRays()

# Send the rays to the focus and plot them
rays.focusX()

cent = analyses.centroid(rays)
fwhm = np.std(rays.x) * 2.355

print('Spectral Resolution: ' + str(-cent[0] / fwhm))


# rays.scatter2d()
import matplotlib.pyplot as plt
plt.figure()
plt.scatter(rays.x,rays.y)
plt.axis('equal')
plt.show()














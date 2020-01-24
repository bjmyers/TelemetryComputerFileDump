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
from prtp.Sources import ConvergingBeam2
import astropy.units as u
import prtp.analyses as analyses
from prtp.CollimatorPlate import CollimatorPlate
from tqdm import tqdm

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

def snappedpfunc(rays):
    '''
    Function snappedpfunc
    This function takes in rays traced to a Grating's surface and calculates
    the groove period that each photon experiences. This method uses snapped
    groove lines and a "scanning" photon, where each photon scans the grating
    in a line and samples one groove period randomly
    '''
    
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
    
    # Make an array of groove positions
    groovePosns = (np.arange(312500) * 160 * u.nm).to('mm')
    
    # Get their slopes
    isl = islope(groovePosns,3300)
    
    for i in tqdm(range(len(rays))):
        
        # Get the y-position of the current photon
        yPos = y[i]
        
        # Find how many y-steps we will need to traverse
        numsteps = (yPos//ystep)
        
        # Propagate the grooves to the ystep before the photon
        prevPosns = propagate(groovePosns,numsteps*ystep,isl)
        
        # Propagate the grooves to the ystep after the photon
        nextPosns = propagate(groovePosns,(numsteps+1)*ystep,isl)
        
        # Snap both sets of positions to the nearest xstep
        prevPosns = fracture(prevPosns,xstep)
        nextPosns = fracture(nextPosns,xstep)
        
        # Find how much the positions have changed
        deltaPosn = nextPosns - prevPosns
        
        # Find how much within the y-step the photon has propagated
        ydisp = (yPos % ystep) / ystep
        
        # Find the groove positions at this photon's exact y-position
        nextPosns = prevPosns + ydisp*deltaPosn
        
        # The preiods are found with np.diff(nextPosns). Now we have to
        # sample these periods to find the one the photon experiences
        d = np.random.choice(np.diff(nextPosns))
        
        # Save it in the period array
        ds.append(d)
    
    # Add units to the groove periods
    ds = np.array(ds) * u.mm
    
    # plt.figure()
    # plt.scatter(yold,ds.to('nm'))
    # plt.plot([-50,50], [160, 155], color='k', linestyle='-', linewidth=2)
    # plt.xlabel('Grating Position [mm]')
    # plt.ylabel('Groove Period [nm]')
    # plt.show()
    
    # These next four lines find the period at the center given the period at the photon's position.
    x,y = grat.getPosns(rays)
    dist = np.sqrt(x**2 + (grat.fdist.value - y)**2)
    ds /= dist
    ds *= grat.fdist.value
    
    return ds.to(u.nm).value


def parallelapprox(rays):
    
    # Define the number of segments used to approximate a radial grating
    n = 10

    dpermm = grat.d / grat.fdist
    
    # Get position of each photon on the Grating's surface
    x,y = grat.getPosns(rays)
    
    # Translate origin so it is now in the bottom center of the grating
    y += (grat.l / 2).value
    # Only consider the right-hand side of the grating (the left side is just mirrored)
    x = np.abs(x)
    
    ## Here, find the period experienced by each photon
    
    # These next four lines find the period at the center given the period at the photon's position.
    x,y = grat.getPosns(rays)
    dist = np.sqrt(x**2 + (grat.fdist.value - y)**2)
    ds /= dist
    ds *= grat.fdist.value
    
    return ds.to(u.nm).value


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
source = Subannulus(1000,rp_back, rp_front,wave=0.83401*u.nm,order=1,dphi=30*u.deg)

# Define Wolter Optic
# The WolterTypeOne object contains primary and secondary mirrors.
# Includes calls to the equivalent of PyXFocus' wolterprimary and 
#   woltersecondary functions
# Currently, we are not considering beckmann scattering
wolter = WolterTypeOne(r0=r0,z0=z0,beckmann_scatter=True,ripple=1.48e-5)

fouralpha = np.arctan(r0/z0)
g = 1.5 * u.deg

# Find Grating Parameters
hub = 3250*u.mm
L =  hub / np.cos(g)
gratz = L * np.cos(fouralpha)
graty = L * np.sin(fouralpha)
gratx = 0*u.mm

yaw = -0.87*u.deg # either 0 or -0.87 deg



snappedradial = True


if (snappedradial):

    grat.periodfunction = snappedpfunc
    
    grat = Grating(x=gratx, y=graty, z=gratz,
                nx=0, ny=1, nz=0,
                sx=0, sy=0, sz=-1,
                d=d, radial=True, fdist=L,
                l = 100*u.mm, w = 100*u.mm)

else:
    
    grat.periodfunction = parallelapprox
    
    grat = Grating(x=gratx, y=graty, z=gratz,
                nx=0, ny=1, nz=0,
                sx=0, sy=0, sz=-1,
                d=d, radial=False, fdist=L,
                l = 100*u.mm, w = 100*u.mm)

grat.pitch(g-fouralpha)
grat.yaw(yaw)

rays = source.generateRays()
wolter.trace(rays)
grat.trace(rays)
rays.focusX()

print('Spectral Resolution: ' + str(rays.spectralResolution()))

rays.scatter2d()


















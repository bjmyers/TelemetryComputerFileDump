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
source = Subannulus(1000,rp_back, rp_front, np.radians(30.)*u.rad,wave=2.16*u.nm,order=0)

# Define Wolter Optic
wolter = WolterTypeOne(r0=r0,z0=z0,beckmann_scatter=True,ripple=1.48e-5)

# Define Grating (values from old code)
grat = Grating(0.*u.mm,151.86466758*u.mm,3247.56521956*u.mm,
            0.,0.99978888,-0.0205472,
            -0.01518378,-0.02054483,-0.99967363,
            l=100*u.mm,w=100*u.mm,d=173.9*u.nm,radial=False,fdist=1.9715*u.m)

def pfunc(rays):
    
    n = grat.steps
    
    d0 = 169.36
    dn = 178.43
    
    # h is the step in grating period per subfield
    h = (dn - d0) / n

    # dy is the step in y-position per subfield
    dy = grat.l.value / n
    
    x,y = grat.getPosns(rays)
    
    y *= -1
    y += (grat.l/2).value
    
    y = y // dy
    
    ds = d0 + y*h

    return ds


# ns = np.logspace(1,8)
# ns = np.linspace(2,1000)
ns = np.arange(1,10)
spects = []

# from tqdm import tqdm
# for n in tqdm(ns):
    
grat.steps = 1000

# grat.periodfunction = pfunc

# Initialize the instrument and add components
i = Instrument(source)
i.addComponent(wolter)
# i.addComponent(grat)

# Simulate the Rays through the instrument
i.simulate()

# Access the final rays
rays = i.getRays()

# grat.trace_to_surf(rays)
# rays.reflect()
grat.trace(rays)

# Send the rays to the focus and plot them
rays.focusX()


rays.scatter2d(c=rays.l)

# cent = analyses.centroid(rays)
# fwhm = np.std(rays.x) * 2.355
# 
# spect = np.abs(cent[0]/fwhm)
# 
# spects.append(spect)
# 
# plt.figure()
# plt.plot(ns,spects)
# plt.xlabel('Number of Subfields')
# plt.ylabel('Spectral Resolution')
# plt.show()















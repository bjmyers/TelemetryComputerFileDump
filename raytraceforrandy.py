import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib
import matplotlib.patches as patches
import sys
import pdb
from copy import deepcopy
from tqdm import tqdm
import astropy.units as u
plt.style.use('seaborn-ticks')
# from prettytable import PrettyTable

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# sys.path.append('C:/python/')
import PyXFocus.sources as sources
import PyXFocus.transformations as trans
import PyXFocus.surfaces as surfaces
import PyXFocus.analyses as analyses
import PyXFocus.conicsolve as conic

# sys.path.append('C:/python/My-Python-Workspace/OGRE/')
import ogre_routines as ogre

# Define Wolter-I parameters.
r0 = 165.  # Radius at the intersection plane of primary/secondary.
z0 = 3500. # 
mirror_len = 100.
mirror_sep = 5.

# Define inner & outer radii to create rays.
rp_front = conic.primrad(z0 + mirror_sep/2 + mirror_len, r0, z0)
rp_back = conic.primrad(z0 + mirror_sep/2, r0, z0)

# Define initial rays in subannulus.
rays = sources.subannulus(rp_back, rp_front, np.radians(30.), 10000)
# rays = sources.subannulus(rp_back, rp_front, np.radians(360.), 100000)

# Transform rays to intersection plane of optic.
trans.transform(rays, 0, 0, -z0, 0, 0, 0)

# Rotate 90 degrees so that diffraction occurs in x-axis.
trans.transform(rays, 0, 0, 0, 0, 0, -np.radians(90.))

# Pass through primary.
surfaces.wolterprimary(rays, r0, z0)
trans.reflect(rays)

# Add some scatter.
ogre.beckmann_scatter(rays, 0, 0, 1.48e-5)

# Pass through secondary.
surfaces.woltersecondary(rays, r0, z0)
trans.reflect(rays)

# Go to optic focus.
surfaces.focusX(rays)

# Define grating parameters.
d = 160.  # Groove period @ 3300 mm [nm]
L = 3250.  # Center of grating [mm]
d *= L  / 3300  # Find groove period at center of grating.
gammy = np.radians(1.5)  # Incidence angle.
wave = 0.83401  # [nm] Al-K wavelength.
yaw = np.radians(0.87)  # Approx. yaw in OGRE geometry.

# Create copy if you need to reference later.
optic_rays = deepcopy(rays)

# Establish starting coordinates.
glob_coords = [trans.tr.identity_matrix()] * 4

# 'Steer' rays to converging beam. Align z-axis with converging beam.
mean_beam = np.mean(rays[5])
trans.transform(rays, 0, 0, 0, -mean_beam, 0, 0, coords=glob_coords)

# Go to center of grating.
conv_length = L / np.cos(gammy)
trans.transform(rays, 0, 0, conv_length, 0, 0, 0, coords=glob_coords)

# Add incidence angle.
trans.transform(rays, 0, 0, 0, gammy, 0, 0, coords=glob_coords)

# Get +y in groove direction.
trans.transform(rays, 0, 0, 0, -np.pi/2, 0, 0, coords=glob_coords)
trans.transform(rays, 0, 0, 0, 0, 0, np.pi, coords=glob_coords)

# Add yaw.
trans.transform(rays, 0, 0, 0, 0, 0, yaw, coords=glob_coords)

# Project rays onto grating surface.
surfaces.flat(rays)

# Move to coordinate system to hub location (required for radgrat).
trans.transform(rays, 0, -L, 0, 0, 0, 0, coords=glob_coords)

good_ind = np.where((rays[2] < 3300.) & (rays[2] > 3200.))[0]
rays = [r[good_ind] for r in deepcopy(rays)]

# Reflect rays.
trans.reflect(rays)

    
#Iterate through every photon in the rays object (PyXFocus does not allow different grating periods for each photon)
trans.radgrat(rays, d/L, -1, wave)
   
#Go back to starting coordinate system.
rays = trans.applyT(rays, glob_coords, inverse=True)
    
#Focus rays.
d = surfaces.focusX(rays)

plt.figure()
plt.scatter(rays[1],rays[2])
plt.show()

# # 
# # y = rays[2]
# # cen = np.mean(y)
# # lower = cen - .000001
# # upper = cen + .000001
# # 
# # while True:
# # 
# #     count = np.where((y > lower) & (y < upper))[0]
# #     
# #     if len(count) > len(rays[0])/2:
# #         print(lower)
# #         print(upper)
# #         print(upper-lower)
# #         break
# #     
# #     lower -= .000001
# #     upper += .000001

# Initialize the list of standard deviations we want to test
# sigmas = np.logspace(2.5,6)
# 
# Initialize the list our spectral resolutions will go into
# res = []
# 
# Create a copy of the rays at this point (right before they get reflected by the grating) so that all our trials start from the same point
# copied_rays = deepcopy(rays)
# 
# for s in sigmas:
#     
#     Retreive the copied rays
#     rays = deepcopy(copied_rays)
#     
#     Iterate through every photon in the rays object (PyXFocus does not allow different grating periods for each photon)
#     for i in range(len(rays[0])):
#         The grating period is normally distributed with mean d/L and standard deviation d/L/s.
#         s varies from ~300 to 100,000. s is defined such that if s was 1000, a random period would be selected from a distribution with a standard deviation that was one part in 1000 the grating period (d/L)
#         trans.radgrat(rays, np.random.normal(d/L,d/L/s), -int(5), wave, ind = [i])
#     
#     Go back to starting coordinate system.
#     rays = trans.applyT(rays, glob_coords, inverse=True)
#     
#     Focus rays.
#     surfaces.focusX(rays)
#     
#     Calculate and store the spectral resolution
#     cent = analyses.centroid(rays)
#     fwhm = np.std(rays[1]) * 2.355
#     The resolution is x/(deltax) ~= mean(x) / fwhm(x)
#     res.append(-cent[0] / fwhm)
# 
# Plot the resolution as a function of standard deviation
# plt.figure()
# plt.plot(160/sigmas,res)
# plt.xlabel("Standard Deviation nm)")
# plt.ylabel("Spectral Resolution")
# plt.show()























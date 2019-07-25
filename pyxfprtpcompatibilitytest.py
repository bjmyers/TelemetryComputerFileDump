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
import pickle
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

import time
t = time.time()

# Define Wolter-I parameters.
r0 = 165.  # Radius at the intersection plane of primary/secondary.
z0 = 3500. # 
mirror_len = 100.
mirror_sep = 5.

# Define inner & outer radii to create rays.
rp_front = conic.primrad(z0 + mirror_sep/2 + mirror_len, r0, z0)
rp_back = conic.primrad(z0 + mirror_sep/2, r0, z0)

# Define initial rays in subannulus.
rays = sources.subannulus(rp_back, rp_front, np.radians(30.), 1000)

# Transform rays to intersection plane of optic.
trans.transform(rays, 0, 0, -z0, 0, 0, 0)

# Rotate 90 degrees so that diffraction occurs in x-axis.
trans.transform(rays, 0, 0, 0, 0, 0, -np.radians(90.))

# Pass through primary.
surfaces.wolterprimary(rays, r0, z0)
trans.reflect(rays)

# # Add some scatter.
# ogre.beckmann_scatter(rays, 0, 0, 1.48e-5)

# Pass through secondary.
surfaces.woltersecondary(rays, r0, z0)
trans.reflect(rays)


# # Go to optic focus.
x = surfaces.focusX(rays)

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


##

# grat = sources.pointsource(0,1)
# grat[8] = np.array([1.])
# grat[2] = np.array([3250.])
# grat = trans.applyT(grat, glob_coords, inverse=True)
# 
# from prtp.Grating import Grating
# from prtp.Rays import Rays
# from prtp.WolterTypeOne import WolterTypeOne
# from prtp.Instrument import Instrument
# from prtp.CollimatorPlate import CollimatorPlate
# 
# # Define initial rays in subannulus.
# prays = Rays.subannulus(rp_back, rp_front, np.radians(30.), 1000)
# # Rotate 90 degrees so that diffraction occurs in x-axis.
# prays.transform(0, 0, 0, 0, 0, np.radians(90.))
# 
# w = WolterTypeOne(r0=r0,z0=z0)
# g1 = Grating(0.,151.86466758,3247.56521956,
#             0.,0.99978888,-0.0205472,
#             0.01518378,0.02054483,0.99967363,
#             l=100,w=1000,d=d,radial=True,fdist=3250)
# g2 = Grating(grat[1][0],grat[2][0],grat[3][0],grat[4][0],grat[5][0],grat[6][0],grat[7][0],grat[8][0],grat[9][0],l=100,w=1000,d=d,radial=True,fdist=3250)
# c2 = CollimatorPlate(grat[1][0],grat[2][0],grat[3][0],grat[4][0],grat[5][0],grat[6][0],grat[7][0],grat[8][0],grat[9][0],l=100,w=1000)
# 
# i = Instrument(prays,waves=wave,orders=0)
# i.addComponent(w)
# i.addComponent(g2)
# i.simulate()
# 
# prays.focusX()

##


good_ind = np.where((rays[2] < 3300.) & (rays[2] > 3200.))[0]
rays = [r[good_ind] for r in deepcopy(rays)]

# Reflect rays.
trans.reflect(rays)

# Diffract rays.
trans.radgrat(rays, d/L, -1, wave)

# Go back to starting coordinate system.
rays = trans.applyT(rays, glob_coords, inverse=True)


# xcen = rays[1].min() + (rays[1].max() - rays[1].min()) / 2
# ycen = rays[2].min() + (rays[2].max() - rays[2].min()) / 2
# zcen = rays[3].min() + (rays[3].max() - rays[3].min()) / 2
# ind = np.argmin(np.abs(rays[1]))
# dir = np.array([xcen - rays[1][ind], ycen - rays[2][ind], zcen - rays[3][ind]])
# 
# ind = np.argmin(np.abs(rays[2] - ycen))
# adir = np.array([xcen - rays[1][ind], ycen - rays[2][ind], zcen - rays[3][ind]])
# 
# norm = np.cross(dir,adir) * -1
# 
# 
# Focus rays.
surfaces.focusX(rays)

# plt.figure()
# plt.scatter(rays[4],rays[5],c='b')
# plt.scatter(prays.l,prays.m,c='g')
# plt.axis('equal')
# plt.show()


# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(prays.x,prays.y,prays.z,c='r')
# # ax.scatter(pc.x,pc.y,pc.z,c='b')
# ax.scatter(rays[1],rays[2],rays[3])
# # ax.quiver(grat[1],grat[2],grat[3],-grat[4],-grat[5],-grat[6],length=1)
# # ax.quiver(grat[1],grat[2],grat[3],grat[7],grat[8],grat[9],length=10)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()


# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(prays.x,prays.y,prays.z,c='r')
# ax.scatter(rays[1],rays[2],rays[3])
# plt.xlabel('l')
# plt.ylabel('m')
# plt.show()







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

def gratPeriod(x,xs,periods):
    indexes = np.searchsorted(xs,x)
    indexes = np.where(indexes == len(indexes), indexes-1,indexes)
    selection = np.random.rand(len(x))
    output = []
    for i in range(len(x)):
        if selection[i] < periods[indexes[i]][0][1]:
            output.append(periods[indexes[i]][0][0])
        else:
            output.append(periods[indexes[i]][1][0])
    return output

def idealGratPeriod(x,xs,periods):
    # Can be more or less than one probability
    initial_sep = 0
    for i in range(len(periods[0])):
        initial_sep += periods[0][i][0] * periods[0][i][1]
    final_sep = 0
    for i in range(len(periods[-1])):
        final_sep = periods[-1][i][0] * periods[-1][i][1]
    sep_slope = (final_sep - initial_sep) / (xs[-1] - xs[0])
    return ((x - xs[0]) * sep_slope) + initial_sep


xs1 = np.load('positions_5.npy')
periods1 = np.load('Probabilities_5.npy')
xs2 = np.load('positions_3.npy')
periods2 = np.load('Probabilities_3.npy')


# Define Wolter-I parameters.
r0 = 165.  # Radius at the intersection plane of primary/secondary.
z0 = 3500. # 
mirror_len = 100.
mirror_sep = 5.

# Define inner & outer radii to create rays.
rp_front = conic.primrad(z0 + mirror_sep/2 + mirror_len, r0, z0)
rp_back = conic.primrad(z0 + mirror_sep/2, r0, z0)

# Define initial rays in subannulus.
rays = sources.subannulus(rp_back, rp_front, np.radians(30.), 100000)

# Transform rays to intersection plane of optic.
trans.transform(rays, 0, 0, -z0, 0, 0, 0)

# Rotate 90 degrees so that diffraction occurs in x-axis.
trans.transform(rays, 0, 0, 0, 0, 0, -np.radians(90.))

# Pass through primary.
surfaces.wolterprimary(rays, r0, z0)
trans.reflect(rays)

# # Add some scatter.
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

ys = np.array(rays[2]).copy()
zero_point = L+50
ys -= zero_point
ys *= -1e6

controlds = gratPeriod(ys,xs1,periods1)
controlds = np.round(controlds,decimals=2)

expds = gratPeriod(ys,xs2,periods2)
expds = np.round(expds,decimals=2)

# Make a copy for the special grating
control_rays = deepcopy(rays)
experiment_rays = deepcopy(rays)

for period in np.unique(controlds):
    ind = np.where(controlds==period)[0]
    trans.grat(control_rays, period, np.full(len(ind), 5, dtype=float), np.full(len(ind), wave),ind=ind)
    
for period in np.unique(expds):
    ind = np.where(expds==period)[0]
    trans.grat(experiment_rays, period, np.full(len(ind), 5, dtype=float), np.full(len(ind), wave),ind=ind)


# # Diffract rays.
trans.radgrat(rays, d/L, -int(5), wave)

# Go back to starting coordinate system.
rays = trans.applyT(rays, glob_coords, inverse=True)
experiment_rays = trans.applyT(experiment_rays, glob_coords, inverse=True)
control_rays = trans.applyT(control_rays, glob_coords, inverse=True)

# Focus rays.
surfaces.focusX(rays)
surfaces.focusX(experiment_rays)
surfaces.focusX(control_rays)


plt.figure()
plt.scatter(experiment_rays[1],experiment_rays[2],s=0.5)
plt.scatter(control_rays[1],control_rays[2],s=0.5)
plt.scatter(rays[1],rays[2],s=0.5)
plt.legend(('0.3nm Snapping','0.5nm Snapping','Radial'))
# plt.axis('equal')
plt.show()


cent = analyses.centroid(rays)
fwhm = np.std(rays[1]) * 2.355
print('Radial Resolution: ' + str(cent[0] / fwhm))

cent = analyses.centroid(control_rays)
fwhm = np.std(control_rays[1]) * 2.355
print('0.1nm Snapped Resolution: ' + str(cent[0] / fwhm))

cent = analyses.centroid(experiment_rays)
fwhm = np.std(experiment_rays[1]) * 2.355
print('0.5nm Snapped Resolution: ' + str(cent[0] / fwhm))


# compare perfect radial versus various step sizes




















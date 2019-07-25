#!/usr/bin/env python
# coding: utf-8

# In[1]:
import time
t = time.time()

# get_ipython().magic(u'matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
# from tabulate import tabulate

# sys.path.append('/home/bdonovan')
import PyXFocus.sources as sources
import PyXFocus.transformations as trans
import PyXFocus.surfaces as surfaces
import PyXFocus.analyses as analyses
import PyXFocus.conicsolve as conic


# We will examine the diffraction arc for the imminent testing that will take place in the long cell. We first define the parameters for the optic and grating which will be used, as well as several parameters for the long cell itself.

# In[14]:


# Optic Parameters
z0 = 3500.  # [mm] Focal length.
r0 = 165.  # [mm] Radius at intersection node.
mirror_len = 100.  # [mm] Axial length of primary / secondary mirror.
mirror_sep = 5.  # [mm] Separation between primary and secondary mirrors.
mirror_width = 60.

# Grating Parameters 
hub = 3250.  # [mm] Hub length.
d = 160.  # [nm] At 3300 mm from grating hub.
d = d * 3300 / hub # [nm] Redefine to value at center of grating.
i = np.radians(1.5)  # [rad.] Graze angle.
blaze = 0.  # [rad.] Blaze angle.
yaw = 0. # [rad.] Yaw of grating.

# Long Cell Parameters
# L = 48000.  # [mm] Distance from source to end of test chamber (from E. Bray).
L = 48300. 
L -= 12750.  # [mm] Accounts for focal length of optic in the finite conjugate.
# wave = 0.98903  # [nm] Mg-K wavelength.
wave = 0.83401  # [nm] Al-K wavelength.


# We begin by defining some rays from a source at the end of the test chamber. 

# In[15]:


# Define inner and outer subannulus radii.
z_in = z0 + mirror_sep/2
z_out = z_in + mirror_len

r_in = conic.primrad(z_in, r0, z0)
r_out = conic.primrad(z_out, r0, z0)

# Define full angular width of subannulus.
dphi = np.radians(30.)

# # Define subannulus of rays.
rays = sources.subannulus(r_in, r_out, dphi, 1000000)

# rays = np.loadtxt('rays.txt')

# Rotate so that dispersion direction is in the x-dimension.
trans.transform(rays, 0, 0, 0, 0, 0, -np.pi/2)

# Find centroid of rays and move rays down to center of beamline.
cen_optic = analyses.centroid(rays)

trans.transform(rays, 0, cen_optic[1], 0, 0, 0, 0)
trans.transform(rays, 0, -30, 0, 0, 0, 0)

# Change ray direction cosines to emanate from source (48 m away).
hyp = np.sqrt(L**2 + rays[1]**2 + rays[2]**2)
l = rays[1]/hyp
m = rays[2]/hyp
n = -np.sqrt(1. - l**2 - m**2)
rays = [rays[0], rays[1], rays[2], rays[3],
        l, m, n, rays[7], rays[8], rays[9]]


# Shift rays back up.
trans.transform(rays, 0, -cen_optic[1], 0, 0, 0, 0)
trans.transform(rays, 0, 30, 0, 0, 0, 0)

# Plot rays.
# plt.figure()
# plt.scatter(rays[1], rays[2], s=0.5, c='k')
# plt.axis('equal')
# plt.xlabel('X [mm]')
# plt.ylabel('Y [mm]')


# This is what the rays look like at the entrance aperture of the primary mirror.

# We will now propagate these rays through the primary and secondary mirrors. Here, we will approximate the PSF to have a Gaussian distribution. This is a simplification, but will be fine for determining the location of orders on the focal plane. Further, we assume both the primary and secondary mirrors have infinite length in the axial dimension.

# In[16]:


# Move to intersection plane.
rays[3] = np.ones(len(rays[1])) * z0
 
# Pass rays through primary mirror.
surfaces.wolterprimary(rays, r0, z0)
trans.reflect(rays)

# Look at distribution of y-hat.
# good_ind = np.where(rays[5] > -0.010168)
# good_ind = np.where(rays[5] > -0.010164)
# rays = [row[good_ind] for row in rays]

# Pass rays through secondary mirror.
# surfaces.wsSecondary(rays,r0,z0,psi=1) # Testing
surfaces.woltersecondary(rays, r0, z0)
trans.reflect(rays)

# Add Gaussian scatter to ray direction cosines.
rays[4] = rays[4] + np.random.normal(scale=2.e-6, size=len(rays[4]))
rays[5] = rays[5] + np.random.normal(scale=1.5e-5, size=len(rays[5]))
# errorl = np.loadtxt('errorl')
# rays[4] += errorl
# errorm = np.loadtxt('errorm')
# rays[5] += errorm
rays[6] = -np.sqrt(1. - rays[5]**2 - rays[4]**2)

# Go to focal plane of optic.
f0 = surfaces.focusX(rays)
eff_f0 = z0 - f0

print('Effective focal length: ' + str(eff_f0/1000) + ' m')
# 
# # Move such that centroid of focus is (0,0).
cen_focus = analyses.centroid(rays)
trans.transform(rays, cen_focus[0], cen_focus[1], 0, 0, 0, 0)
# 
# # Plot optic focus.
# # plt.figure()
# # plt.scatter(rays[1], rays[2], s=0.5, c='k')
# # plt.axis('equal')
# # plt.xlabel('X [mm]')
# # plt.ylabel('Y [mm]')
# 
# # Create copy of rays to use when plotting diffraction arc.
# copied_rays = trans.copy_rays(rays)
# 
# 
# # Let's now look to see what the beam looks like at the grating.
# 
# # In[17]:
# 
# 
trans.transform(rays, 0, 0, 3250, 0, 0, 0)
surfaces.flat(rays)

print("Runtime: " + str(time.time() - t) + " seconds")

# Plot optic focus.
plt.figure()
plt.scatter(rays[1], rays[2], s=0.5, c='k')
plt.axis('equal')
plt.xlabel('X [mm]')
plt.ylabel('Y [mm]')
plt.title('Rays at Grating')
plt.show()

# In[20]:


np.max(rays[1]) - np.min(rays[1])


# In[ ]:





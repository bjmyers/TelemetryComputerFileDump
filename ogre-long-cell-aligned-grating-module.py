
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys
import pdb
from copy import deepcopy
from tqdm import tqdm
import time
plt.style.use('seaborn')

sys.path.append('/bdonovan/python/')
import PyXFocus.sources as sources
import PyXFocus.transformations as trans
import PyXFocus.surfaces as surfaces
import PyXFocus.analyses as analyses
import PyXFocus.conicsolve as conic

import PyXFocus.transformMod as tr

# sys.path.append('/home/bdonovan/python/My-Python-Workspace/OGRE/')
# import ogre_routines as ogre


# In this notebook, we will simulate the testing of three (3) aligned gratings with a stack of three (3) optic shells. 
# 
# The three shells are assumed to be the inner three shells of the OGRE optic module:
# * $R_{int}$ = 165.0, 167.5502851426756, 170.1193298483941 mm
# * P-H separation: 5 mm
# * Azimuthal extent: 30 deg.
# * $Z_{int}$ = 3500 mm
# * Length of mirrors along optical axis: 100 mm
# 
# The three gratings have the following specifications:
# * d = 160 nm (@ 3300 mm)
# * Stacked vertically
# * Pitch = 1.5 deg.
# * Yaw = 0.87 deg.
# * L = 3250 mm (at center of grating)
# 
# We first will simulate the optic.
# 
# ## Optic Simulation

t = time.time()
# In[2]:
num = 100000

# Optic Parameters
z0 = 3500.  # [mm] Focal length.
r_int = [165.0, 167.550, 170.119]  # [mm] Radius at intersection node.
mirror_length = 100.  # [mm] Axial length of primary / secondary mirror.
mirror_sep = 5.  # [mm] Separation between primary and secondary mirrors.

# Long Cell Parameters
# L = 48000.  # [mm] Distance from source to end of test chamber (from E. Bray).
L = 48300. 
L -= 4500.  # [mm] Accounts for focal length of optic in the finite conjugate.
# wave = 0.98903  # [nm] Mg-K wavelength.
wave = 0.83401  # [nm] Al-K wavelength.


# In[3]:


# Define inner and outer subannulus radii.
z_in = z0 + mirror_sep/2
z_out = z_in + mirror_length

r_in = conic.primrad(z_in, r_int[0], z0)
r_out = conic.primrad(z_out, r_int[-1], z0)

# Define full angular width of subannulus.
dphi = np.radians(30.)

# Define subannulus of rays.
rays = sources.subannulus(r_in, r_out, dphi, num)

np.save('originalpyxfrays.npy',rays)

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

# Move to intersection node.
trans.transform(rays, 0, 0, -z0, 0, 0, 0)
# 
# # Plot rays.
# plt.figure()
# plt.scatter(rays[1], rays[2], s=0.5, c='k')
# plt.axis('equal')
# plt.xlabel('X [mm]')
# plt.ylabel('Y [mm]')
# plt.title('Rays @ Optic Front [Z=3602.5 mm]')
# 
# # print(rays)
# 
# 
# # We first propagate these rays through the three-shell optic.
# 
# # In[4]:
# 
# 
# Define blank PyXFocus ray object.
new_rays = sources.annulus(0,0,0)

# Define mirror parameters.
zp_back = z0 + mirror_sep/2  # Axial position of parabola front.
zp_front = zp_back + mirror_length  # axial position of parabola back.
rp_front = [conic.primrad(zp_front, r, z0) for r in r_int]
rp_back = [conic.primrad(zp_back, r, z0) for r in r_int]

# Loop through each mirror in mirror stack.
for i in range(len(r_int)):
    # Find photons which will hit this mirror shell.
    r = np.sqrt(rays[1]**2 + rays[2]**2)
    ind = np.where((r > rp_back[i]) & (r < rp_front[i]))[0]
    
    # Create new ray object with only these rays.
    ind_rays = [r[ind] for r in rays]
    
    # Propagate these photons to primary mirror.
    surfaces.wolterprimary(ind_rays, r_int[i], z0)
    
    # Find which photons interact with this mirror.
    ind = np.where((ind_rays[3] > (z0 + mirror_sep/2)) &
                   (ind_rays[3] < (z0 + mirror_sep/2 + mirror_length)))[0]

    # Keep only the photons which interact with the actual size
    # of the mirror
    ind_rays = [r[ind] for r in ind_rays]
    
    # Reflect photons off of primary.
    trans.reflect(ind_rays)
    
    # # Add Beckmann + Gaussian scatter.
    # # ogre.beckmann_scatter(ind_rays, 0, 0, 2e-5)
    # ind_rays[4] = ind_rays[4] + np.random.normal(scale=1e-7, size=len(ind_rays[4]))
    # ind_rays[5] = ind_rays[5] + np.random.normal(scale=1e-6, size=len(ind_rays[5]))
    # ind_rays[6] = -np.sqrt(1. - ind_rays[5]**2 - ind_rays[4]**2)
    
    # Propagate photons to the secondary mirror.
    surfaces.woltersecondary(ind_rays, r_int[i], z0)

    # Find which photons will interact with hyperboloid.
    ind = np.where((ind_rays[3] < (z0 - mirror_sep/2)) &
                   (ind_rays[3] > (z0 - mirror_sep/2 - mirror_length)))[0]
    
    # keep only photons which interact with mirror.
    ind_rays = [r[ind] for r in ind_rays]
    
    # Reflect the photons off of the secondary.
    trans.reflect(ind_rays)
    
    # Add to 'master' PyXFocus ray object.
    new_rays = [np.append(new_rays[i], ind_rays[i])
                for i in range(len(new_rays))]
    
# Copy rays. 
rays = trans.copy_rays(new_rays)


# With rays propagated through optic, we can now go to the focus.

# In[5]:


# Go to the X-ray focus.
f0 = surfaces.focusX(rays)

# Put optic focus at y=0.
cen_y = np.mean(rays[2])
trans.transform(rays, 0, cen_y, 0, 0, 0, 0)

# Copy rays to reference later.
of_rays = deepcopy(rays)

# # Plot rays.
# plt.figure()
# plt.scatter(rays[1], rays[2], s=0.5, c='k')
# plt.axis('equal')
# plt.xlabel('X [mm]')
# plt.ylabel('Y [mm]')
# plt.title('Rays @ Optic Focus [Z=' + str(round(f0, 2)) + ' mm]')
# 
# fwhm_disp = np.std(rays[1]) * 2.355 / 3800 * 206265
# hpd_xdisp = analyses.hpdY(rays) / 3800 * 206265
# plt.text(0.15, 0.2, 'Disp. FWHM = ' + str(round(fwhm_disp, 2)) + '"')
# plt.text(0.15, 0.15, 'X-Disp. HPD = ' + str(round(hpd_xdisp, 2)) + '"')
# print fwhm_disp, hpd_xdisp
# 
# 
# With rays now passed through the optic stack, we can focus on the grating stack. We first define the grating parameters.

# In[6]:


# Grating Parameters 
hub = 3250.  # [mm] Hub length.
d = 160.  # [nm] At 3300 mm from grating hub.
d = d * 3300 / hub # [nm] Redefine to value at center of grating.
gammy = np.radians(1.5)  # [rad.] Graze angle.
blaze = 0.  # [rad.] Blaze angle.
yaw = 0.87 # [rad.] Yaw of grating.
throw = hub / np.cos(gammy)

grat_length = 75.


# Now we figure out what the optic looks like were the gratings should be. This will tell us how many gratings we need to test. 

# In[17]:


# Calculate mean convergence angle of rays.
r = np.mean(r_int) - cen_y
z = z0 - f0

conv_ang = np.arctan(r/z)


# Find out how far we need to move from focal plane.
z_grat = np.cos(conv_ang) * throw

# Copy optic focus rays to use.
rays = deepcopy(of_rays)

# Transform rays to grating location.
trans.transform(rays, 0, 0, z_grat, 0, 0, 0)
surfaces.flat(rays)

# # Plot rays at this position.
# plt.figure()
# plt.scatter(rays[1], rays[2], s=0.5, c='k')
# plt.axis('equal')
# plt.xlabel('X [mm]')
# plt.ylabel('Y [mm]')
# plt.title('Rays @ Grating [Z=' + str(round(f0 + z_grat, 2)) + ' mm]')
# plt.text(10, 165, 'Vertical Extent: ' + str(round(np.max(rays[2])-np.min(rays[2]),2)) + ' mm')
# 
# 
# The vertical extent of the optic is 10.14 mm, which means we would need five (5) 2.3mm thick gratings to fully sample the three-shell optic. 
# 
# **Note: It's not fully clear what optic we will receive (if any). I'll proceed with this design, but it might not be what we receive.**
# 
# We will start from the bottom and then place gratings until we can sample the entire optic.
# 
# In[18]:
# 

# Find y-positions of each grating.
r_cen = [np.min(rays[2]) + grat_length/2 * np.sin(gammy)]
while (r_cen[-1] + grat_length/2 * np.sin(gammy)) < np.max(rays[2]):
    r_cen.append(r_cen[-1]+2.3)


# We require five gratings to sample the three-shell optic. The gratings will be aligned to be at the same z-position at a pitch of zero degrees. This means that when they are pitched, some gratings will be slightly closer to and some will be slightly further from the focal plane. The middle grating will be placed at the exact distance away from the focal plane as determined above.

# In[19]:


# Find z-positions of each grating.
z_cen = np.arange(-2, 3) * 2.3 * np.sin(gammy) + z_grat

# We now will propagate the photons through the grating array.

# In[20]:


# Define array to store diffracted rays.
diff_inds = np.array([])
waves = np.ones(len(rays[0])) * wave

# Copy arrays to reference later.
rays = deepcopy(of_rays)

for i in range(len(z_cen)):
    # Define coordinate system so that we can return to it later.
    glob_coords = [np.identity(4)] * 4
    # Move to grating location.
    trans.transform(rays, 0, r_cen[i], z_cen[i], 0, 0, 0, coords=glob_coords)
    # Rotate to angle of the beam.
    trans.transform(rays, 0, 0, 0, -np.pi/2, 0, 0, coords=glob_coords)
    trans.transform(rays, 0, 0, 0, -conv_ang, 0, 0, coords=glob_coords)
    # Put incidence angle onto grating.
    trans.transform(rays, 0, 0, 0, gammy, 0, 0, coords=glob_coords)
    # Get +y to point towards grating surface for radgrat function.
    trans.transform(rays, 0, 0, 0, 0, 0, np.pi, coords=glob_coords)
    # Add yaw.
    trans.transform(rays, 0, 0, 0, 0, 0, yaw, coords=glob_coords)
    # Go to hub location.
    trans.transform(rays, 0, -hub, 0, 0, 0, 0, coords=glob_coords)
    # Project photons onto x-y plane (the grating surface).
    surfaces.flat(rays)
    # Find which photons will interact with this grating.
    ind = np.where((rays[2] > (hub - grat_length/2)) &
                   (rays[2] < (hub + grat_length/2)))[0]
    # Remove indices which have already been reflected/diffracted.
    ind = ind[np.isin(ind, diff_inds, invert=True)]
    ind = np.array(ind, dtype=int)
    print(ind)
    diff_inds = np.concatenate((ind, diff_inds))
    # If there are no rays that interact with the grating, continue.
    
    if len(ind) < 1:
        rays = trans.applyT(rays, glob_coords, inverse=True)
        continue
    # Reflect photons which fall onto this grating.
    trans.reflect(rays, ind=ind)
    # Diffract photons.
    trans.radgrat(rays, 1e6/d/hub, 1, waves, ind=ind)
    # trans.radgrat(rays,1e6/d/hub,1,waves)
    # # Return back to original coordinate system.
    rays = trans.applyT(rays, glob_coords, inverse=True)
    
# Only keep rays which have been diffracted.
diff_inds = np.array(diff_inds, dtype=int)
print(len(rays[0]), len(diff_inds))
rays = [r[diff_inds] for r in rays]


# With the rays propagated through the grating module, we will now go to the spectral focus.

# In[22]:


# Find optimal focus position.
grat_focus = surfaces.focusX(rays)
print(grat_focus)



print('PyXFocus Time: ' + str(time.time() - t) + ' seconds')
print()
# Plot rays at this position.
plt.figure()
plt.scatter(rays[1], rays[2], s=0.5, c='k')
plt.axis('equal')
plt.xlabel('X [mm]')
plt.ylabel('Y [mm]')
plt.show()

np.save('pyxfocusrays.npy',rays)

# plt.title('Rays @ Grating [Z=' + str(round(grat_focus, 2)) + ' mm]')


# According to the plot above, if we had an optic with a PSF with a 3" HPD in the cross-dispersion direction, we should be able to easily separate the five components of the grating array. This also requires that the gratings perform perfectly. However, we know that this is not the case.
# 
# This test will use the "direct-write" grating from the 2018 PANTER test campaign. If the SCIL replicated the gratings perfectly, we would have a grating LSF on the focal plane that is ~4 mm in the cross-dispersion direction. Therefore, we would not be able to separate out each of the five grating LSFs. Further, yaw errors from dicing and alignment will cause the LSFs to move in the cross-dispersion direction, further complicating matters. If there are roll errors though, we might be able to see these. 

import numpy as np
import matplotlib.pyplot as plt
import time
import os
import astropy.units as u
import prtp.transformationsf as transE
from prtp.Rays import Rays
from prtp.WolterOptic import WolterPrimary
from prtp.WolterOptic import WolterSecondary
from prtp.WolterOptic import WolterTypeOne
from prtp.WolterOptic import WolterModule
from prtp.Grating import Grating
from prtp.CollimatorPlate import CollimatorPlate
from prtp.Detector import Detector
from prtp.Combination import Combination
from prtp.GratingStack import GratingStack
from prtp.Instrument import Instrument
from prtp.Sources import CircularBeam
from prtp.Sources import Subannulus
from prtp.Sources import PointSource
from prtp.Modification import Modification
from prtp.FlatComponent import FlatComponent
import ogre_routines as ogre
import PyXFocus.surfaces as surf
import PyXFocus.analyses as analyses
from tqdm import tqdm

## Detector Noise Testing

import matplotlib.pyplot as plt
from prtp.Detector import Detector
from prtp.Sources import CircularBeam
import astropy.units as u

s = CircularBeam(num=100,rad=8*u.mm,wave=100*u.eV)
rays = s.generateRays()

d = Detector(x=0*u.mm,y=0*u.mm,z=2*u.mm,
    nx=0,ny=0,nz=1,sx=0,sy=1,sz=0,q=1.,
    l=10*u.mm,w=10*u.mm,xpix=100,ypix=100,
    considersplits=True,steps=50)

d.trace(rays)
arr = d.view(rays)

plt.figure()
plt.imshow(arr)
plt.show()

# rays = CircularBeam(rad=5*u.mm).generateRays()
# 
# rays.addTag('xaaaa',rays.x < 0)
# rays.remove(tags='!x')
# 
# rays.scatter2d()


## OGRE Sim

# rays = ogre.create_rays(n=100000)
# rays = ogre.ogre_mirror_module(rays,scatter='b')
# rays = ogre.ogre_grating_module(rays)
# 
# # # # 3D Plot
# # # from mpl_toolkits.mplot3d import Axes3D
# # # fig = plt.figure()
# # # ax = fig.add_subplot(111, projection='3d')
# # # # c = ax.scatter(rays[4],rays[5],rays[6],c=rays[1])
# # # c = ax.scatter(rays[1],rays[2],rays[3],c=rays[4])
# # # plt.xlabel('l')
# # # plt.ylabel('m')
# # # plt.colorbar(c)
# # # plt.show()
# 
# x = surf.focusX(rays)
# 
# print(np.std(rays[1] * 2.355))
# print(analyses.hpdY(rays))
# 
# plt.figure()
# plt.scatter(rays[1],rays[2])
# plt.show()



# from prtp.Missions.OGRE import OgreSource
# from prtp.Missions.OGRE import OgreMirrorModule
# from prtp.Missions.OGRE import OgreGratingModule
# 
# s = OgreSource(num=10000)
# s.order = -1
# s.wave = 4.76 * u.nm
# mm = OgreMirrorModule()
# gm = OgreGratingModule()
# lstack = gm.componentlist[0]
# rstack = gm.componentlist[1]
# 
# for g in lstack.componentlist:
#     norm = g.Normal()
#     g.rotate(theta=0.022*u.deg,ux=norm[0],uy=norm[1],uz=norm[2])
# 
# for g in rstack.componentlist:
#     norm = g.Normal()
#     g.rotate(theta=-0.022*u.deg,ux=norm[0],uy=norm[1],uz=norm[2])
# 
# i = Instrument(s)
# i.addComponent(mm)
# i.addComponent(gm)
# i.addFocus()
# i.simulate()
# 
# rays = i.getRays()
# 
# # print(rays.fwhm())
# # print(rays.hpdY())
# 
# rays.scatter2d()

# rays.scatter3d('dir',c=rays.x)
# rays.scatter3d(c=rays.l)


# 
# # gm.defineRotationPoint()
# 
# leftfwhms = []
# lefthpdys = []
# 
# for j in range(len(lstack.componentlist)):
# 
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(lstack.componentlist[j])
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     leftfwhms.append(rays.fwhm())
#     lefthpdys.append(rays.hpdY())
# 
# rightfwhms = []
# righthpdys = []
# 
# for j in range(len(rstack.componentlist)):
# 
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(rstack.componentlist[j])
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     rightfwhms.append(rays.fwhm())
#     righthpdys.append(rays.hpdY())
#     
# leftfwhms = np.array(leftfwhms)
# lefthpdys = np.array(lefthpdys)
# rightfwhms = np.array(rightfwhms)
# righthpdys = np.array(righthpdys)
# 
# arr = np.vstack((leftfwhms,lefthpdys,rightfwhms,righthpdys)).transpose()
# 
# np.savetxt('Grating_Contributions',arr,header='Left-FWHM Left-HPDY Right-FWHM Right-HPDY',delimiter=' ')



## Documentation Plot Generation
# 
# from prtp.FlatComponent import FlatComponent
# from prtp.Sources import RectBeam
# import astropy.units as u
# 
# f = CollimatorPlate(x=0*u.mm,y=0*u.mm,z=0*u.mm,
#     nx=0,ny=0,nz=1,sx=1,sy=0,sz=0)
# 
# s = RectBeam(num=1000,xhalfwidth=2*u.mm,yhalfwidth=2*u.mm)
# rays = s.generateRays()
# 
# # f.unitrotate(2*u.deg,axis=1)
# 
# f.roll(theta=2*u.deg)
# 
# f.trace(rays)
# rays.scatter3d() 



## Save This: (3D plotting Script)
# # 3D Plot
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(orays.x,orays.y,orays.z,c='r')
# ax.scatter(rays[1],rays[2],rays[3])
# # ax.scatter(optic_rays[1],optic_rays[2],optic_rays[3],c='r')
# # ax.quiver(xcen,ycen,zcen,dir[0],dir[1],dir[2],arrow_length_ratio = .01)
# # ax.quiver(xcen,ycen,zcen,adir[0],adir[1],adir[2],arrow_length_ratio = .01)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

# # Put this block before a 3D plt.show() to make sure all 3 axes are equal
# max_range = np.array([r.x.max()-r.x.min(), r.y.max()-r.y.min(), r.z.max()-r.z.min()]).max()
# Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(r.x.max()+r.x.min())
# Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(r.y.max()+r.y.min())
# Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(r.z.max()+r.z.min())
# # Comment or uncomment following both lines to test the fake bounding box:
# for xb, yb, zb in zip(Xb, Yb, Zb):
#    ax.plot([xb], [yb], [zb], 'w')

## Runtime Profiling:
# import cProfile, pstats
# pr = cProfile.Profile()
# pr.enable()
# 
# # code to profile goes here
# 
# pr.disable()
# sortby = 'cumulative'
# ps = pstats.Stats(pr).sort_stats(sortby)
# ps.print_stats()






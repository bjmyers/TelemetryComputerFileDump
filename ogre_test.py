import numpy as np
import matplotlib.pyplot as plt
from prtp.Missions.OGRE import OgreSource
from prtp.Missions.OGRE import OgreMirrorModule
from prtp.Missions.OGRE import OgreGratingModule
from prtp.Instrument import Instrument
import astropy.units as u

## Basic Script

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()
# gm = OgreGratingModule()
# 
# 
# # for g in gm.getSubComponents():
# #     g.pitch(2*u.deg)
# 
# i = Instrument(s)
# i.addComponent(mm)
# i.addComponent(gm)
# # i.addFocus()
# i.simulate()
# 
# rays = i.getRays()
# 
# x = rays.focusX()
# 
# print(rays.fwhm())
# print(rays.hpdY())
# 
# rays.scatter2d()


# ##############################
# ----------------------------##
# ------BEGIN-TESTING---------##
# ----------------------------##
# ##############################

# Done:
# - Optic Assembly Alignment
# - Module Alignment
# - Grating-Grating Alignment

from tqdm import tqdm

def f(x):
    '''
    This function is used to print out the values in a format that is easy to
        copy and paste into a spreadsheet
    '''
    for i in x:
        print(i)


## Optic X Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# gm = OgreGratingModule()
# 
# mis = np.linspace(-1.27,1.27,21) * u.mm
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
# 
#     mm = OgreMirrorModule()
#     
#     mm.translate(dx=m)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Optic Z Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# gm = OgreGratingModule()
# 
# mis = np.linspace(-1.27,1.27,21) * u.mm
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
# 
#     mm = OgreMirrorModule()
#     
#     mm.translate(dz=m)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Optic Pitch Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# gm = OgreGratingModule()
# 
# mis = np.linspace(-300,300,21) * u.arcsec
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
# 
#     mm = OgreMirrorModule()
#     
#     mm.unitrotate(theta=m,axis=1)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())

## Optic Yaw Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# gm = OgreGratingModule()
# 
# mis = np.linspace(-300,300,21) * u.arcsec
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
# 
#     mm = OgreMirrorModule()
#     
#     mm.unitrotate(theta=m,axis=2)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Module Pitch Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# # Find average position
# gxs = []
# gys = []
# gzs = []
# for g in gm.getSubComponents():
#     gxs.append(g.x.value)
#     gys.append(g.y.value)
#     gzs.append(g.z.value)
# 
# mis = np.linspace(-300,300,21) * u.arcsec
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     gm.defineRotationPoint(x=np.mean(gxs)*u.mm,y=np.mean(gys)*u.mm,z=np.mean(gzs)*u.mm)
#     
#     gm.unitrotate(theta=m,axis=1)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Module Yaw Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# # Find average position
# gxs = []
# gys = []
# gzs = []
# for g in gm.getSubComponents():
#     gxs.append(g.x.value)
#     gys.append(g.y.value)
#     gzs.append(g.z.value)
# 
# mis = np.linspace(-300,300,21) * u.arcsec
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     gm.defineRotationPoint(x=np.mean(gxs)*u.mm,y=np.mean(gys)*u.mm,z=np.mean(gzs)*u.mm)
#     
#     gm.unitrotate(theta=m,axis=2)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Module Roll Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# # Find average position
# gxs = []
# gys = []
# gzs = []
# for g in gm.getSubComponents():
#     gxs.append(g.x.value)
#     gys.append(g.y.value)
#     gzs.append(g.z.value)
# 
# mis = np.linspace(-300,300,21) * u.arcsec
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     gm.defineRotationPoint(x=np.mean(gxs)*u.mm,y=np.mean(gys)*u.mm,z=np.mean(gzs)*u.mm)
#     
#     gm.unitrotate(theta=m,axis=3)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())

## Module X Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(-1.27,1.27,21) * u.mm
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     gm.translate(dx=m)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())
    

## Module Y Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(-1.27,1.27,21) * u.mm
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     gm.translate(dy=m)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Module Z Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(-1.27,1.27,21) * u.mm
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     gm.translate(dz=m)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Grating-Grating Pitch Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(0,60,31)
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     for g in gm.getSubComponents():
#         g.pitch(theta=np.random.normal(0,m)*u.arcsec)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Grating-Grating Roll Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(0,30,20)
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     for g in gm.getSubComponents():
#         g.roll(theta=np.random.normal(0,m)*u.arcsec)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Grating-Grating Yaw Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(0,60,31)
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     for g in gm.getSubComponents():
#         g.yaw(theta=np.random.normal(0,m)*u.arcsec)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Grating-Grating X Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(0,1.27,21)
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     for g in gm.getSubComponents():
#         g.translate(dx=np.random.normal(0,m)*u.mm)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Grating-Grating Y Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(0,1.27,21)
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     for g in gm.getSubComponents():
#         g.translate(dy=np.random.normal(0,m)*u.mm)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Grating-Grating Z Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(0,1.27,21)
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     for g in gm.getSubComponents():
#         g.translate(dz=np.random.normal(0,m)*u.mm)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())


## Stack-Stack Pitch Misalignments
    
# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# # Find average left stack position
# lgxs = []
# lgys = []
# lgzs = []
# for g in gm.componentlist[0].componentlist:
#     lgxs.append(g.x.value)
#     lgys.append(g.y.value)
#     lgzs.append(g.z.value)
# 
# # Find average right stack position
# rgxs = []
# rgys = []
# rgzs = []
# for g in gm.componentlist[1].componentlist:
#     rgxs.append(g.x.value)
#     rgys.append(g.y.value)
#     rgzs.append(g.z.value)
# 
# 
# mis = np.linspace(-100,100,21)*u.arcsec
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     gm.componentlist[0].defineRotationPoint(x=np.mean(lgxs)*u.mm,y=np.mean(lgys)*u.mm,z=np.mean(lgzs)*u.mm)
#     gm.componentlist[1].defineRotationPoint(x=np.mean(rgxs)*u.mm,y=np.mean(rgys)*u.mm,z=np.mean(rgzs)*u.mm)
# 
#     gm.componentlist[0].unitrotate(theta=m,axis=1)
#     gm.componentlist[1].unitrotate(theta=-1*m,axis=1)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm())  


## Stack-Stack Yaw Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# # Find average left stack position
# lgxs = []
# lgys = []
# lgzs = []
# for g in gm.componentlist[0].componentlist:
#     lgxs.append(g.x.value)
#     lgys.append(g.y.value)
#     lgzs.append(g.z.value)
# 
# # Find average right stack position
# rgxs = []
# rgys = []
# rgzs = []
# for g in gm.componentlist[1].componentlist:
#     rgxs.append(g.x.value)
#     rgys.append(g.y.value)
#     rgzs.append(g.z.value)
# 
# 
# mis = np.linspace(-100,100,21)*u.arcsec
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     gm.componentlist[0].defineRotationPoint(x=np.mean(lgxs)*u.mm,y=np.mean(lgys)*u.mm,z=np.mean(lgzs)*u.mm)
#     gm.componentlist[1].defineRotationPoint(x=np.mean(rgxs)*u.mm,y=np.mean(rgys)*u.mm,z=np.mean(rgzs)*u.mm)
# 
#     gm.componentlist[0].unitrotate(theta=m,axis=2)
#     gm.componentlist[1].unitrotate(theta=-1*m,axis=2)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm()) 


## Stack-Stack Roll Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# # Find average left stack position
# lgxs = []
# lgys = []
# lgzs = []
# for g in gm.componentlist[0].componentlist:
#     lgxs.append(g.x.value)
#     lgys.append(g.y.value)
#     lgzs.append(g.z.value)
# 
# # Find average right stack position
# rgxs = []
# rgys = []
# rgzs = []
# for g in gm.componentlist[1].componentlist:
#     rgxs.append(g.x.value)
#     rgys.append(g.y.value)
#     rgzs.append(g.z.value)
# 
# 
# mis = np.linspace(-100,100,21)*u.arcsec
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
#     
#     gm.componentlist[0].defineRotationPoint(x=np.mean(lgxs)*u.mm,y=np.mean(lgys)*u.mm,z=np.mean(lgzs)*u.mm)
#     gm.componentlist[1].defineRotationPoint(x=np.mean(rgxs)*u.mm,y=np.mean(rgys)*u.mm,z=np.mean(rgzs)*u.mm)
# 
#     gm.componentlist[0].unitrotate(theta=m,axis=3)
#     gm.componentlist[1].unitrotate(theta=-1*m,axis=3)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm()) 


## Stack-Stack X Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(-1.27,1.27,21)*u.mm
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
# 
#     gm.componentlist[0].translate(dx=m)
#     gm.componentlist[1].translate(dx=-1*m)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm()) 


## Stack-Stack Y Misalignments

# s = OgreSource(num=100000)
# s.order = -1
# s.wave = 4.762969 * u.nm
# mm = OgreMirrorModule()    
# gm = OgreGratingModule()
# 
# mis = np.linspace(-1.27,1.27,21)*u.mm
# 
# xs = []
# ys = []
# hpds = []
# fwhms = []
# 
# for m in tqdm(mis):
#     
#     gm = OgreGratingModule()
# 
#     gm.componentlist[0].translate(dy=m)
#     gm.componentlist[1].translate(dy=-1*m)
#     
#     i = Instrument(s)
#     i.addComponent(mm)
#     i.addComponent(gm)
#     i.addFocus()
#     i.simulate()
#     
#     rays = i.getRays()
#     
#     x,y = rays.centroid()
#     xs.append(x)
#     ys.append(y)
#     hpds.append(rays.hpdY())
#     fwhms.append(rays.fwhm()) 


## Stack-Stack Z Misalignments

s = OgreSource(num=100000)
s.order = -1
s.wave = 4.762969 * u.nm
mm = OgreMirrorModule()    
gm = OgreGratingModule()

mis = np.linspace(-1.27,1.27,21)*u.mm

xs = []
ys = []
hpds = []
fwhms = []

for m in tqdm(mis):
    
    gm = OgreGratingModule()

    gm.componentlist[0].translate(dz=m)
    gm.componentlist[1].translate(dz=-1*m)
    
    i = Instrument(s)
    i.addComponent(mm)
    i.addComponent(gm)
    i.addFocus()
    i.simulate()
    
    rays = i.getRays()
    
    x,y = rays.centroid()
    xs.append(x)
    ys.append(y)
    hpds.append(rays.hpdY())
    fwhms.append(rays.fwhm()) 
    
    
    
    
    
    
    
    
    
    
    
    
    

# TESTED & CONFIRMED: (using point source, rad=0.5)
# surfacesf: flat, sphere, cyl, cylconic, conic, paraxial, torus
# transformationsf: rotatevector, rotateaxis, reflect, refract, transform, itransform, radgrat, radgratw
# woltsurf: wolterprimary, woltersecondary, woltersine

# NOTES:
# cylconic must have nonzero Y to have any effect
# transformationsf-rotateaxis now outputs ux,uy, and uz as well as x,y, and z

# ISSUES:
# cylconic: behaves oddly, sometimes the new version will successfully calculate a photon but the old version will produce nan's
# conic: old version uses 0s while new version uses nan's
# Surfacesf: prtp code can be refactored, elimination method from specialfunctions should be used if possible
# transformationsf.f95: rotatevector only seems to rotate the first photon in the array, all others are left alone
# grat: I couldn't get the old version to work, got error: 'l' not compatible to 'd'. New version gave all nans

# FIXED:
# conic: changed (+) to (-) in denom calculation
# torus: moved more calculations to within the errstate invalid='ignore' block
# rotatevector: changed "axiz" to "axis" in header
# rotatevector: added temporary variables to ensure important values were not overwritten
# prtp (in general): removed "[:]" in several places to make code look cleaner
# transformationsf: removed some unnecessary code
# itransform: Added a negative sign to all of the angles
# woltprimary: fixed typo, 'Fn' should be 'n'
# woltsurf: changed some threshs to 1e-9 (from 1e-10) to reduce number of nans, error will not be mostly below 1e-10 even with 1000 iterations
# Added focus functions to surfacesf.py

import numpy as np
import matplotlib.pyplot as plt
import time
from PyXFocus import sources
import PyXFocus.surfaces as oldsurf
import prtp.surfacesf as newsurf
import PyXFocus.transformations as oldtrans
import prtp.transformationsf as newtrans
import PyXFocus.woltsurf as oldwolt
import prtp.woltsurf as newwolt
import PyXFocus.conicsolve as conic
import PyXFocus.analyses as analyses
import prtp.Rays as Rays

def printrays(rays, stop = 1):
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    print("opd: " + str(opd[:stop]))
    print("x:   " + str(x[:stop]))
    print("y:   " + str(y[:stop]))
    print("z:   " + str(z[:stop]))
    print("l:   " + str(l[:stop]))
    print("m:   " + str(m[:stop]))
    print("n:   " + str(n[:stop]))
    print("ux:  " + str(ux[:stop]))
    print("uy:  " + str(uy[:stop]))
    print("uz:  " + str(uz[:stop]))

def shallowcopy(rays):
    num = len(rays[0])
    output = np.zeros((10,num))
    for i in range(10):
        for j in range(num):
            output[i][j] = rays[i][j]
    return output

def checkeq(rays1,rays2):
    num = len(rays[0])
    for i in range(10):
        for j in range(num):
            if (rays1[i][j] != rays2[i][j]):
                print(i)
                print(j)
                return False
    return True

def getfirstray(rays):
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    return [opd[0],x[0],y[0],z[0],l[0],m[0],n[0],ux[0],uy[0],uz[0]]

def checkerror(ray1,ray2):
    errors = []
    for i in range(len(ray1)):
        errors.append((ray1[i] - ray2[i]) / ray1[i])
    for e in errors:
        print(e)

num = 1

rays = sources.pointsource(.5,num)
rays[2] += 0
rays[3] += 0
rays[7] += 0
rays[8] += 0
rays[9] += 0

opd,x,y,z,l,m,n,ux,uy,uz = rays


## Testing below this line
amp = 10
freq = 1


printrays(rays,3)
print('-------------')

rays1 = shallowcopy(rays)
opd1,x1,y1,z1,l1,m1,n1,ux1,uy1,uz1 = rays1
rays1 = Rays.Rays(x1,y1,z1,l1,m1,n1,ux1,uy1,uz1)
rays = [opd,x,y,z,l,m,n,ux,uy,uz]
# Modify rays here:
oldsurf.conic(rays,1,2)

printrays(rays,3)
print('-----------')

# Modify rays1 here:
rays1.conic(1,2)


# r1 = getfirstray(rays)
# r2 = getfirstray(rays1)
# 
# print(checkeq(rays,rays1))


# rays1 = shallowcopy(rays)
# opd,x,y,z,l,m,n,ux,uy,uz = rays
# 
# 
# printrays(rays,3)
# rays1 = newwolt.wolterprimary(rays1,r0,z0,psi)
# 
# print('----New:----')
# printrays(rays1,3)
# oldwolt.wolterprimary(x,y,l,m,n,ux,uy,uz,r0,z0,psi,num)
# print('----Old:----')
# printrays(rays,3)
# 
# print('------------')
# print(checkeq(rays,rays1))

















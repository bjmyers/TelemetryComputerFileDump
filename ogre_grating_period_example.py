import astropy.units as u
import numpy as np

# Define points on either end of the grating array.
x3300 = np.arange(-100, 101) * 160. * u.nm
x3200 = np.arange(-100, 101) * 155.151515 * u.nm
y3300 = 3300. * u.mm
y3200 = 3200 * u.mm

# Find run/rise.
islope = (x3300 - x3200) / (y3300 - y3200)

# Define middle of grating.
y3250 = 3250 * u.mm

ys = np.linspace(3200,3300,10000) * u.mm

conserved = True
for y in ys:
    xs = islope * y
    diffs = np.diff(xs)
    if len(np.unique(np.round(diffs.value,decimals=4))) > 1:
        print('Period not Conserved at y=' + str(y))
        conserved = False
        break

if conserved:
    print('WOOOO! Period was always conserved! YEAHHH!')

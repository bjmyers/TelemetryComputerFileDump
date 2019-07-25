import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

with open('attenuation_lens') as f:
    c = f.readlines()

def func(ev):
    return 5.86902987e-06*ev**2 + -4.12807400e-03*ev + 9.79687932e-01

evs = []
lens = []

for i in c:
    e,l = i[3:].split()
    evs.append(float(e))
    lens.append(float(l))

evs = np.array(evs)
lens = np.array(lens)

# p = np.polyfit(evs,lens,2)
# p,temp = curve_fit(func,evs,lens,p0=[.1,.000001,1])

plt.figure()
plt.scatter(evs,lens,s=5,marker='x',c='r')
plt.plot(evs,func(evs))
plt.show()
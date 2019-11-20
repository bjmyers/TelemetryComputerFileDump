from numpy import *
from scipy import io 
from scipy.ndimage.measurements import center_of_mass
import numpy as np
from scipy.optimize import curve_fit,leastsq
from astropy.modeling import models,fitting
import lmfit
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from photutils import datasets, DAOStarFinder


def CenterOfMass1(im):
    
    vec = im
    r,c = shape(im)
    
    """break into a bunch of columns"""
    columns = hsplit(vec,c)
    """sum each column"""
    summedx = array([sum(i) for i in columns])
    print(summedx)
    """find the max"""
    m = max(summedx)
    """WHERE IT AT THO"""
    x = where(summedx==m)
    
    """break into a bunch of rows"""
    rows = vsplit(vec,r)
    """sum each column"""
    summedy = array([sum(i) for i in rows])
    print(summedy)
    """find the max"""
    m = max(summedy)
    """WHERE IT AT THO"""
    y = where(summedy==m)
    
    return y,x


def CenterOfMass2(im):
    
    vec = im
    rows,cols = shape(im)
    
    y = 0
    x = 0
    
    per = vec/sum(vec)
    
    y,x = center_of_mass(vec)
    
    return y,x

def gauss(x,a,mu,sigma):
    return (a/np.sqrt(2*np.pi*sigma*sigma))*np.exp(-(x-mu)**2/(2*sigma**2))

def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (sqrt(2*pi) * wid)) * exp(-(x-cen)**2 / (2*wid**2))
    
def line(x, slope, intercept):
    """a line"""
    return slope*x + intercept

def gaussian2(x, amp, cen, wid, a, b):
    return amp * np.exp(-(x-cen)**2 / wid) + a * x + b

def CenterOfMass3(im):
    
    arr = im+1
    
    r,c = arr.shape
    
    x = np.array(np.hsplit(arr,c))
    y = np.array(np.vsplit(arr,r))
    
    xra = np.arange(0,c)
    yra = np.arange(0,r)
    
    sum_x = np.array([np.sum(i) for i in x])
    sum_y = np.array([np.sum(j) for j in y])
    
    max_locx = np.where(sum_x==sum_x.max())[0]
    max_locy = np.where(sum_y==sum_y.max())[0]
    
    sigx = sum_x.std()
    sigy = sum_y.std()
    
    wx = 1./np.sqrt(sum_x)
    wx[np.isinf(wx)] = 0
    wx[np.isnan(wx)] = 0
    wy = 1./np.sqrt(sum_y)
    wy[np.isinf(wy)] = 0
    wy[np.isnan(wy)] = 0
    
    # gmod = lmfit.models.GaussianModel()
    # # lmod = lmfit.models.LinearModel()
    # 
    # # gmod -= lmod
    # 
    # pargx = gmod.guess(sum_x, x=xra)
    # # pargx = gmod.make_params(amplitude=1000,center=max_locx,sigma=sigx,slope=5,intercept=500)
    # gx = gmod.fit(sum_x, pargx, x=xra)
    # 
    # pargy = gmod.guess(sum_y, x=yra)
    # # pargy = gmod.make_params(amplitude=1000,center=max_locy,sigma=sigy,slope=5,intercept=500)
    # gy = gmod.fit(sum_y, pargy, x=yra)
    # 
    # # gmod += lmod

    mod = lmfit.Model(gaussian) + lmfit.Model(line) + lmfit.models.ConstantModel()
    
    pars = mod.make_params(amp=np.max(sum_x),cen=max_locx,wid=sigx,slope=0,intercept=0,c=np.median(sum_x))
    result = mod.fit(sum_x, pars, x=xra)
    
    comps = result.eval_components()
    
    plt.figure()
    plt.plot(xra,sum_x,'b*',label='data x')
    plt.plot(xra,result.best_fit,'r--')
    plt.plot(xra,comps['gaussian'],'g--')
    plt.plot(xra,comps['line'],'g--')
    plt.show()
    
    # plt.figure()
    # plt.plot(xra,sum_x,'b*',label='data x')
    # plt.plot(yra,sum_y,'r*',label='data y')
    # plt.plot(xra,gx.best_fit,'b-', label='lmfit x')
    # plt.plot(yra,gy.best_fit,'r--',label='lmfit y')
    # plt.xlabel('Position')
    # plt.ylabel('Counts')
    # plt.legend()
    # plt.show()
    
    x = np.mean(np.where(gx.best_fit==np.max(gx.best_fit)))
    y = np.mean(np.where(gy.best_fit==np.max(gy.best_fit)))
    
    return y,x


def CenterOfMass4(im):
    
    data = im 
    
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)    
    
    print((mean, median, std))    
    
    daofind = DAOStarFinder(fwhm=50, threshold=5.*std)    
    
    sources = daofind(data - median)   
    
    sources = np.array(sources)
    
    flux = np.array([sources[i][9] for i in range(len(sources))])
    f_m = np.max(flux)
    ind = np.where(flux==f_m)
    
    return sources[ind][0][2],sources[ind][0][1]

def CenterOfMass5(im):
    
    arr = im+1
    
    r,c = arr.shape
    
    x = np.array(np.hsplit(arr,c))
    y = np.array(np.vsplit(arr,r))
    
    xra = np.arange(0,c)
    yra = np.arange(0,r)
    
    sum_x = np.array([np.sum(i) for i in x])
    sum_y = np.array([np.sum(j) for j in y])
    
    max_locx = np.where(sum_x==sum_x.max())[0][0]
    max_locy = np.where(sum_y==sum_y.max())[0][0]
    
    sigx = sum_x.std()
    sigy = sum_y.std()
    
    gmodel = lmfit.Model(gaussian2)
    #pargx = gmodel.guess(sum_x, x=xra)
    print(np.max(sum_x),max_locx,sigx)
    gx = gmodel.fit(sum_x, x=xra,amp=np.max(sum_x),cen=max_locx,wid = sigx,a=0.,b=0.)
    #pargy = gmodel.guess(sum_y, x=yra)
    gy = gmodel.fit(sum_y, x=yra,amp=np.max(sum_y),cen=max_locy,wid = sigy,a=0.,b=0.)
    
    plt.figure()
    plt.plot(xra,sum_x,'b*',label='data x')
    plt.plot(yra,sum_y,'r*',label='data y')
    plt.plot(xra,gx.best_fit,'b-', label='lmfit x')
    plt.plot(yra,gy.best_fit,'r--',label='lmfit y')
    plt.xlabel('Position')
    plt.ylabel('Counts')
    plt.legend()
    plt.show()
    
    x = np.mean(np.where(gx.best_fit==np.max(gx.best_fit)))
    y = np.mean(np.where(gy.best_fit==np.max(gy.best_fit)))
    
    return y,x

if __name__ == '__main__':
    # abc = io.loadmat('C:/Users/PSU_Telemetry/Documents/WORKPLACE(CHRIS)/Python/Test Images/Test Images/gaussPixelTest')
    # 
    # abc = abc['gFilter']
    
    abc = np.loadtxt('C:/Users/PSU_Telemetry/Documents/WORKPLACE(CHRIS)/Python/Test Images/Test Images/plot_1_22.txt')
    
    # plt.figure()
    # plt.imshow(abc)
    # plt.show()
    
    y,x = CenterOfMass5(abc)
    
    print('x =',x)
    print('y =',y)
    
    
    
    
    
    
    
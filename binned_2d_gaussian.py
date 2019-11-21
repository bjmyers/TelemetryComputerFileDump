import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from tqdm import tqdm
import lmfit 

@np.vectorize
def weightedgauss(r):
    '''
    Weighted Gauss is a function that needs to be integrated by the split() 
        function. If a Gaussian function is defined as G(r), the weighted
        gaussian is W(r) = r*G(r).
    '''
    # Assumes mean of 0 and a standard deviation of zero
    return r * np.exp(-.5 * r**2) / np.sqrt(2 * np.pi)


def split(pix, x, y, sigma, steps = 50):
    '''
    Function split:
    Calculates how the charge cloud produced by each photon spill out into
        adjacent pixels. The algorithm works by considering small slices
        of the charge cloud and figuring out when the slice moves into 
        adjacent pixels. Using the area of the slice and a Gaussian integral,
        each slice will deposit charge into nearby pixels.
    
    Inputs:
    rays - The photons you want to simulate
    steps - The number of slices you want to simulate. A smaller number is
        faster but less accurate. Accuracy is most important for large 
        charge clouds.
    
    Outputs:
    arr - An array the same size as the Detector's pixel array, containing
        the charge cloud binned into pixels. This array is added to the 
        noise array by the view() function.
    '''

    # Initialize the list of slices we will be using in the future
    thetas = np.linspace(1e-9,2*np.pi,steps)
    
    # Find the angle that once slice subtends
    arc = np.mean(np.diff(thetas) / (2*np.pi))
    
    for i in range(len(thetas)):
        
        theta = thetas[i]
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        
        # Find the coordinates of the pixel that each photon starts at
        xpos = int(x)
        ypos = int(y)
        
        # Find the offset of each photon's initial position from the center
        # of the pixel
        dx = x % 1
        dy = y % 1
        
        # We need to keep track of how many vertical and horizontal steps
        # that each photon has undergone so far
        num_horiz_steps = 0
        num_vert_steps = 0
        
        # This will keep track of the distance between the cloud's center 
        # and the previous pixel boundary
        prev_r = 0
        
        # X and Y step are either 1 or -1 that depends on the direction of theta
        ystep = int(costheta / np.abs(costheta))
        xstep = int(sintheta / np.abs(sintheta))
        
        # This will keep track of the distance between the cloud's center 
        # and the next pixel boundary
        r = 0
        
        while True:
            
            # calculate the two potential step sizes
            next_horiz_step = (((2*num_horiz_steps + 1) * ystep / 2) - dx) / np.abs(costheta)
            next_vert_step = (((2*num_vert_steps + 1) * xstep / 2) - dy) / np.abs(sintheta)
            
            # Find the photons which need to take a horizontal step
            # --------------------------------------------------------
            needhoriz = np.abs(next_horiz_step) <= np.abs(next_vert_step)
            
            # update the positions of the photons
            if needhoriz:
                r = next_horiz_step
                
                pix[xpos,ypos] += quad(lambda x: weightedgauss(x), prev_r/sigma, r/sigma)[0] / steps
                
                # Update ypositions
                ypos += ystep
                num_horiz_steps += 1
            
            # Find the photons which need to take a vertical step
            # --------------------------------------------------------
            needvert = np.abs(next_vert_step) < np.abs(next_horiz_step)
            
            if needvert:
                
                # update the positions of the photons
                r = next_vert_step
                
                pix[xpos,ypos] += quad(lambda x: weightedgauss(x), prev_r/sigma, r/sigma)[0] / steps
                
                # Update ypositions
                xpos += xstep
                num_vert_steps += 1
            
            prev_r = r
            
            if (abs(r) > sigma*2.355):
                break
    
    return pix

def fast_split(pix, x, y, sigma, num = 50):
    
    xs = np.random.normal(loc=x,scale=sigma,size=num).astype(int)
    ys = np.random.normal(loc=y,scale=sigma,size=num).astype(int)
    
    for i in range(num):
        if (xs[i] > 0 and xs[i] < len(pix) and ys[i] > 0 and ys[i] < len(pix[0])):
            pix[xs[i],ys[i]] += 1
    
    return pix


def gaussian(x, amp, cen, wid, const):
    """1-d gaussian: gaussian(x, amp, cen, wid, const)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2)) + const


## Generating Basic Data

# pix = np.zeros(shape=(1600,1600))
# 
# pix = split(pix,800.5,700.5,sigma=10,steps=100)
# 
# pix += np.random.normal(loc=0.0005,scale=0.00005,size=(1600,1600))
# 
# plt.figure()
# plt.imshow(pix)
# plt.show()

## Generating Realistic Data

pix = np.zeros(shape=(1600,1600))

# Make Real Event:
pix = fast_split(pix,np.random.rand()*1600,np.random.rand()*1600,sigma=30,num=2000)

# Add between 1 and 6 false events
for i in range(np.random.random_integers(1,6)):
    pix = fast_split(pix,np.random.rand()*1600,np.random.rand()*1600,sigma=np.random.rand()*12+3,num=2000)

# Add Hot Pixels
for i in range(np.random.random_integers(1,3)):
    pix[np.random.random_integers(1,1600),np.random.random_integers(1,1600)] += 15

# Add Noise
pix += np.random.normal(loc=.2,scale=0.1,size=(1600,1600))

plt.figure()
plt.imshow(pix)
plt.show()

cols = pix.sum(axis=0)
rows = pix.sum(axis=1)

plt.figure()
plt.plot(cols)
plt.show()

## Fitting Data


# cols = pix.sum(axis=0)
# rows = pix.sum(axis=1)
# 
# r = np.arange(1600)
# 
# 
# gmodel = lmfit.Model(gaussian)
# 
# 
# gx = gmodel.fit(cols,x=r,amp=np.max(cols),cen=710,wid=2,const = 0.01)
# 
# plt.figure()
# plt.plot(r,cols,label='Data')
# plt.plot(r,gx.best_fit,label='Model')
# plt.legend()
# plt.show()




















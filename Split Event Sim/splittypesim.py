import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm



def split(d = 10,w = 16,steps = 150,dx = 0,dy = 0):
    sigma = d / (2 * np.sqrt(2*np.log(2)))
    thetas = np.linspace(0.01,2*np.pi,steps)

    gridsize = 7 # Must be Odd
    
    arr = np.zeros((gridsize,gridsize))
    

    for i in range(len(thetas)):
        theta = thetas[i]
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        
        maxstep = gridsize//2 + 1
        # NOTE: When accessing numpy arrays, x and y are flipped from their cartesian orders
        xpos = gridsize//2
        ypos = gridsize//2
        
        num_horiz_steps = 0
        num_vert_steps = 0
        prev_r = 0
        
        # X and Y step are either 1 or -1 that depends on the direction of theta
        ystep = int(costheta / np.abs(costheta))
        xstep = int(sintheta / np.abs(sintheta))
        
        while True:
            
            # calculate the two potential step sizes
            next_horiz_step = (((2*num_horiz_steps + 1) * w * ystep / 2) - dx) / np.abs(costheta)
            next_vert_step = (((2*num_vert_steps + 1) * w * xstep / 2) - dy) / np.abs(sintheta)
            
            if np.abs(next_horiz_step) > d and np.abs(next_vert_step) > d:
                # Final Step
                r = d
                arr[xpos,ypos] += (stats.norm.cdf(r/sigma) - stats.norm.cdf(prev_r/sigma)) * (r**2 - prev_r**2)
                break
            
            if np.abs(next_horiz_step) < np.abs(next_vert_step):
                # Take a horizontal step

                r = np.abs(next_horiz_step)
                arr[xpos,ypos] += (stats.norm.cdf(r/sigma) - stats.norm.cdf(prev_r/sigma)) * (r**2 - prev_r**2)
                
                # Move Y-Position by 1 or -1 (without a conditional)
                ypos += ystep
                
                num_horiz_steps += 1
                
            else:
                # Take a vertical step

                r = np.abs(next_vert_step)
                arr[xpos,ypos] += (stats.norm.cdf(r/sigma) - stats.norm.cdf(prev_r/sigma)) * (r**2 - prev_r**2)
                
                # Move X-Position by 1 or -1 (without a conditional)
                xpos += xstep
                
                num_vert_steps += 1
            
            prev_r = r
    
    # Normalize
    arr /= np.sum(arr)
    
    return arr

## Does the same thing as splittypesim.py but for each call to split, finds how many pixels the charge
## cloud is distributed over, does this for many events to find how common each event type is

num = 1000
d = 10 #FWHM in microns
w = 16 #Pixel width in microns (square with dimensions wxw)
xs = (np.random.rand(num) - .5) * (w/2)
ys = (np.random.rand(num) - .5) * (w/2)

counts = np.zeros(10)

for i in tqdm(range(len(xs))):
    arr = split(d=d,w=w,dx=xs[i],dy=ys[i],steps=100)
    x = (arr > .05) # The .05 indicates that a pixel will be included if it contains at least 5% of the total charge cloud
    c = np.count_nonzero(x)
    counts[c] += 1

for i in range(1,len(counts)):
    print(str(i) + '-pixel events: ' + str(counts[i]) + ' --- ' + str(counts[i] / np.sum(counts) * 100) + "%")






















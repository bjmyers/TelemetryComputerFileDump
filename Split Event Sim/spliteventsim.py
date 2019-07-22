import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


def split(d = 10,w = 16,steps = 200,dx = 0,dy = 0):
    '''
    Complicated Algorithm to distribute a charge cloud over a square array of pixels.
    Basically generates a number of slices of the circle and finds when the slices go into
    adjacent pixels to know how much charge to add to each pixel. If you're curious about the
    algorithm ask Bailey, he has the math in his lab notebook
    '''
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

    
num = 100

d = 10 #FWHM in microns
w = 16 #Pixel width in microns (square with dimensions wxw)

# Random offsets used to change the center of the charge randomly throughout the center pixel
xs = (np.random.rand(num) - .5) * (w/2)
ys = (np.random.rand(num) - .5) * (w/2)

# Pixel array
arr = np.zeros((7,7))

for i in tqdm(range(len(xs))):
    arr += split(d=d,w=w,dx=xs[i],dy=ys[i],steps=100)

# Normalize
arr /= np.sum(arr)

plt.figure()
plt.imshow(arr,origin='lower')
plt.title('Sum of Offset Rays')
plt.show()

























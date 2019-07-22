import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


# Create a num*num grid which will take up a 3x3 pixel area, should be odd
# Changing num changes the resolution of the final image
num = 75
arr = np.zeros((num,num))

# Pixel and charge cloud parameters
w = 16  # Width of each pixel
d = 10  # Radius of Charge cloud
sigma = d / (2 * np.sqrt(2*np.log(2))) # Standard deviation of gaussian given that d=FWHM

# Scale is the width of each block in the num*num grid
scale = (3*w) / num

##-----------------------------------------------------------------------------
# Generates a charge cloud for a single event
centx = num//2 
centy = num//2

for i in range(num):
    for j in range(num):
        r = np.sqrt((i*scale - centx*scale)**2 + (j*scale - centy*scale)**2)
        arr[i,j] += 1 - stats.norm.cdf(r/sigma)


plt.figure()
plt.imshow(arr)
plt.title('Singular Central Event')
plt.axvline(x=num//3, color='r')
plt.axvline(x=2*num//3,color='r')
plt.axhline(y=num//3, color='r')
plt.axhline(y=2*num//3,color='r')
plt.colorbar()
plt.show()
##-----------------------------------------------------------------------------



##-----------------------------------------------------------------------------
# Sums up several events

# Number of events which will be summed up
num_events = 250

for i in tqdm(range(num_events)):
    centx = (num//2) + ((np.random.rand() - .5) * num//3 * scale)
    centy = (num//2) + ((np.random.rand() - .5) * num//3 * scale)
    
    for i in range(num):
        for j in range(num):
            r = np.sqrt((i*scale - centx*scale)**2 + (j*scale - centy*scale)**2)
            arr[i,j] += 1 - stats.norm.cdf(r/sigma)

plt.figure()
plt.imshow(arr)
plt.title('Sum of ' + str(num_events) + ' events')
plt.axvline(x=num//3, color='r')
plt.axvline(x=2*num//3,color='r')
plt.axhline(y=num//3, color='r')
plt.axhline(y=2*num//3,color='r')
plt.colorbar()
plt.show()
##-----------------------------------------------------------------------------







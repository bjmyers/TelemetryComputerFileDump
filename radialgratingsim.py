import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import astropy.units as u

## Functions

def islope(pos,focus):
    '''
    pos - initial position of the groove (could eventually be calculated using the index, pos = init_sep * index
    focus - distance to the focus from bottom of wafer
    
    returns: inverse slope (i.e; run/rise). Multiplying this by the change in height will give you the change in x
    '''
    return pos / focus

def propagate(posns,dist,islopes=[]):
    '''
    posns - initial position of the groove positions
    dist - the distance by which you want to move the gr ooves
    
    returns - a list of the new groove positions
    '''
    return posns - islopes*dist

def fracture(posns,nearest=0.3):
    '''
    posns - The list of current groove positions
    nearest - What you are rounding to. By default, rounds to the nearest 3 Angstroms (0.3 nm)
    
    returns - a new array, with each value rounded to the nearest specified value
    '''
    return np.round(posns / nearest) * nearest

def findprobs(arr):
    '''
    Inputs an array, finds the abundance of each unique element. For example,
    if arr contains the values 500 and 501, this finds the percentage of indexes
    that contain 500 and the percentage of indexes that contain 501.
    Output - formatted array of the form:
    [(element1,prob),(element2,prob)...(elementn,prob)]
    '''
    l = len(arr)
    elements = np.unique(arr)
    output = []

    for x in elements:
        prob = np.count_nonzero(arr==x)/l
        output.append((x,prob))
    
    return output
    
    
'''
To calculate the separation between grooves, use np.diff(posns)
'''


## Constants
num = 10000     # Number of grooves
init_sep = 160   # Initial separation of grooves in nm
focus = .33e10    # Distance to focus in nm
inner_rad = 165
outer_rad = 189.7




## Code:


## Plot average separation as a function of position
# means = []
# devs = []
# indexes = np.arange(num)
# posns = ((init_sep) * indexes).astype(np.float64)
# isl = islope(posns,focus)
# 
# xs = np.arange(0,int(1e8),int(1e5))
# for i in xs:
#     newposns = propagate(posns,i*0.3,isl)
#     newposns = fracture(newposns)
#     means.append(np.mean(np.diff(newposns)))
#     devs.append(np.std(np.diff(newposns)))
# 
# plt.figure()
# plt.scatter(xs,means,c=devs)
# plt.xlabel('Displacement (nm)')
# plt.ylabel('Mean Separation (nm)')
# plt.colorbar()
# plt.show()


## Finds Period Distribution at each Position
# indexes = np.arange(num)
# posns = ((init_sep) * indexes).astype(np.float64)
# isl = islope(posns,focus)
# periods = []
# 
# xs = np.arange(0,int(1.1e8),int(1e5))
# for i in xs:
#     newposns = propagate(posns,i,isl)
#     newposns = fracture(newposns,.1)
#     p = fracture(np.diff(newposns),.1)
#     periods.append(findprobs(p))
# 
#     
# # Saves the two arrays to files:
# np.save('positions_01',xs)
# np.save('Probabilities_01',periods)


## Given a Position Array, Returns a Grating Period
# xs = np.load('positions_3.npy')
# periods = np.load('Probabilities_3.npy')
# 
# def gratPeriod(x):
#     indexes = np.searchsorted(xs,x)
#     indexes = np.where(indexes == len(indexes), indexes-1,indexes)
#     selection = np.random.rand(len(x))
#     output = []
#     for i in range(len(x)):
#         if selection[i] < periods[indexes[i]][0][1]:
#             output.append(periods[indexes[i]][0][0])
#         else:
#             output.append(periods[indexes[i]][1][0])
#     return output
# 
# def idealGratPeriod(x):
#     initial_sep = periods[0][0][0] * periods[0][0][1] + periods[0][1][0] * periods[0][1][1]
#     final_sep = periods[-1][0][0] * periods[-1][0][1] + periods[-1][1][0] * periods[-1][1][1]
#     sep_slope = (final_sep - initial_sep) / (xs[-1] - xs[0])
#     return ((x - xs[0]) * sep_slope) + initial_sep
# 
# 
# x = np.random.rand(1000) * 1e8
# ys = gratPeriod(x)
# 
# plt.figure()
# plt.scatter(x,ys)
# plt.show()


## New Code (E-Beam Method)

## E-Beam Method - Plot average separation as a function of position
# means = []
# devs = []
# indexes = np.arange(num)
# posns = ((init_sep) * indexes).astype(np.float64)
# isl = islope(posns,focus)
# 
# xs = np.arange(0,int(1e8),int(1e5))
# for i in xs:
#     newposns = propagate(posns,i*0.3,isl)
#     newposns = fracture(newposns,0.1)
#     means.append(np.mean(np.diff(newposns)))
#     devs.append(np.std(np.diff(newposns)))
# 
# plt.figure()
# plt.scatter(xs,means,c=devs)
# plt.xlabel('Displacement (nm)')
# plt.ylabel('Mean Separation (nm)')
# plt.colorbar()
# plt.show()

## Finding distribution of Periods as a function of grid size and height over a large y-disp

# indexes = np.arange(num)
# posns = ((init_sep) * indexes).astype(np.float64)
# isl = islope(posns,focus)
# 
# xstep = .5  # Grid size in the x-direction in nm
# ystep = .5  # Grid size in the y-direction (towards focus) in nm
# 
# # Starting Position
# y = 0
# 
# # We need the groove positions along two horizontal steps, then we can compare their separations in between the steps
# 
# # Initialize the lists that will store our means and stds
# ys = []
# tot_means = []
# tot_stds = []
# 
# from tqdm import tqdm
# 
# # This loop dictates how many vertical steps we want to simulate
# for i in tqdm(range(1000)):
#     
#     # Define Current Positions
#     currentxs = propagate(posns,y,isl)
#     currentxs = fracture(currentxs,xstep)
# 
#     # Get the Next Positions
#     nextxs = propagate(posns,y+ystep,isl)
#     nextxs = fracture(nextxs,xstep)
#     
#     # Find if the grooves have "jumped" an x-step
#     deltax = nextxs - currentxs
#     
#     # Initialize the mean and std period for this step
#     means = []
#     stds = []
#     
#     # This loop dictates how many random slices will be taken
#     for j in range(100):
#     
#         # Find what fraction of the y-step you want to move (0 to 1)
#         partialstep = np.random.rand()
#         
#         # Find how far the grooves have travelled (if at all) 
#         newposns = currentxs + partialstep*deltax/ystep
#         
#         # Find the period of each groove
#         periods = np.diff(newposns)
#         
#         # Store the mean and standard deviations of this sample
#         means.append(np.mean(periods))
#         stds.append(np.std(periods))
#     
#     # Store combined means and standard deviations:
#     tot_means.append(np.mean(means))
#     tot_stds.append(np.mean(stds))
#     # stds = np.array(stds)
#     # tot_stds.append(np.sqrt(np.sum(stds**2)))
#     ys.append(y)
#         
#         
#     # Now we're done with this step, go to the next y-step:
#     y += ystep * 1000
#     currentxs = nextxs
# 
# ys = np.array(ys)
# ys *= 1e-6
# 
# 
# fig,ax1 = plt.subplots()
# color = 'tab:red'
# ax1.set_xlabel('Y-Displacement (mm)')
# ax1.set_ylabel('Mean Groove Period (nm)',color=color)
# ax1.plot(ys,tot_means,color=color)
# ax1.tick_params(axis='y', labelcolor=color)
# 
# ax2 = ax1.twinx() 
# color = 'tab:blue'
# ax2.set_ylabel('Standard Deviation Groove Period (nm)', color=color)
# ax2.plot(ys, tot_stds, color=color)
# ax2.tick_params(axis='y', labelcolor=color)
# 
# fig.tight_layout()
# plt.show()



## Find Period distribution on one horizontal slice
# indexes = np.arange(num)
# # Move the indexes farther over to see more extreme slopes
# indexes += 1000000
# posns = ((init_sep) * indexes).astype(np.float64)
# isl = islope(posns,focus)
# 
# xstep = .1  # Grid size in the x-direction in nm
# ystep = .1  # Grid size in the y-direction (towards focus) in nm
# 
# # Starting Position
# y = (10*u.mm).to(u.nm).value
# 
# # We need the groove positions along two horizontal steps, then we can compare their separations in between the steps
# 
# # Define Current Positions
# currentxs = propagate(posns,y,isl)
# # currentxs = fracture(currentxs,xstep)
# 
# # Get the Next Positions
# nextxs = propagate(posns,y+ystep,isl)
# # nextxs = fracture(nextxs,xstep)
# 
# # Find if the grooves have "jumped" an x-step
# deltax = nextxs - currentxs
# 
# # Initialize the mean and std period for this step
# periods = []
# 
# # This loop dictates how many random slices will be taken
# for j in range(1000):
# 
#     # Find what fraction of the y-step you want to move (0 to 1)
#     partialstep = np.random.rand()
#     
#     # Find how far the grooves have travelled (if at all) 
#     newposns = currentxs + partialstep*deltax/ystep
#     
#     # Find the period of each groove
#     ds = np.diff(newposns)
#     
#     periods.extend(ds)
# 
# 
# plt.figure()
# plt.hist(periods)
# plt.show()


import astropy.units as u
from tqdm import tqdm

indexes = np.arange(1e3)
indexes *= 125
indexes2 = indexes + 1

posns = ((init_sep) * indexes).astype(np.float64)
posns2 = ((init_sep) * indexes2).astype(np.float64)
# x-Positions go from 0 - 20mm
isl = islope(posns,focus)
isl2 = islope(posns2,focus)

xstep = .5  # Grid size in the x-direction in nm
ystep = (1.63 * u.um).to(u.nm).value  # Grid size in the y-direction (towards focus) in nm

# Starting Position
y = (00*u.mm).to(u.nm).value

# Initialize bins
bins = np.linspace(157.0,160.5,num=1000) 
myhist = np.zeros(999, dtype='int32')

yend = (50*u.mm).to(u.nm)
m = np.floor(((yend-y) / ystep).value).astype('int32')

for j in tqdm(range(m)):
    
    y += ystep
    
    # Define Current Positions
    currentxs = propagate(posns,y,isl)
    currentxs2 = propagate(posns2,y,isl2)
    # currentxs = fracture(currentxs,xstep)
    # currentxs2 = fracture(currentxs2,xstep)
    
    # Get the Next Positions
    nextxs = propagate(posns,y+ystep,isl)
    nextxs2 = propagate(posns2,y+ystep,isl2)
    # nextxs = fracture(nextxs,xstep)
    # nextxs2 = fracture(nextxs2,xstep)
    
    # Find if the grooves have "jumped" an x-step
    deltax = nextxs - currentxs
    deltax2 = nextxs2 - currentxs2
    
    # n - Number of samples within eachsubsection
    n = 50
    for i in range(n):
        # Sample within a subsection
    
        # Find what fraction of the y-step you want to move (0 to 1)
        # partialstep = np.random.rand()
        
        # Find how far the grooves have travelled (if at all) 
        # newposns = currentxs + partialstep*deltax
        # newposns2 = currentxs2 + partialstep*deltax2
        newposns = currentxs + (i/n)*deltax
        newposns2 = currentxs2 + (i/n)*deltax2
        
        # Find the period of each groove
        temp, junk = np.histogram(newposns2 - newposns,bins) 
        myhist += temp



fig, ax = plt.subplots()
ax.bar(bins[:-1],myhist,width = (160.5-157)/1000)
plt.show()











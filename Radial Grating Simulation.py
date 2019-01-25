import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

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
    dist - the distance by which you want to move the grooves
    
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


'''
To calculate the separation between grooves, use np.diff(posns)
'''


## Constants
num = 1000     # Number of grooves
init_sep = 160   # Initial separation of grooves in nm
focus = .33e10    # Distance to focus in nm
inner_rad = 165
outer_rad = 189.7




## Code:

## Plots the number of unique separations at different propogation values
# indexes = np.arange(num)
# posns = ((init_sep) * indexes).astype(np.float64)
# isl = islope(posns,focus)
# 
# seps = []
# n = 5000
# for i in tqdm(range(n),desc='Propogating Grooves'):
#     
#     newposns = (fracture(posns))
#     
#     newposns = propagate(newposns, 0.3 * i, isl)
#     
#     newposns = fracture(newposns)
#     sep = np.diff(newposns)
#     
#     sep = fracture(sep,0.1)
#     
#     seps.append(len(np.unique(sep)))
# 
# 
# plt.figure()
# plt.scatter(np.arange(n),np.array(seps))
# plt.show()

## Plots the positions of the grooves over several propogation values
# num = 10
# focus = 1
# indexes = np.arange(num)
# posns = ((init_sep) * indexes).astype(np.float64)
# isl = islope(posns,focus)
# xs,ys = [],[]
# n = 10
# 
# for i in range(n):
#     newposns = (fracture(posns))
#     newposns = propagate(newposns, 0.3 * i, isl)
#     
#     newposns = fracture(newposns)
#     
#     xs.append(newposns)
#     ys.append([0.3 * i]*num)
# 
# 
# plt.figure()
# plt.scatter(xs,ys)
# plt.show()
    
## Propogate the entire grating
indexes = np.arange(num)
posns = ((init_sep) * indexes).astype(np.float64)
isl = islope(posns,focus)

newposns = posns

newposns = propagate(newposns, 1e8, isl)

error = np.random.normal(30,5,(num))
newposns += error

newposns = fracture(newposns)
sep = np.diff(newposns)

# sep = fracture(sep,0.1)

plt.figure()
plt.hist(sep)
plt.show()



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

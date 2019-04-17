import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
num = 1000     # Number of grooves
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
indexes = np.arange(num)
posns = ((init_sep) * indexes).astype(np.float64)
isl = islope(posns,focus)
periods = []

xs = np.arange(0,int(1.1e8),int(1e5))
for i in xs:
    newposns = propagate(posns,i,isl)
    newposns = fracture(newposns,.1)
    p = fracture(np.diff(newposns),.1)
    periods.append(findprobs(p))

    
# Saves the two arrays to files:
np.save('positions_01',xs)
np.save('Probabilities_01',periods)


## Given a Position Array, Returns a Grating Period
xs = np.load('positions_01.npy')
periods = np.load('Probabilities_01.npy')

def gratPeriod(x):
    indexes = np.searchsorted(xs,x)
    indexes = np.where(indexes == len(indexes), indexes-1,indexes)
    selection = np.random.rand(len(x))
    output = []
    for i in range(len(x)):
        if selection[i] < periods[indexes[i]][0][1]:
            output.append(periods[indexes[i]][0][0])
        else:
            output.append(periods[indexes[i]][1][0])
    return output

def idealGratPeriod(x):
    initial_sep = periods[0][0][0] * periods[0][0][1] + periods[0][1][0] * periods[0][1][1]
    final_sep = periods[-1][0][0] * periods[-1][0][1] + periods[-1][1][0] * periods[-1][1][1]
    sep_slope = (final_sep - initial_sep) / (xs[-1] - xs[0])
    return ((x - xs[0]) * sep_slope) + initial_sep


x = np.random.rand(1000) * 1e8
# x = np.linspace(0,1e8,1000)
ys = gratPeriod(x)

plt.figure()
plt.scatter(x,ys)
plt.show()


















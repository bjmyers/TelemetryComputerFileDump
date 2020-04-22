import numpy as np

def rotatevector(vec,theta,axis):
    '''
    Function rotatevector:
    Takes in a vector, an angle to rotate, and a unit vector about which to
        rotate
    
    Inputs:
    vec - The vector you want to rotate as a tuple (x,y,z)
    theta - The angle through which the vector is to be rotated, must be in 
        radians
    axis - Specifies the axis about which to rotate:
        1 - Rotation about x-axis
        2 - Rotation about y-axis
        Other - Rotation about z-axis
    
    Outputs:
    (x, y, z) - The components of the rotated vector
    
    Notes:
    - Rotations follow the right-hand rule
    '''
    
    x, y, z = vec
    
    if (axis==1):
        tempy = np.cos(theta)*y - np.sin(theta)*z
        z = np.sin(theta)*y + np.cos(theta)*z
        return x,tempy,z
    elif (axis==2):
        tempx = np.cos(theta)*x + np.sin(theta)*z
        z = -np.sin(theta)*x + np.cos(theta)*z
        return tempx,y,z
    else:
        tempx = np.cos(theta)*x - np.sin(theta)*y
        y = np.sin(theta)*x + np.cos(theta)*y
        return tempx,y,z


def rotateaxis(vec,theta,axis):
    '''
    Function rotateaxis:
    Takes in a vector, an angle to rotate, and a vector about which to rotate
    
    Inputs:
    vec - The vector you want to rotate as a tuple (x,y,z)
    theta - The angle through which the vector is to be rotated, must be in 
        radians
    axis - The axis about which you want to rotate as a tuple (x,y,z)
    
    Outputs:
    (x, y, z) - The components of the rotated vector
    
    Notes:
    - Rotations follow the right-hand rule
    '''
    x, y, z = vec
    ux, uy, uz = axis

    # Normalize the axis
    mag = np.sqrt(ux**2 + uy**2 + uz**2)
    with np.errstate(invalid='ignore'):
        ux /= mag
        uy /= mag
        uz /= mag
    
    c = np.cos(theta)
    s = np.sin(theta)
    tempx = (c+ux**2*(1-c))*x + (ux*uy*(1-c)-uz*s)*y + (ux*uz*(1-c)+uy*s)*z
    tempy = (uy*ux*(1-c)+uz*s)*x + (c+uy**2*(1-c))*y + (uy*uz*(1-c)-ux*s)*z
    z = (uz*ux*(1-c)-uy*s)*x + (uz*uy*(1-c)+ux*s)*y + (c+uz**2*(1-c))*z
    return tempx,tempy,z



def findLED(init_vec, x, y, d):
    '''
    Function findLED:
    Given the position of the LED on the image and an initial vector pointing
        towards the central calibrated point, find the coordinates of the LED
        in x, y, and z.
    
    Inputs:
    init_vec - Tuple of the form (x, y, z), gives the components of the initial
        vector which points from the center of the detector to the central
        calibrated point. This should be in Hexapod Coordinates
    x - The x-coordinate of the LED, with the calibrated center at x=0
    y - The y-coordinate of the LED, with the calibrated center at y=0
    d - The distance to the LED, calculated with the distance function
    
    Output:
    (x, y, z) - The coordinates of the LED in Hexapod coordinates
    '''
    
    # Find the angles phi and theta
    phi = findPhi(x,y)
    theta = findTheta(x,y)
    
    # To find the LEDs coordinates, we need two transformations. First rotate
    # by theta about the x-axis, then rotate by phi about init_vec
    
    c_prime = rotatevector(init_vec,theta,1)
    
    return rotateaxis(c_prime,phi,init_vec)


# *** IMPORTANT: This function depends on the directions of the x and y axes in
# the image. When someone can go in and look at the image analysis script, they
# should test this function.
def findPhi(x, y):
    '''
    Function findPhi:
    Given the coordinates of the LED on the image, finds the angle Phi as 
        described in the metrology summary
    
    Inputs:
    x - The x-coordinate of the LED, with the calibrated center at x=0
    y - The y-coordinate of the LED, with the calibrated center at y=0
    
    Outputs:
    phi - The angle phi (in radians)
    '''
    return np.arctan(x/y)

def findTheta(x, y):
    '''
    Function findTheta:
    Given the coordinates of the LED on the image, finds the angle Theta as 
        described in the metrology summary
    
    Inputs:
    x - The x-coordinate of the LED, with the calibrated center at x=0
    y - The y-coordinate of the LED, with the calibrated center at y=0
    
    Outputs:
    theta - The angle theta (in radians)
    
    Notes:
    - This function needs to be calibrated with proper values for a and b, 
        it follows the general form described in the metrology summary
    '''
    # Initialize constants
    a = 1
    b = 1
    
    # Find the distance between the central calibrated point and the LED
    d = np.sqrt(x**2 + y**2)
    
    # Find theta
    return a*np.arctan(b*d)




# Define calibrated central vector (a vector pointing from the detector to
# the calibrated central point on the plate). This vector must be in Hexapod
# coordinates
init_vec = (0,1,0)

# Define the coordinates of the LED on the image (in pixels)
x = 50
y = 50

# Define the distance to the LED, in the future this will be done with a
# distance function
d = 1

ledx, ledy, ledz = findLED(init_vec, x, y, d)










#
# Reimplementing /Users/phillicl/HoomdWork/smac/src/shpdescfourier.cpp
# So not reliant on SMAC
# See /Users/phillicl/HoomdWork/smac/unit_tests/fourierdescriptor2d_test.cpp for test
#   based on passing polygons to descriptor

import numpy as np
import coordinate_helper

# Calculating the Invariant 2D Fourier Descriptor (weights are 1.0)
def fourier2D(coordinates,L,frequencies):


    # Note, that I do not perform the calculation below on the central particle
    # Messes with results!
    
    N = coordinates.shape[0]-1;
    
    # Note.  The first coordinates is the center particle.
    #Translate all coordinates so the center is the origin
    coordinates = coordinate_helper.pbc(coordinates[1:,:]-coordinates[0,:],L)

    #RSQ is only calculated in SMAC to partition particles by shells.
    #We are not partitioning particles by shells
    #rsq = coordinate_helper.rsq(coordinates)
    
    # calculate theta of all particles (including center)
    theta = np.arctan2(coordinates[:,1],coordinates[:,0]) + np.pi
    
    shape_descriptor = []
    
    # For each frequency.. calculate a complex number which
    # is sum [np.cos(k*theta), -np.sin(k*theta)] for each theta
    # Applying weight of 1.0 to all
    for k in frequencies:
    
        components = np.zeros(shape=(coordinates.shape[0],2))

        #components += [np.cos(k*theta), -np.sin(k*theta)]
        # weights
        w=1.0
        components[:,0] += w*np.cos(k*theta)
        components[:,1] += -w*np.sin(k*theta)
    
        # c is a complex number
        c= np.sum(components,axis=0)
        
        # scale the components
        c /=float(N)
        
        
        #print c
        # because invariant (return magnitude of complex number)
        c = np.sqrt(np.sum(c*c))
    
        shape_descriptor.append(c)

    return shape_descriptor
    
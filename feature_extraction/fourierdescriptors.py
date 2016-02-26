#
# Reimplementing /Users/phillicl/HoomdWork/smac/src/shpdescfourier.cpp
# So not reliant on SMAC
# See /Users/phillicl/HoomdWork/smac/unit_tests/fourierdescriptor2d_test.cpp for test
#   based on passing polygons to descriptor

import numpy as np
import scipy.special as sp
import coordinate_helper


#
# Pad the shape descriptor out to a new length with zeros
#
def pad_sd(sd, new_length):

    if new_length > len(sd):
        tmp = np.zeros(new_length)
        tmp[:len(sd)]=sd
        sd = tmp
    
    return sd

#
# Calculating a normalized distance as is done in SMAC
#
def distance_func(sd1, sd2):

    # Just in case the two shape descriptors are not the same length, pad with zeros
    if len(sd1) != len(sd2):
        new_length = np.max([len(sd1),len(sd2)])
        sd1 = pad_sd(sd1,new_length)
        sd2 = pad_sd(sd2,new_length)

    max_difference = np.linalg.norm(sd1) + np.linalg.norm(sd2)
    residual = np.linalg.norm(sd1-sd2)

    rnorm = residual/max_difference

    return rnorm

#
# Takes the normalized distance and returns the similarity instead
# This is what is returned by matchfundist in SMAC.
#
def similarity_func(sd1,sd2,min=0,max=1.0):
    rnorm = distance_func(sd1,sd2)
    return min + (max-min)*(1-rnorm)


#
# Calculating the Invariant 2D Fourier Descriptor (weights are 1.0)
#
def fourier2D(coordinates,L,frequencies):


    # Note, that I do not perform the calculation below on the central particle
    # Messes with results!
    
    N = coordinates.shape[0]-1;
    
    if N == 0:
        return [0]
    
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


#
# Calculating the Invariant 3D Fourier Descriptor (weights are 1.0)
#
def fourier3D(coordinates,L,frequencies):


    # Note, that I do not perform the calculation below on the central particle
    # Messes with results!
    
    N = coordinates.shape[0]-1;
    
    if N == 0:
        return [0]
    
    # Note.  The first coordinates is the center particle.
    #Translate all coordinates so the center is the origin
    coordinates = coordinate_helper.pbc(coordinates[1:,:]-coordinates[0,:],L)

    #-- STEP 1 - Convert to spherical  coordinates
    # calculate the radial distance of all the coordinates
    r = np.linalg.norm(coordinates)                            #Radial
    phi = np.arctan2(coordinates[:,1],coordinates[:,0]) +np.pi #Azimuthal
    theta = np.arccos(coordinates[:,2]/r)                      #Inclination

    #-- STEP 2: compute Ylm(theta, phi) for each theta, phi

    # Using spherical harmonic function provided by scipy
    #http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.special.sph_harm.html

    # scipy.special.sph_harm(m, l, azimuthal, polar) azimuthal in [0,2pi], polar in [0,pi]  (n == l)

    # NEED MORE WORK BELOW

    # In the original smac code, for a given frequency l, loops over the n particles, and over m values from -l to +l, and then applies a "q" invarinace that combines all the m values generated for all the particles.

    shape_descriptor= []
    for l in frequencies:

        components = np.zeros(shape=(coordinates.shape[0],2*l+1,2))
        
        for m in range(-l,l+1):
            w = 1.0
            ylm = sp.sph_harm(m,l,phi,theta)
            components[:,m,0] = w*ylm.real
            components[:,m,1] = w*ylm.imag


        # Summing over all the particles
        ylm_average = np.sum(components,axis=0)/components.shape[0]

        # Applying invariance
        # computes the "q" invariant for a 3d fourier coefficient
        sd = np.sqrt(4.0*np.pi/ylm.shape[1]*np.sum(ylm_average**2))
        shape_descriptor.append(sd)


    return shape_descriptor



    
#
#
#
import numpy as np

def pbc(coordinates,L):

    # If a single coordinate gets handed, expand it so slicing works
    if len(coordinates.shape) ==1:
        coordinates=np.expand_dims(coordinates,axis=0)
    
    for i in range(3):
        mask = coordinates[:,i] < -L[i]/2.
        coordinates[mask,i] +=L[i]
    
        mask = coordinates[:,i] > L[i]/2.
        coordinates[mask,i] -=L[i]

    return coordinates


def rsq(coordinates):
    return np.sum(coordinates*coordinates,axis=1)
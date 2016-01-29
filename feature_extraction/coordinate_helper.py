#
#
#


# Helper Function to map back into box  (coordinates must be numpy array)
def pbc(coordinates,L):

    for i in range(3):
        mask = coordinates[:,i] < -L[i]/2.
        coordinates[mask,i] +=L[i]/2.
    
        mask = coordinates[:,i] > L[i]/2.
        coordinates[mask,i] -=L[i]/2.


    return coordinates


def rsq(coordinates):
    return np.sum(coordinates*coordinates,axis=1)
# 01/28/2016 - Carolyn L. Phillips
# Very simple minded cluster extraction
# No optimization (e.g. cell list)
# But sufficient to get the job done for testing
import numpy as np
import coordinate_helper
import copy

class Cluster:

    def __init__(self,coordinates,L,N):
        self.coordinates = coordinates
        self.L = L
        self.N = N

    # Returns all clusters defined by a central particle and all particles
    # within a distance rmax
    def getClusters(self,rmax):

        clusters = []

        for i in range(self.N):
            
            mask_notme = [True]*(self.N)
            mask_notme[i] = False
            
            dd = self.coordinates - self.coordinates[i]
            dd = coordinate_helper.pbc(dd,self.L)
            rij = np.linalg.norm(dd, axis=1)

            mask = (rij < rmax) & mask_notme
            
            myparticles = self.coordinates[mask]
            
            tmp = np.zeros(shape=(myparticles.shape[0]+1, 3))
            tmp[0,:] =copy.copy(self.coordinates[i])
            tmp[1:,:]= copy.copy(myparticles)
            
            clusters.append(tmp)

        return clusters


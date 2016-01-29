import numpy as np
import matplotlib.pyplot as plt
import extractdata
import rdf
import extractClusters
import fourierdescriptors as fd

X=extractdata.parseHOOMDxml()

X.parse('ljg_2.10_25_5.0.xml')

X.printStats()


rdf = rdf.RDF(rmax = 5., dr = 0.01)

x_i, hist = rdf.calculate(coordinates = X.getCoordinates(),
                          L = X.getBox(),
                          N = X.getNumberParticles())

# Show a plot of the histogram that I just created
plt.figure("Radial Distribution Function")
plt.plot(x_i,hist, 'b-')
plt.xlabel('r')
plt.ylabel('RDF')


# peakline is dependent on the scaling of the histogram!  So either need to get this right or pick a new peak line
#  Same is true for the other parameters!!!
valley,peaks = rdf.getPeaks(15., num_peaks=2)


# Extract Clusters
C = extractClusters.Cluster(coordinates = X.getCoordinates(),
                          L = X.getBox(),
                          N = X.getNumberParticles())


print "Using First Valley as cutoff r_cutoff:",valley[0].radial_position
clusters = C.getClusters(valley[0].radial_position)

# Generate Cluster Size Fingerprint
cluster_sizes = [len(x) for x in clusters]
max_cluster_size = 12
fingerprint =np.zeros(max_cluster_size+1)
for i in range(max_cluster_size+1):
    fingerprint[i] = cluster_sizes.count(i)

fingerprint /= sum(fingerprint)

print fingerprint

frequencies= [3,4,5,6,7,8,9,10,11,12]

# Generating Invariant Frequencies using a center Particle
sd = fd.fourier2D(clusters[0],L = X.getBox(), frequencies=frequencies)

square = np.zeros(shape=(5,3))
square[0,:] = np.array([0,0,0])
square[1,:] = np.array([1,0,0])
square[2,:] = np.array([-1,0,0])
square[3,:] = np.array([0,1,0])
square[4,:] = np.array([0,-1,0])

print
print
sd = fd.fourier2D(square,L = X.getBox(), frequencies=frequencies)


print sd

plt.show()



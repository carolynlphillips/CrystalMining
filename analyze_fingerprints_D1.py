import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
import matplotlib.pyplot as plt
import matplotlib.cm as cm

test_eps_values = False
min_pts = 7
metric = "manhattan"

# Reading file of Neighbor size histogram
with open("fingerprints/D1_fingerprint.dat",'r') as f:
    data = np.loadtxt(f)

print "Fingerprints read:", data.shape[0]

r = data[:,0];
epsilon = data[:,1]
features = data[:,2:]

#  My Features are bins of a histogram.  This removes the mean and scales each feature
#  to unit variance.
features = StandardScaler().fit_transform(features)

# This line below is used without the Standaridzed Scaler and
# reflects how clustering was performed in the original paper
# Note that the HON and DEC RT regions don't separate unless eps is very low.
#db = DBSCAN(eps=0.027, min_samples=1).fit(features)

#
if test_eps_values :
    eps_range = np.linspace(0.05,1,30)
    silhouette_score = np.zeros_like(eps_range)
    num_clusters = np.zeros_like(eps_range)
    fraction_noise = np.zeros_like(eps_range)

    for i,eps in enumerate(eps_range):
        print eps
        db = DBSCAN(eps=eps, min_samples=min_pts, metric=metric).fit(features)
        labels = db.labels_
        silhouette_score[i] = metrics.silhouette_score(features, labels)
        num_clusters[i] = len(set(labels)) - (1 if -1 in labels else 0)
        fraction_noise[i] = sum(labels==-1)/float(len(labels))

    plt.figure("Silhouette Score")
    plt.plot(eps_range,silhouette_score, 'bo-')
    plt.xlabel('eps')
    plt.ylabel('Silhouette Score')
    plt.figure("Number of Clusters")
    plt.plot(eps_range,num_clusters, 'ro-')
    plt.xlabel('eps')
    plt.ylabel('Number of Clusters')
    plt.figure("Fraction Noise")
    plt.plot(eps_range,fraction_noise, 'go-')
    plt.xlabel('eps')
    plt.ylabel('Fraction Noise')

# This line is an experiment used with the Standardized Scaler
# The value of eps = 0.57 was selected because it was a local peak in the Silhouette Score
db = DBSCAN(eps=0.57, min_samples=min_pts, metric=metric).fit(features)



#np.zeros_like returns an array of zeros the same size as that provided, but in this case. booleans (False) instead of zeros)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

print('Estimated number of clusters: %d' % n_clusters_)


# There are 50 values of epsilon  (0.1 to 5)
ncolumns = 50;
# There are 110 values of r (1.01 to 2.01)
nrows = 110;
grid = labels.reshape((nrows,ncolumns))

#extent=(epsilon.min(), epsilon.max(), r.max(), r.min())

##############################################################################
# Plot result

plt.figure("Phase Diagram")
# Using a Qualitative Color Map from http://matplotlib.org/examples/color/colormaps_reference.html
plt.imshow(np.flipud(grid.T), extent=(r.min(), r.max(), epsilon.min(), epsilon.max()),
           interpolation='nearest', cmap=cm.Paired,  aspect=0.1)

plt.ylabel('epsilon')
plt.xlabel('r0')
plt.show()
############################

# Should check the fingerprint being generated for the rhombus data. eps = 4.8 r = 1.9
# Its being grouped with square
# Compare to eps = 5 r = 1.2


# Most of the measures from the example require having ground truth (i.e. labels_true) to compute
'''
print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
print("Adjusted Rand Index: %0.3f"
      % metrics.adjusted_rand_score(labels_true, labels))
print("Adjusted Mutual Information: %0.3f"
      % metrics.adjusted_mutual_info_score(labels_true, labels))
'''

# The silhoette coefficient is the mean intra-cluster distance minus the mean nearest cluster distance (scaled by the max of both).  A value of zero means overlapping clusters.
print("Silhouette Coefficient: %0.3f"
      % metrics.silhouette_score(features, labels))


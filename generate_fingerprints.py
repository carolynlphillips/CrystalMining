import numpy as np
import matplotlib.pyplot as plt
import feature_extraction.extractdata as extractdata
import feature_extraction.rdf as rdf
import feature_extraction.extractClusters as extractClusters
import feature_extraction.fourierdescriptors as fd
import glob

X=extractdata.parseHOOMDxml()
rdf = rdf.RDF(rmax = 5., dr = 0.01)

max_cluster_size_d1 = 8
max_cluster_size_d2 = 16

simulation_files = glob.glob("Engel_dataset/*.xml")


for nfile, file in enumerate(simulation_files):
    
    # Extract parameters from file names
    split_name = file[:-4].split('_')
    r = split_name[2]
    epsilon = split_name[4]
    

    X.parse(file)
    X.printStats()

    # Calculate RDF
    x_i, hist = rdf.calculate(coordinates = X.getCoordinates(),
                              L = X.getBox(),
                              N = X.getNumberParticles())

    # Get the valleys and peaks of the rdf
    # Using Parameters of Soft Matter Paper (Phillips,2013) - Section 3.1.1
    valley,peaks = rdf.getPeaks(peak_line=15., integrated_peak_fraction =0.1, fraction_valley= 0.5, num_peaks=2)

    # Extract particle clusters based on the first valley
    C = extractClusters.Cluster(coordinates = X.getCoordinates(),
                          L = X.getBox(),
                          N = X.getNumberParticles())

    ####  Extract Clusters using First Valley #####

    clusters = C.getClusters(valley[0].radial_position)

    # Generate cluster size fingerprint
    cluster_sizes = [len(x) for x in clusters]
    
    # Just a check that we set the max_cluster_size appropriately
    if np.max(cluster_sizes) > max_cluster_size_d1:
        print "Warning - Found a cluster of size larger than ", max_cluster_size_d1

    # Create fingerprint
    fingerprint =np.zeros(max_cluster_size_d1+1)
    for i in range(max_cluster_size_d1+1):
        fingerprint[i] = cluster_sizes.count(i)
    fingerprint /= sum(fingerprint)


    # Write to file (don't want to lose our work if it crashes)
    outline = np.zeros(shape=max_cluster_size_d1+3)
    outline[0] = r
    outline[1] = epsilon
    outline[2:] = fingerprint
    outline = np.expand_dims(outline,axis=0)
    with open("D1_fingerprint.dat",'a') as f:
        np.savetxt(f,outline)


    ####  Repeat using Second Valley #####
    clusters = C.getClusters(valley[1].radial_position)

    # Generate cluster size fingerprint
    cluster_sizes = [len(x) for x in clusters]
    
    # Just a check that we set the max_cluster_size appropriately
    if np.max(cluster_sizes) > max_cluster_size_d2:
        print "Warning - Found a cluster of size", np.max(cluster_sizes)," which is larger than ", max_cluster_size_d2

    # Create fingerprint
    fingerprint =np.zeros(max_cluster_size_d2+1)
    for i in range(max_cluster_size_d2+1):
        fingerprint[i] = cluster_sizes.count(i)
    fingerprint /= sum(fingerprint)


    # Write to file (don't want to lose our work if it crashes)
    outline = np.zeros(shape=max_cluster_size_d2+3)
    outline[0] = r
    outline[1] = epsilon
    outline[2:] = fingerprint
    outline = np.expand_dims(outline,axis=0)
    with open("D2_fingerprint.dat",'a') as f:
        np.savetxt(f,outline)


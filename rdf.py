# 01/28/2016 - Carolyn L. Phillips
# This is an very simple implementation of a radial distribution function calculator
# It uses no optimization
# But it is sufficient to support the calculation
from operator import sub
import numpy as np
import copy
import operator

import coordinate_helper

class rdf_position:

    def  __init__(self,index,radial_position,height):
    
        self.index = index
        self.radial_position = radial_position
        self.height = height
        return

    def setValues(self,index,radial_position,height):
    
        self.index = index
        self.radial_position = radial_position
        self.height = height
        return

class RDF:

    def __init__(self,rmax,dr):
        self.rmax = rmax
        self.dr = dr
        self.maxbin = int(rmax/dr)
    
        self.x_i = []
        self.hist = []

    # Calculate the pair distribution function histogram
    # Determine Normalization Factor for each Bin from
    #~/HoomdWork/smac/src/diagnostics.cpp
    
    def calculate(self, coordinates, L, N):

        # Calculate RDF  - Doing this Brute Force for now - no cell list
        
        # radial position of the bin
        x_i = np.array(range(self.maxbin))*self.dr
        hist = np.zeros(self.maxbin)

        # Convienent Lambdas
        makeint = lambda x: int(x)
        
        # Note that currently all distances are being calculated twice.
        # If I change this, need to double the weight of the count in the histogram
        for rollval in range(1,N):
        
            dd = np.array(map(sub, coordinates, np.roll(coordinates,rollval,axis=0)))  # Box Periodicity is not used right now
            
            dd = coordinate_helper.pbc(dd,L)

            rij = np.linalg.norm(dd, axis=1)
            bin = map(makeint, rij/self.dr)

            for b in bin:
                if b < self.maxbin :
                    hist[b] += 1
    
    
        # Calculate the array of normalization factors
        delta = self.dr/2.
        # density
        number_density = float(N)/np.product(L)
        
        volume_shell = 4./3.*np.pi*np.power(x_i+delta, 3) - 4./3.*np.pi*np.power(x_i-delta, 3);
        
        # Correction for first shell
        volume_shell[0] += 4./3.*np.pi*np.power(-delta, 3);
        
        np_ideal_gas = number_density * volume_shell;
        normalization =  np_ideal_gas * N;

        hist /= normalization;

        self.x_i = x_i
        self.hist = hist
        return x_i,hist


    # Based on gofr_helper.h of original code
    # Returns the index and radial position of the valleys and peaks of the histogram
    # If integratedpeak_fraction is provided, then will not include peaks that are too small
    
    def getPeaks(self, peak_line, integratedpeak_fraction=None, num_peaks=None):

        zero=peak_line
    
        # Could make valley and peak structs
        valley = [];
        for i in range(1,self.maxbin):
            if self.hist[i] < zero and not(self.hist[i-1] < zero) :
                #valley.append((i,self.x_i[i]))
                valley.append(rdf_position(i,self.x_i[i],self.hist[i]))
            

        if len(valley)==1:
            print "Warning Simulation may be too hot."
        
        Num_valleys = len(valley)

        # Check the volume under the peak
        # Check the peak height
        integrated_peak =np.zeros(Num_valleys)
        #peak_height = np.zeros(Num_valleys)
        
        peak = [];
        j = 0
        tmp = copy.copy(valley[0])
        
        for i in range(self.maxbin):
            if j < Num_valleys:
                
                if self.hist[i] > zero : integrated_peak[j] += self.hist[i]
                
                #if self.hist[i] > peak_height[j]:
                if self.hist[i] > tmp.height :
                    tmp.setValues(i,self.x_i[i],self.hist[i])

                    #peak[j] = (i,self.x_i[i])
                    
                if i == valley[j].index:
                    peak.append(tmp)
                    j += 1
                    if j < Num_valleys: tmp = copy.copy(valley[j])
                    

        # Print a few statistics to the screen
        peak_height = [x.height for x in peak]
        index, max_peak = max(enumerate(peak_height), key=operator.itemgetter(1))
        print "Peak", index+1, "has the heighest peak of value", max_peak

        index, max_integrated_peak = max(enumerate(integrated_peak), key=operator.itemgetter(1))
        print "Peak", index+1, "has the largest volume", max_integrated_peak

        if integratedpeak_fraction:
            # Remove peaks that are too small
            for j in range(Num_valleys):
                if integrated_peak[j] < integratedpeak_fraction * max_integrated_peak:
                    valley[j]= []
                    peak[j] = []
                    integrated_peak[j] = []

        if num_peaks:
            valley = valley[:num_peaks]
            peak = peak[:num_peaks]

        return valley,peak











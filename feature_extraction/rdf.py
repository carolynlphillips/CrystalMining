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
    # Note.  Produces same g(r) as VMD analysis
    
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
            
            # DEBUGGING
            mask = rij < 0.5
            if sum(mask) :
                print
                print coordinates[mask], np.roll(coordinates,rollval,axis=0)[mask], dd[mask]
            
            
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
    # If integrated_peak_fraction is provided, then will not include peaks that are too small
    
    # Note - An improvement to this method may be to fit a curve to the rdf and extract peaks and valleys from it.  May be more resilient to noise!
    
    def getPeaks(self, peak_line, integrated_peak_fraction=None, num_peaks=None, fraction_valley=None, verbose = False):

        zero=peak_line
    
    
        valley = [];
        maxes = []
        mins = []
        for i in range(1,self.maxbin-1):
        
            # We collect local maximums and local minimums over the range until the rdf drops below the zero line
            
            if self.hist[i] > zero and fraction_valley:
                # (Very) Local Maximum
                if self.hist[i-1] < self.hist[i] and self.hist[i+1] < self.hist[i]:
                    maxes.append(rdf_position(i,self.x_i[i],self.hist[i]))
                
        
                # (Very) Local Minimum
                if self.hist[i-1] > self.hist[i] and self.hist[i+1] > self.hist[i]:
                    mins.append(rdf_position(i,self.x_i[i],self.hist[i]))

            #Valley Criteria #1 - GOFR Drops to Zero
            if self.hist[i] < zero and not(self.hist[i-1] < zero) :
                
                # Valley Criteria #2 - Large Max to Min Drop
                if len(maxes) > 1 :  # (must be at least two maxes)
                    for ind in range(len(maxes)-1):
                    
                        drop = maxes[ind].height - mins[ind].height
                        rise = maxes[ind+1].height - mins[ind].height
                        
                        print rise, drop
                        
                        if drop > maxes[ind].height * fraction_valley and rise > maxes[ind+1].height*fraction_valley :
            
                            valley.append(mins[ind])
            
                            if verbose:
                                print "Criteria 2: Adding Valley at ",valley[-1].radial_position
            
            
                valley.append(rdf_position(i,self.x_i[i],self.hist[i]))
                maxes = []
                mins = []
                if verbose:
                        print "Criteria 1: Adding Valley at ",valley[-1].radial_position
            

        if len(valley)==1:
            print "Warning Simulation may be too hot."
        
        Num_valleys = len(valley)
        
        

        # Check the volume under the peak
        integrated_peak = []
        peak = [];
        j = 0
        tmp = copy.copy(valley[0])
        tmp_integrate = 0


        # Valley Criteria #3 - Must have sufficient integrated volume under prior peak
        for i in range(self.maxbin):
            if j < Num_valleys:
                
                # Only integrating values over zero line
                if self.hist[i] > zero : tmp_integrate += self.hist[i]
                
                # Looking for the peak over the interval between the valleys
                if self.hist[i] > tmp.height :
                    tmp.setValues(i,self.x_i[i],self.hist[i])
                
                if i == valley[j].index:
                    peak.append(tmp)
                    integrated_peak.append(tmp_integrate)
                    tmp_integrate = 0;
                    j += 1
                    if j < Num_valleys: tmp = copy.copy(valley[j])
                    

        # Print a few statistics to the screen
        if verbose:
            peak_height = [x.height for x in peak]
            index, max_peak = max(enumerate(peak_height), key=operator.itemgetter(1))
            print "Peak", index+1, "has the highest peak of value", max_peak, "at radial position", peak[index].radial_position
        
        index, max_integrated_peak = max(enumerate(integrated_peak), key=operator.itemgetter(1))
        if verbose:
            print "Peak", index+1, "has the largest volume", max_integrated_peak, "at radial position", peak[index].radial_position

        if integrated_peak_fraction:
            # Remove peaks that are too small
            for j in range(Num_valleys-1,-1,-1):
                if integrated_peak[j] < integrated_peak_fraction * max_integrated_peak:
                    if verbose: print "Dropping Peak ", j
                    valley.pop(j)
                    peak.pop(j)
                    integrated_peak.pop(j)

        if num_peaks:
            valley = valley[:num_peaks]
            peak = peak[:num_peaks]
    
        if verbose:
            for i, v in enumerate(valley):
                print "Valley Cutoff", i, ":", v.radial_position

        return valley,peak











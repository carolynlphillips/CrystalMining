# 01/28/2016 - Carolyn L. Phillips
#This is an extremely parser that ignores much of the information
# that a HOOMD xml file may contain. It is designed to only get what is needed to support this calculation

import xml.etree.ElementTree as ET
import numpy as np


# Reconfigure this to a class and definition

class parseHOOMDxml:

    def __init__(self):
        self.filename = None
        self.coordinates = []
        self.L = []
        self.N = None
        
        return

    def parse(self,filename):
        self.filename = filename
        tree = ET.parse(filename)
        root = tree.getroot()

        # Extract Coordinates
        text_coordinates = root[0].findall('position')[0].text.splitlines()
        coordinates = []

        convert = lambda x: float(x)
        for i,c in enumerate(text_coordinates):
            c=c.split(" ")
            if len(c) == 3:  # First line can just be empty
                coordinates.append(map(convert,c))

        self.coordinates=np.array(coordinates)
        self.N = len(coordinates)
        
        # THERE ARE MORE EFFICIENT WAYS TO REPACK A LIST INTO A NUMPY ARRAY
        tmp = np.zeros(shape=(self.N,3))
        for i in range(self.N):
            tmp[i,:] = coordinates[i]
        self.coordinates = tmp

        
        #print "Number of Particle Coordinates Found:", self.N

        # Extract Box
        box=root[0].findall('box')
        Lx = float(box[0].attrib['lx'])
        Ly = float(box[0].attrib['ly'])
        Lz = float(box[0].attrib['lz'])
        self.L = [Lx,Ly,Lz]
        #print "Box: ", Lx, Ly, Lz

    def getCoordinates(self):
        return self.coordinates

    def getBox(self):
        return self.L

    def getNumberParticles(self):
        return self.N

    def printStats(self):
        if self.filename:
            print "File:",self.filename
            print "Particles:",self.N
            print "Box:", self.L
        else :
            print "No File Read"
        return



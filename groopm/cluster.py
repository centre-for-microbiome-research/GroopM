#!/usr/bin/env python
###############################################################################
#                                                                             #
#    cluster.py                                                               #
#                                                                             #
#    A collection of classes / methods used when clustering contigs           #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#          .d8888b.                                    888b     d888          #
#         d88P  Y88b                                   8888b   d8888          #
#         888    888                                   88888b.d88888          #
#         888        888d888 .d88b.   .d88b.  88888b.  888Y88888P888          #
#         888  88888 888P"  d88""88b d88""88b 888 "88b 888 Y888P 888          #
#         888    888 888    888  888 888  888 888  888 888  Y8P  888          #
#         Y88b  d88P 888    Y88..88P Y88..88P 888 d88P 888   "   888          #
#          "Y8888P88 888     "Y88P"   "Y88P"  88888P"  888       888          #
#                                             888                             #
#                                             888                             #
#                                             888                             #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################
import sys
import math
import colorsys
import random

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

import numpy as np
import scipy.ndimage as ndi
import scipy.spatial.distance as ssdist
from scipy.stats import kstest

import time

# GroopM imports
import PCA
import mstore

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class DataBlob:
    """Interacts with the groopm datamanager and local data fields
    
    Simple a wrapper around a group of numpy arrays
    """
    def __init__(self, dbFileName, force=False, scaleFactor=1000):
        # data
        self.indicies = np.array([])        # indicies into the data structure based on condition
        self.covProfiles = np.array([])
        self.contigNames = np.array([])
        self.contigLengths = np.array([])
        self.contigColours = np.array([])
        self.kmerSigs = np.array([])
        self.bins = np.array([])
        self.cores = np.array([])
        self.transformedData = np.array([]) # the munged data points
        
        self.condition = ""                 # condition will be supplied at loading time
        self.numContigs = 0                 # this depends on the condition given
        self.numStoits = 0                  # this depends on the data which was parsed

        # misc
        self.dataManager = mstore.GMDataManager()       # most data is saved to hdf
        self.dbFileName = dbFileName        # db containing all the data we'd like to use
        self.forceWriting = force           # overwrite existng values silently?
        self.scaleFactor = scaleFactor      # scale every thing in the transformed data to this dimension

    def loadData(self,
                 condition="",
                 silent=False,
                 loadBins=False,
                 loadCores=False):
        """Load pre-parsed data"""
        try:
            self.numStoits = self.getNumStoits()
            self.condition = condition
            print "\tLoading indicies (", condition,")"
            self.indicies = self.dataManager.getConditionalIndicies(self.dbFileName, condition=condition)
            self.numContigs = len(self.indicies)
            print "\tWorking with:",self.numContigs,"contigs"

            print "\tLoading coverage profiles"
            self.covProfiles = self.dataManager.getCoverageProfiles(self.dbFileName, indicies=self.indicies)

            print "\tLoading kmer sigs"
            self.kmerSigs = self.dataManager.getKmerSigs(self.dbFileName, indicies=self.indicies)

            print "\tCreating colour profiles"
            colourProfile = self.makeColourProfile()
            # use HSV to RGB to generate colours
            S = 1       # SAT and VAL remain fixed at 1. Reduce to make
            V = 1       # Pastels if that's your preference...
            for val in colourProfile:
                self.contigColours = np.append(self.contigColours, [colorsys.hsv_to_rgb(val, S, V)])
            self.contigColours = np.reshape(self.contigColours, (self.numContigs, 3))            

            print "\tLoading contig names"
            self.contigNames = self.dataManager.getContigNames(self.dbFileName, indicies=self.indicies)
            
            print "\tLoading contig lengths"
            self.contigLengths = self.dataManager.getContigLengths(self.dbFileName, indicies=self.indicies)
            
            if(loadBins):
                print "\tLoading bins"
                self.bins = self.dataManager.getBins(self.dbFileName, indicies=self.indicies)

            if(loadCores):
                print "\tLoading core info"
                self.cores = self.dataManager.getCores(self.dbFileName, indicies=self.indicies)
            
        except:
            print "Error loading DB:", self.dbFileName, sys.exc_info()[0]
            raise

#------------------------------------------------------------------------------
# GET / SET 

    def getNumStoits(self):
        """return the value of numStoits in the metadata tables"""
        return self.dataManager.getNumStoits(self.dbFileName)
            
    def getMerColNames(self):
        """return the value of merColNames in the metadata tables"""
        return self.dataManager.getMerColNames(self.dbFileName)
            
    def getMerSize(self):
        """return the value of merSize in the metadata tables"""
        return self.dataManager.getMerSize(self.dbFileName)

    def getNumMers(self):
        """return the value of numMers in the metadata tables"""
        return self.dataManager.getNumMers(self.dbFileName)

### USE the member vars instead!
#    def getNumCons(self):
#        """return the value of numCons in the metadata tables"""
#        return self.dataManager.getNumCons(self.dbFileName)

    def getNumBins(self):
        """return the value of numBins in the metadata tables"""
        return self.dataManager.getNumBins(self.dbFileName)
        
    def setNumBins(self, numBins):
        """set the number of bins"""
        self.dataManager.setNumBins(self.dbFileName, numBins)
        
    def getStoitColNames(self):
        """return the value of stoitColNames in the metadata tables"""
        return self.dataManager.getStoitColNames(self.dbFileName)
    
    def isClustered(self):
        """Has the data been clustered already"""
        return self.dataManager.isClustered(self.dbFileName)
    
    def setClustered(self):
        """Save that the db has been clustered"""
        self.dataManager.setClustered(self.dbFileName, True)
    
    def isComplete(self):
        """Has the data been *completely* clustered already"""
        return self.dataManager.isComplete(self.dbFileName)
    
    def setComplete(self):
        """Save that the db has been completely clustered"""
        self.dataManager.setComplete(self.dbFileName, True)
    
    def saveBins(self, update):
        """Save our bins into the DB"""
        self.dataManager.setBins(self.dbFileName, update)
    
    def saveCores(self, update):
        """Save our core flags into the DB"""
        self.dataManager.setCores(self.dbFileName, update)

#------------------------------------------------------------------------------
# DATA TRANSFORMATIONS 

    def transformData(self, kf=10, sf=1):
        """Perform all the necessary data transformations"""
        # Update this guy now we know how big he has to be
        # do it this way because we may apply successive transforms to this
        # guy and this is a neat way of clearing the data 
        s = (self.numContigs,3)
        self.transformedData = np.zeros(s)
        tmp_data = np.array([])

        print "\tRadial mapping"
        # first we shift the edge values accordingly and then 
        # map each point onto the surface of a hyper-sphere
        # the vector we wish to move closer to...
        radialVals = np.array([])        
        ax = np.zeros_like(self.covProfiles[0])
        ax[0] = 1
        las = self.getAngBetween(ax, np.ones_like(self.covProfiles[0]))
        for point in self.covProfiles:
            norm = np.linalg.norm(point)
            radialVals = np.append(radialVals, norm)
            point /= np.abs(np.log(norm+1)) # make sure we're always taking a log of something greater than 1
            tmp_data = np.append(tmp_data, self.rotateVectorAndScale(point, las, phi_max=8))

        # it's nice to think that we can divide through by the min
        # but we need to make sure that it's not at 0!
        min_r = np.amin(radialVals)
        if(0 == min_r):
            min_r = 1
        # reshape this guy
        tmp_data = np.reshape(tmp_data, (self.numContigs,self.numStoits))
    
        # now we use PCA to map the surface points back onto a 
        # 2 dimensional plane, thus making the data usefuller
        index = 0
        if(self.numStoits == 2):
            print "Skip dimensionality reduction (dim < 3)"
            for point in self.covProfiles:
                self.transformedData[index,0] = tmp_data[index,0]
                self.transformedData[index,1] = tmp_data[index,1]
                self.transformedData[index,2] = math.log(radialVals[index]/min_r)
                index += 1
        else:    
            # Project the points onto a 2d plane which is orthonormal
            # to the Z axis
            print "\tDimensionality reduction"
            PCA.Center(tmp_data,verbose=0)
            p = PCA.PCA(tmp_data)
            components = p.pc()
            for point in components:
                self.transformedData[index,0] = components[index,0]
                self.transformedData[index,1] = components[index,1]
                if(0 > radialVals[index]):
                    self.transformedData[index,2] = 0
                else:
                    self.transformedData[index,2] = math.log(radialVals[index]/min_r)
                index += 1

        # finally scale the matrix to make it equal in all dimensions                
        min = np.amin(self.transformedData, axis=0)
        max = np.amax(self.transformedData, axis=0)
        max = max - min
        max = max / (self.scaleFactor-1)
        for i in range(0,3):
            self.transformedData[:,i] = (self.transformedData[:,i] -  min[i])/max[i]

    def makeColourProfile(self):
        """Make a colour profile based on ksig information"""
        ret_array = np.array([0.0]*np.size(self.indicies))
        working_data = np.array(self.kmerSigs, copy=True) 
        PCA.Center(working_data,verbose=0)
        p = PCA.PCA(working_data)
        components = p.pc()
        
        # now make the colour profile based on PC1
        index = 0
        for point in components:
            ret_array[index] = float(components[index,0])
            index += 1
        
        # normalise to fit between 0 and 1
        ret_array -= np.min(ret_array)
        ret_array /= np.max(ret_array)
        if(False):
            print ret_array
            plt.figure(1)
            plt.subplot(111)
            plt.plot(components[:,0], components[:,1], 'r.')
            plt.show()
        return ret_array
    
    def rotateVectorAndScale(self, point, las, phi_max=6):
        """
        Move a vector closer to the center of the positive quadrant
        
        Find the co-ordinates of its projection
        onto the surface of a hypersphere with radius R
        
        What?...  ...First some definitions:
       
        For starters, think in 3 dimensions, then take it out to N.
        Imagine all points (x,y,z) on the surface of a sphere
        such that all of x,y,z > 0. ie trapped within the positive
        quadrant.
       
        Consider the line x = y = z which passes through the origin
        and the point on the surface at the "center" of this quadrant.
        Call this line the "main mapping axis". Let the unit vector 
        coincident with this line be called A.
       
        Now think of any other vector V also located in the positive
        quadrant. The goal of this function is to move this vector
        closer to the MMA. Specifically, if we think about the plane
        which contains both V and A, we'd like to rotate V within this
        plane about the origin through phi degrees in the direction of
        A.
        
        Once this has been done, we'd like to project the rotated co-ords 
        onto the surface of a hypersphere with radius R. This is a simple
        scaling operation.
       
        The idea is that vectors closer to the corners should be pertubed
        more than those closer to the center.
        
        Set phi max as the divisor in a radial fraction.
        Ie set to '12' for pi/12 = 15 deg; 6 = pi/6 = 30 deg etc
        """
        # the vector we wish to move closer to and unitise
        center_vector = np.ones_like(point)
        center_vector /= np.linalg.norm(center_vector)

        # find the existing angle between them -> theta
        theta = self.getAngBetween(point, center_vector)
        V_p = point
        # theta == 0 means we are already on the center!
        if(0 != theta):
            # at the boundary we want max rotation
            # at the center we like 0 rotation. For simplicity, let's use the logistic function!
            # The denominator produces a value of ~0 if theta == 0 and ~1 is theta == las  
            phi = (np.pi/phi_max) / (2*(1 + np.exp(-1*(2*np.pi)*(theta/las) + np.pi)))  
            
            # now we can find a vector which approximates the rotation of unit(V)
            # by phi. It's norm will be a bit wonky but we're going to scale it anywho...
            V_p = ((point / np.linalg.norm(point)) * ( theta - phi ) + center_vector * phi ) / theta
            norm = np.linalg.norm(V_p) * np.linalg.norm(point)
            if(0 == norm):
                return np.zeros_like(point)
            else:
                return (V_p / norm) 
        return V_p
        
    def getAngBetween(self, P1, P2):
        """Return the angle between two points (in radians)"""
        # find the existing angle between them theta
        c = np.dot(P1,P2)/np.linalg.norm(P1)/np.linalg.norm(P2) 
        # rounding errors hurt everyone...
        if(c > 1):
            c = 1
        elif(c < -1):
            c = -1
        return np.arccos(c) # in radians

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def plotTransViews(self, tag="fordens"):
        """Plot top, side and front views of the transformed data"""
        self.renderTransData(tag+"_top.png",azim = 0, elev = 90)
        self.renderTransData(tag+"_front.png",azim = 0, elev = 0)
        self.renderTransData(tag+"_side.png",azim = 90, elev = 0)

    def renderTransData(self, fileName="", show=True, elev=45, azim=45, all=False):
        """Plot transformed data in 3D"""
        fig = plt.figure()
        if(all):
            myAXINFO = {
                'x': {'i': 0, 'tickdir': 1, 'juggled': (1, 0, 2),
                'color': (0, 0, 0, 0, 0)},
                'y': {'i': 1, 'tickdir': 0, 'juggled': (0, 1, 2),
                'color': (0, 0, 0, 0, 0)},
                'z': {'i': 2, 'tickdir': 0, 'juggled': (0, 2, 1),
                'color': (0, 0, 0, 0, 0)},
            }

            ax = fig.add_subplot(131, projection='3d')
            ax.scatter(self.transformedData[:,0], self.transformedData[:,1], self.transformedData[:,2], edgecolors=self.contigColours, c=self.contigColours, marker='.')
            ax.azim = 0
            ax.elev = 0
            for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
                for elt in axis.get_ticklines() + axis.get_ticklabels():
                    elt.set_visible(False)
            ax.w_xaxis._AXINFO = myAXINFO
            ax.w_yaxis._AXINFO = myAXINFO
            ax.w_zaxis._AXINFO = myAXINFO
            
            ax = fig.add_subplot(132, projection='3d')
            ax.scatter(self.transformedData[:,0], self.transformedData[:,1], self.transformedData[:,2], edgecolors=self.contigColours, c=self.contigColours, marker='.')
            ax.azim = 90
            ax.elev = 0
            for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
                for elt in axis.get_ticklines() + axis.get_ticklabels():
                    elt.set_visible(False)
            ax.w_xaxis._AXINFO = myAXINFO
            ax.w_yaxis._AXINFO = myAXINFO
            ax.w_zaxis._AXINFO = myAXINFO
            
            ax = fig.add_subplot(133, projection='3d')
            ax.scatter(self.transformedData[:,0], self.transformedData[:,1], self.transformedData[:,2], edgecolors=self.contigColours, c=self.contigColours, marker='.')
            ax.azim = 0
            ax.elev = 90
            for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
                for elt in axis.get_ticklines() + axis.get_ticklabels():
                    elt.set_visible(False)
            ax.w_xaxis._AXINFO = myAXINFO
            ax.w_yaxis._AXINFO = myAXINFO
            ax.w_zaxis._AXINFO = myAXINFO
        else:
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.transformedData[:,0], self.transformedData[:,1], self.transformedData[:,2], edgecolors=self.contigColours, c=self.contigColours, marker='.')
            ax.azim = azim
            ax.elev = elev
            ax.set_axis_off()

        if(fileName != ""):
            try:
                if(all):
                    fig.set_size_inches(42,12)
                else:
                    fig.set_size_inches(12,12)            
                plt.savefig(fileName,dpi=300)
                plt.close(fig)
            except:
                print "Error saving image",fileName, sys.exc_info()[0]
                raise
        elif(show):
            try:
                plt.show()
                plt.close(fig)
            except:
                print "Error showing image", sys.exc_info()[0]
                raise
        del fig
        
###############################################################################
###############################################################################
###############################################################################
###############################################################################
class ClusterEngine:
    """Top level interface for clustering contigs"""
    def __init__(self, dbFileName, plot=False, outFile="", force=False):
        # worker classes
        self.dataBlob = DataBlob(dbFileName)
        self.clusterBlob = ClusterBlob(self.dataBlob, debugPlots=plot)
    
        # misc
        self.plot = plot
        self.outFile = outFile
        self.forceWriting = force
        
    def makeCores(self, coreCut, minSize, minVol):
        """Cluster the contigs to make bin cores"""
        # check that the user is OK with nuking stuff...
        if(not self.promptForOverwrite()):
            return False

        # get some data
        t0 = time.time()
        print "Load data"
        self.dataBlob.loadData(condition="length >= "+str(coreCut))
        t1 = time.time()
        print "\tTHIS: [",self.secondsToStr(t1-t0),"]\tTOTAL: [",self.secondsToStr(t1-t0),"]"
        
        # transform the data
        print "Apply data transformations"
        self.dataBlob.transformData()
        # plot the transformed space (if we've been asked to...)
        if(self.plot):
            self.dataBlob.renderTransData()
        t2 = time.time()
        print "\tTHIS: [",self.secondsToStr(t2-t1),"]\tTOTAL: [",self.secondsToStr(t2-t0),"]"
        
        # cluster and bin!
        print "Create cores"
        cum_contigs_used_good = self.clusterBlob.initialiseCores(minVol)
        t3 = time.time()
        print "\tTHIS: [",self.secondsToStr(t3-t2),"]\tTOTAL: [",self.secondsToStr(t3-t0),"]"
        
        # now we assume that some true bins may be separated across two cores
        # try to condense things a little
        #print "Condense cores"
        #self.clusterBlob.condenseCores(cum_contigs_used_good, minSize)
        #t4 = time.time()
        #print "\tTHIS: [",self.secondsToStr(t4-t3),"]\tTOTAL: [",self.secondsToStr(t4-t0),"]"
        
        # Now save all the stuff to disk!
        print "Saving cores"
        (bin_update, core_update) = self.clusterBlob.getCoreBinUpdates()
        self.dataBlob.saveBins(bin_update)
        self.dataBlob.saveCores(core_update)
        self.dataBlob.setClustered()
        self.dataBlob.setNumBins(np.size(self.dataBlob.bins))
        t4 = time.time()
        print "\tTHIS: [",self.secondsToStr(t4-t3),"]\tTOTAL: [",self.secondsToStr(t4-t0),"]"
        
    def expandBins(self):
        """Load cores and expand bins"""
        # check that the user is OK with nuking stuff...
        if(not self.promptForOverwrite()):
            return False
        
        # get some data
        t0 = time.time()
        print "Load data"
        self.dataBlob.loadData(condition="length >= 10000")#(length >= 4000 ) & (length <= 4300)")
        t1 = time.time()
        print "\tTHIS: [",self.secondsToStr(t1-t0),"]\tTOTAL: [",self.secondsToStr(t1-t0),"]"

        # Now we use SOMs to classify the remaininfg contigs
        print "Start SOM classification"
        t2 = time.time()
        print "\tTHIS: [",self.secondsToStr(t2-t1),"]\tTOTAL: [",self.secondsToStr(t2-t0),"]"

    def promptForOverwrite(self):
        """Check that the user is ok with possibly overwriting the DB"""
        if(not self.forceWriting):
            if(self.dataBlob.isClustered()):
                option = raw_input(" ****WARNING**** Database: '"+self.dataBlob.dbFileName+"' has already been clustered.\n" \
                                   " If you continue you *MAY* overwrite existing bins!\n" \
                                   " Overwrite? (y,n) : ")
                print "****************************************************************"
                if(option.upper() != "Y"):
                    print "Operation cancelled"
                    return False
                else:
                    print "Overwriting database",self.dataBlob.dbFileName
        return True
    
    def secondsToStr(self, t):
        rediv = lambda ll,b : list(divmod(ll[0],b)) + ll[1:]
        return "%d:%02d:%02d.%03d" % tuple(reduce(rediv,[[t*1000,],1000,60,60]))
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################
class ClusterBlob:
    """Main class for performing contig clustering
    
    All the bits and bobs you'll need to cluster and bin out 
    pre-transformed primary data
    """    
    def __init__(self, dataBlob, debugPlots=False):
        # See DataTransformer for details about these variables
        self.dataBlob = dataBlob
        self.scaleFactor = self.dataBlob.scaleFactor
        
        # get enough memory for three heat maps
        self.imageMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        self.blurredMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        #self.maxMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        
        # we need a way to reference from the imageMaps back onto the transformed data
        self.mappedIndicies = {}
        self.binnedIndicies = {}

        # store our bins
        self.numBins = 0
        self.bins = {}                  # bins we kinda trust
        self.badBins = {}               # bins we don't kinda trust
        
        # housekeeping / parameters
        self.roundNumber = 0            # how many times have we done this?
        self.span = 30                  # amount we can travel about when determining "hot spots"

        # When blurring the raw image maps I chose a radius to suit my data, you can vary this as you like
        self.blurRadius = 12
        
        self.debugPlots = debugPlots
        self.imageCounter = 1           # when we print many images

    def initialiseCores(self, minVol):
        """Process contigs and form CORE bins"""
        small_bin_cutoff = 10           # anything with less than this many contigs is not a real bin
        num_below_cutoff = 0            # how many consecutive attempts have produced small bins
        breakout_point = 50             # how many will we allow before we stop this loop
        tmp_num_bins = 10000            # gotta keep count, high numbered bins ar baaad!
        
        cum_contigs_used_good = 0
        cum_contigs_used_bad = 0
        # First we need to find the centers of each blob.
        # We can make a heat map and look for hot spots
        self.populateImageMaps()
        print "\t",
        nl_counter = 0
        while(num_below_cutoff < breakout_point):
            #if(self.numBins > 1):
            #    break
            sys.stdout.flush()
            # apply a gaussian blur to each image map to make hot spots
            # stand out more from the background 
            self.blurMaps()
    
            # now search for the "hottest" spots on the blurred map
            # this is a new bin centroid
            (center_indicies, max_blur_value) = self.findNewClusterCenter()
            if(np.size(center_indicies) == 0):
                break
            else:
                # make sure we got something
                self.roundNumber += 1
                
                # time to make a bin
                tmp_num_bins += 1
                bin = Bin(center_indicies, self.dataBlob.kmerSigs, tmp_num_bins)
                
                # work out the distribution in points in this bin
                bin.makeBinDist(self.dataBlob.transformedData, self.dataBlob.kmerSigs)     
                
                # Plot?
                if(self.debugPlots):          
                    bin.plotBin(self.dataBlob.transformedData, self.dataBlob.contigColours, fileName="Image_"+str(self.imageCounter), tag="Initial")
                    self.imageCounter += 1
        
                # make the bin more gooder
                is_good_bin = True
                bin_size = bin.recruit(self.dataBlob.transformedData, self.dataBlob.kmerSigs, self.mappedIndicies, self.binnedIndicies)
                if(bin.calcTotalSize(self.dataBlob.contigLengths) < minVol):    # less than the good volume
                    if(bin_size < small_bin_cutoff):
                        is_good_bin = False
                        cum_contigs_used_bad += bin_size
                        self.badBins[tmp_num_bins] = bin
                        num_below_cutoff += 1
                        print "-",
                    # else bin is large enough!
                if(is_good_bin):
                    # make this bin legit!
                    cum_contigs_used_good += bin_size
                    self.numBins += 1
                    bin.id = self.numBins 
                    self.bins[self.numBins] = bin
                    # Plot?
                    if(self.debugPlots):          
                        bin.plotBin(self.dataBlob.transformedData, self.dataBlob.contigColours, fileName="Image_"+str(self.imageCounter), tag="CORE")
                        self.imageCounter += 1
                    num_below_cutoff = 0
                    print "+",

                # make the printing prettier
                nl_counter += 1
                if(nl_counter > 9):
                    nl_counter = 0
                    print "\n\t",
                    
                if(self.debugPlots):
                    self.plotHeat("3X3_"+str(self.roundNumber)+".png", max=max_blur_value)

                # append this bins list of mapped indicies to the main list
                self.updatePostBin(bin)
        print ""
        perc = "%.2f" % round((float(cum_contigs_used_good)/float(self.dataBlob.numContigs))*100,2)
        print "\t",cum_contigs_used_good,"contigs are distributed across",self.numBins,"cores (",perc,"% )"
        print "\t",cum_contigs_used_bad,"contigs are distributed across",len(self.badBins),"pseudo cores"
        
        return cum_contigs_used_good

    def getCoreBinUpdates(self):
        """Merge the bin information with the raw DB indexes so we can save to disk"""
        core_update = dict(zip(self.dataBlob.indicies, [False]*np.size(self.dataBlob.indicies)))
        bin_update = dict(zip(self.dataBlob.indicies, [0]*np.size(self.dataBlob.indicies)))

        # we need a mapping from cid (or local index) to binID
        c2b = dict(zip(range(0,np.size(self.dataBlob.indicies)), [0]*np.size(self.dataBlob.indicies)))
        for bid in self.bins:
            for index in self.bins[bid].indicies:
                c2b[index] = bid
        
        # at this stage, all bins are cores
        for index in range(0, self.dataBlob.numContigs):
            if index in self.binnedIndicies:
                bin_update[self.dataBlob.indicies[index]] = c2b[index]
                core_update[self.dataBlob.indicies[index]] = True

        return (bin_update, core_update)
            
    def populateImageMaps(self):
        """Load the transformed data into the main image maps"""
        # reset these guys... JIC
        self.imageMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        self.mappedIndicies = {}
        
        # add to the grid wherever we find a contig
        index = -1
        for point in np.around(self.dataBlob.transformedData):
            index += 1

            # can only bin things once!
            if index not in self.binnedIndicies:
                # readability
                px = point[0]
                py = point[1]
                pz = point[2]
                
                # add to the index dict so we can relate the 
                # map back to individual points later
                if (px,py,pz) in self.mappedIndicies:
                    self.mappedIndicies[(px,py,pz)].append(index)
                else:
                    self.mappedIndicies[(px,py,pz)] = [index]
                
                # now increment in the grid
                # for each point we encounter we incrmement
                # it's position + the positions to each side
                # and touching each corner
                self.incrementAboutPoint(0, px, py)
                self.incrementAboutPoint(1, self.scaleFactor - pz - 1, py)
                self.incrementAboutPoint(2, self.scaleFactor - pz - 1, self.scaleFactor - px - 1)

    def updatePostBin(self, bin):
        """Update data structures after assigning contigs to a new bin"""
        for index in bin.indicies:
            self.binnedIndicies[index] = True
            
            # now update the image map, decrement
            point = np.around(self.dataBlob.transformedData[index])
            # readability
            px = point[0]
            py = point[1]
            pz = point[2]
            self.decrementAboutPoint(0, px, py)
            self.decrementAboutPoint(1, self.scaleFactor - pz - 1, py)
            self.decrementAboutPoint(2, self.scaleFactor - pz - 1, self.scaleFactor - px - 1)

    def incrementAboutPoint(self, index, px, py, valP=1, valS=0.6, valC=0.2 ):
        """Increment value at a point in the 2D image maps
        
        Increment point by valP, increment neighbouring points at the
        sides and corners of the target point by valS and valC
        """
        if px > 0:
            if py > 0:
                self.imageMaps[index,px-1,py-1] += valC      # Top left corner
            self.imageMaps[index,px-1,py] += valS            # Top
            if py < self.scaleFactor-1:             
                self.imageMaps[index,px-1,py+1] += valC      # Top right corner

        if py > 0:
            self.imageMaps[index,px,py-1] += valS            # Left side
        self.imageMaps[index,px,py] += valP                  # Point
        if py < self.scaleFactor-1:             
            self.imageMaps[index,px,py+1] += valS            # Right side

        if px < self.scaleFactor-1:
            if py > 0:
                self.imageMaps[index,px+1,py-1] += valC      # Bottom left corner
            self.imageMaps[index,px+1,py] += valS            # Bottom
            if py < self.scaleFactor-1:             
                self.imageMaps[index,px+1,py+1] += valC      # Bottom right corner

    def decrementAboutPoint(self, index, px, py, valP=1, valS=0.6, valC=0.2 ):
        """Decrement value at a point in the 2D image maps"""
        if px > 0:
            if py > 0:
                self.imageMaps[index,px-1,py-1] -= valC      # Top left corner
                if self.imageMaps[index,px-1,py-1] < np.finfo(float).eps:
                    self.imageMaps[index,px-1,py-1] = 0
                
            self.imageMaps[index,px-1,py] -= valS            # Top
            if self.imageMaps[index,px-1,py] < np.finfo(float).eps:
                self.imageMaps[index,px-1,py] = 0
            if py < self.scaleFactor-1:             
                self.imageMaps[index,px-1,py+1] -= valC      # Top right corner
                if self.imageMaps[index,px-1,py+1] < np.finfo(float).eps:
                    self.imageMaps[index,px-1,py+1] = 0

        if py > 0:
            self.imageMaps[index,px,py-1] -= valS            # Left side
            if self.imageMaps[index,px,py-1] < np.finfo(float).eps:
                self.imageMaps[index,px,py-1] = 0
            
        self.imageMaps[index,px,py] -= valP                  # Point
        if self.imageMaps[index,px,py] < np.finfo(float).eps:
            self.imageMaps[index,px,py] = 0
        if py < self.scaleFactor-1:             
            self.imageMaps[index,px,py+1] -= valS            # Right side
            if self.imageMaps[index,px,py+1] < np.finfo(float).eps:
                self.imageMaps[index,px,py+1] = 0

        if px < self.scaleFactor-1:
            if py > 0:
                self.imageMaps[index,px+1,py-1] -= valC      # Bottom left corner
                if self.imageMaps[index,px+1,py-1] < np.finfo(float).eps:
                    self.imageMaps[index,px+1,py-1] = 0
            self.imageMaps[index,px+1,py] -= valS            # Bottom
            if self.imageMaps[index,px+1,py] < np.finfo(float).eps:
                self.imageMaps[index,px+1,py] = 0
            if py < self.scaleFactor-1:             
                self.imageMaps[index,px+1,py+1] -= valC      # Bottom right corner
                if self.imageMaps[index,px+1,py+1] < np.finfo(float).eps:
                    self.imageMaps[index,px+1,py+1] = 0

    def incrementAboutPoint3D(self, workingBlock, px, py, pz, vals=(6.4,4.9,2.5,1.6)):
        """Increment a point found in a 3D column
        
        used when finding the centroid of a hot area
        update the 26 points which surround the centre point
        z spans the height of the entire column, x and y have been offset to
        match the column subspace
        """
        
        # top slice
        if pz < self.scaleFactor-1:
            self.subIncrement3D(workingBlock, px, py, pz+1, vals, 1)
        
        # center slice
        self.subIncrement3D(workingBlock, px, py, pz, vals, 0)
        
        # bottom slice
        if pz > 0:
            self.subIncrement3D(workingBlock, px, py, pz-1, vals, 1)
        
    def subIncrement3D(self, workingBlock, px, py, pz, vals, offset):
        """AUX: Called from incrementAboutPoint3D does but one slice"""       
        # get the size of the working block
        shape = np.shape(workingBlock)
        if px > 0:
            if py > 0:
                workingBlock[px-1,py-1,pz] += vals[offset + 2]      # Top left corner
            workingBlock[px-1,py,pz] += vals[offset + 1]            # Top
            if py < shape[1]-1:             
                workingBlock[px-1,py+1,pz] += vals[offset + 2]      # Top right corner

        if py > 0:
            workingBlock[px,py-1,pz] += vals[offset + 1]            # Left side
        workingBlock[px,py,pz] += vals[offset]                      # Point
        if py < shape[1]-1:             
            workingBlock[px,py+1,pz] += vals[offset + 1]            # Right side

        if px < shape[0]-1:
            if py > 0:
                workingBlock[px+1,py-1,pz] += vals[offset + 2]      # Bottom left corner
            workingBlock[px+1,py,pz] += vals[offset + 1]            # Bottom
            if py < shape[1]-1:             
                workingBlock[px+1,py+1,pz] += vals[offset + 2]      # Bottom right corner
    
    def blurMaps(self):
        """Blur the 2D image maps"""
        self.blurredMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        #self.maxMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        
        for i in range (0,3): # top, front and side
            self.blurredMaps[i,:,:] = ndi.gaussian_filter(self.imageMaps[i,:,:]**0.5, (self.blurRadius,self.blurRadius)) 

        # there's still a lot of background signal to remove
        # we wish to remove 90% of the data, this will leave just the really hot spots
        # Make a histogram of the data (use the top face)
        [vals,points] = np.histogram(np.reshape(self.blurredMaps[0,:,:], (self.scaleFactor, self.scaleFactor,1)), 50)
        total = np.sum(vals)*0.80
        lop_index = 1       # where we lop off the low values
        for val in vals:
            total -= val
            if total <= 0:
                break
            lop_index += 1
        lop_val = points[lop_index]

        # remove these low values and down normalise so that the largest value is equal to exactly 1
        for i in range (0,3): # top, front and side
            self.blurredMaps[i,:,:] = np.where(self.blurredMaps[i,:,:] >= lop_val, self.blurredMaps[i,:,:], 0)/lop_val

    def makeCoordRanges(self, pos, span):
        """Make search ranges which won't go out of bounds"""
        lower = pos-span
        upper = pos+span+1
        if(lower < 0):
            lower = 0
        if(upper >= self.scaleFactor):
            upper = self.scaleFactor - 1
        return (lower, upper)

    def findNewClusterCenter(self):
        """Find a putative cluster"""
        # we work from the top view as this has the base clustering
        max_index = np.argmax(self.blurredMaps[0])
        max_value = self.blurredMaps[0].ravel()[max_index]

        max_x = int(max_index/self.scaleFactor)
        max_y = max_index - self.scaleFactor*max_x
        max_z = -1

        this_span = int(1.5 * self.span)
        span_len = 2*this_span+1
        
        # work out the region this max value lives in
        x_density = np.zeros(span_len)
        x_offset = max_x - this_span
        
        y_density = np.zeros(span_len)
        y_offset = max_y - this_span
        
        if(self.debugPlots):
            self.plotRegion(max_x,max_y,max_z, fileName="Image_"+str(self.imageCounter), tag="column", column=True)
            self.imageCounter += 1

        # make a 3d grid to hold the values
        working_block = np.zeros((span_len, span_len, self.scaleFactor))
        
        # go through the entire column
        (x_lower, x_upper) = self.makeCoordRanges(max_x, this_span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, this_span)
        for z in range(0, self.scaleFactor):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    # check that the point is real and that it has not yet been binned
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if index not in self.binnedIndicies:
                                # this is an unassigned point. 
                                self.incrementAboutPoint3D(working_block, x-x_lower, y-y_lower, z)

        # blur and find the highest value
        bwb = ndi.gaussian_filter(working_block, self.blurRadius)
        
        densest_index = np.unravel_index(np.argmax(bwb), (np.shape(bwb)))
        max_x = densest_index[0] + x_lower
        max_y = densest_index[1] + y_lower
        max_z = densest_index[2]
       
        if(self.debugPlots):
            self.plotRegion(max_x,max_y,max_z, fileName="Image_"+str(self.imageCounter), tag="first approx")
            self.imageCounter += 1

        # now get the basic color of this dense point
        (x_lower, x_upper) = self.makeCoordRanges(max_x, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, self.span)
        (z_lower, z_upper) = self.makeCoordRanges(max_z, self.span)
        center_values = np.array([])
        cv_colours = np.array([])
        c_inc = 0
        for z in range(z_lower, z_upper):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if index not in self.binnedIndicies:
                                center_values = np.append(center_values, self.dataBlob.kmerSigs[index])
                                cv_colours = np.append(cv_colours, self.dataBlob.contigColours[index])
                                c_inc += 1

        # make sure we have something to go on here
        if(np.size(center_values) == 0):
            return (np.array([]), -1)

        # reshape these guys!
        center_values = np.reshape(center_values, (c_inc, np.size(self.dataBlob.kmerSigs[0])))
        cv_colours = np.reshape(cv_colours, (c_inc, 3))
        
        # transform them into one dimensional points
        oneD_center_values = np.zeros(c_inc)
        working_centres = np.array(center_values, copy=True) 
        PCA.Center(working_centres,verbose=0)
        p = PCA.PCA(working_centres)
        components = p.pc()
        index = 0
        for point in components:
            oneD_center_values[index] = components[index,0]
            index += 1

        if(False):
            plt.figure(1)
            plt.subplot(111)
            cm = mpl.colors.LinearSegmentedColormap('my_colormap', cv_colours, 1024)
            plt.scatter(components[:,0], components[:,1], edgecolors=cv_colours, c=cv_colours, cmap=cm, marker='.')
            plt.show()

        # find the sig which lies at the center
        cf = CenterFinder()
        oneD_center_values -= np.min(oneD_center_values)
        oneD_center_values /= np.max(oneD_center_values)
        centroid_sig = center_values[cf.findArrayCenter(oneD_center_values)]
        
        # now we need to work out how close to this sig we need to be...
        dists = np.array([])
        for sig in center_values:
            dists = np.append(dists, np.linalg.norm(sig - centroid_sig))
            
        # get a first approximation on the upper cutoff
        tol = 2
        dists = np.sort(dists)
        sub_dists = dists[0:int(np.size(dists)/4):1]
        upper_dist = np.mean(sub_dists) + tol * np.std(sub_dists)

        # now refine (expand) based on the first approximation
        sub_dists = np.array([x for x in dists if x < upper_dist])
        upper_dist = np.mean(sub_dists) + tol * np.std(sub_dists)

        # now scoot out around this point and soak up similar points
        # get all the real indicies of these points so we can use them in
        # the primary data map
        center_indicies = np.array([])
        for z in range(z_lower, z_upper):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    # check that the point is real and that it has not yet been binned
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if(index not in center_indicies) and (index not in self.binnedIndicies):
                                # make sure the kmer sig is close enough
                                dist = np.linalg.norm(self.dataBlob.kmerSigs[index] - centroid_sig)
                                if(dist < upper_dist):
                                    center_indicies = np.append(center_indicies, index)
        if(np.size(center_indicies) > 0):
            return (center_indicies, max_value)
        
        return (np.array([]), -1)

    def Ablur(self, blur, density, incAtPoint, index, offset, size):
        """AUX: Used when finding the densest point in a small block"""
        point = index + offset;
        if(point >= 0 and point < size):
            blur[point] += incAtPoint[abs(offset)] * density[index]
    
#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 
    def plotRegion(self, px, py, pz, fileName="", tag="", column=False):
        """Plot the region surrounding a point """
        disp_vals = np.array([])
        disp_cols = np.array([])
        num_points = 0
        # plot all points within span
        (z_lower, z_upper) = self.makeCoordRanges(pz, self.span)
        if(column):
            z_lower = 0
            z_upper = self.scaleFactor - 1

        (x_lower, x_upper) = self.makeCoordRanges(px, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(py, self.span)
        for z in range(z_lower, z_upper):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if index not in self.binnedIndicies:
                                num_points += 1
                                disp_vals = np.append(disp_vals, self.dataBlob.transformedData[index])
                                disp_cols = np.append(disp_cols, self.dataBlob.contigColours[index])
        
        # make a black mark at the max values
        small_span = self.span/2
        (x_lower, x_upper) = self.makeCoordRanges(px, small_span)
        (y_lower, y_upper) = self.makeCoordRanges(py, small_span)
        (z_lower, z_upper) = self.makeCoordRanges(pz, small_span)
        for z in range(z_lower, z_upper):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if index not in self.binnedIndicies:
                                num_points += 1
                                disp_vals = np.append(disp_vals, self.dataBlob.transformedData[index])
                                disp_cols = np.append(disp_cols, colorsys.hsv_to_rgb(0,0,0))
        # reshape
        disp_vals = np.reshape(disp_vals, (num_points, 3))
        disp_cols = np.reshape(disp_cols, (num_points, 3))
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        cm = mpl.colors.LinearSegmentedColormap('my_colormap', disp_cols, 1024)
        result = ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors=disp_cols, c=disp_cols, cmap=cm, marker='.')
        title = str.join(" ", ["Focus at: (",str(px), str(py), str(self.scaleFactor - pz - 1),")\n",tag])
        plt.title(title)
      
        if(fileName != ""):
            fig.set_size_inches(6,6)
            plt.savefig(fileName,dpi=300)
        elif(show):
            plt.show()
        
        plt.close(fig)
        del fig
    
    def plotHeat(self, fileName = "", max=-1):
        """Print the main heat maps
        
        Useful for debugging
        """
        fig = plt.figure()
        images = []

        ax = fig.add_subplot(231)
        images.append(ax.imshow(self.blurredMaps[0,:,:]**0.5))
        if(max > 0):
            title = str.join(" ", ["Max value:",str(max)])
            plt.title(title)
        ax = fig.add_subplot(232)
        images.append(ax.imshow(self.blurredMaps[1,:,:]**0.5))
        ax = fig.add_subplot(233)
        images.append(ax.imshow(self.blurredMaps[2,:,:]**0.5))

        ax = fig.add_subplot(234)
        images.append(ax.imshow(self.imageMaps[0,:,:]**0.5))
        ax = fig.add_subplot(235)
        images.append(ax.imshow(self.imageMaps[1,:,:]**0.5))
        ax = fig.add_subplot(236)
        images.append(ax.imshow(self.imageMaps[2,:,:]**0.5))
        
        if(fileName != ""):
            fig.set_size_inches(18,18)
            plt.savefig(fileName,dpi=300)
        elif(show):
            plt.show()

        plt.close(fig)
        del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CoreValidator:
    """Plot images of transformed data and bin cores for user validation"""
    def __init__(self, dbFileName):
        self.dataBlob = DataBlob(dbFileName)        # based on user specified length
        self.dataBlobC = DataBlob(dbFileName)       # cores!

        # core stats
        self.cores = {}
        self.ccData = np.array([])
        self.ccColour = np.array([])
        
        # munged data 
        self.numCores = 0
        self.coreMembers = {}

    def validate(self, coreCut):
        """Main wrapper for core validation"""
        self.loadData(coreCut)
        self.findCoreCentres()
        self.measureBinKVariance()
        #self.renderValImg()
        
    def loadData(self, coreCut):
        """Load data from the DB and transform"""
        self.dataBlob.loadData(condition="length >= "+str(coreCut))
        self.dataBlobC.loadData(loadBins=True, loadKSigs=True, condition="core >= True")
        print "\t( length >=",coreCut,") -> ",np.size(self.dataBlob.kmerSigs)
        print "\tCores -> ",np.size(self.dataBlobC.kmerSigs)
        self.numCores = self.dataBlobC.getNumBins()
        self.dataBlob.transformData()
        self.dataBlobC.transformData()

    def findCoreCentres(self):
        """Find the point representing the centre of each core"""
        self.ccData = np.zeros((self.numCores,3))
        for index in range(0,self.numCores):
            self.coreMembers[index+1] = []
        
        # fill them up
        for index in range(0, np.size(self.dataBlobC.indicies)):
            self.coreMembers[self.dataBlobC.bins[index]].append(index)

        # remake the cores and populate the centres
        S = 1       # SAT and VAL remain fixed at 1. Reduce to make
        V = 1       # Pastels if that's your preference...
        for index in range(0,self.numCores):
            # add 1 to the index as 0 is the null core!
            self.cores[index+1] = Bin(np.array(self.coreMembers[index+1]), self.dataBlobC.kmerSigs, index+1)
            self.cores[index+1].makeBinDist(self.dataBlobC.transformedData, self.dataBlobC.kmerSigs)
            for i in range (0,3):
                self.ccData[index][i] = self.cores[index+1].mean[i]
            self.ccColour = np.append(self.ccColour, [colorsys.hsv_to_rgb(self.cores[index+1].mean[3], S, V)])
            
        self.ccColour = np.reshape(self.ccColour, (self.numCores, 3))            

    def measureBinKVariance(self):
        """Measure within and between bin variance of kmer sigs"""
        means = np.array([])
        stdevs = np.array([])
        bids = np.array([])
        
        # work out the mean and stdev for the kmer sigs for each bin
        for bid in self.cores:
            bkworking = np.array([])
            for index in self.cores[bid].indicies:
                bkworking = np.append(bkworking, self.dataBlobC.kmerSigs[index])
            bkworking = np.reshape(bkworking, (self.cores[bid].binSize, np.size(self.dataBlobC.kmerSigs[0])))
            bids = np.append(bids, [bid])
            means = np.append(means, np.mean(bkworking, axis=0))
            stdevs = np.append(stdevs, np.std(bkworking, axis=0))
            
        means = np.reshape(means, (self.numCores, np.size(self.dataBlobC.kmerSigs[0])))
        stdevs = np.reshape(stdevs, (self.numCores, np.size(self.dataBlobC.kmerSigs[0])))
        
        # now work out the between and within core variances
        between = np.std(means, axis=0)
        within = np.median(stdevs, axis=0)

        B = np.arange(0, np.size(self.dataBlobC.kmerSigs[0]), 1)
        names = self.dataBlobC.getMerColNames().split(',')
        
        # we'd like to find the indicies of the worst 10% for each type so we can ignore them
        # specifically, we'd like to remove the least variable between core kms and the 
        # most variable within core kms.
        sort_between_indicies = np.argsort(between)
        sort_within_indicies = np.argsort(within)[::-1]
        number_to_trim = int(0.1* float(np.size(self.dataBlobC.kmerSigs[0])))
        print "BETWEEN"
        for i in range(0,number_to_trim):
            print names[sort_between_indicies[i]]
        print "WITHIN" 
        for i in range(0,number_to_trim):
            print names[sort_within_indicies[i]] 
        
        plt.figure(1)
        plt.subplot(211)
        plt.plot(B, between, 'r--', B, within, 'b--')
        plt.xticks(B, names, rotation=90)
        plt.grid()
        plt.subplot(212)
        ratio = between/within
        plt.plot(B, ratio, 'r--')
        plt.xticks(B, names, rotation=90)
        plt.grid()
        plt.show()
        
        
#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def renderValImg(self):
        """Render the image for validating cores"""
        fig = plt.figure()
        ax1 = fig.add_subplot(121, projection='3d')
        ax1.scatter(self.dataBlob.transformedData[:,0], self.dataBlob.transformedData[:,1], self.dataBlob.transformedData[:,2], edgecolors=self.dataBlob.contigColours, c=self.dataBlob.contigColours, marker='.')

        ax2 = fig.add_subplot(122, projection='3d')
        ax2.scatter(self.ccData[:,0], self.ccData[:,1], self.ccData[:,2], edgecolors=self.ccColour, c=self.ccColour, marker='.')

        try:
            plt.show()
            plt.close(fig)
        except:
            print "Error showing image", sys.exc_info()[0]
            raise

        del fig


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Bin:
    """Class for managing collections of contigs
    
    To (perhaps) simplify things think of a "bin" as an index into the
    column names array. The ClusterBlob has a list of bins which it can
    update etc...
    """
    def __init__(self, indicies, kmerSigs, id, covtol=3, mertol=2):
        self.id = id
        self.indicies = indicies           # all the indicies belonging to this bin
        self.binSize = self.indicies.shape[0]
        self.totalBP = 0
        
        # we need some objects to manage the distribution of contig proerties
        self.covMeans = np.zeros((3))
        self.covStdevs = np.zeros((3))
        self.covTolerance = covtol
        self.covLowerLimits = np.zeros((3)) # lower and upper limits based on tolerance
        self.covUpperLimits = np.zeros((3))
        
        self.merMeans = np.array([])
        self.merStdevs = np.array([])
        self.merCentroid = np.array([])
        self.kDistMean = 0
        self.kDiststdev = 0
        self.kDistTolerance = mertol
        self.kDistUpperLimit = 0

#------------------------------------------------------------------------------
# Tools used for condensing 
    
    def __cmp__(self, alien):
        """Sort bins based on aux values"""
        if self.mean[3] < alien.mean[3]:
            return -1
        elif self.mean[3] == alien.mean[3]:
            return 0
        else:
            return 1

    def removeFromBin(self, transformedData, kmerSigs, contigLengths, numTigs):
        """Remove contigs from a bin; returns the global indicies of those removed"""
        # get the first numTigs
        if(numTigs > self.binSize):
            numTigs = self.binSize
            
        ret_indicies = np.array([])
        while(numTigs > 0):
            index = random.randint(0, len(self.indicies)-1)
            ret_indicies = np.append(ret_indicies, self.indicies[index])
            self.indicies = np.delete(self.indicies, index)
            numTigs -= 1
            self.binSize -= 1
            
        # fix the stats on our bin
        self.makeBinDist(transformedData, kmerSigs)
        self.calcTotalSize(contigLengths)

        return ret_indicies
        
    def isSimilar(self, compBin, stdevs=1):
        """Check whether two bins are similar"""
        this_lowers = self.mean - stdevs * self.stdev
        this_uppers = self.mean + stdevs * self.stdev
        that_lowers = compBin.mean - stdevs * compBin.stdev
        that_uppers = compBin.mean + stdevs * compBin.stdev
        #print "\n\n",this_uppers,"\n",compBin.mean,"\n",this_lowers,"\n---\n",that_uppers,"\n",self.mean,"\n",that_lowers
        # reciprocial test!
        for index in range(0,4):
            if(self.mean[index] < that_lowers[index] or self.mean[index] > that_uppers[index]):
                return False
            if(compBin.mean[index] < this_lowers[index] or compBin.mean[index] > this_uppers[index]):
                return False
        # got here? Must be similar!
        return True
    
    def consume(self, transformedData, kmerSigs, contigLengths, deadBin):
        """Combine the contigs of another bin with this one"""
        # consume all the other bins indicies
        self.indicies = np.concatenate([self.indicies, deadBin.indicies])
        self.binSize += deadBin.binSize
        
        # fix the stats on our bin
        self.makeBinDist(transformedData, kmerSigs)
        self.calcTotalSize(contigLengths)
        
#------------------------------------------------------------------------------
# Stats and properties 

    def clearBinDist(self):
        """Clear any set distribution statistics"""
        self.covMeans = np.zeros((3))
        self.covStdevs = np.zeros((3))
        self.covLowerLimits = np.zeros((3))
        self.covUpperLimits = np.zeros((3))
        
        self.merMeans = np.array([])
        self.merStdevs = np.array([])
        self.merCentroid = np.array([])
        self.kDistMean = 0
        self.kDiststdev = 0
        self.kDistUpperLimit = 0
        
        
    def makeBinDist(self, transformedData, kmerSigs):
        """Determine the distribution of the points in this bin
        
        The distribution is largely normal, except at the boundaries.
        """
        self.clearBinDist()
        self.merCentroid = np.zeros((np.size(kmerSigs[0])))
        
        # Get some data!
        cov_working_array = np.zeros((self.binSize,3))
        mer_working_array = np.zeros((self.binSize,np.size(kmerSigs[0])))
        outer_index = 0
        for index in self.indicies:
            for i in range(0,3):
                cov_working_array[outer_index][i] = transformedData[index][i]
            mer_working_array[outer_index] = kmerSigs[index]
            outer_index += 1
        
        # calculate the coverage mean and stdev 
        self.covMeans = np.mean(cov_working_array,axis=0)
        self.covStdevs = np.std(cov_working_array,axis=0)

        # now do the kmerSigs
        # z-normalise each column in the working array
        self.merMeans = np.mean(mer_working_array, axis=0)
        tmpMerStdevs = np.std(mer_working_array, axis=0)
        # no zeros!
        self.merStdevs = np.array([x if x !=0 else 1.0 for x in tmpMerStdevs])
        for index in range(0,np.size(self.indicies)):
            mer_working_array[index] = (mer_working_array[index]-self.merMeans)/self.merStdevs
        
        # work out the distribution of distances from z-normed sigs to the centroid
        k_dists = np.array([])
        for sig in mer_working_array:
            k_dists = np.append(k_dists, np.linalg.norm(sig-self.merCentroid))
        self.kDistMean = np.mean(k_dists)
        self.kDiststdev = np.std(k_dists)
        # set the acceptance ranges
        self.makeLimits()
        
    def makeLimits(self, pt=-1, st=-1):
        """Set inclusion limits based on mean, variance and tolerance settings"""
        if(-1 == pt):
            pt=self.covTolerance
        if(-1 == st):
            st=self.kDistTolerance
        for i in range(0,3):
            self.covLowerLimits[i] = int(self.covMeans[i] - pt * self.covStdevs[i])
            self.covUpperLimits[i] = int(self.covMeans[i] + pt * self.covStdevs[i]) + 1  # so range will look neater!
        self.kDistUpperLimit = self.kDistMean + st * self.kDiststdev
        
    def getKDist(self, sig):
        """Get the distance of this sig from the centroid"""
        # z-norm and then distance!
        return np.linalg.norm((sig-self.merMeans)/self.merStdevs - self.merCentroid)
    
#------------------------------------------------------------------------------
# Grow the bin 
    
    def recruit(self, transformedData, kmerSigs, mappedIndicies, binnedIndicies):
        """Iteratively grow the bin"""
        self.makeBinDist(transformedData, kmerSigs)

        # save these
        pt = self.covTolerance
        st = self.kDistTolerance

        self.binSize = self.indicies.shape[0]
        num_recruited = self.recruitRound(transformedData, kmerSigs, mappedIndicies, binnedIndicies) 
        while(num_recruited > 0):
            # reduce these to force some kind of convergence
            self.covTolerance *= 0.8
            self.kDistTolerance *= 0.8
            # fix these
            self.binSize = self.indicies.shape[0]
            self.makeBinDist(transformedData, kmerSigs)
            # go again
            num_recruited = self.recruitRound(transformedData, kmerSigs, mappedIndicies, binnedIndicies)
        
        self.covTolerance = pt
        self.kDistTolerance = st
        
        # finally, fix this guy
        return self.binSize
        
    def recruitRound(self, transformedData, kmerSigs, mappedIndicies, binnedIndicies):
        """Recruit more points in from outside the current blob boundaries"""
        num_recruited = 0
        for x in range(int(self.covLowerLimits[0]), int(self.covUpperLimits[0])):
            for y in range(int(self.covLowerLimits[1]), int(self.covUpperLimits[1])):
                for z in range(int(self.covLowerLimits[2]), int(self.covUpperLimits[2])):
                    if((x,y,z) in mappedIndicies):
                        for index in mappedIndicies[(x,y,z)]:
                            if (index not in binnedIndicies) and (index not in self.indicies):
                                k_dist = self.getKDist(kmerSigs[index])
                                if(k_dist <= self.kDistUpperLimit):
                                    self.indicies = np.append(self.indicies,index)
                                    num_recruited += 1
        return num_recruited

    def calcTotalSize(self, contigLengths):
        """Work out the total size of this bin in BP"""
        totalBP = 0
        for index in self.indicies:
            totalBP += contigLengths[index]
        self.totalBP = totalBP

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 
#
    def plotBin(self, transformedData, contigColours, fileName="", tag=""):
        """Plot a bin"""
        disp_vals = np.array([])
        disp_cols = np.array([])
        num_points = 0
        for index in self.indicies:
            num_points += 1
            disp_vals = np.append(disp_vals, transformedData[index])
            disp_cols = np.append(disp_cols, contigColours[index])

        # make a black mark at the max values
        self.makeLimits(pt=1, st=1)
        px = int(self.covMeans[0])
        py = int(self.covMeans[1])
        pz = int(self.covMeans[2])
        num_points += 1
        disp_vals = np.append(disp_vals, [px,py,pz])
        disp_cols = np.append(disp_cols, colorsys.hsv_to_rgb(0,0,0))
        
        # fix these
        self.makeLimits()
        
        # reshape
        disp_vals = np.reshape(disp_vals, (num_points, 3))
        disp_cols = np.reshape(disp_cols, (num_points, 3))

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors=disp_cols, c=disp_cols, marker='.')
        from locale import format, setlocale, LC_ALL # purdy commas
        setlocale(LC_ALL, "")
        title = str.join(" ", ["Bin:",str(self.id),"--",tag,"\n",
                               "Focus at: (",str(px), str(py), str(pz),")\n",
                               "Contains:",str(self.binSize),"contigs\n",
                               "Total:",format('%d', self.totalBP, True),"BP"
                               ])
        plt.title(title)
        
        if(fileName != ""):
            try:
                fig.set_size_inches(6,6)
                plt.savefig(fileName,dpi=300)
            except:
                print "Error saving image:", fileName, sys.exc_info()[0]
                raise
        elif(show):
            try:
                plt.show()
            except:
                print "Error showing image:", sys.exc_info()[0]
                raise
        plt.close(fig)
        del fig
    
    def printContents(self):
        """Dump the contents of the object"""
        print "--------------------------------------"
        print "Bin:", self.id
        print "Bin size:", self.binSize
        print "Total BP:", self.totalBP
        print "--------------------------------------"
    
    def dumpContigIDs(self, contigNames):
        """Print out the contigIDs"""
        from cStringIO import StringIO
        file_str = StringIO()
        for index in self.indicies:
            file_str.write(contigNames[index]+"\t")
        return file_str.getvalue()

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class CenterFinder:
    """When a plain old mean won't cut it

    Uses a bouncing ball algorithm. Imagine walking along a "path",
    (through the array) hitting a ball into the air each time you
    come across a value. Gravity is bringing the ball down. If we plot
    the height of the ball vs array index then the highest the ball
    reaches is the index in the center of the densest part of the array 
    
    NOTE: Assumes the data is normalised between 0 and 1!
    """
    def __init__(self): pass
    
    def findArrayCenter(self, vals):
        """Find the center of the numpy array vals, return the index of the center"""
        # parameters
        current_val_max = -1
        delta = 0
        bounce_amount = 0.1
        height = 0
        last_val= 0

        working = np.array([])
        final_index = -1
        # run through in one direction
        vals_sorted = np.sort(vals)
        for val in vals_sorted:
            # calculate delta
            delta = val - last_val
            # reduce the current value according to the delta value
            height = self.reduceViaDelta(height, bounce_amount, delta)
            # bounce the ball up
            height += bounce_amount
            
            # store the height
            working = np.append(working, height)
            final_index += 1

            # save the last val            
            last_val = val

        current_val_max = -1
        height = 0
        last_val = 0
        # run through in the reverse direction
        vals_sorted = vals_sorted[::-1]
        for val in vals_sorted:
            if last_val == 0:
                delta = 0
            else:
                delta = last_val - val
            height = self.reduceViaDelta(height, bounce_amount, delta)
            height += bounce_amount
            # add to the old heights
            working[final_index] += height
            final_index -= 1
            last_val = val
        
        max_index = np.argmax(working)
        vals_sorted = np.sort(vals_sorted)
        max_value = vals_sorted[max_index]
        
        # find the original index!
        index = 0
        for val in vals:
            if(val == max_value):
                return index
            index += 1
        return -1
    
    def reduceViaDelta(self, height, bounce_amount, delta):
        """Reduce the height of the 'ball'"""
        perc = (delta / bounce_amount)**0.5
        if(perc > 1):
            #print height, delta, 1, " H: ", 0
            return 0
        #print height, delta, (1-perc), " H: ", (height * (1-perc)) 
        return height * (1-perc)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

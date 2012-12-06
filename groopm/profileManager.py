#!/usr/bin/env python
###############################################################################
#                                                                             #
#    profileManager.py                                                          #
#                                                                             #
#    GroopM - High level data management                                      #
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
__version__ = "0.2.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Alpha"

###############################################################################

from sys import exc_info, exit, stdout
from operator import itemgetter

from colorsys import hsv_to_rgb as htr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

from numpy import abs as np_abs, amax as np_amax, amin as np_amin, append as np_append, arccos as np_arccos, argmin as np_argmin, argsort as np_argsort, array as np_array, ceil as np_ceil, concatenate as np_concatenate, delete as np_delete, log10 as np_log10, max as np_max, mean as np_mean, median as np_median, min as np_min, pi as np_pi, reshape as np_reshape, seterr as np_seterr, size as np_size, sort as np_sort, sqrt as np_sqrt, std as np_std, where as np_where, zeros as np_zeros, cos as np_cos, sin as np_sin
from numpy.linalg import norm as np_norm 
import scipy.ndimage as ndi
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree as kdt
from scipy.stats import f_oneway, distributions

# GroopM imports
from PCA import PCA, Center
from mstore import GMDataManager
from bin import Bin
import groopmExceptions as ge

np_seterr(all='raise')     

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ProfileManager:
    """Interacts with the groopm DataManager and local data fields
    
    Mostly a wrapper around a group of numpy arrays and a pytables quagmire
    """
    def __init__(self, dbFileName, force=False, scaleFactor=1000):
        # data
        self.dataManager = GMDataManager()  # most data is saved to hdf
        self.dbFileName = dbFileName        # db containing all the data we'd like to use
        self.condition = ""                 # condition will be supplied at loading time
        # --> NOTE: ALL of the arrays in this section are in sync
        # --> each one holds information for an individual contig 
        self.indices = np_array([])        # indices into the data structure based on condition
        self.covProfiles = np_array([])     # coverage based coordinates
        self.transformedCP = np_array([])   # the munged data points
        self.averageCoverages = np_array([]) # average coverage across all stoits
        self.kmerSigs = np_array([])        # raw kmer signatures
        self.kmerVals = np_array([])        # PCA'd kmer sigs

        self.contigNames = np_array([])
        self.contigLengths = np_array([])
        self.contigColours = np_array([])   # calculated from kmerVals
        
        self.binIds = np_array([])          # list of bin IDs
        # --> end section

        # meta                
        self.validBinIds = {}               # valid bin ids -> numMembers
        self.binnedRowIndicies = {}         # dictionary of those indices which belong to some bin
        self.restrictedRowIndicies = {}     # dictionary of those indices which can not be binned yet
        self.numContigs = 0                 # this depends on the condition given
        self.numStoits = 0                  # this depends on the data which was parsed

        # contig links
        self.links = {}
        
        # misc
        self.forceWriting = force           # overwrite existng values silently?
        self.scaleFactor = scaleFactor      # scale every thing in the transformed data to this dimension

    def loadData(self,
                 condition="",              # condition as set by another function
                 bids=[],                   # if this is set then only load those contigs with these bin ids
                 verbose=True,              # many to some output messages
                 silent=False,              # some to no output messages
                 loadCovProfiles=True,
                 loadKmerSigs=True,
                 makeColours=True,
                 loadContigNames=True,
                 loadContigLengths=True,
                 loadBins=False,
                 loadLinks=False):
        """Load pre-parsed data"""
        if(verbose):
            print "Loading data from:", self.dbFileName
        
        # check to see if we need to override the condition
        if(len(bids) != 0):
            condition = "((bid == "+str(bids[0])+")"
            for index in range (1,len(bids)):
                condition += " | (bid == "+str(bids[index])+")"
            condition += ")"
        if(silent):
            verbose=False
        try:
            self.numStoits = self.getNumStoits()
            self.condition = condition
            if(verbose):
                print "    Loading indices (", condition,")"
            self.indices = self.dataManager.getConditionalIndicies(self.dbFileName, condition=condition)
            self.numContigs = len(self.indices)
            
            if(not silent):
                print "    Working with: %d contigs" % self.numContigs

            if(loadCovProfiles):
                if(verbose):
                    print "    Loading coverage profiles"
                self.covProfiles = self.dataManager.getCoverageProfiles(self.dbFileName, indices=self.indices)

                # work out average coverages
                self.averageCoverages = np_array([sum(i)/self.numStoits for i in self.covProfiles])

            if(loadKmerSigs):
                if(verbose):
                    print "    Loading kmer sigs"
                self.kmerSigs = self.dataManager.getKmerSigs(self.dbFileName, indices=self.indices)

                if(makeColours):
                    if(verbose):
                        print "    Creating colour profiles"
                    self.makeColourProfile()
                    # use HSV to RGB to generate colours
                    S = 1       # SAT and VAL remain fixed at 1. Reduce to make
                    V = 1       # Pastels if that's your preference...
                    self.contigColours = np_array([htr(val, S, V) for val in self.kmerVals])

            if(loadContigNames):
                if(verbose):
                    print "    Loading contig names"
                self.contigNames = self.dataManager.getContigNames(self.dbFileName, indices=self.indices)
            
            if(loadContigLengths):
                if(verbose):
                    print "    Loading contig lengths"
                self.contigLengths = self.dataManager.getContigLengths(self.dbFileName, indices=self.indices)
                print "    Contigs contain %d BP" % ( sum(self.contigLengths) )
            
            if(loadBins):
                if(verbose):
                    print "    Loading bins"
                self.binIds = self.dataManager.getBins(self.dbFileName, indices=self.indices)
                if len(bids) != 0: # need to make sure we're not restricted in terms of bins
                    tmp_bids = self.getBinStats()
                    for bid in bids:
                        self.validBinIds[bid] = tmp_bids[bid]
                else:
                    self.validBinIds = self.getBinStats()

                # fix the binned indices
                self.binnedRowIndicies = {}
                for i in range(len(self.indices)):
                    if(self.binIds[i] != 0):
                        self.binnedRowIndicies[i] = True 
            else:
                # we need zeros as bin indicies then...
                self.binIds = np_zeros(len(self.indices))
                
            if(loadLinks):
                self.loadLinks()
            
        except:
            print "Error loading DB:", self.dbFileName, exc_info()[0]
            raise

    def reduceIndicies(self, deadRowIndicies):
        """purge indices from the data structures
        
        Be sure that deadRowIndicies are sorted ascending
        """
        # strip out the other values        
        self.indices = np_delete(self.indices, deadRowIndicies, axis=0)
        self.covProfiles = np_delete(self.covProfiles, deadRowIndicies, axis=0)
        self.transformedCP = np_delete(self.transformedCP, deadRowIndicies, axis=0)
        self.contigNames = np_delete(self.contigNames, deadRowIndicies, axis=0)
        self.contigLengths = np_delete(self.contigLengths, deadRowIndicies, axis=0)
        self.contigColours = np_delete(self.contigColours, deadRowIndicies, axis=0)
        self.kmerSigs = np_delete(self.kmerSigs, deadRowIndicies, axis=0)
        self.kmerVals = np_delete(self.kmerVals, deadRowIndicies, axis=0)
        self.binIds = np_delete(self.binIds, deadRowIndicies, axis=0)
        
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

    def getBinStats(self):
        """Go through all the "bins" array and make a list of unique bin ids vs number of contigs"""
        return self.dataManager.getBinStats(self.dbFileName)
    
    def setBinStats(self, binStats):
        """Store the valid bin Ids and number of members
                
        binStats is a dictionary which looks like:
        { tableRow : [bid , numMembers] }
        """
        self.dataManager.setBinStats(self.dbFileName, binStats)
        self.setNumBins(len(binStats.keys()))

    def setBinAssignments(self, assignments):
        """Save our bins into the DB"""
        self.dataManager.setBinAssignments(self.dbFileName, assignments)

    def loadLinks(self):
        """Extra wrapper 'cause I am dumb"""
        self.links = self.getLinks()
        
    def getLinks(self):
        """Get contig links"""
        # first we get the absolute links
        absolute_links = self.dataManager.restoreLinks(self.dbFileName, self.indices)
        # now convert this into plain old row_indices
        reverse_index_lookup = {} 
        for i in range(len(self.indices)):
            reverse_index_lookup[self.indices[i]] = i

        # now convert the absolute links to local ones
        relative_links = {}
        for cid in self.indices:
            local_cid = reverse_index_lookup[cid]
            relative_links[local_cid] = []
            try:
                for link in absolute_links[cid]:
                    relative_links[local_cid].append([reverse_index_lookup[link[0]], link[1], link[2], link[3]])
            except KeyError: # not everyone is linked
                pass

        return relative_links
                 
#------------------------------------------------------------------------------
# DATA TRANSFORMATIONS 

    def getAverageCoverage(self, rowIndex):
        """Return the average coverage for this contig across all stoits"""
        return sum(self.transformedCP[rowIndex])/self.numStoits

    def transformCP(self, silent=False, nolog=False, min=None, max=None):
        """Do the main ransformation on the coverage profile data"""
        shrinkFn = np_log10
        if(nolog):
            shrinkFn = lambda x:x
         
        s = (self.numContigs,3)
        self.transformedCP = np_zeros(s)

        if(not silent):
            print "    Dimensionality reduction"

        # get the median distance from the origin
        unit_vectors = [(np_cos(i*2*np_pi/self.numStoits),np_sin(i*2*np_pi/self.numStoits)) for i in range(self.numStoits)]
        for i in range(len(self.indices)):
            norm = np_norm(self.covProfiles[i])
            if(norm != 0):
                radial = shrinkFn(norm)
            else:
                radial = norm
            shifted_vector = np_array([0.0,0.0])
            flat_vector = (self.covProfiles[i] / sum(self.covProfiles[i]))
            
            for j in range(self.numStoits):
                shifted_vector[0] += unit_vectors[j][0] * flat_vector[j]
                shifted_vector[1] += unit_vectors[j][1] * flat_vector[j]

            # log scale it towards the centre
            scaling_vector = shifted_vector * self.scaleFactor
            sv_size = np_norm(scaling_vector)
            if(sv_size > 1):
                shifted_vector /= shrinkFn(sv_size)

            self.transformedCP[i,0] = shifted_vector[0]
            self.transformedCP[i,1] = shifted_vector[1]
            self.transformedCP[i,2] = radial

        if(not silent):
            print "    Reticulating splines"
            
        # finally scale the matrix to make it equal in all dimensions
        if(min is None):                
            min = np_amin(self.transformedCP, axis=0)
            max = np_amax(self.transformedCP, axis=0)
            max = max - min
            max = max / (self.scaleFactor-1)

        for i in range(0,3):
            self.transformedCP[:,i] = (self.transformedCP[:,i] -  min[i])/max[i]
            
        return(min,max)

    def makeColourProfile(self):
        """Make a colour profile based on ksig information"""
        working_data = np_array(self.kmerSigs, copy=True) 
        Center(working_data,verbose=0)
        p = PCA(working_data)
        components = p.pc()
        
        # now make the colour profile based on PC1
        self.kmerVals = np_array([float(i) for i in components[:,0]])
        
        # normalise to fit between 0 and 1
        self.kmerVals -= np_min(self.kmerVals)
        self.kmerVals /= np_max(self.kmerVals)
        if(False):
            plt.figure(1)
            plt.subplot(111)
            plt.plot(components[:,0], components[:,1], 'r.')
            plt.show()
        
    def rotateVectorAndScale(self, point, las, centerVector, delta_max=0.25):
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
        
        Set delta_max as the max percentage of the existing angle to be removed
        """
        theta = self.getAngBetween(point, centerVector)
        A = delta_max/((las)**2)
        B = delta_max/las
        delta = 2*B*theta - A *(theta**2) # the amount to shift
        V_p = point*(1-delta) + centerVector*delta
        return V_p/np_norm(V_p)
    
    def rad2deg(self, anglein):
        return 180*anglein/np_pi

    def getAngBetween(self, P1, P2):
        """Return the angle between two points (in radians)"""
        # find the existing angle between them theta
        c = np_dot(P1,P2)/np_norm(P1)/np_norm(P2) 
        # rounding errors hurt everyone...
        if(c > 1):
            c = 1
        elif(c < -1):
            c = -1
        return np_arccos(c) # in radians

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def plotUnbinned(self, coreCut):
        """Plot all contigs over a certain length which are unbinned"""
        self.loadData(condition="((length >= "+str(coreCut)+") & (bid == 0))")
        self.transformCP()
        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.scatter(self.transformedCP[:,0], self.transformedCP[:,1], self.transformedCP[:,2], edgecolors=self.contigColours, c=self.contigColours, marker='.')
        try:
            plt.show()
            plt.close(fig)
        except:
            print "Error showing image", exc_info()[0]
            raise
        del fig


    def plotTransViews(self, tag="fordens"):
        """Plot top, side and front views of the transformed data"""
        self.renderTransData(tag+"_top.png",azim = 0, elev = 90)
        self.renderTransData(tag+"_front.png",azim = 0, elev = 0)
        self.renderTransData(tag+"_side.png",azim = 90, elev = 0)

    def renderTransCPData(self, fileName="", show=True, elev=45, azim=45, all=False, showAxis=False, primaryWidth=12, primarySpace=3, dpi=300, format='png', fig=None):
        """Plot transformed data in 3D"""
        del_fig = False
        if(fig is None):
            fig = plt.figure()
            del_fig = True
        else:
            plt.clf()
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
            ax.scatter(self.transformedCP[:,0], self.transformedCP[:,1], self.transformedCP[:,2], edgecolors=self.contigColours, c=self.contigColours, marker='.')
            ax.azim = 0
            ax.elev = 0
            ax.set_xlim3d(0,self.scaleFactor)
            ax.set_ylim3d(0,self.scaleFactor)
            ax.set_zlim3d(0,self.scaleFactor)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
                for elt in axis.get_ticklines() + axis.get_ticklabels():
                    elt.set_visible(False)
            ax.w_xaxis._AXINFO = myAXINFO
            ax.w_yaxis._AXINFO = myAXINFO
            ax.w_zaxis._AXINFO = myAXINFO
            
            ax = fig.add_subplot(132, projection='3d')
            ax.scatter(self.transformedCP[:,0], self.transformedCP[:,1], self.transformedCP[:,2], edgecolors=self.contigColours, c=self.contigColours, marker='.')
            ax.azim = 90
            ax.elev = 0
            ax.set_xlim3d(0,self.scaleFactor)
            ax.set_ylim3d(0,self.scaleFactor)
            ax.set_zlim3d(0,self.scaleFactor)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
                for elt in axis.get_ticklines() + axis.get_ticklabels():
                    elt.set_visible(False)
            ax.w_xaxis._AXINFO = myAXINFO
            ax.w_yaxis._AXINFO = myAXINFO
            ax.w_zaxis._AXINFO = myAXINFO
            
            ax = fig.add_subplot(133, projection='3d')
            ax.scatter(self.transformedCP[:,0], self.transformedCP[:,1], self.transformedCP[:,2], edgecolors=self.contigColours, c=self.contigColours, marker='.')
            ax.azim = 0
            ax.elev = 90
            ax.set_xlim3d(0,self.scaleFactor)
            ax.set_ylim3d(0,self.scaleFactor)
            ax.set_zlim3d(0,self.scaleFactor)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            for axis in ax.w_xaxis, ax.w_yaxis, ax.w_zaxis:
                for elt in axis.get_ticklines() + axis.get_ticklabels():
                    elt.set_visible(False)
            ax.w_xaxis._AXINFO = myAXINFO
            ax.w_yaxis._AXINFO = myAXINFO
            ax.w_zaxis._AXINFO = myAXINFO
        else:
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.transformedCP[:,0], self.transformedCP[:,1], self.transformedCP[:,2], edgecolors='none', c=self.contigColours, s=2, marker='.')
            ax.azim = azim
            ax.elev = elev
            ax.set_xlim3d(0,self.scaleFactor)
            ax.set_ylim3d(0,self.scaleFactor)
            ax.set_zlim3d(0,self.scaleFactor)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            if(not showAxis):
                ax.set_axis_off()

        if(fileName != ""):
            try:
                if(all):
                    fig.set_size_inches(3*primaryWidth+2*primarySpace,primaryWidth)
                else:
                    fig.set_size_inches(primaryWidth,primaryWidth)            
                plt.savefig(fileName,dpi=dpi,format=format)
            except:
                print "Error saving image",fileName, exc_info()[0]
                raise
        elif(show):
            try:
                plt.show()
            except:
                print "Error showing image", exc_info()[0]
                raise
        if del_fig:
            plt.close(fig)
            del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

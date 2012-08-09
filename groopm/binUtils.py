#!/usr/bin/env python
###############################################################################
#                                                                             #
#    binUtils.py                                                              #
#                                                                             #
#    Bins Bins Bins                                                           #
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
import mstore

np.seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BinManager:
    """Class used for manipulating bins"""
    def __init__(self, dbFileName):
        self.DB = mstore.DataBlob(dbFileName)
        
        # all about bins
        self.bins = {}                              # bid -> Bin
        self.binCentroidPoints = np.array([])       # array of bin centers
        self.binCentroidColours = np.array([])      # average colour of each bin 
        self.binSizes = {}                          # total size in bp of each bin 
        self.binMembers = {}                        # number of contigs in each bin

#------------------------------------------------------------------------------
# LOADING / SAVING
        
    def loadBins(self,
                 getUnbinned=False,
                 bids=[],
                 makeBins=False,
                 silent=True,
                 loadKmerSigs=False,
                 loadCovProfiles=True):
        """Load data and make bin objects"""
        # fix the condition
        condition=""
        if(len(bids) == 0):
            condition='bid != 0'
        # if we're going to make bins then we'll need kmer sigs
        if(makeBins):
            loadKmerSigs=True
            loadCovProfiles=True
        
        self.DB.loadData(bids=bids,
                         condition=condition,
                         silent=silent,
                         loadCovProfiles=loadCovProfiles,
                         loadKmerSigs=loadKmerSigs,
                         makeColours=True,
                         loadContigNames=True,
                         loadContigLengths=True,
                         loadBins=True,
                         loadCores=False
                        )
        self.initialiseContainers()
        if(makeBins):
            self.DB.transformCP()
            self.makeBins()

    def initialiseContainers(self):
        """Munge the raw data into something more usable"""
        # initialise these containers
        self.binMembers[0] = []
        self.binSizes[0] = 0
        for bid in self.DB.validBinIds.keys():
            self.binSizes[bid] = 0;
            self.binMembers[bid] = []
        # fill them up
        for index in range(0, np.size(self.DB.indicies)):
            self.binMembers[self.DB.binIds[index]].append(index)
            self.binSizes[self.DB.binIds[index]] += self.DB.contigLengths[index]

    def makeBins(self):
        """Make bin objects from loaded data"""
        for bid in self.DB.validBinIds.keys():
            self.bins[bid] = Bin(np.array(self.binMembers[bid]), self.DB.kmerSigs, bid)
            self.bins[bid].makeBinDist(self.DB.transformedCP, self.DB.kmerSigs)       

#------------------------------------------------------------------------------
# BIN UTILITIES 
    
    def split(self, bid, parts):
        """split a bin into n parts"""
        pass

    def merge(self, bids, force=False, newBid=0):
        """Merge two or more bins"""
        pass

#------------------------------------------------------------------------------
# BIN STATS 

    def findCoreCentres(self):
        """Find the point representing the centre of each core"""
        print "\tFinding bin centers"
        self.binCentroidPoints = np.zeros((len(self.bins),3))
        # remake the cores and populate the centres
        S = 1       # SAT and VAL remain fixed at 1. Reduce to make
        V = 1       # Pastels if that's your preference...
        outer_index = 0
        for bid in self.bins.keys():
            self.binCentroidPoints[outer_index] = self.bins[bid].covMeans
            cum_colour = np.array([])
            for index in self.bins[bid].indicies:
                cum_colour = np.append(cum_colour, self.DB.contigColours[index])
            cum_colour = np.reshape(cum_colour, (self.bins[bid].binSize, 3))
            ave_colour = np.mean(cum_colour, axis=0)
            self.binCentroidColours = np.append(self.binCentroidColours, ave_colour)
            outer_index += 1
            
        self.binCentroidColours = np.reshape(self.binCentroidColours, (len(self.bins), 3))            

    def measureBinKVariance(self, outlierTrim=0.1, plot=False):
        """Measure within and between bin variance of kmer sigs
        
        return a list of potentially confounding kmer indicies
        """
        print "\tMeasuring kmer type variances"        
        means = np.array([])
        stdevs = np.array([])
        bids = np.array([])
        
        # work out the mean and stdev for the kmer sigs for each bin
        for bid in self.bins:
            bkworking = np.array([])
            for index in self.bins[bid].indicies:
                bkworking = np.append(bkworking, self.DB.kmerSigs[index])
            bkworking = np.reshape(bkworking, (self.bins[bid].binSize, np.size(self.DB.kmerSigs[0])))
            bids = np.append(bids, [bid])
            means = np.append(means, np.mean(bkworking, axis=0))
            stdevs = np.append(stdevs, np.std(bkworking, axis=0))
            
        means = np.reshape(means, (len(self.bins), np.size(self.DB.kmerSigs[0])))
        stdevs = np.reshape(stdevs, (len(self.bins), np.size(self.DB.kmerSigs[0])))
        
        # now work out the between and within core variances
        between = np.std(means, axis=0)
        within = np.median(stdevs, axis=0)

        B = np.arange(0, np.size(self.DB.kmerSigs[0]), 1)
        names = self.DB.getMerColNames().split(',')
        
        # we'd like to find the indicies of the worst 10% for each type so we can ignore them
        # specifically, we'd like to remove the least variable between core kms and the 
        # most variable within core kms.
        sort_between_indicies = np.argsort(between)
        sort_within_indicies = np.argsort(within)[::-1]
        number_to_trim = int(outlierTrim* float(np.size(self.DB.kmerSigs[0])))
        
        return_indicies =[]
        for i in range(0,number_to_trim):
            if(sort_between_indicies[i] not in return_indicies):
                return_indicies.append(sort_between_indicies[i])
            if(sort_within_indicies[i] not in return_indicies):
                return_indicies.append(sort_within_indicies[i]) 
        
        if(plot):
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

        return return_indicies

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def printBins(self, outFormat, fileName=""):
        """Wrapper for print handles piping to file or stdout"""
        if("" != fileName):
            try:
                # redirect stdout to a file
                sys.stdout = open(fileName, 'w')
                self.printInner()
            except:
                print "Error diverting stout to file:", fileName, sys.exc_info()[0]
                raise
        else:
            self.printInner(outFormat)           
        
    def printInner(self, outFormat):
        """Print bin information to STDOUT"""
        if(outFormat == 'summary'):
            print "#\"bid\"\t\"totalBP\"\t\"numCons\""
            for bid in self.binMembers:
                if(np.size(self.binMembers[bid]) > 0):
                    print str(bid)+"\t"+str(self.binSizes[bid])+"\t"+str(self.DB.validBinIds[bid])
        elif(outFormat == 'full'):
            for bid in self.binMembers:
                if(np.size(self.binMembers[bid]) > 0):
                    print "#bid_"+str(bid)+"_totalBP_"+str(self.binSizes[bid])+"_numCons_"+str(self.DB.validBinIds[bid])
                    print "#\"bid\"\t\"cid\"\t\"length\""            
                    for member in self.binMembers[bid]:
                        print bid, self.DB.contigNames[member], self.DB.contigLengths[member]
        elif(outFormat == 'minimal'):
            print "#\"bid\"\t\"cid\"\t\"length\""            
            for bid in self.binMembers:
                if(np.size(self.binMembers[bid]) > 0):
                    for member in self.binMembers[bid]:
                        print bid, self.DB.contigNames[member], self.DB.contigLengths[member]
            pass
        else:
            print "Error: Unrecognised format:", outFormat

    def plotProfileDistributions(self):
        """Plot the coverage and kmer distributions for each bin"""
        for bid in self.bins:
            self.bins[bid].plotProfileDistributions(self.DB.transformedCP, self.DB.kmerSigs, fileName="PROFILE_"+str(bid))

    def plotBins(self, FNPrefix="BIN"):
        """Make plots of all the bins"""
        for bid in self.bins:
            self.bins[bid].plotBin(self.DB.transformedCP, self.DB.contigColours, fileName=FNPrefix+"_"+str(bid),)

    def plotSideBySide(self, bids, fileName="", tag=""):
        """Plot two bins side by side in 3d"""
        fig = plt.figure()
        front_nums = 100 + len(bids) * 10
        end_num = 1
        for bid in bids:
            spn = front_nums+end_num
            end_num +=1 
            title = self.bins[bid].plotOnAx(fig, self.DB.transformedCP, self.DB.contigColours, subPlotNum=spn, fileName=fileName, tag=tag)
            plt.title(title)
        if(fileName != ""):
            try:
                fig.set_size_inches(12,6)
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BinExplorer:
    """Inspect bins, used for validation"""
    def __init__(self, dbFileName, bids=[]):
        self.DB = mstore.DataBlob(dbFileName)   # based on user specified length
        self.BM = BinManager(dbFileName)        # bins
        if bids is None:
            self.bids = []
        else:
            self.bids = bids

    def plotBinProfiles(self):
        """Plot the distributions of kmer and coverage signatures"""
        self.BM.loadBins(makeBins=True,silent=False,bids=self.bids)
        print "Plotting bin profiles"
        self.BM.plotProfileDistributions()
    
    def plotSideBySide(self, coreCut):
        """Plot cores side by side with their contigs"""
        self.DB.loadData(condition="length >= "+str(coreCut))
        self.DB.transformCP()
        self.BM.loadBins(makeBins=True,bids=self.bids)
        print "Creating side by side plots"
        self.BM.findCoreCentres()
        self.BM.measureBinKVariance()
        self.plotCoresVsContigs()

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def plotCoresVsContigs(self):
        """Render the image for validating cores"""
        fig = plt.figure()
        ax1 = fig.add_subplot(121, projection='3d')
        ax1.scatter(self.DB.transformedCP[:,0], self.DB.transformedCP[:,1], self.DB.transformedCP[:,2], edgecolors=self.DB.contigColours, c=self.DB.contigColours, marker='.')
        ax2 = fig.add_subplot(122, projection='3d')
        ax2.scatter(self.BM.binCentroidPoints[:,0], self.BM.binCentroidPoints[:,1], self.BM.binCentroidPoints[:,2], edgecolors=self.BM.binCentroidColours, c=self.BM.binCentroidColours)
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
        self.indicies = indicies             # all the indicies belonging to this bin
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
        self.merZeros = np.array([])
        self.kDistMean = 0
        self.kDistStdev = 0
        self.kDistTolerance = mertol
        self.kDistUpperLimit = 0

#------------------------------------------------------------------------------
# Tools used for condensing 
    
    def __cmp__(self, alien):
        """Sort bins based on aux values"""
        if self.kDistMean < alien.kDistMean:
            return -1
        elif self.kDistMean == alien.kDistMean:
            return 0
        else:
            return 1

    def removeFromBin(self, transformedCP, kmerSigs, contigLengths, numTigs):
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
        self.makeBinDist(transformedCP, kmerSigs)
        self.calcTotalSize(contigLengths)

        return ret_indicies
        
    def isSimilar(self, compBin, stdevs=1):
        """Check whether two bins are similar"""
        this_lowers = self.covMeans - stdevs * self.covStdevs
        this_uppers = self.covMeans + stdevs * self.covStdevs
        that_lowers = compBin.covMeans - stdevs * compBin.covStdevs
        that_uppers = compBin.covMeans + stdevs * compBin.covStdevs
        # reciprocial test on x and y co-ords only
        for index in range(0,2):
            if(self.covMeans[index] < that_lowers[index] or self.covMeans[index] > that_uppers[index]):
                return False
            if(compBin.covMeans[index] < this_lowers[index] or compBin.covMeans[index] > this_uppers[index]):
                return False
        
        # now test for overlaps in the z dimension (this dimension varys widely and is representitive of gross coverage)
        if((this_uppers[2] > that_lowers[2]) and (this_uppers[2] < that_uppers[2])) or ((that_uppers[2] > this_lowers[2]) and (that_uppers[2] < this_uppers[2])):
            return True

        return False
    
    def consume(self, transformedCP, kmerSigs, contigLengths, deadBin):
        """Combine the contigs of another bin with this one"""
        # consume all the other bins indicies
        self.indicies = np.concatenate([self.indicies, deadBin.indicies])
        self.binSize += deadBin.binSize
        
        # fix the stats on our bin
        self.makeBinDist(transformedCP, kmerSigs)
        self.calcTotalSize(contigLengths)
        
#------------------------------------------------------------------------------
# Stats and properties 

    def clearBinDist(self, kmerSigs):
        """Clear any set distribution statistics"""
        self.covMeans = np.zeros((3))
        self.covStdevs = np.zeros((3))
        self.covLowerLimits = np.zeros((3))
        self.covUpperLimits = np.zeros((3))
        
        self.merMeans = np.array([])
        self.merStdevs = np.array([])
        self.merCentroid = np.zeros((np.size(kmerSigs[0])))
        self.merZeros = np.zeros((np.size(kmerSigs[0])))
        self.kDistMean = 0
        self.kDistStdev = 0
        self.kDistUpperLimit = 0
        
    def makeBinDist(self, transformedCP, kmerSigs):
        """Determine the distribution of the points in this bin
        
        The distribution is largely normal, except at the boundaries.
        """
        self.clearBinDist(kmerSigs)
        if(0 == np.size(self.indicies)):
            return
        
        # Get some data!
        cov_working_array = np.zeros((self.binSize,3))
        mer_working_array = np.zeros((self.binSize,np.size(kmerSigs[0])))
        outer_index = 0
        for index in self.indicies:
            cov_working_array[outer_index] = transformedCP[index]
            #print transformedCP[index]
            mer_working_array[outer_index] = kmerSigs[index]
            self.merCentroid += kmerSigs[index]
            outer_index += 1
        self.merCentroid /= float(np.size(self.indicies))
        
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
            k_dists = np.append(k_dists, np.linalg.norm(sig-self.merZeros))
        self.kDistMean = np.mean(k_dists)
        self.kDistStdev = np.std(k_dists)
        
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
        self.kDistUpperLimit = self.kDistMean + st * self.kDistStdev
        
    def getKDist(self, sig, centroid=None):
        """Get the distance of this sig from the centroid"""
        # z-norm and then distance!
        if centroid is None:
            centroid = self.merZeros
        return np.linalg.norm((sig-self.merMeans)/self.merStdevs - centroid)
    
#------------------------------------------------------------------------------
# Grow the bin 
    
    def recruit(self, transformedCP, kmerSigs, mappedIndicies, binnedIndicies):
        """Iteratively grow the bin"""
        # save these
        pt = self.covTolerance
        st = self.kDistTolerance

        self.binSize = self.indicies.shape[0]
        num_recruited = self.recruitRound(transformedCP, kmerSigs, mappedIndicies, binnedIndicies) 
        while(num_recruited > 0):
            # reduce these to force some kind of convergence
            self.covTolerance *= 0.8
            self.kDistTolerance *= 0.8
            # fix these
            self.binSize = self.indicies.shape[0]
            self.makeBinDist(transformedCP, kmerSigs)
            # go again
            num_recruited = self.recruitRound(transformedCP, kmerSigs, mappedIndicies, binnedIndicies)
        
        # put everything back where we found it...
        self.binSize = self.indicies.shape[0]
        self.covTolerance = pt
        self.kDistTolerance = st
        self.makeBinDist(transformedCP, kmerSigs)
        
        # finally, fix this guy
        return self.binSize
        
    def recruitRound(self, transformedCP, kmerSigs, mappedIndicies, binnedIndicies):
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
    def plotProfileDistributions(self, transformedCP, kmerSigs, fileName=""):
        """plot the profile distibutions for this bin"""
        cov_working_array = np.zeros((self.binSize,3))
        mer_working_array = np.zeros((self.binSize,np.size(kmerSigs[0])))
        outer_index = 0
        for index in self.indicies:
            for i in range(0,3):
                cov_working_array[outer_index][i] = transformedCP[index][i]
            mer_working_array[outer_index] = kmerSigs[index]
            outer_index += 1
        
        # calculate the mean and stdev 
        covMeans = np.mean(cov_working_array,axis=0)
        covStdevs = np.std(cov_working_array,axis=0)
        merMeans = np.mean(mer_working_array, axis=0)
        merStdevs = np.std(mer_working_array, axis=0)

        # z-normalise each column in each working array
        for index in range(0,np.size(self.indicies)):
            mer_working_array[index] = (mer_working_array[index]-merMeans)/merStdevs
            cov_working_array[index] = (cov_working_array[index]-covMeans)/covStdevs
        
        # work out the distribution of distances from z-normed sigs to the centroid
        k_dists = np.array([])
        c_dists = np.array([])
        merZeros = np.zeros((np.size(kmerSigs[0])))
        covZeros = np.zeros((3))
        
        for i in range(0,self.binSize):
            k_dists = np.append(k_dists, np.linalg.norm(mer_working_array[i]-merZeros))
            c_dists = np.append(c_dists, np.linalg.norm(cov_working_array[i]-covZeros))

        k_dists = np.sort(k_dists)
        c_dists = np.sort(c_dists)

        kDistMean = np.mean(k_dists)
        kDistStdev = np.std(k_dists)
        cDistMean = np.mean(c_dists)
        cDistStdev = np.std(c_dists)

        for i in range(0,self.binSize):
            k_dists[i] = (k_dists[i] - kDistMean)/kDistStdev
            c_dists[i] = (c_dists[i] - cDistMean)/cDistStdev
        
        B = np.arange(0, self.binSize, 1)
        
        fig = plt.figure()
        plt.subplot(211)
        plt.plot(B, k_dists, 'r-')
        plt.xlabel("kmer distribution")
        plt.subplot(212)
        plt.plot(B, c_dists, 'b-')
        plt.xlabel("coverage distribution")
        if(fileName != ""):
            try:
                fig.set_size_inches(10,4)
                plt.savefig(fileName,dpi=300)
            except:
                print "Error saving image:", fileName, sys.exc_info()[0]
                raise
        else:
            try:
                plt.show()
            except:
                print "Error showing image:", sys.exc_info()[0]
                raise
        del fig
            
        
    def plotBin(self, transformedCP, contigColours, fileName="", tag=""):
        """Plot a single bin"""
        fig = plt.figure()
        title = self.plotOnAx(fig, transformedCP, contigColours, fileName=fileName, tag=tag)
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

    def plotOnAx(self, fig, transformedCP, contigColours, subPlotNum=111, fileName="", tag=""):
        """Plot a bin in a given subplot"""
        disp_vals = np.array([])
        disp_cols = np.array([])
        num_points = 0
        for index in self.indicies:
            num_points += 1
            disp_vals = np.append(disp_vals, transformedCP[index])
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

        ax = fig.add_subplot(subPlotNum, projection='3d')
        ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors=disp_cols, c=disp_cols, marker='.')
        from locale import format, setlocale, LC_ALL # purdy commas
        setlocale(LC_ALL, "")
        title = str.join(" ", ["Bin:",str(self.id),"--",tag,"\n",
                               "Focus at: (",str(px), str(py), str(pz),")\n",
                               "Contains:",str(self.binSize),"contigs\n",
                               "Total:",format('%d', self.totalBP, True),"BP"
                               ])
        return title
    
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

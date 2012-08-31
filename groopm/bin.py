#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bin.py                                                                   #
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

import networkx as nx

import time

# GroopM imports
import dataManagers

np.seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Bin:
    """Class for managing collections of contigs
    
    To (perhaps) simplify things think of a "bin" as an row_index into the
    column names array. The ClusterBlob has a list of bins which it can
    update etc...
    """
    def __init__(self, rowIndicies, kmerSigs, id, upperCov, covtol=2, mertol=2):
        self.id = id
        self.rowIndicies = rowIndicies             # all the indicies belonging to this bin
        self.binSize = self.rowIndicies.shape[0]
        self.upperCov = upperCov
        self.totalBP = 0

        self.covTolerance = covtol
        self.kDistTolerance = mertol
        
        # we need some objects to manage the distribution of contig proerties
        self.covMeans = np.zeros((3))
        self.covStdevs = np.zeros((3))
        self.covLowerLimits = np.zeros((3)) # lower and upper limits based on tolerance
        self.covUpperLimits = np.zeros((3))
        
        self.merMeans = np.array([])
        self.merStdevs = np.array([])
        self.merZeros = np.zeros((np.size(kmerSigs[0])))
        self.kDistMean = 0.0
        self.kDistStdev = 0.0
        self.kValMean = 0.0
        self.kValStdev = 0.0
        self.kDistUpperLimit = 0.0

#------------------------------------------------------------------------------
# Tools used for comparing / condensing 
    
    def __cmp__(self, alien):
        """Sort bins based on PC1 of kmersig values"""
        if self.kValMean < alien.kValMean:
            return -1
        elif self.kValMean == alien.kValMean:
            return 0
        else:
            return 1

#------------------------------------------------------------------------------
# Grow and shrink 
    
    def consume(self, transformedCP, kmerSigs, contigLengths, contigColours, deadBin, verbose=False):
        """Combine the contigs of another bin with this one"""
        # consume all the other bins rowIndicies
        if(verbose):
            print "    BIN:",deadBin.id,"will be consumed by BIN:",self.id
        self.rowIndicies = np.concatenate([self.rowIndicies, deadBin.rowIndicies])
        self.binSize  = self.rowIndicies.shape[0]
        
        # fix the stats on our bin
        self.makeBinDist(transformedCP, kmerSigs)
        self.calcTotalSize(contigLengths)
        self.getKmerColourStats(contigColours)

    def purge(self, deadIndicies, transformedCP, kmerSigs, contigLengths, contigColours):
        """Delete some rowIndicies and remake stats"""
        old_ri = self.rowIndicies
        self.rowIndicies = np.array([])
        for i in old_ri:
            if i not in deadIndicies:
                self.rowIndicies = np.append(self.rowIndicies, i)
            
        # fix the stats on our bin
        self.makeBinDist(transformedCP, kmerSigs)
        self.calcTotalSize(contigLengths)
        self.getKmerColourStats(contigColours)
        
#------------------------------------------------------------------------------
# Stats and properties 

    def clearBinDist(self, kmerSigs):
        """Clear any set distribution statistics"""
        self.totalBP = 0
        
        self.covMeans = np.zeros((3))
        self.covStdevs = np.zeros((3))
        self.covLowerLimits = np.zeros((3)) # lower and upper limits based on tolerance
        self.covUpperLimits = np.zeros((3))
        
        self.merMeans = np.zeros((np.size(kmerSigs[0])))
        self.merStdevs = np.zeros((np.size(kmerSigs[0])))
        self.kDistMean = 0.0
        self.kDistStdev = 0.0
        self.kDistUpperLimit = 0.0
        self.kValMean = 0.0
        self.kValStdev = 0.0
        
    def makeBinDist(self, transformedCP, kmerSigs):
        """Determine the distribution of the points in this bin
        
        The distribution is largely normal, except at the boundaries.
        """
        #print "MBD", self.id, self.binSize 
        self.binSize = self.rowIndicies.shape[0]
        self.kValMean = 0.0
        self.kValStdev = 0.0
        if(0 == np.size(self.rowIndicies)):
            return

        # get the centroids
        (self.covMeans, self.covStdevs) = self.getCentroidStats(transformedCP)
        (self.merMeans, self.merStdevs) = self.getCentroidStats(kmerSigs)
        
        # work out the distribution of distances of kmersigs
        (self.kDistMean, self.kDistStdev, range) = self.getInnerVariance(kmerSigs, mode="kmer")
        
        # set the acceptance ranges
        self.makeLimits()
        
    def makeLimits(self, covTol=-1, merTol=-1):
        """Set inclusion limits based on mean, variance and tolerance settings"""
        if(-1 == covTol):
            covTol=self.covTolerance
        if(-1 == merTol):
            merTol=self.kDistTolerance
        for i in range(0,3):
            self.covLowerLimits[i] = int(self.covMeans[i] - covTol * self.covStdevs[i])
            if(self.covLowerLimits[i] < 0):
                self.covLowerLimits[i] = 0.0
            self.covUpperLimits[i] = int(self.covMeans[i] + covTol * self.covStdevs[i])
            if(self.covUpperLimits[i] > self.upperCov):
                self.covUpperLimits[i] = self.upperCov            
        self.kDistUpperLimit = self.kDistMean + merTol * self.kDistStdev

    def calcTotalSize(self, contigLengths):
        """Work out the total size of this bin in BP"""
        totalBP = 0
        for row_index in self.rowIndicies:
            totalBP += contigLengths[row_index]
        self.totalBP = totalBP

    def getCentroidStats(self, profile):
        """Calculate the centroids of the profile"""
        working_list = np.zeros((self.binSize, np.size(profile[0])))
        outer_index = 0
        for row_index in self.rowIndicies:
            working_list[outer_index] = profile[row_index]
            outer_index += 1
        # return the mean and stdev 
        return (np.mean(working_list,axis=0), np.std(working_list,axis=0))
        
    def getInnerVariance(self, profile, mode="kmer"):
        """Work out the variance for the coverage/kmer profile"""
        dists = []
        if(mode == "kmer"):
            for row_index in self.rowIndicies:
                dist = self.getKDist(profile[row_index])
                dists.append(dist)
        elif(mode =="cov"):
            for row_index in self.rowIndicies:
                dist = self.getCDist(profile[row_index])
                dists.append(dist)
        else:
            raise ModeNotAppropriateException("Mode",mode,"unknown")
        range = np.max(np.array(dists)) - np.min(np.array(dists))
        return (np.mean(np.array(dists)), np.std(np.array(dists)), range)
        
    def getKDist(self, Ksig, centroid=None):
        """Get the distance of this kmer sig from the centroid"""
        # z-norm and then distance!
        if centroid is None:
            centroid = self.merMeans
        return np.linalg.norm(Ksig-centroid)

    def getCDist(self, Csig, centroid=None):
        """Get the distance of this contig from the coverage centroid"""
        # z-norm and then distance!
        if centroid is None:
            centroid = self.covMeans
        return np.linalg.norm(Csig-centroid)
    
    def getKmerColourStats(self, contigColours):
        """Determine the mean and stdev of the kmer profile colours"""
        kmer_vals = np.array([])
        for row_index in self.rowIndicies:
            kmer_vals = np.append(kmer_vals, colorsys.rgb_to_hsv(contigColours[row_index][0],
                                                                 contigColours[row_index][1],
                                                                 contigColours[row_index][2]
                                                                 )[0]
                                  )

        self.kValMean = np.mean(kmer_vals)
        self.kValStdev = np.std(kmer_vals)
  
    def findOutliers(self, transformedCP, kmerSigs, percent=0.1, mode="kmer"):
        """Return the list of row indicies which least match the profile of the bin"""

        # check we're not trying to do something stupid
        num_to_purge = int(self.binSize * percent)
        if(num_to_purge == self.binSize):
            return []

        # make a list of all the profile distances
        dists = []
        if(mode == "kmer"):
            for row_index in self.rowIndicies:
                dists.append(self.getKDist(kmerSigs[row_index]))
        elif(mode =="cov"):
            for row_index in self.rowIndicies:
                dists.append(self.getCDist(transformedCP[row_index]))
        else:
            raise ModeNotAppropriateException("Mode",mode,"unknown")
        
        # find the bottom x
        sorted_dists = np.argsort(dists)[::-1]
        ret_list = []
        for i in range(num_to_purge):
            ret_list.append(self.rowIndicies[sorted_dists[i]])
        return ret_list
        
#------------------------------------------------------------------------------
# Grow the bin 
    
    def recruit(self, transformedCP, kmerSigs, im2RowIndicies, binnedRowIndicies):
        """Iteratively grow the bin"""
        # save these
        pt = self.covTolerance
        st = self.kDistTolerance

        self.binSize = self.rowIndicies.shape[0]
        num_recruited = self.recruitRound(transformedCP, kmerSigs, im2RowIndicies, binnedRowIndicies) 
        while(num_recruited > 0):
            # reduce these to force some kind of convergence
            self.covTolerance *= 0.8
            self.kDistTolerance *= 0.8
            # fix these
            self.binSize = self.rowIndicies.shape[0]
            self.makeBinDist(transformedCP, kmerSigs)
            # go again
            num_recruited = self.recruitRound(transformedCP, kmerSigs, im2RowIndicies, binnedRowIndicies)
        
        # put everything back where we found it...
        self.binSize = self.rowIndicies.shape[0]
        self.covTolerance = pt
        self.kDistTolerance = st
        self.makeBinDist(transformedCP, kmerSigs)
        
        # finally, fix this guy
        return self.binSize
        
    def recruitRound(self, transformedCP, kmerSigs, im2RowIndicies, binnedRowIndicies):
        """Recruit more points in from outside the current blob boundaries"""
        #print "RRR",self.id,self.binSize,"(",self.covLowerLimits[0],self.covUpperLimits[0],")","(",self.covLowerLimits[1],self.covUpperLimits[1],")","(",self.covLowerLimits[2],self.covUpperLimits[2],")" 
        num_recruited = 0
        for x in range(int(self.covLowerLimits[0]), int(self.covUpperLimits[0])):
            for y in range(int(self.covLowerLimits[1]), int(self.covUpperLimits[1])):
                for z in range(int(self.covLowerLimits[2]), int(self.covUpperLimits[2])):
                    if((x,y,z) in im2RowIndicies):
                        for row_index in im2RowIndicies[(x,y,z)]:
                            if (row_index not in binnedRowIndicies) and (row_index not in self.rowIndicies):
                                k_dist = self.getKDist(kmerSigs[row_index])
                                if(k_dist <= self.kDistUpperLimit):
                                    self.rowIndicies = np.append(self.rowIndicies,row_index)
                                    num_recruited += 1
        return num_recruited

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 
#
    def plotProfileDistributions(self, transformedCP, kmerSigs, fileName=""):
        """plot the profile distibutions for this bin"""
        cov_working_array = np.zeros((self.binSize,3))
        mer_working_array = np.zeros((self.binSize,np.size(kmerSigs[0])))
        outer_index = 0
        for row_index in self.rowIndicies:
            for i in range(0,3):
                cov_working_array[outer_index][i] = transformedCP[row_index][i]
            mer_working_array[outer_index] = kmerSigs[row_index]
            outer_index += 1
        
        # calculate the mean and stdev 
        covMeans = np.mean(cov_working_array,axis=0)
        covStdevs = np.std(cov_working_array,axis=0)
        merMeans = np.mean(mer_working_array, axis=0)
        merStdevs = np.std(mer_working_array, axis=0)

        # z-normalise each column in each working array
        for index in range(0,np.size(self.rowIndicies)):
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
            
        
    def plotBin(self, transformedCP, contigColours, fileName=""):
        """Plot a single bin"""
        fig = plt.figure()
        title = self.plotOnAx(fig, 1, 1, 1, transformedCP, contigColours, fileName=fileName)
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

    def plotOnAx(self, fig, plot_rows, plot_cols, plot_num, transformedCP, contigColours, fileName=""):
        """Plot a bin in a given subplot"""
        disp_vals = np.array([])
        disp_cols = np.array([])
        num_points = 0
        for row_index in self.rowIndicies:
            num_points += 1
            disp_vals = np.append(disp_vals, transformedCP[row_index])
            disp_cols = np.append(disp_cols, contigColours[row_index])

        # make a black mark at the max values
        self.makeLimits()
        px = int(self.covMeans[0])
        py = int(self.covMeans[1])
        pz = int(self.covMeans[2])
        num_points += 1
        disp_vals = np.append(disp_vals, [px,py,pz])
        disp_cols = np.append(disp_cols, colorsys.hsv_to_rgb(0,0,0))
        
        # fix these
        self.makeLimits()
        self.getKmerColourStats(contigColours)

        # reshape
        disp_vals = np.reshape(disp_vals, (num_points, 3))
        disp_cols = np.reshape(disp_cols, (num_points, 3))

        ax = fig.add_subplot(plot_rows, plot_cols, plot_num, projection='3d')
        ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors=disp_cols, c=disp_cols, marker='.')
        from locale import format, setlocale, LC_ALL # purdy commas
        setlocale(LC_ALL, "")
        title = str.join(" ", ["Bin:",str(self.id),":",str(self.binSize),"contigs : ",format('%d', self.totalBP, True),"BP\n",
                               "Coverage centroid: (",str(px), str(py), "[",str(self.covLowerLimits[2]),"-",str(self.covUpperLimits[2]),"])\n",
                               "Kmers: mean: %.4f stdev: %.4f" % (self.kValMean, self.kValStdev),"\n",
                               
                               ])
        return title

    def printBin(self, contigNames, contigLengths, outFormat="summary", separator="\t"):
        """print this bin info in csvformat"""
        kvm_str = "%.4f" % self.kValMean
        kvs_str = "%.4f" % self.kValStdev
        if(outFormat == 'summary'):
            #print separator.join(["#\"bid\"","\"totalBP\"","\"numCons\"","\"kMean\"","\"kStdev\""]) 
            print separator.join([str(self.id), str(self.totalBP), str(self.binSize), kvm_str, kvs_str])
        elif(outFormat == 'full'):
            print("#bid_"+str(self.id)+
                  "_totalBP_"+str(self.totalBP)+
                  "_numCons_"+str(self.binSize)+
                  "_kMean_"+kvm_str+
                  "_kStdev_"+kvs_str
                  )
            print separator.join(["#\"bid\"","\"cid\"","\"length\""])
            for row_index in self.rowIndicies:
                print separator.join([str(self.id), contigNames[row_index], str(contigLengths[row_index])])
        elif(outFormat == 'minimal'):
            #print separator.join(["#\"bid\"","\"cid\"","\"length\""])            
            for row_index in self.rowIndicies:
                print separator.join([str(self.id), contigNames[row_index], str(contigLengths[row_index])])
        else:
            print "--------------------------------------"
            print "Bin:", self.id
            print "Bin size:", self.binSize
            print "Total BP:", self.totalBP
            print "--------------------------------------"

###############################################################################
###############################################################################
###############################################################################
###############################################################################

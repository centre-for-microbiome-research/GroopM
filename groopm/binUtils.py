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
class BinNotFoundException(BaseException):
    pass

import traceback
class Tracer:
    def __init__(self, oldstream):
        self.oldstream = oldstream
        self.count = 0
        self.lastStack = None
   
    def write(self, s):
        newStack = traceback.format_stack()
        if newStack != self.lastStack:
            self.oldstream.write("".join(newStack))
            self.lastStack = newStack
        self.oldstream.write(s)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BinManager:
    """Class used for manipulating bins"""
    def __init__(self, dbFileName="", pm=None):
        # data storage
        if(dbFileName != ""):
            self.PM = mstore.ProfileManager(dbFileName)
        elif(pm is not None):
            self.PM = pm
        
        # all about bins
        self.nextFreeBinId = 0                      # increment before use!
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
        
        self.PM.loadData(bids=bids,
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
            self.PM.transformCP(silent=silent)
            self.makeBins()

    def initialiseContainers(self):
        """Munge the raw data into something more usable"""
        # initialise these containers
        self.binMembers[0] = []
        self.binSizes[0] = 0
        for bid in self.PM.validBinIds.keys():
            self.binSizes[bid] = 0;
            self.binMembers[bid] = []
        
        # fill them up
        for row_index in range(0, np.size(self.PM.indicies)):
            self.binMembers[self.PM.binIds[row_index]].append(row_index)
            self.binSizes[self.PM.binIds[row_index]] += self.PM.contigLengths[row_index]

        # we need to get the largest BinId in use
        bids = self.PM.getBinStats().keys()
        if(len(bids) > 0):
            self.nextFreeBinId = np.max(bids)

    def makeBins(self):
        """Make bin objects from loaded data"""
        for bid in self.PM.validBinIds.keys():
            self.bins[bid] = Bin(np.array(self.binMembers[bid]), self.PM.kmerSigs, bid)
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)      
            self.bins[bid].calcTotalSize(self.PM.contigLengths)

    def saveBins(self, doCores=True, saveBinStats=True, updateBinStats=True):
        """Save binning results"""
        c2b_update = {}
        core_update = {}
        if doCores:
            (c2b_update, core_update) = self.getCoreBinUpdates()
            self.PM.saveCores(core_update)
        else:
            c2b_update = self.getBinUpdates()
        self.PM.saveBinIds(c2b_update)
        self.PM.setClustered()
        if saveBinStats:
            self.saveBinStats()
        elif updateBinStats:
            # we need to make a list of updates
            updates = {}
            for bid in self.bins.keys():
                updates[bid] = self.bins[bid].binSize
            self.updateBinStats(updates)
    
    def saveBinStats(self):
        """Update / overwrite the table holding the bin stats
        
        Note that this call effectively nukes the existing table
        and should only be used during initial coring. BID must be somewhere!
        """
        bin_updates = {}
        for bid in self.bins:
            bin_updates[bid] = np.size(self.bins[bid].rowIndicies)
        self.PM.saveValidBinIds(bin_updates)

    def updateBinStats(self, updates):
        """Update the table holding the bin stats
        
        This is an update in the true sense of the word
        updates looks like { bid : numMembers }
        if numMember == 0 then this bid is removed
        from the table
        """
        self.PM.updateValidBinIds(updates)

    def getCoreBinUpdates(self):
        """Merge the bids, raw DB indexes and core information so we can save to disk"""
        # at this stage, all bins are cores
        core_update = {}
        for row_index in range(len(self.PM.indicies)):
            if row_index in self.PM.binnedRowIndicies:
                core_update[self.PM.indicies[row_index]] = True

        bin_update = self.getBinUpdates()

        return (bin_update, core_update)

    def getBinUpdates(self):
        """Merge the bids, raw DB indexes and core information so we can save to disk"""
        # we need a mapping from cid (or local index) to binID
        bin_update = {}
        c2b = {}
        for bid in self.bins:
            for row_index in self.bins[bid].rowIndicies:
                c2b[row_index] = bid

        for row_index in range(len(self.PM.indicies)):
            if row_index in self.PM.binnedRowIndicies and row_index in c2b:
                bin_update[self.PM.indicies[row_index]] = c2b[row_index]

        return bin_update

#------------------------------------------------------------------------------
# BIN UTILITIES 

    def condenseWrapper(self,
                      lowerKmerStdev,      # lower limits mean we can just munch away
                      lowerCoverageStdev,
                      upperKmerStdev,      # upper limits mean we can consume but we need to ask for permission first
                      upperCoverageStdev, 
                      auto=False,          # do we need to ask permission at all?
                      save=False
                      ):
        """Iterative wrapper for the condense function"""
        self.printMergeInstructions()
        total_num_bins_condensed = 0
        num_bins_condensed = 0
        while True: # do while loop anyone?
            num_bins_condensed = self.condenseBins(lowerKmerStdev,
                                                      lowerCoverageStdev,
                                                      upperKmerStdev,
                                                      upperCoverageStdev,
                                                      auto=auto,      
                                                      save=save
                                                      )
            total_num_bins_condensed += num_bins_condensed 
            if(num_bins_condensed == 0):
                break
        print "\t",total_num_bins_condensed,"bins condensed.",len(self.bins),"bins remain"
                        
    def condenseBins(self, 
                      lowerKmerStdev,      # lower limits mean we can just munch away
                      lowerCoverageStdev,
                      upperKmerStdev,      # upper limits mean we can consume but we need to ask for permission first
                      upperCoverageStdev, 
                      verbose=False,
                      auto=False,          # do we need to ask permission at all?
                      save=False           # should we save the merges as we go?
                      ):
        """combine similar bins"""
        any_merged = False
        merged = []       # who is getting merged by who?
        merged_order = [] # and in what order to merge!

        # go through all the bins, sorted according to kmer profile
        num_cores_merged = 0
        
        # we need all the bins
        ordered_bids = self.bins.keys()
        num_ordered_bids = len(ordered_bids)  

        # now do an all versus all search
        for subject_index in range(num_ordered_bids):
            if(subject_index not in merged):
                subject_bin = self.bins[ordered_bids[subject_index]]
                subject_upper_kmer_limit = subject_bin.kDistMean + upperKmerStdev * subject_bin.kDistStdev
                subject_lower_kmer_limit = subject_bin.kDistMean + lowerKmerStdev * subject_bin.kDistStdev
                # need only search the triangle!  
                query_index = subject_index+1
                while(query_index < num_ordered_bids):
                    if(query_index not in merged):
                        query_bin = self.bins[ordered_bids[query_index]]
                        query_upper_kmer_limit = query_bin.kDistMean + upperKmerStdev * query_bin.kDistStdev 
                        query_lower_kmer_limit = query_bin.kDistMean + lowerKmerStdev * query_bin.kDistStdev 
                        
                        # only worth comparing if their kmerSigs are similar
                        continue_merge = False
                        ask_kmer = True           # in the case where we may need to ask permission
                        ask_permission = True     # we'll check these two variables
                        k_dists = np.array([])
                        median_k_dist = 1000000     # something large
                        # pick the bin with the highest limit
                        if(query_upper_kmer_limit > subject_upper_kmer_limit):
                            for row_index in subject_bin.rowIndicies:
                                k_dists = np.append(k_dists, query_bin.getKDist(self.PM.kmerSigs[row_index]))
                            median_k_dist = np.median(k_dists) 
                            if median_k_dist <= query_upper_kmer_limit:
                                continue_merge = True
                                if median_k_dist < query_lower_kmer_limit:
                                    ask_kmer = False
                        else:
                            for row_index in query_bin.rowIndicies:
                                k_dists = np.append(k_dists, subject_bin.getKDist(self.PM.kmerSigs[row_index]))
                            median_k_dist = np.median(k_dists) 
                            if median_k_dist <= subject_upper_kmer_limit:
                                continue_merge = True
                                if median_k_dist < subject_lower_kmer_limit:
                                    ask_kmer = False

                        # the coverage test is symmetrical, so it doesn't matter
                        # who compares with who
                        if(continue_merge):
                            # this is a recursive function which tests both limits
                            similar = subject_bin.isSimilar(query_bin, stdevs=[upperCoverageStdev, lowerCoverageStdev])
                            if(0 == similar):
                                # no need to ask permission, unless kmers said so!
                                ask_permission = ask_kmer
                            elif(1 == similar):
                                # we need to ask on the bounds of coverage profile, who cares about kmers
                                ask_permission = True
                            else:
                                # coverage test failed
                                continue_merge = False

                        if(auto):       # override if we've been told to
                            ask_permission = False

                        # we ask permission at the time of merging, just save te flag for now
                        if(continue_merge):
                            merged.append(query_index) # subject consumes query
                            merged_order.append(([subject_bin.id,query_bin.id],ask_permission))   
                    query_index += 1
            subject_index += 1
        
        # Merge them in the order they were seen in
        for merge_job in merged_order:
            num_cores_merged += 1
            self.merge(merge_job[0], auto=(merge_job[1]^True), saveBins=save, verbose=verbose, printInstructions=False)
            any_merged = True

        return num_cores_merged

    def chimeraWrapper(self, auto=False, save=False, verbose=False, printInstructions=False):
        """Automatically search for and separate chimeric bins"""
        b_means = np.array([])
        b_stdevs = np.array([])
        names = []
        nums = []
        
        for bid in self.bins:
            self.bins[bid].getKmerColourStats(self.PM.contigColours)
            b_means = np.append(b_means, self.bins[bid].kValMean)
            b_stdevs = np.append(b_stdevs, self.bins[bid].kValStdev)
            names.append("BIN: %d" % bid)
            nums.append(np.log10(self.bins[bid].binSize)/10)
            
        sorted_stdevs = np.argsort(b_stdevs)
        sorted_means = np.argsort(b_means)
        
        sm = []
        nm = []
        cm = []
        ss = []
        ns = []
        cs = []
        stdev_median = np.median(b_stdevs)
        for i in range(len(self.bins)):
            if(b_stdevs[sorted_stdevs[i]] > stdev_median):
                ss.append(b_stdevs[sorted_stdevs[i]])
                ns.append(names[sorted_stdevs[i]])
                cs.append(nums[sorted_stdevs[i]])
            if(b_stdevs[sorted_means[i]] > stdev_median):
                sm.append(b_stdevs[sorted_means[i]])
                nm.append(names[sorted_means[i]])
                cm.append(nums[sorted_means[i]])

        bs = np.arange(0, len(ns), 1)
        bm = np.arange(0, len(nm), 1)
        
        plt.figure(1)
        
        plt.subplot(211)
        plt.plot(bs, ss, 'b--', bs, cs, 'g--')
        plt.xticks(bs, ns, rotation=90)
        plt.grid()

        plt.subplot(212)
        plt.plot(bm, sm, 'b--', bm, cm, 'g--')
        plt.xticks(bm, nm, rotation=90)
        plt.grid()
        
        plt.show()

    def split(self, bid, n, mode='kmer', auto=False, saveBins=False, verbose=False, printInstructions=True):
        """split a bin into n parts"""
        # we need to work out which profile to cluster on
        if(printInstructions):
            self.printSplitInstructions()
       
        bin = self.getBin(bid)
        obs = np.array([])
        if(mode=='kmer'):
            for row_index in bin.rowIndicies:
                obs = np.append(obs, self.PM.kmerSigs[row_index])
            obs = np.reshape(obs, (len(bin.rowIndicies),len(self.PM.kmerSigs[0])))
        elif(mode=='cov'):
            for row_index in bin.rowIndicies:
                obs = np.append(obs, self.PM.covProfiles[row_index])
            obs = np.reshape(obs, (len(bin.rowIndicies),len(self.PM.covProfiles[0])))
        
        # do the clustering
        from scipy.cluster.vq import kmeans,vq
        centroids,_ = kmeans(obs,n)
        idx,_ = vq(obs,centroids)
        
        # plot some stuff
        idx_sorted = np.argsort(np.array(idx))
        current_group = -1
        bids = [bid]
        holding_array = np.array([])
        for i in idx_sorted:
            if(idx[i] != current_group):
                if(current_group != -1):
                    # bin is full!
                    split_bin = self.makeNewBin(holding_array)
                    split_bin.makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
                    split_bin.calcTotalSize(self.PM.contigLengths)
                    bids.append(split_bin.id)
                    holding_array = np.array([])
                current_group = idx[i]
            holding_array = np.append(holding_array, bin.rowIndicies[i])
        if(np.size(holding_array) != 0):
            split_bin = self.makeNewBin(holding_array)
            split_bin.makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
            split_bin.calcTotalSize(self.PM.contigLengths)
            bids.append(split_bin.id)
            
        self.plotSideBySide(bids)
        continue_split = self.promptOnSplit(n)
        if(continue_split == 'Y'):
            # go ahead
            pass
        elif(continue_split == 'N'):
            return
        elif(continue_split == 'C'):
            if(mode == "cov"):
                print "Already doing split based on coverage profile"
            else:
                self.split(bid, n, mode='cov', auto=auto, saveBins=saveBins, verbose=verbose, printInstructions=False)
        elif(continue_split == 'K'):
            if(mode == "kmer"):
                print "Already doing split based on kmer profile"
            else:
                self.split(bid, n, mode='kmer', auto=auto, saveBins=saveBins, verbose=verbose, printInstructions=False)
        

    def merge(self, bids, auto=False, newBid=False, saveBins=False, verbose=False, printInstructions=True):
        """Merge two or more bins"""
        parent_bin = None
        bin_stats = {}
        if(printInstructions):
            self.printMergeInstructions()

        if(newBid):
            # we need to make this into a new bin
            parent_bin = makeNewBin()
            # now merge it with the first in the new list
            dead_bin = self.getBin(bids[0])
            parent_bin.consume(self.PM.transformedCP, self.PM.kmerSigs, self.PM.contigLengths, dead_bin, verbose=verbose)
            self.deleteBin(bids[0])
            bin_stats[bids[0]] = 0
            bin_stats[parent_bin.id] = parent_bin.binSize
        else:
            # just use the first given as the parent
            parent_bin = self.getBin(bids[0])
        
        # let this guy consume all the other guys
        for i in range(1,len(bids)):
            continue_merge = False
            dead_bin = self.getBin(bids[i])
            if(auto):
                continue_merge = True
            else:
                self.plotSideBySide([parent_bin.id,dead_bin.id])
                if(self.promptOnMerge(bids=[parent_bin.id,dead_bin.id])):
                    continue_merge=True
            if(continue_merge):
                parent_bin.consume(self.PM.transformedCP, self.PM.kmerSigs, self.PM.contigLengths, dead_bin, verbose=verbose)
                self.deleteBin(bids[i])
                bin_stats[bids[i]] = 0
                bin_stats[parent_bin.id] = parent_bin.binSize

        if(saveBins):
            self.updateBinStats(bin_stats)
            self.saveBins(doCores=False, saveBinStats=False)
            
    def printMergeInstructions(self):
        raw_input( "****************************************************************\n" +
                   " MERGING INSTRUCTIONS - PLEASE READ CAREFULLY\n"+
                   "****************************************************************\n" +
                   " The computer cannot always be trusted to perform bin mergers\n"
                   " automatically, so during merging you may be shown a 3D plot\n"
                   " which should help YOU determine whether or not the bins should\n" +
                   " be merged. Look carefully at each plot and then close the plot\n" +
                   " to continue with the merging operation.\n\n" +
                   " Press any key to produce plots...")
        print "****************************************************************"        

    def promptOnMerge(self, bids=[], minimal=False):
        """Check that the user is ok with this merge"""
        input_not_ok = True
        valid_responses = ['Y','N']
        bin_str = ""
        if(len(bids) != 0):
            bin_str = ": "+str(bids[0])
            for i in range(1, len(bids)):
                bin_str += " and "+str(bids[i])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" Merge? (y,n) : ")
            else: 
                option = raw_input(" ****WARNING**** About to merge bins"+bin_str+"\n" \
                                   " If you continue you *WILL* overwrite existing bins!\n" \
                                   " You have been shown a 3d plot of the bins to be merged.\n" \
                                   " Continue only if you're sure this is what you want to do!\n" \
                                   " Merge? (y,n) : ")
            if(option.upper() in valid_responses):
                if(option.upper() != "Y"):
                    print "Merge skipped"
                    print "****************************************************************"
                    return False
                print "****************************************************************"
                return True
            else:
                print "Error, unrecognised choice '"+option.upper()+"'"
                minimal = True

    def printSplitInstructions(self):
        raw_input( "****************************************************************\n" +
                   " SPLITTING INSTRUCTIONS - PLEASE READ CAREFULLY\n"+
                   "****************************************************************\n" +
                   " The computer cannot always be trusted to perform bin splits\n"
                   " automatically, so during splitting you may be shown a 3D plot\n"
                   " which should help YOU determine whether or not the bin should\n" +
                   " be split. Look carefully at each plot and then close the plot\n" +
                   " to continue with the splitting operation.\n\n" +
                   " Press any key to produce plots...")
        print "****************************************************************"        

    def promptOnSplit(self, parts, minimal=False):
        """Check that the user is ok with this split"""
        input_not_ok = True
        valid_responses = ['Y','N','C','K']
        while(input_not_ok):
            if(minimal):
                option = raw_input(" Split? (y,n,c,k) : ")
            else: 
                option = raw_input(" ****WARNING**** About to split bin into "+str(parts)+" parts\n" \
                                   " If you continue you *WILL* overwrite existing bins!\n" \
                                   " You have been shown a 3d plot of the bin after splitting.\n" \
                                   " Continue only if you're sure this is what you want to do!\n" \
                                   " y = yes, n = no, c = redo but use coverage profile,\n" \
                                   " k = redo but use kmer profile\n" \
                                   " Split? (y,n,c,k) : ")
            if(option.upper() in valid_responses):
                print "****************************************************************"
                return option.upper()
            else:
                print "Error, unrecognised choice '"+option.upper()+"'"
                minimal = True

    def getBin(self, bid):
        """get a bin or raise an error"""
        if bid in self.bins:
            return self.bins[bid]
        else:
            raise BinNotFoundException("Cannot find: "+str(bid)+" in bins dicts")
            
    def deleteBin(self, bid):
        """Purge a bin from our lists"""
        if bid in self.bins:
            del self.bins[bid]
        else:
            raise BinNotFoundException("Cannot find: "+str(bid)+" in bins dicts")
        
    def makeNewBin(self, rowIndicies=np.array([]), bid=None):
        """Make a new bin and add to the list of existing bins"""
        if bid is None:
            self.nextFreeBinId +=1
            bid = self.nextFreeBinId
        self.bins[bid] = Bin(rowIndicies, self.PM.kmerSigs, bid)        
        return self.bins[bid]
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
            for row_index in self.bins[bid].rowIndicies:
                cum_colour = np.append(cum_colour, self.PM.contigColours[row_index])
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
            for row_index in self.bins[bid].rowIndicies:
                bkworking = np.append(bkworking, self.PM.kmerSigs[row_index])
            bkworking = np.reshape(bkworking, (self.bins[bid].binSize, np.size(self.PM.kmerSigs[0])))
            bids = np.append(bids, [bid])
            means = np.append(means, np.mean(bkworking, axis=0))
            stdevs = np.append(stdevs, np.std(bkworking, axis=0))
            
        means = np.reshape(means, (len(self.bins), np.size(self.PM.kmerSigs[0])))
        stdevs = np.reshape(stdevs, (len(self.bins), np.size(self.PM.kmerSigs[0])))
        
        # now work out the between and within core variances
        between = np.std(means, axis=0)
        within = np.median(stdevs, axis=0)

        B = np.arange(0, np.size(self.PM.kmerSigs[0]), 1)
        names = self.PM.getMerColNames().split(',')
        
        # we'd like to find the indicies of the worst 10% for each type so we can ignore them
        # specifically, we'd like to remove the least variable between core kms and the 
        # most variable within core kms.
        sort_between_indicies = np.argsort(between)
        sort_within_indicies = np.argsort(within)[::-1]
        number_to_trim = int(outlierTrim* float(np.size(self.PM.kmerSigs[0])))
        
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
                    print str(bid)+"\t"+str(self.binSizes[bid])+"\t"+str(self.PM.validBinIds[bid])
        elif(outFormat == 'full'):
            for bid in self.binMembers:
                if(np.size(self.binMembers[bid]) > 0):
                    print "#bid_"+str(bid)+"_totalBP_"+str(self.binSizes[bid])+"_numCons_"+str(self.PM.validBinIds[bid])
                    print "#\"bid\"\t\"cid\"\t\"length\""            
                    for member in self.binMembers[bid]:
                        print bid, self.PM.contigNames[member], self.PM.contigLengths[member]
        elif(outFormat == 'minimal'):
            print "#\"bid\"\t\"cid\"\t\"length\""            
            for bid in self.binMembers:
                if(np.size(self.binMembers[bid]) > 0):
                    for member in self.binMembers[bid]:
                        print bid, self.PM.contigNames[member], self.PM.contigLengths[member]
            pass
        else:
            print "Error: Unrecognised format:", outFormat

    def plotProfileDistributions(self):
        """Plot the coverage and kmer distributions for each bin"""
        for bid in self.bins:
            self.bins[bid].plotProfileDistributions(self.PM.transformedCP, self.PM.kmerSigs, fileName="PROFILE_"+str(bid))

    def plotBins(self, FNPrefix="BIN", sideBySide=False):
        """Make plots of all the bins"""
        for bid in self.bins:
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
        if(sideBySide):
            self.plotSideBySide(self.bins.keys(), tag=FNPrefix)
        else:
            for bid in self.bins:
                self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigColours, fileName=FNPrefix+"_"+str(bid))

    def plotSideBySide(self, bids, fileName="", tag=""):
        """Plot two bins side by side in 3d"""
        fig = plt.figure()
        # we need to work out how to shape the plots
        num_plots = len(bids)
        plot_rows = float(int(np.sqrt(num_plots)))
        plot_cols = np.ceil(float(num_plots)/plot_rows)
        plot_num = 1
        for bid in bids:
            title = self.bins[bid].plotOnAx(fig, plot_rows, plot_cols, plot_num, self.PM.transformedCP, self.PM.contigColours, fileName=fileName)
            plot_num += 1
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
        self.PM = mstore.ProfileManager(dbFileName)   # based on user specified length
        self.BM = BinManager(dbFileName=dbFileName)   # bins
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
        self.PM.loadData(condition="length >= "+str(coreCut))
        self.PM.transformCP()
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
        ax1.scatter(self.PM.transformedCP[:,0], self.PM.transformedCP[:,1], self.PM.transformedCP[:,2], edgecolors=self.PM.contigColours, c=self.PM.contigColours, marker='.')
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
    
    To (perhaps) simplify things think of a "bin" as an row_index into the
    column names array. The ClusterBlob has a list of bins which it can
    update etc...
    """
    def __init__(self, rowIndicies, kmerSigs, id, covtol=3, mertol=2):
        self.id = id
        self.rowIndicies = rowIndicies             # all the indicies belonging to this bin
        self.binSize = self.rowIndicies.shape[0]
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
        
    def isSimilar(self, compBin, stdevs=[2,2], failValue=2, index=0):
        """Check whether two bins are similar
        
        When calling DO NOT set failValue or index. 
        stdevs = [upperStdev,lowerStdev]
        Function returns:
            0 if lower limit is met (stringent)
            1 if only upper limit is met (less stringent)
            2 otherwise (not similar)
        return failValue returns 2 on the first call and 1 on the recursive call
        """
        this_lowers = self.covMeans - stdevs[index] * self.covStdevs
        this_uppers = self.covMeans + stdevs[index] * self.covStdevs
        that_lowers = compBin.covMeans - stdevs[index] * compBin.covStdevs
        that_uppers = compBin.covMeans + stdevs[index] * compBin.covStdevs
        # reciprocial test on x and y co-ords only
        for index in range(0,2):
            if(self.covMeans[index] < that_lowers[index] or self.covMeans[index] > that_uppers[index]):
                return failValue 
            if(compBin.covMeans[index] < this_lowers[index] or compBin.covMeans[index] > this_uppers[index]):
                return failValue
        
        # now test for overlaps in the z dimension (this dimension varys widely and is representitive of gross coverage)
        if((this_uppers[2] > that_lowers[2]) and (this_uppers[2] < that_uppers[2])) or ((that_uppers[2] > this_lowers[2]) and (that_uppers[2] < this_uppers[2])):
            # we know we can return at least 1
            if(index == 1):
                # this is the recursive call, so return 0
                return 0
            elif(stdevs[0] == stdevs[1]):
                # no point testing the same thing twice
                return 0
            else:
                # this is the call for the upper limit, call recursively
                return self.isSimilar(compBin, stdevs=stdevs, failValue=(failValue-1), index=index+1)

        return failValue
    
    def consume(self, transformedCP, kmerSigs, contigLengths, deadBin, verbose=False):
        """Combine the contigs of another bin with this one"""
        # consume all the other bins rowIndicies
        if(verbose):
            print "\tBIN:",deadBin.id,"will be consumed by BIN:",self.id
        self.rowIndicies = np.concatenate([self.rowIndicies, deadBin.rowIndicies])
        self.binSize  = self.rowIndicies.shape[0]
        
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
        self.kValMean = 0
        self.kValStdev = 0
        self.kDistUpperLimit = 0
        
    def makeBinDist(self, transformedCP, kmerSigs):
        """Determine the distribution of the points in this bin
        
        The distribution is largely normal, except at the boundaries.
        """
        self.clearBinDist(kmerSigs)
        if(0 == np.size(self.rowIndicies)):
            return
        
        # Get some data!
        cov_working_array = np.zeros((self.binSize,3))
        mer_working_array = np.zeros((self.binSize,np.size(kmerSigs[0])))
        outer_index = 0
        for row_index in self.rowIndicies:
            cov_working_array[outer_index] = transformedCP[row_index]
            #print transformedCP[row_index]
            mer_working_array[outer_index] = kmerSigs[row_index]
            self.merCentroid += kmerSigs[row_index]
            outer_index += 1
        self.merCentroid /= float(np.size(self.rowIndicies))
        
        # calculate the coverage mean and stdev 
        self.covMeans = np.mean(cov_working_array,axis=0)
        self.covStdevs = np.std(cov_working_array,axis=0)

        # now do the kmerSigs
        # z-normalise each column in the working array
        self.merMeans = np.mean(mer_working_array, axis=0)
        tmpMerStdevs = np.std(mer_working_array, axis=0)
        # no zeros!
        self.merStdevs = np.array([x if x !=0 else 1.0 for x in tmpMerStdevs])
        for index in range(0,np.size(self.rowIndicies)):
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
        
    def getKDist(self, sig, centroid=None):
        """Get the distance of this sig from the centroid"""
        # z-norm and then distance!
        if centroid is None:
            centroid = self.merZeros
        return np.linalg.norm((sig-self.merMeans)/self.merStdevs - centroid)
    
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

    def calcTotalSize(self, contigLengths):
        """Work out the total size of this bin in BP"""
        totalBP = 0
        for row_index in self.rowIndicies:
            totalBP += contigLengths[row_index]
        self.totalBP = totalBP

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
        for row_index in self.rowIndicies:
            file_str.write(contigNames[row_index]+"\t")
        return file_str.getvalue()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

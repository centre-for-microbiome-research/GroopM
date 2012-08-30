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

import networkx as nx

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

class ModeNotAppropriateException(BaseException):
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
        
        bin_members = self.initialiseContainers()
        if(makeBins):
            self.PM.transformCP(silent=silent)
            self.makeBins(bin_members)

    def initialiseContainers(self):
        """Munge the raw data into something more usable"""
        # initialise these containers
        bin_members = {0:[]}
        bin_sizes = {0:0}
        for bid in self.PM.validBinIds.keys():
            bin_sizes[bid] = 0;
            bin_members[bid] = []
        
        # fill them up
        for row_index in range(0, np.size(self.PM.indicies)):
            bin_members[self.PM.binIds[row_index]].append(row_index)
            bin_sizes[self.PM.binIds[row_index]] += self.PM.contigLengths[row_index]

        # we need to get the largest BinId in use
        bids = self.PM.getBinStats().keys()
        if(len(bids) > 0):
            self.nextFreeBinId = np.max(bids)
        
        return bin_members

    def makeBins(self, binMembers):
        """Make bin objects from loaded data"""
        for bid in self.PM.validBinIds.keys():
            self.bins[bid] = Bin(np.array(binMembers[bid]), self.PM.kmerSigs, bid, self.PM.scaleFactor-1)
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


    def removeBinAndIndicies(self, bid):
        """Remove indicies from the PM based on bin identity
        
        "unload" some data
        """
        # get some info
        rem_bin = self.getBin(bid)
        original_length = len(self.PM.indicies)
        rem_list = np.sort(rem_bin.rowIndicies)
        
        # affect the raw data in the PM
        self.PM.reduceIndicies(rem_list)
        del self.PM.validBinIds[bid]
        
        # remove the bin here
        del self.bins[bid]
        
        # now fix all the rowIndicies in all the other bins
        for bid in self.bins:
            self.bins[bid].rowIndicies = self.fixRowIndexLists(original_length, np.sort(self.bins[bid].rowIndicies), rem_list)


    def fixRowIndexLists(self, originalLength, oldList, remList):
        """Fix up row index lists which reference into the
        data structure after a call to reduceIndicies
        
        originalLength is the length of all possible row indicies
        before the removal (ie self.indicies)
        oldList is the old list of row indicies
        remList is the list of indicies to be removed
        
        BOTH OLD AND REM LIST MUST BE SORTED ASCENDING!
        """
        shift_down = 0;
        old_list_index = 0
        new_list = np.array([])
        for i in range(originalLength):
            if(i in remList):
                shift_down+=1
            elif(i in oldList):
                new_list = np.append(new_list, oldList[old_list_index]-shift_down)
                old_list_index += 1
        return new_list

#------------------------------------------------------------------------------
# BIN UTILITIES 

    def getCentroidProfiles(self, mode="mer"):
        """Return an array containing the centroid stats for each bin"""
        if(mode == "mer"):
            ret_vecs = np.zeros((len(self.bins), len(self.PM.kmerSigs[0])))
            outer_index = 0
            for bid in self.bins:
                ret_vecs[outer_index] = self.bins[bid].merMeans
                outer_index += 1
            return ret_vecs
        elif(mode == "cov"):
            ret_vecs = np.zeros((len(self.bins), len(self.PM.transformedCP[0])))
            outer_index = 0
            for bid in self.bins:
                ret_vecs[outer_index] = self.bins[bid].covMeans
                outer_index += 1
            return ret_vecs
        else:
            raise ModeNotAppropriateException("Mode",mode,"unknown")            

    def removeChimeras(self):
        """identify and remove chimeric bins"""
        (kill_list, M_cut) = self.measureBinVariance(makeKillList=True, verbose=True)
        print "    Removing chimeras"
        for bid in kill_list:
            # make these guys here
            if(not self.split(bid, 2, M_cut, auto=True, printInstructions=False)):
                # the bin could not be split, delete the parent
                self.deleteBins([bid], force=True, freeBinnedRowIndicies=True, saveBins=False)
        for bid in self.bins:
            # make these guys here
            self.bins[bid].getKmerColourStats(self.PM.contigColours)

    def split(self, bid, n, mode='kmer', MCut=0.0, auto=False, test=False, saveBins=False, verbose=False, printInstructions=True):
        """split a bin into n parts
        
        if auto == True, then just railroad the split
        if test == True, then test via merging
        if savebins == True, save the split (if you will do it)
        if MCut != 0, carry the split through only if both daughter bins have an M
          less than MCut
        """
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
        try:
            centroids,_ = kmeans(obs,n)
        except ValueError:
            if(verbose):
                print "Error splitting"
            return False
        idx,_ = vq(obs,centroids)

        # build some temp bins        
        idx_sorted = np.argsort(np.array(idx))
        current_group = -1
        bids = [bid]
        bin_stats = {} # bin id to bin size
        bin_stats[bid]=0
        bin_update = {} # row index to bin id
        holding_array = np.array([])
        split_bin = None
        for i in idx_sorted:
            if(idx[i] != current_group):
                if(current_group != -1):
                    # bin is full!
                    split_bin = self.makeNewBin(holding_array)
                    for row_index in holding_array:
                        bin_update[self.PM.indicies[row_index]] = split_bin.id
                    bin_stats[split_bin.id] = split_bin.binSize  
                    split_bin.makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
                    split_bin.calcTotalSize(self.PM.contigLengths)
                    bids.append(split_bin.id)
                    holding_array = np.array([])
                current_group = idx[i]
            holding_array = np.append(holding_array, bin.rowIndicies[i])
        # do the last one
        if(np.size(holding_array) != 0):
            split_bin = self.makeNewBin(holding_array)
            for row_index in holding_array:
                bin_update[self.PM.indicies[row_index]] = split_bin.id  
            bin_stats[split_bin.id] = split_bin.binSize  
            split_bin.makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
            split_bin.calcTotalSize(self.PM.contigLengths)
            bids.append(split_bin.id)

        if(auto and saveBins):
            # charge on through
            self.updateBinStats(bin_stats)
            self.PM.saveBinIds(bin_update)
            return
        
        if(test and saveBins):
            # implementation of autosplitting
            del bids[0]
            bin1 = self.getBin(bids[0])
            bin2 = self.getBin(bids[1])
            if(not self.shouldMerge(bin1, bin2)):
                # should not merge the split, so the split can stay
                print "merge failed"
                self.updateBinStats(bin_stats)
                self.PM.saveBinIds(bin_update)
            else:
                print "merge passed"
            return
        
        if(MCut != 0):
            # see if the split bins have a tighted distribution
            del bids[0]
            bin1 = self.getBin(bids[0])
            bin2 = self.getBin(bids[1])
            (b1_kM, b1_kS, b1_kR)= bin1.getInnerVariance(self.PM.kmerSigs)
            (b2_kM, b2_kS, b2_kR)= bin2.getInnerVariance(self.PM.kmerSigs)
            print "Split", b1_kM,b1_kR,b2_kM,b2_kR
            if((b1_kM < MCut) and (b2_kM < MCut)):
                # ok!, delete the parent bin
                self.deleteBins([bid], force=True, saveBins=False)
                return True
            else:
                # delete the temp bins
                self.deleteBins(bids, force=True, saveBins=False)
                return False

        # we will need to confer with the user
        # plot some stuff
        self.plotSideBySide(bids)
        # remove this query from the list so we don't delete him
        del bids[0]
        user_option = self.promptOnSplit(n)
        if(user_option == 'Y'):
            if(saveBins):
                # save the temp bins
                self.updateBinStats(bin_stats)
                self.PM.saveBinIds(bin_update)
            return
        # delete the temp bins
        self.deleteBins(bids, force=True)
        
        # see what the user wants to do
        if(user_option == 'N'):
            return
        elif(user_option == 'C'):
            if(mode == "cov"):
                print "Already doing split based on coverage profile"
            else:
                self.split(bid, n, mode='cov', auto=auto, saveBins=saveBins, verbose=verbose, printInstructions=False)
        elif(user_option == 'K'):
            if(mode == "kmer"):
                print "Already doing split based on kmer profile"
            else:
                self.split(bid, n, mode='kmer', auto=auto, saveBins=saveBins, verbose=verbose, printInstructions=False)
        elif(user_option == 'P'):
            try:
                parts = int(raw_input("Enter new number of parts:"))
            except ValueError:
                print "You need to enter an integer value!"
                parts = int(raw_input("Enter new number of parts:"))
            self.split(bid, parts, mode=mode, auto=auto, saveBins=saveBins, verbose=verbose, printInstructions=False)   
        
    def condenseWrapper(self,
                      manual=False,          # do we need to ask permission every time?
                      save=False,
                      plotter=False
                      ):
        """Iterative wrapper for the condense function"""
        if(plotter):
            self.plotterCondenseBins()
            return
        if(manual):
            self.printMergeInstructions()
        print "    Calculating preliminary stats"
        (meanMs, stdMs, medianSs, stdSs) = self.measureBinVariance(verbose=False)
        total_num_bins_condensed = 0
        num_bins_condensed = 0
        while True: # do while loop anyone?
            (num_bins_condensed,continue_merging) = self.autoCondenseBins(manual=manual,     
                                                                      save=save,
                                                                      medianVariance=medianSs
                                                                      )
            print "********************************************************************************"
            print "    End of round:",num_bins_condensed,"condensed"
            print "********************************************************************************"
            total_num_bins_condensed += num_bins_condensed 
            if(num_bins_condensed == 0 or continue_merging == False):
                break
        print "    ",total_num_bins_condensed,"bins condensed.",len(self.bins),"bins remain"
                        
    def plotterCondenseBins(self):
        """combine similar bins using 3d plots"""
        self.printCondensePlotterInstructions()
        self.plotBinIds()
        continue_merge = True
        while(continue_merge):
            user_option = self.promptOnPlotterCondense()
            if(user_option == 'R'):
                self.plotBinIds()
            elif(user_option == 'M'):
                merge_bids = self.getPlotterMergeIds()
                if(not len(merge_bids) == 0):
                    self.merge(merge_bids, auto=False, manual=True, newBid=False, saveBins=True, verbose=False, printInstructions=False)
            else:
                return
        return
                    
    def autoCondenseBins(self, 
                      verbose=True,
                      manual=False,          # do we need to ask permission every merge?
                      save=False,            # should we save the merges as we go?
                      medianVariance=-1
                      ):
        """combine similar bins using super auto algorithm"""
        any_merged = False
        merged = {}       # who is getting merged by who?
        merged_order = [] # and in what order to merge!
        num_cores_merged = 0

        # we'd like to know beforehand, whom can merge with who
        # the following function returns a graph of such relationships 
        (merge_allowed, scores) = self.calculateAllowableMergers(self.bins.keys())
        # we should join large cliques first
        all_cliques = []
        for clique in nx.find_cliques(merge_allowed):
            all_cliques.append(clique)

        ordered_cliques = sorted(all_cliques, key=len)[::-1]
        for clique in ordered_cliques:
            cl = len(clique)
            if(cl > 1):
                for i in range(cl):
                    subject_bid = clique[i] # subject may already be merged. It's ok, we'll fix it up below
                    subject_bin = self.getBin(subject_bid)
                    (subject_cM, subject_cS, subject_cR) = subject_bin.getInnerVariance(self.PM.transformedCP, mode="cov")
                    (subject_kM, subject_kS, subject_kR) = subject_bin.getInnerVariance(self.PM.kmerSigs)
                    j = i + 1
                    while(j < cl):
                        query_bid = clique[j]
                        if(query_bid not in merged):                    
                            query_bin = self.getBin(query_bid)
                            # always test against the unchanged subject
                            (continue_merge,m_auto) = self.shouldMerge(subject_bin, query_bin, bin1_cM=subject_cM, bin1_kM=subject_kM, medianVariance=medianVariance)
                            if(continue_merge):
                                # work out if this merge can be done on the fly or if we need someone to look it over
                                auto = True
                                if(manual):
                                    auto = False
                                else:
                                    if(not m_auto):
                                        if(scores[self.makeBidKey(subject_bid,query_bid)] < 2):
                                            auto = False
                                        elif(np.abs(subject_bin.kValMean - query_bin.kValMean) > 0.05):
                                            auto = False
                                
                                # change the subject now if it's been merged
                                if(subject_bid in merged):
                                    subject_bid = merged[subject_bid]   
                                if(subject_bid != query_bid):          # make sure there are no loops 
                                    merged[query_bid] = subject_bid    # subject consumes query, enforce the tree structure
                                    merged_order.append([subject_bid,query_bid,auto])
                        j += 1
        
        # Merge them in the order they were seen in
        for merge_job in merged_order:
            m_bids = [merge_job[0],merge_job[1]]
            merge_status = self.merge(m_bids, auto=merge_job[2], manual=manual, saveBins=save, verbose=verbose, printInstructions=False)
            if(merge_status == 2):
                num_cores_merged += 1
            elif(merge_status == 0):
                return (num_cores_merged,False)
        return (num_cores_merged,True)

    def calculateMergeLimits(self, bid, tolerance):
        """Calculate allowable upper and lower limits"""
        bin = self.getBin(bid)
        return [[bin.covMeans[0] - tolerance * bin.covStdevs[0],
                 bin.covMeans[1] - tolerance * bin.covStdevs[1],
                 bin.covMeans[2] - tolerance * bin.covStdevs[2],
                 bin.kValMean - tolerance * bin.kValStdev],
                [bin.covMeans[0] + tolerance * bin.covStdevs[0],
                 bin.covMeans[1] + tolerance * bin.covStdevs[1],
                 bin.covMeans[2] + tolerance * bin.covStdevs[2],
                 bin.kValMean + tolerance * bin.kValStdev]]

    def calculateAllowableMergers(self, bids, tolerance=3.0):
        """Create a list of allowable but not necessary mergers"""
        merge_allowed = nx.Graph()
        scores = {}
        lower_limits = {}
        upper_limits = {}
        bins = []

        # make sure everyone has stats!
        for bid in bids:
            bin = self.getBin(bid)
            bins.append(bin)
            bin.makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
            bin.getKmerColourStats(self.PM.contigColours)
            [lower_limits[bid],upper_limits[bid]] = self.calculateMergeLimits(bid, tolerance)
            merge_allowed.add_node(bid)
            
        # sort bins according to kmer profile (kValMean)
        bins_sorted = sorted(bins)
        for i in range(len(bids)):
            subject_bin = bins_sorted[i]
            #print "SUB",subject_bin.id,lower_limits[subject_bin.id],upper_limits[subject_bin.id]
            for j in range(i+1,len(bids)):
                query_bin = bins_sorted[j]
                #print "QRY",query_bin.id,query_bin.covMeans[0],query_bin.covMeans[1],query_bin.covMeans[2],query_bin.kValMean
                k_test_1 = (query_bin.kValMean >= lower_limits[subject_bin.id][3] and
                            query_bin.kValMean <= upper_limits[subject_bin.id][3])
                k_test_2 = (subject_bin.kValMean >= lower_limits[query_bin.id][3] and
                            subject_bin.kValMean <= upper_limits[query_bin.id][3])
                still_allowed = False
                score = 0
                if(k_test_1 and k_test_2):
                    still_allowed = True
                    score = 2
                elif(k_test_1 or k_test_2):
                    still_allowed = True
                    score= 1
                if(still_allowed):
                    #print "t1",
                    for k in range(3):
                        still_allowed &= ((query_bin.covMeans[k] >= lower_limits[subject_bin.id][k] and
                                           query_bin.covMeans[k] <= upper_limits[subject_bin.id][k]) or
                                          (subject_bin.covMeans[k] >= lower_limits[query_bin.id][k] and
                                           subject_bin.covMeans[k] <= upper_limits[query_bin.id][k]))
                        #print still_allowed, 
                    if(still_allowed):
                        #print "\nt2"
                        merge_allowed.add_edge(query_bin.id, subject_bin.id)
                        scores[self.makeBidKey(query_bin.id, subject_bin.id)] = score
                    #else:
                    #    print "\nf2"
                else:
                    #print "f1"
                    break

        return (merge_allowed, scores)

    def shouldMerge(self, bin1, bin2, kDistWobble=1.1, cDistWobble=1.1, bin1_cM=0.0, bin1_kM=0.0, medianVariance=-1):
        """Should two bins be merged?"""
        
        # make the bin that would be if it should be
        tmp_bin = self.makeNewBin(np.concatenate([bin1.rowIndicies,bin2.rowIndicies]))
        tmp_bin.makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
        tmp_bin.calcTotalSize(self.PM.contigLengths)
        
        # coverage is cheaper to test first
        (bin2_cM, bin2_cS, bin2_cR) = bin2.getInnerVariance(self.PM.transformedCP, mode="cov")
        if(bin1_cM==0):
            (bin1_cM, bin1_cS, bin1_cR) = bin1.getInnerVariance(self.PM.transformedCP, mode="cov")
        c_cut = cDistWobble * np.max([bin1_cM, bin2_cM])
        (tmp_cM, tmp_cS, tmp_cR) = tmp_bin.getInnerVariance(self.PM.transformedCP, mode="cov")

        #print "------------------"
        #bin1.printBin(self.PM.contigNames, self.PM.contigLengths, outFormat="summary", separator="\t")
        #bin2.printBin(self.PM.contigNames, self.PM.contigLengths, outFormat="summary", separator="\t")

        if(tmp_cM < c_cut):
            # coverage passed, test kmerness
            (bin2_kM, bin2_kS, bin2_kR) = bin2.getInnerVariance(self.PM.kmerSigs)
            if(bin1_kM==0):
                (bin1_kM, bin1_kS, bin1_kR) = bin1.getInnerVariance(self.PM.kmerSigs) 
            k_cut = kDistWobble * np.max([bin1_kM, bin2_kM])
            (tmp_kM, tmp_kS, tmp_kR) = tmp_bin.getInnerVariance(self.PM.kmerSigs)
            if(tmp_kM < k_cut):
                self.deleteBins([tmp_bin.id], force=True)
                if(medianVariance != -1):        # check whether the variance of the new bin is less than the median 
                    if(tmp_kM < medianVariance): # varience observed across all bins
                        return (True,True)
                    else:
                        return (True, False)
                else:
                    #print bin1.id, bin2.id, bin1_cM, bin2_cM, c_cut, tmp_cM, bin1_kM, bin2_kM, k_cut, tmp_kM, True
                    #print "------------------"
                    return True
        self.deleteBins([tmp_bin.id], force=True)
        #print bin1.id, bin2.id, bin1_cM, bin2_cM, c_cut, tmp_cM, False
        #print "------------------"
        return False

    def merge(self, bids, auto=False, manual=False, newBid=False, saveBins=False, verbose=False, printInstructions=True):
        """Merge two or more bins
        
        It's a bit strange to have both manual and auto params
        NOTE: manual ALWAYS overrides auto. In the condensing code, auto is
        set programmaticaly, manual is always set by the user. So we listen
        to manual first
        """
        parent_bin = None

        bin_stats = {}
        if(printInstructions and not auto):
            self.printMergeInstructions()

        if(newBid):
            # we need to make this into a new bin
            parent_bin = makeNewBin()
            # now merge it with the first in the new list
            dead_bin = self.getBin(bids[0])
            parent_bin.consume(self.PM.transformedCP, self.PM.kmerSigs, self.PM.contigLengths, self.PM.contigColours, dead_bin, verbose=verbose)
            self.deleteBins([bids[0]], force=True)
            bin_stats[bids[0]] = 0
            bin_stats[parent_bin.id] = parent_bin.binSize
        else:
            # just use the first given as the parent
            parent_bin = self.getBin(bids[0])
        
        # let this guy consume all the other guys
        ret_val = 0
        for i in range(1,len(bids)):
            continue_merge = False
            dead_bin = self.getBin(bids[i])
            if(auto and not manual):
                ret_val = 2
                continue_merge = True
            else:
                tmp_bin = self.makeNewBin(np.concatenate([parent_bin.rowIndicies,dead_bin.rowIndicies]))
                tmp_bin.makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
                tmp_bin.calcTotalSize(self.PM.contigLengths)
                tmp_bin.getKmerColourStats(self.PM.contigColours)
                self.plotSideBySide([parent_bin.id,dead_bin.id,tmp_bin.id])
                self.deleteBins([tmp_bin.id], force=True)
                user_option = self.promptOnMerge(bids=[parent_bin.id,dead_bin.id]) 
                if(user_option == "N"):
                    print "Merge skipped"
                    ret_val = 1                    
                    continue_merge=False
                elif(user_option == "Q"):
                    print "All mergers skipped"
                    return 0
                else:
                    ret_val = 2
                    continue_merge=True
            if(continue_merge):
                parent_bin.consume(self.PM.transformedCP, self.PM.kmerSigs, self.PM.contigLengths, self.PM.contigColours, dead_bin, verbose=verbose)
                self.deleteBins([bids[i]], force=True)
                bin_stats[bids[i]] = 0
                bin_stats[parent_bin.id] = parent_bin.binSize

        if(saveBins):
            self.updateBinStats(bin_stats)
            self.saveBins(doCores=False, saveBinStats=False)
            
        return ret_val

    def makeBidKey(self, bid1, bid2):
        """Make a unique key from two bids"""
        if(bid1 < bid2):
            return bid1 * 1000000 + bid2
        return bid2 * 1000000 + bid1

    def getBin(self, bid):
        """get a bin or raise an error"""
        if bid in self.bins:
            return self.bins[bid]
        else:
            raise BinNotFoundException("Cannot find: "+str(bid)+" in bins dicts")
            
    def deleteBins(self, bids, force=False, freeBinnedRowIndicies=False, saveBins=False):
        """Purge a bin from our lists"""
        if(not force):
            user_option = self.promptOnDelete(bids)
            if(user_option != 'Y'):
                return False
        bin_stats = {}
        bin_update = {}
        for bid in bids:
            if bid in self.bins:
                if(freeBinnedRowIndicies):
                    for row_index in self.bins[bid].rowIndicies:
                        if row_index in self.PM.binnedRowIndicies:
                            del self.PM.binnedRowIndicies[row_index]
                        else:
                            print bid, row_index, "FUNG"
                        bin_update[self.PM.indicies[row_index]] = 0 
                bin_stats[bid] = 0
                del self.bins[bid]
            else:
                raise BinNotFoundException("Cannot find: "+str(bid)+" in bins dicts")
            
        if(saveBins):
            self.updateBinStats(bin_stats)
            self.PM.saveBinIds(bin_update)
        return True
        
    def makeNewBin(self, rowIndicies=np.array([]), bid=None):
        """Make a new bin and add to the list of existing bins"""
        if bid is None:
            self.nextFreeBinId +=1
            bid = self.nextFreeBinId
        self.bins[bid] = Bin(rowIndicies, self.PM.kmerSigs, bid, self.PM.scaleFactor-1)        
        return self.bins[bid]

#------------------------------------------------------------------------------
# UI 
    
    def printCondensePlotterInstructions(self):
        raw_input( "****************************************************************\n"
                   " CONDENSING INSTRUCTIONS - PLEASE READ CAREFULLY\n"+
                   "****************************************************************\n"
                   " You have chosen to refine in plotter mode. Congratulations!\n"
                   " You will be shown a 3d plot of all the bins, coloured by kmer\n"
                   " profile. Bin Ids in close proximity and similar colour may need\n"
                   " to be merged. Follow the instructions to merge these bins\n"
                   " Good Luck!\n"
                   " Press any key to produce plots...")
        print "****************************************************************"        
    
    def printMergeInstructions(self):
        raw_input( "****************************************************************\n"
                   " MERGING INSTRUCTIONS - PLEASE READ CAREFULLY\n"
                   "****************************************************************\n"
                   " The computer cannot always be trusted to perform bin mergers\n"
                   " automatically, so during merging you may be shown a 3D plot\n"
                   " which should help YOU determine whether or not the bins should\n"
                   " be merged. Look carefully at each plot and then close the plot\n"
                   " to continue with the merging operation.\n"
                   " The image on the far right shows the bins after merging\n"
                   " Press any key to produce plots...")
        print "****************************************************************"        

    def printSplitInstructions(self):
        raw_input( "****************************************************************\n"
                   " SPLITTING INSTRUCTIONS - PLEASE READ CAREFULLY\n"
                   "****************************************************************\n"
                   " The computer cannot always be trusted to perform bin splits\n"
                   " automatically, so during splitting you may be shown a 3D plot\n"
                   " which should help YOU determine whether or not the bin should\n"
                   " be split. Look carefully at each plot and then close the plot\n"
                   " to continue with the splitting operation.\n\n"
                   " Press any key to produce plots...")
        print "****************************************************************"        

    def promptOnPlotterCondense(self, minimal=False):
        """Find out what the user wishes to do next when refining bins"""
        input_not_ok = True
        valid_responses = ['R','M','Q']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" What next? ("+vrs+") : ")
            else:
                option = raw_input(" You have been shown a 3D plot of the bins\n" \
                                   " How do you want to continue?\n" \
                                   " r = replot, m = merge, q = quit\n" \
                                   " What next? ("+vrs+") : ")
            if(option.upper() in valid_responses):
                return option.upper()
            else:
                print "Error, unrecognised choice '"+option+"'"

    def getPlotterMergeIds(self):
        """Prompt the user for ids to be merged and check that it's all good"""
        input_not_ok = True
        ret_bids = []
        while(input_not_ok):
            ret_bids = []
            option = raw_input("Please enter 'space' separated bin Ids or 'q' to quit: ")
            if(option.upper() == 'Q'):
                return []
            bids = option.split(" ")
            for bid in bids:
                try:
                    # check that it's an int
                    i_bid = int(bid)
                    # check that it's in the bins list
                    if(i_bid not in self.bins):
                        print "**Error: bin",bid,"not found"
                        input_not_ok = True
                        break
                    input_not_ok = False
                    ret_bids.append(i_bid)
                except ValueError:
                    print "**Error: invalid value:", bid
                    input_not_ok = True
                    break
        return ret_bids

    def promptOnMerge(self, bids=[], minimal=False):
        """Check that the user is ok with this merge"""
        input_not_ok = True
        valid_responses = ['Y','N','Q']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        bin_str = ""
        if(len(bids) != 0):
            bin_str = ": "+str(bids[0])
            for i in range(1, len(bids)):
                bin_str += " and "+str(bids[i])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" Merge? ("+vrs+") : ")
            else: 
                option = raw_input(" ****WARNING**** About to merge bins"+bin_str+"\n" \
                                   " If you continue you *WILL* overwrite existing bins!\n" \
                                   " You have been shown a 3d plot of the bins to be merged.\n" \
                                   " Continue only if you're sure this is what you want to do!\n" \
                                   " y = yes, n = no, q = no and quit merging\n" \
                                   " Merge? ("+vrs+") : ")
            if(option.upper() in valid_responses):
                print "****************************************************************"
                return option.upper()
            else:
                print "Error, unrecognised choice '"+option.upper()+"'"
                minimal = True

    def promptOnSplit(self, parts, minimal=False):
        """Check that the user is ok with this split"""
        input_not_ok = True
        valid_responses = ['Y','N','C','K','P']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" Split? ("+vrs+") : ")
            else: 
                option = raw_input(" ****WARNING**** About to split bin into "+str(parts)+" parts\n" \
                                   " If you continue you *WILL* overwrite existing bins!\n" \
                                   " You have been shown a 3d plot of the bin after splitting.\n" \
                                   " Continue only if you're sure this is what you want to do!\n" \
                                   " y = yes, n = no, c = redo but use coverage profile,\n" \
                                   " k = redo but use kmer profile, p = choose new number of parts\n" \
                                   " Split? ("+vrs+") : ")
            if(option.upper() in valid_responses):
                print "****************************************************************"
                return option.upper()
            else:
                print "Error, unrecognised choice '"+option.upper()+"'"
                minimal = True

    def promptOnDelete(self, bids, minimal=False):
        """Check that the user is ok with this split"""
        input_not_ok = True
        valid_responses = ['Y','N']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        bids_str = ",".join([str.lower(str(x)) for x in bids])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" Delete? ("+vrs+") : ")
            else: 
                option = raw_input(" ****WARNING**** About to delete bin(s):\n" \
                                   " "+bids_str+"\n" \
                                   " If you continue you *WILL* overwrite existing bins!\n" \
                                   " Continue only if you're sure this is what you want to do!\n" \
                                   " y = yes, n = no\n"\
                                   " Delete? ("+vrs+") : ")
            if(option.upper() in valid_responses):
                print "****************************************************************"
                return option.upper()
            else:
                print "Error, unrecognised choice '"+option.upper()+"'"
                minimal = True

#------------------------------------------------------------------------------
# BIN STATS 

    def measureBinVariance(self, mode='kmer', makeKillList=False, tolerance=1.0, verbose=False):
        """Get the stats on M's across all bins
        
        If specified, will return a list of all bins which
        fall outside of the average M profile
        """
        Ms = {}
        Ss = {}
        Rs = {}
        for bid in self.bins:
            if(mode == 'kmer'):
                (Ms[bid], Ss[bid], Rs[bid]) = self.bins[bid].getInnerVariance(self.PM.kmerSigs)
            elif(mode == 'cov'):
                (Ms[bid], Ss[bid], Rs[bid]) = self.bins[bid].getInnerVariance(self.PM.transformedCP, mode="cov")
        
        # find the mean and stdev 
        if(not makeKillList):
            return (np.mean(np.array(Ms.values())), np.std(np.array(Ms.values())), np.median(np.array(Ss.values())), np.std(np.array(Ss.values())))
        
        else:
            cutoff = np.mean(np.array(Ms.values())) + tolerance * np.std(np.array(Ms.values()))  
            kill_list = []
            for bid in Ms:
                if(Ms[bid] > cutoff):
                    kill_list.append(bid)
            return (kill_list, cutoff)

    def findCoreCentres(self):
        """Find the point representing the centre of each core"""
        print "    Finding bin centers"
        bin_centroid_points = np.zeros((len(self.bins),3))
        bin_centroid_colours = np.zeros((len(self.bins),3))
        # remake the cores and populate the centres
        S = 1       # SAT and VAL remain fixed at 1. Reduce to make
        V = 1       # Pastels if that's your preference...
        outer_index = 0
        for bid in self.bins.keys():
            cum_colour = np.array([])
            for row_index in self.bins[bid].rowIndicies:
                cum_colour = np.append(cum_colour, self.PM.contigColours[row_index])
            cum_colour = np.reshape(cum_colour, (self.bins[bid].binSize, 3))
            ave_colour = np.mean(cum_colour, axis=0)

            bin_centroid_points[outer_index] = self.bins[bid].covMeans
            bin_centroid_colours[outer_index] = ave_colour
            outer_index += 1
            
        return (bin_centroid_points, bin_centroid_colours)

    def analyseBinKVariance(self, outlierTrim=0.1, plot=False):
        """Measure within and between bin variance of kmer sigs
        
        return a list of potentially confounding kmer indicies
        """
        print "    Measuring kmer type variances"        
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
                self.printInner(outFormat)
            except:
                print "Error diverting stout to file:", fileName, sys.exc_info()[0]
                raise
        else:
            self.printInner(outFormat)           

    def printInner(self, outFormat):
        """Print bin information to STDOUT"""
        # handle the headers first
        separator = "\t"
        if(outFormat == 'summary'):
            print separator.join(["#\"bid\"","\"totalBP\"","\"numCons\"","\"kMean\"","\"kStdev\""]) 
        elif(outFormat == 'minimal'):
            print separator.join(["#\"bid\"","\"cid\"","\"length\""])            
        elif(outFormat == 'full'):
            pass
        else:
            print "Error: Unrecognised format:", outFormat
            return

        for bid in self.bins:
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)
            self.bins[bid].getKmerColourStats(self.PM.contigColours)
            self.bins[bid].calcTotalSize(self.PM.contigLengths)
            self.bins[bid].printBin(self.PM.contigNames, self.PM.contigLengths, outFormat=outFormat, separator=separator)

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

    def plotBinIds(self):
        """Render 3d image of core ids"""
        (bin_centroid_points, bin_centroid_colours) = self.findCoreCentres()
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        outer_index = 0
        for bid in self.bins.keys():
            ax.text(bin_centroid_points[outer_index,0], 
                    bin_centroid_points[outer_index,1], 
                    bin_centroid_points[outer_index,2], 
                    str(bid), 
                    color=bin_centroid_colours[outer_index]
                    )
            outer_index += 1
        
        ax.set_xlim3d(0, 1000)
        ax.set_ylim3d(0, 1000)
        ax.set_zlim3d(0, 1000)
        try:
            plt.show()
            plt.close(fig)
        except:
            print "Error showing image", sys.exc_info()[0]
            raise
        del fig

    def plotBinPoints(self):
        """Render the image for validating cores"""
        (bin_centroid_points, bin_centroid_colours) = self.findCoreCentres()
        fig = plt.figure()
        ax2 = fig.add_subplot(111, projection='3d')
        ax2.scatter(bin_centroid_points[:,0], bin_centroid_points[:,1], bin_centroid_points[:,2], edgecolors=bin_centroid_colours, c=bin_centroid_colours)
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

class BinExplorer:
    """Inspect bins, used for validation"""
    def __init__(self, dbFileName, bids=[]):
        self.PM = mstore.ProfileManager(dbFileName)   # based on user specified length
        self.BM = BinManager(dbFileName=dbFileName)   # bins
        if bids is None:
            self.bids = []
        else:
            self.bids = bids

    def plotFlyOver(self, fps=10.0, totalTime=120.0):
        """Plot a flyover of the data with bins being removed"""
        self.loadBins(makeBins=True,silent=True,bids=self.bids)
        all_bids = self.bins.keys()

        # control image form and output
        current_azim = 45.0
        current_elev = 0.0
        current_frame = 0.0
        total_frames = fps * totalTime
        total_azim_shift = 720.0
        total_elev_shift = 360.0
        azim_increment = total_azim_shift / total_frames
        elev_increment = total_elev_shift / total_frames
        
        print "Need",total_frames,"frames:"
        # we need to know when to remove each bin
        bid_remove_rate = total_frames / float(len(all_bids))
        bid_remove_indexer = 1.0
        bid_remove_counter = 0.0
        current_bid_index = 0
        current_bid = all_bids[current_bid_index]
        
        while(current_frame < total_frames):
            print "Frame",int(current_frame)
            file_name = "%04d" % current_frame +".jpg"
            self.PM.renderTransCPData(fileName=file_name,
                                         elev=current_elev,
                                         azim=current_azim,
                                         primaryWidth=6,
                                         dpi=200,
                                         showAxis=True,
                                         format='jpeg'
                                         )
            current_frame += 1
            current_azim += azim_increment
            current_elev += elev_increment

            bid_remove_counter += 1.0
            if(bid_remove_counter >= (bid_remove_rate*bid_remove_indexer)):
                # time to remove a bin!
                self.removeBinAndIndicies(current_bid)
                bid_remove_indexer+=1
                current_bid_index += 1
                if(current_bid_index < len(all_bids)):
                    current_bid = all_bids[current_bid_index]
                else:
                    return

    def plotBinProfiles(self):
        """Plot the distributions of kmer and coverage signatures"""
        self.loadBins(makeBins=True,silent=False,bids=self.bids)
        print "Plotting bin profiles"
        self.plotProfileDistributions()
    
    def plotPoints(self):
        """plot points"""
        self.loadBins(makeBins=True,silent=False,bids=self.bids)
        self.plotBinPoints()
    
    def plotSideBySide(self, coreCut):
        """Plot cores side by side with their contigs"""
        self.PM.loadData(condition="length >= "+str(coreCut))
        self.PM.transformCP()
        self.loadBins(makeBins=True,bids=self.bids)
        print "Creating side by side plots"
        (bin_centroid_points, bin_centroid_colours) = self.findCoreCentres()
        self.analyseBinKVariance()
        self.plotCoresVsContigs(bin_centroid_points, bin_centroid_colours)

    def plotIds(self):
        """Make a 3d plot of the bins but use IDs instead of points
        
        This function will help users know which bins to merge
        """
        self.loadBins(makeBins=True,silent=False,bids=self.bids)
        self.plotBinIds()

    def plotUnbinned(self, coreCut):
        """Plot all contigs over a certain length which are unbinned"""
        self.PM.loadData(condition="((length >= "+str(coreCut)+") & (bid == 0))")
        self.PM.transformCP()
        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')
        ax1.scatter(self.PM.transformedCP[:,0], self.PM.transformedCP[:,1], self.PM.transformedCP[:,2], edgecolors=self.PM.contigColours, c=self.PM.contigColours, marker='.')
        try:
            plt.show()
            plt.close(fig)
        except:
            print "Error showing image", sys.exc_info()[0]
            raise
        del fig
            
#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def plotCoresVsContigs(self, binCentroidPoints, binCentroidColours):
        """Render the image for validating cores"""
        fig = plt.figure()
        ax1 = fig.add_subplot(121, projection='3d')
        ax1.scatter(self.PM.transformedCP[:,0], self.PM.transformedCP[:,1], self.PM.transformedCP[:,2], edgecolors=self.PM.contigColours, c=self.PM.contigColours, marker='.')
        ax2 = fig.add_subplot(122, projection='3d')
        ax2.scatter(binCentroidPoints[:,0], binCentroidPoints[:,1], binCentroidPoints[:,2], edgecolors=binCentroidColours, c=binCentroidColours)
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

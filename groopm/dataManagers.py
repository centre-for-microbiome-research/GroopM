#!/usr/bin/env python
###############################################################################
#                                                                             #
#    dataManagers.py                                                          #
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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import sys
import math
import random
import os
import string

import colorsys
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

import numpy as np
import scipy.ndimage as ndi
import scipy.spatial.distance as ssdist
from scipy.stats import f_oneway, distributions

import tables

# GroopM imports
import PCA
import mstore
import bin
import groopmExceptions as ge
import som

np.seterr(all='raise')     

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BinManager:
    """Class used for manipulating bins"""
    def __init__(self, dbFileName="", pm=None):
        # data storage
        if(dbFileName != ""):
            self.PM = ProfileManager(dbFileName)
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
                 loadCovProfiles=True,
                 min=None,
                 max=None):
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
            self.PM.transformCP(silent=silent, min=min, max=max)
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
            self.bins[bid] = bin.Bin(np.array(binMembers[bid]), bid, self.PM.scaleFactor-1)
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)      

    def saveBins(self, doCores=True, saveBinStats=True, updateBinStats=True, unbinned={}):
        """Save binning results"""
        c2b_update = {}
        core_update = {}
        if doCores:
            (c2b_update, core_update) = self.getCoreBinUpdates()
            self.PM.saveCores(core_update)
        else:
            c2b_update = self.getBinUpdates(unbinned)
            
        self.PM.saveBinIds(c2b_update)
        self.PM.setClustered()
        
        if saveBinStats:
            self.saveBinStats()
        elif updateBinStats:
            # we need to make a list of updates
            updates = {}
            for bid in self.getBids():
                updates[bid] = self.bins[bid].binSize
            self.updateBinStats(updates)
    
    def saveBinStats(self):
        """Update / overwrite the table holding the bin stats
        
        Note that this call effectively nukes the existing table
        and should only be used during initial coring. BID must be somewhere!
        """
        bin_updates = {}
        for bid in self.getBids():
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

    def getBinUpdates(self, c2b={}):
        """Merge the bids, raw DB indexes and core information so we can save to disk"""
        # we need a mapping from cid (or local index) to binID
        bin_update = {}
        for row_index in c2b: # load any now-unbinned critters
            bin_update[self.PM.indicies[row_index]] = c2b[row_index]
            
        for bid in self.getBids():
            for row_index in self.bins[bid].rowIndicies:
                bin_update[self.PM.indicies[row_index]] = bid
            
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
        for bid in self.getBids():
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
# BIN REFINEMENT 

    def condenseWrapper(self,
                      manual=False,          # do we need to ask permission every time?
                      save=False,
                      plotter=False,
                      chimeraCheck=False,
                      outlierCheck=False
                      ):
        """Iterative wrapper for the condense function"""
        if(plotter):
            self.plotterCondenseBins()
        elif(chimeraCheck):
            self.identifyChimeras()
        elif(outlierCheck):
            self.removeOutliers(saveBins=True, remove=True)
                        
    def plotterCondenseBins(self):
        """combine similar bins using 3d plots"""
        self.printCondensePlotterInstructions()
        self.plotBinIds()
        continue_merge = True
        while(continue_merge):
            user_option = self.promptOnPlotterCondense()
            if(user_option == 'R'):
                self.plotBinIds()
            
            elif(user_option == 'P'):
                self.plotBinPoints()
            
            elif(user_option == 'M'):
                # merge bins
                merge_bids = self.getPlotterMergeIds()
                if(not len(merge_bids) == 0):
                    self.merge(merge_bids, auto=False, manual=True, newBid=False, saveBins=True, verbose=False, printInstructions=False)
            
            elif(user_option == 'K'):
                # display a subset only!
                have_range = False
                krange=0
                while(not have_range):
                    try:
                        krange = int(raw_input(" Enter kmer range (0-9):"))
                        have_range = True
                    except ValueError:
                        print "You need to enter an integer value!"
                self.plotBinIds(krange=krange)
                
            elif(user_option == 'B'):
                # print single bin
                have_bid = False
                while(not have_bid):
                    try:
                        bid = int(raw_input(" Enter bid to plot:"))
                        if bid not in self.bins:
                            print "ERROR: Bin %d not found!" % bid
                        else:
                            have_bid = True
                    except ValueError:
                        print "You need to enter an integer value!"
                self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals)
                
            elif(user_option == 'S'):
                # split bins
                have_bid = False
                have_parts = False
                while(not have_bid):
                    try:
                        bid = int(raw_input(" Enter bid to split:"))
                        if bid not in self.bins:
                            print "ERROR: Bin %d not found!" % bid
                        else:
                            have_bid = True
                    except ValueError:
                        print "You need to enter an integer value!"
                while(not have_parts):
                    try:
                        parts = int(raw_input(" Enter number of parts to split into:"))
                        if(parts < 2):
                            print "ERROR: Need to choose 2 or more parts"
                        else:
                            have_parts = True
                    except ValueError:
                        print "You need to enter an integer value!"
                self.split(bid, parts, mode='kmer', auto=False, saveBins=True, printInstructions=False)   
            else:
                return
        return

    def autoCondenseBins(self, searchSpan):
        """Automagically condense bins"""
        print "Condense cores [begin: %d]" % len(self.bins)
        all_dists = self.findAllNeighbours(searchSpan)
        condense_map = {}   # keep track of where bin ids go to
        for key in all_dists.keys():
            bin1_id = key[0]
            bin2_id = key[1]
            if(bin1_id in condense_map):
                bin1_id = condense_map[bin1_id]
            if(bin2_id in condense_map):
                bin2_id = condense_map[bin2_id]

            if(bin1_id != bin2_id):
                if(self.shouldMerge(self.getBin(bin1_id), self.getBin(bin2_id))):
                    # merge!
                    self.merge([bin1_id,bin2_id], auto=True, manual=False, newBid=False, saveBins=False, verbose=False, printInstructions=False)
                    condense_map[bin2_id] = bin1_id # add this one in so we know how to map it next time
                    for key in condense_map:
                        if condense_map[key] == bin2_id:
                            # bin2 consumed this key, so we'll need to update it
                            condense_map[key] = bin1_id

    def removeOutliers(self, saveBins=False, remove=False):
        """Identify and remove outlying contigs in each bin"""
        print "Removing outliers"
        bins_affected = 0
        contigs_removed = 0
        updates = {}
        bids = self.getBids()
        for bid in bids:
            bin = self.getBin(bid)
            outliers = bin.identifyOutliers(self.PM.averageCoverages, self.PM.kmerVals)
            original_size = bin.binSize
            lo = len(outliers)
            if(remove and lo != 0):
                # we delete these guys!
                contigs_removed += lo
                bins_affected += 1
                nri = np.array([i for i in bin.rowIndicies if i not in outliers])
                bin.rowIndicies = nri
                bin.binSize = original_size - lo
                for i in outliers:
                    updates[i] = 0  # this will set the bid to 0 when we save
                    del self.PM.binnedRowIndicies[i]
        print "    Removed: %d contigs across %d bins" % (contigs_removed, bins_affected)
        print "    Saving updated bins"
        # save what we did
        self.saveBins(doCores=False, saveBinStats=False, unbinned=updates)

    def identifyChimeras(self, remove=False):
        """identify and remove chimeric bins"""
        bids = self.getBids()
        for bid in bids:
            self.isChimera(bid, 'kmer')
            self.isChimera(bid, 'cov')

    def isChimera(self, bid, mode='kmer', saveBins=False):
        """Is this bin chimeric?
        
        Change mode to test for coverage of kmer
        """
        (bin_stats, bin_update, bids) = self.getSplitties(bid, 2, mode)
        del bids[0]
        bin1 = self.getBin(bids[0])
        bin2 = self.getBin(bids[1])
        bins_saved = False
        if mode == 'kmer':
            ik = False
            ic = True
        else:
            ik = True
            ic = False
        if(self.shouldMerge(bin1, bin2, confidence=0.85, ignoreCov=ic, ignoreMer=ik, verbose=True)):
            # split means nothing
            print "OK: %d, %s" % (bid, mode)
        else:
            # should not merge the split, so the split can stay
            print "Chimera: %d, %s" % (bid, mode)
            if(saveBins):
                bins_saved = True
                self.deleteBins([bids[0]], force=True)  # delete the combined bin
                self.updateBinStats(bin_stats)
                self.PM.saveBinIds(bin_update)
        
        # did we change stuff?
        if(not bins_saved):
            # If we're here than we don't need the temp bins        
            # remove this query from the list so we don't delete him
            del bids[0]
            self.deleteBins(bids, force=True)
            
        return

#------------------------------------------------------------------------------
# BIN UTILITIES 

    def getBids(self):
        """Return a sorted list of bin ids"""
        return sorted(self.bins.keys())

    def getCentroidProfiles(self, mode="mer"):
        """Return an array containing the centroid stats for each bin"""
        if(mode == "mer"):
            ret_vecs = np.zeros((len(self.bins)))
            outer_index = 0
            for bid in self.getBids():
                ret_vecs[outer_index] = self.bins[bid].kValMean
                outer_index += 1
            return ret_vecs
        elif(mode == "cov"):
            ret_vecs = np.zeros((len(self.bins), len(self.PM.transformedCP[0])))
            outer_index = 0
            for bid in self.getBids():
                ret_vecs[outer_index] = self.bins[bid].covMeans
                outer_index += 1
            return ret_vecs
        else:
            raise ge.ModeNotAppropriateException("Mode",mode,"unknown")            

    def split(self, bid, n, mode='kmer', auto=False, saveBins=False, verbose=False, printInstructions=True):
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

        # make some split bins
        (bin_stats, bin_update, bids) = self.getSplitties(bid, n, mode)
        
        if(auto and saveBins):
            # charge on through
            self.deleteBins([bids[0]], force=True)  # delete the combined bin
            self.updateBinStats(bin_stats)
            self.PM.saveBinIds(bin_update)
            return

        # we will need to confer with the user
        # plot some stuff
        self.plotSideBySide(bids)
        user_option = self.promptOnSplit(n,mode)
        if(user_option == 'Y'):
            if(saveBins):
                # save the temp bins
                self.deleteBins([bids[0]], force=True)  # delete the combined bin
                self.updateBinStats(bin_stats)
                self.PM.saveBinIds(bin_update)
            return

        # If we're here than we don't need the temp bins        
        # remove this query from the list so we don't delete him
        del bids[0]
        self.deleteBins(bids, force=True)
        
        # see what the user wants to do
        if(user_option == 'N'):
            return
        elif(user_option == 'C'):
            self.split(bid, n, mode='cov', auto=auto, saveBins=saveBins, verbose=verbose, printInstructions=False)
        elif(user_option == 'K'):
            self.split(bid, n, mode='kmer', auto=auto, saveBins=saveBins, verbose=verbose, printInstructions=False)
        elif(user_option == 'P'):
            not_got_parts = True
            parts = 0
            while(not_got_parts):
                try:
                    parts = int(raw_input("Enter new number of parts:"))
                except ValueError:
                    print "You need to enter an integer value!"
                    parts = 0
                if(1 == parts):
                    print "Don't be a silly sausage!"
                elif(0 != parts):
                    not_got_parts = False
                    self.split(bid, parts, mode=mode, auto=auto, saveBins=saveBins, verbose=verbose, printInstructions=False)   

    def getSplitties(self, bid, n, mode):
        """Return a set of split bins"""
        bin = self.getBin(bid)
        obs = np.array([])
        if(mode=='kmer'):
            obs = np.array([self.PM.kmerVals[i] for i in bin.rowIndicies])
        elif(mode=='cov'):
            obs = np.array([self.PM.covProfiles[i] for i in bin.rowIndicies])
        
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
        # this way we can show the user what the split will look like       
        idx_sorted = np.argsort(np.array(idx))
        current_group = 0
        bids = [bid]
        bin_stats = {} # bin id to bin size
        bin_stats[bid]=0 # this will ensure that the old bin id will be deleted!
        bin_update = {} # row index to bin id
        holding_array = np.array([])
        split_bin = None
        for i in idx_sorted:
            if(idx[i] != current_group):
                # bin is full!
                split_bin = self.makeNewBin(holding_array)
                for row_index in holding_array:
                    bin_update[self.PM.indicies[row_index]] = split_bin.id
                bin_stats[split_bin.id] = split_bin.binSize  
                split_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
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
            split_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
            bids.append(split_bin.id)

        return (bin_stats, bin_update, bids)
            
    def findAllNeighbours(self, maxDist):
        """Create a lookup of all the bins and their nearest neighbours"""
        # first we make three sorted lists for X, Y and Z
        # this will save us an all Vs all comparison
        bids = self.getBids()
        max_index = len(bids)
        Xs = []
        Ys = []
        Zs = []
        for bid in bids:
            CM = self.getBin(bid).covMeans
            Xs.append(CM[0])
            Ys.append(CM[1])
            Zs.append(CM[2])
        sorted_Xs = np.argsort(Xs)
        sorted_Ys = np.argsort(Ys)
        sorted_Zs = np.argsort(Zs)

        x_dists = {}
        for i in range(max_index):
            for j in range(i + 1, max_index):
                dist = np.abs(Xs[sorted_Xs[i]] - Xs[sorted_Xs[j]])
                if(dist <= maxDist):
                    key = self.makeBidKey(bids[sorted_Xs[i]], bids[sorted_Xs[j]])
                    x_dists[key] = dist
                else:
                    break 
        y_dists = {}
        for i in range(max_index):
            for j in range(i + 1, max_index):
                dist = np.abs(Ys[sorted_Ys[i]] - Ys[sorted_Ys[j]])
                if(dist <= maxDist):
                    key = self.makeBidKey(bids[sorted_Ys[i]], bids[sorted_Ys[j]])
                    if(key in x_dists):
                        y_dists[key] = dist
                else:
                    break 
        actual_dists = {}
        for i in range(max_index):
            for j in range(i + 1, max_index):
                dist = np.abs(Zs[sorted_Zs[i]] - Zs[sorted_Zs[j]])
                if(dist <= maxDist):
                    key = self.makeBidKey(bids[sorted_Zs[i]], bids[sorted_Zs[j]])
                    if(key in y_dists):
                        dist = np.sqrt(y_dists[key]**2 + x_dists[key]**2 + dist**2)
                        if(dist <= maxDist):
                            actual_dists[key] = dist
                else:
                    break 

        return actual_dists

    def shouldMerge(self, bin1, bin2, ignoreCov=False, ignoreMer=False, confidence=0.95, verbose=False):
        """Determine whether its wise to merge two bins
        
        Perfoms a one-way anova to determine if the larger bin would be 
        significantly changed if it consumed the smaller
        """
        if(bin1.id != bin2.id):
            if not ignoreCov: # work out coverage distributions
                b1_c_dist = bin1.getAverageCoverageDist(self.PM.averageCoverages) 
                b2_c_dist = bin2.getAverageCoverageDist(self.PM.averageCoverages)
                c_dist_1 = b1_c_dist
                if(bin1.binSize < bin2.binSize):
                    c_dist_1 = b2_c_dist
                c_dist_2 = np.append(b2_c_dist, b1_c_dist)
                cov_match = self.isSameVariance(c_dist_1, c_dist_2, confidence=confidence, verbose=verbose) 
            else:
                cov_match = True
                
            if not ignoreMer: # work out kmer distributions   
                b1_k_dist = bin1.getkmerValDist(self.PM.kmerVals) 
                b2_k_dist = bin2.getkmerValDist(self.PM.kmerVals)
                k_dist_1 = b1_k_dist
                if(bin1.binSize < bin2.binSize):
                    k_dist_1 = b2_k_dist
                k_dist_2 = np.append(b2_k_dist, b1_k_dist)
                mer_match = self.isSameVariance(k_dist_1, k_dist_2, confidence=confidence, verbose=verbose)
            else:
                mer_match = True
                
            return cov_match and mer_match
         
        return False

    def isSameVariance(self, dist1, dist2, confidence=0.95, verbose=False):
        """Test to see if the kmerValues for two bins are the same"""
        F_cutoff =  distributions.f.ppf(confidence, 2, len(dist1)+len(dist2)-2)
        F_value = f_oneway(dist1,dist2)[0]
        if verbose:
           print "[V: %f, C: %f]" % (F_value, F_cutoff)
        return F_value < F_cutoff

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
            parent_bin.consume(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths, dead_bin, verbose=verbose)
            self.deleteBins([bids[0]], force=True)
            bin_stats[bids[0]] = 0
            bin_stats[parent_bin.id] = parent_bin.binSize
        else:
            # just use the first given as the parent
            parent_bin = self.getBin(bids[0])
        
        # let this guy consume all the other guys
        ret_val = 0
        some_merged = False
        for i in range(1,len(bids)):
            continue_merge = False
            dead_bin = self.getBin(bids[i])
            if(auto and not manual):
                ret_val = 2
                continue_merge = True
            else:
                tmp_bin = self.makeNewBin(np.concatenate([parent_bin.rowIndicies,dead_bin.rowIndicies]))
                tmp_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
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
                parent_bin.consume(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths, dead_bin, verbose=verbose)
                self.deleteBins([bids[i]], force=True)
                bin_stats[bids[i]] = 0
                bin_stats[parent_bin.id] = parent_bin.binSize
                some_merged = True

        if(saveBins and some_merged):
            self.updateBinStats(bin_stats)
            self.saveBins(doCores=False, saveBinStats=False)
            
        return ret_val

    def makeBidKey(self, bid1, bid2):
        """Make a unique key from two bids"""
        if(bid1 < bid2):
            return (bid1, bid2)
        return (bid2, bid1)

    def getBin(self, bid):
        """get a bin or raise an error"""
        if bid in self.bins:
            return self.bins[bid]
        else:
            raise ge.BinNotFoundException("Cannot find: "+str(bid)+" in bins dicts")
            
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
                raise ge.BinNotFoundException("Cannot find: "+str(bid)+" in bins dicts")
            
        if(saveBins):
            self.updateBinStats(bin_stats)
            self.PM.saveBinIds(bin_update)
        return True
        
    def makeNewBin(self, rowIndicies=np.array([]), bid=None):
        """Make a new bin and add to the list of existing bins"""
        if bid is None:
            self.nextFreeBinId +=1
            bid = self.nextFreeBinId
        self.bins[bid] = bin.Bin(rowIndicies, bid, self.PM.scaleFactor-1)        
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
        valid_responses = ['R','P','B','M','S','K','Q']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" What next? ("+vrs+") : ")
            else:
                option = raw_input(" How do you want to continue?\n" \
                                   " r = replot ids, p = replot points, b = plot single bin," \
                                   " m = merge, s = split, k = set kmer range, q = quit\n" \
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

    def promptOnSplit(self, parts, mode, minimal=False):
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
                if(option.upper() == 'K' and mode.upper() == 'KMER' or option.upper() == 'C' and mode.upper() == 'COV'):
                    print "Error, you are already using that profile to split!"
                    minimal=True
                else:
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

    def classify(self, rowIndex, bids):
        """Classify a contig based on similarity to a set of bins"""
        min_score = 100000000
        info = [rowIndex, len(bids)]
        classification = 0
        for bid in sorted(bids):
            (score, scores) = self.scoreContig(rowIndex, bid)
            info.append((bid, score, scores))
            if(score < min_score):
                classification = bid
                min_score = score
        return (classification, "".join(str(info)))

    def scoreContig(self, rowIndex, bid):
        """Determine how well a particular contig fits with a bin"""
        return self.getBin(bid).scoreProfile(self.PM.kmerVals[rowIndex], self.PM.transformedCP[rowIndex])

    def measureBinVariance(self, mode='kmer', makeKillList=False, tolerance=1.0, verbose=False):
        """Get the stats on M's across all bins
        
        If specified, will return a list of all bins which
        fall outside of the average M profile
        """
        Ms = {}
        Ss = {}
        Rs = {}
        for bid in self.getBids():
            if(mode == 'kmer'):
                (Ms[bid], Ss[bid], Rs[bid]) = self.bins[bid].getInnerVariance(self.PM.kmerVals)
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

    def findCoreCentres(self, krange=None):
        """Find the point representing the centre of each core"""
        print "    Finding bin centers"
        bin_centroid_points = np.array([])
        bin_centroid_colours = np.array([])
        bids = np.array([])
        k_low = 0.0
        k_high = 0.0
        if krange is not None:
            # we only want to plot a subset of these guys
            k_low = float((krange - 1.5)/10.0)
            k_high = float((krange + 1.5)/10.0)
        num_added = 0
        for bid in self.getBids():
            add_bin = True
            if krange is not None:
                ave_kval = np.mean([self.PM.kmerVals[row_index] for row_index in self.bins[bid].rowIndicies])
                if ave_kval < k_low or ave_kval > k_high:
                    add_bin = False
            if add_bin:
                bin_centroid_points = np.append(bin_centroid_points,
                                                self.bins[bid].covMeans)
                bin_centroid_colours = np.append(bin_centroid_colours, 
                                                 np.mean([
                                                          self.PM.contigColours[row_index] for row_index in 
                                                          self.bins[bid].rowIndicies
                                                          ],
                                                         axis=0)
                                                 )
                bids = np.append(bids, bid)
                num_added += 1
        
        if num_added != 0:
            bin_centroid_points = np.reshape(bin_centroid_points, (num_added, 3))
            bin_centroid_colours = np.reshape(bin_centroid_colours, (num_added, 3))
        
        return (bin_centroid_points, bin_centroid_colours, bids)

    def analyseBinKVariance(self, outlierTrim=0.1, plot=False):
        """Measure within and between bin variance of kmer sigs
        
        return a list of potentially confounding kmer indicies
        """
        print "    Measuring kmer type variances"        
        means = np.array([])
        stdevs = np.array([])
        bids = np.array([])
        
        # work out the mean and stdev for the kmer sigs for each bin
        for bid in self.getBids():
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

    def makeCentroidPalette(self):
        """Return a hash of bin ids to colours"""
        (bin_centroid_points, bin_centroid_colours, bids) = self.findCoreCentres()
        pal = {}
        for i in range(len(bids)):
            pal[bids[i]] = (int(bin_centroid_colours[i][0]*255),
                            int(bin_centroid_colours[i][1]*255),
                            int(bin_centroid_colours[i][2]*255)
                           )
        return pal

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

        for bid in self.getBids():
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
            self.bins[bid].printBin(self.PM.contigNames, self.PM.contigLengths, outFormat=outFormat, separator=separator)

    def plotProfileDistributions(self):
        """Plot the coverage and kmer distributions for each bin"""
        for bid in self.getBids():
            self.bins[bid].plotProfileDistributions(self.PM.transformedCP, self.PM.kmerSigs, fileName="PROFILE_"+str(bid))

    def plotBins(self, FNPrefix="BIN", sideBySide=False):
        """Make plots of all the bins"""
        for bid in self.getBids():
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
        if(sideBySide):
            self.plotSideBySide(self.bins.keys(), tag=FNPrefix)
        else:
            for bid in self.getBids():
                self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, fileName=FNPrefix+"_"+str(bid))

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

    def plotBinIds(self, krange=None):
        """Render 3d image of core ids"""
        (bin_centroid_points, bin_centroid_colours, bids) = self.findCoreCentres(krange=krange)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        outer_index = 0
        for bid in bids:
            ax.text(bin_centroid_points[outer_index,0], 
                    bin_centroid_points[outer_index,1], 
                    bin_centroid_points[outer_index,2], 
                    str(int(bid)), 
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
        (bin_centroid_points, bin_centroid_colours, bids) = self.findCoreCentres()
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
class SOMManager:
    """Manage multiple SOMs"""
    def __init__(self,
                 binManager,
                 numSoms=3,
                 somSide=0,
                 somIterations=1000,
                 makeBins=True,
                 load=False
                 ):
        # raw data storage
        self.BM = binManager
        self.BM.loadBins(makeBins=makeBins,silent=False)
        self.PM = self.BM.PM
        self.DM = self.PM.dataManager

        # pointers to various torus maps
        self.covSoms = {}
        self.merSoms = {}

        # normalisation / training vectors
        self.cVecs = None
        self.cMeans = None
        self.cStdevs = None
        self.cMins = None
        self.cMaxs = None

        self.kVecs = None
        self.kMeans = None
        self.kStdevs = None
        self.kMins = None
        self.kMaxs = None

        # used for automerged bins
        self.collapsedMappings = {}

        # misc
        self.numSoms = 3
        self.somIterations = somIterations
        self.covDim = 0
        self.merDim = 0
        
        if(load):
            # load metadata
            meta = self.DM.getSOMMetaFields(self.PM.dbFileName)
            self.somSide = meta['side']
            self.covDim = meta['covDimension']
            self.merDim = meta['merDimension']            
            
            # load the actual data
            self.loadSoms(self.DM.getSOMDataInfo(self.PM.dbFileName), meta=meta)
            
        else:
            if(somSide == 0):
                self.somSide = somSide
            else:
                self.somSide = int(np.sqrt(75*len(self.BM.bins)))
            self.covDim = len(self.PM.transformedCP[0])
            self.merDim = len(self.PM.kmerSigs[0])
            
#------------------------------------------------------------------------------
# SAVING LOADING

    def loadSoms(self, idsInUse, meta=None):
        """load a bunch of SOM data in one go"""
        print "Loading saved SOM data"
        if(meta is None):
            meta = getSOMMetaFields(self.PM.dbFileName)
            self.covDim = meta['covDimension']
            self.merDim = meta['merDimension']
        for flavour in ["mer","cov"]:
            for type in ["weights","regions"]:
                for index in idsInUse[type][flavour]:
                    self.loadSomData(index, type=type, flavour=flavour, meta=meta)

    def loadSomData(self, index, type="weights", flavour="mer", meta=None):
        """Load a saved SOM"""
        if(meta is None):
            meta = self.DM.getSOMMetaFields(self.PM.dbFileName)
            self.covDim = meta['covDimension']
            self.merDim = meta['merDimension']

        print "    Loading",flavour,type,index
        data = self.DM.getSOMData(self.PM.dbFileName, index, type, flavour)
        map = None
        if(flavour == "mer"):
            if(index not in self.merSoms):  # make the map if it's not already in the hash
                map = som.SOM(self.somSide,self.merDim)
                self.merSoms[index] = map
            map = self.merSoms[index] 
        elif(flavour == "cov"):
            if(index not in self.covSoms):
                map = som.SOM(self.somSide,self.covDim)
                self.covSoms[index] = map
            map = self.covSoms[index]
            
        else:
            raise ge.SOMFlavourException("Unknown SOM flavour: "+flavour)
        
        if(type=="weights"):    # make sure we store the data in the right place
            map.loadWeights(data)
        elif(type=="regions"):
            map.loadRegions(data)
        else:
            raise ge.SOMTypeException("Unknown SOM type: "+type)      

    def promptOnOverwrite(self, minimal=False):
        """Check that the user is ok with overwriting the db"""
        input_not_ok = True
        valid_responses = ['Y','N']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" Overwrite? ("+vrs+") : ")
            else: 
                
                option = raw_input(" ****WARNING**** SOMS for database: '"+self.PM.dbFileName+"' exist.\n" \
                                   " If you continue you *WILL* delete any previous matricies!\n" \
                                   " Overwrite? ("+vrs+") : ")
            if(option.upper() in valid_responses):
                print "****************************************************************"
                return option.upper()
            else:
                print "Error, unrecognised choice '"+option.upper()+"'"
                minimal = True

    def saveCovWeights(self, index):
        """Save some coverage weights"""
        self.PM.dataManager.updateSOMTables(self.PM.dbFileName,
                                            self.somSide,
                                            self.covDim,
                                            self.merDim,
                                            covWeights={index:self.covSoms[index].getWeights()})
        
    def saveMerWeights(self, index):
        """Save some kmer weights"""
        self.PM.dataManager.updateSOMTables(self.PM.dbFileName,
                                            self.somSide,
                                            self.covDim,
                                            self.merDim,
                                            merWeights={index:self.merSoms[index].getWeights()})
        
    def saveCovRegions(self, index):
        """Save some coverage regions"""
        self.PM.dataManager.updateSOMTables(self.PM.dbFileName,
                                            self.somSide,
                                            self.covDim,
                                            self.merDim,
                                            covRegions={index:self.covSoms[index].getRegions()})
        
    def saveMerRegions(self, index):
        """Save some coverage regions"""
        self.PM.dataManager.updateSOMTables(self.PM.dbFileName,
                                            self.somSide,
                                            self.covDim,
                                            self.merDim,
                                            merRegions={index:self.merSoms[index].getRegions()})

    def loadTrainingVectors(self):
        """Load and whiten training vectors"""
        (self.cVecs, self.cMeans, self.cStdevs, self.cMins, self.cMaxs) = self.whiten(self.BM.getCentroidProfiles(mode="cov"))
        (self.kVecs, self.kMeans, self.kStdevs, self.kMins, self.kMaxs) = self.whiten(self.BM.getCentroidProfiles(mode="mer"))

#------------------------------------------------------------------------------
# PIPELINING
    
    def DoSOMPipeline(self, merge=True, force=False, tag=""):
        """Wrap the various tasks needed to produce SOMs"""
        if(not self.buildSomWeights(force=force)):
            return
        self.regionalise(force=True)
        self.findRegionNeighbours(merge=merge)
        self.validateRegions()
        if(tag != ""):
            self.renderWeights(tag)
            self.renderRegions(tag)


#------------------------------------------------------------------------------
# CLASSIFICATION

    def remapCollapsed(self, bid):
        """Replace this bid with it's match in the collaped mappings"""
        while(bid in self.collapsedMappings):
            bid = self.collapsedMappings[bid]
        return bid

    def classify(self, rowIndex):
        """Classify a contig (rowIndex) against the SOMS et al
        
        If the soms retuen a clear majority for both profiles 
        then we just go with that. If there are any problems at all
        we check the bin stats of all the regions hit by the soms and
        all of the neighbouring regions. 
        """
        # mers first!
        mc_bids = [{}, {}]
        choices = [0,0]
        whiteVectors = [self.whitenKVector(self.PM.kmerSigs[rowIndex]),
                        self.whitenCVector(self.PM.transformedCP[rowIndex])]
        soms = (self.merSoms, self.covSoms)
        for i in range(2):
            for j in soms[i]:
                tmp_bid = self.remapCollapsed(soms[i][j].classify(whiteVectors[i]))
                if(tmp_bid not in mc_bids[i]):
                    mc_bids[i][tmp_bid] = 1
                else:
                    mc_bids[i][tmp_bid] += 1
            for bid in mc_bids[i]:
                if(mc_bids[i][bid] > 1):    # choose based on consensus of mer classifications (assumes 3 soms)
                    choices[i] = self.remapCollapsed(bid)
                    break

        # if they agree, we are done        
        if(choices[0] == choices[1] and choices[0] != 0):
            return (choices[0], "=")
        
        # m_choice and c_choice disagree
        # we will need to make
        n_query_bids =  []
        n_query_bids.extend([self.remapCollapsed(i) for i in mc_bids[0].keys() if i not in n_query_bids])
        n_query_bids.extend([self.remapCollapsed(i) for i in mc_bids[1].keys() if i not in n_query_bids])
        
        (classification, info) = self.BM.classify(rowIndex, self.getNeighbours(n_query_bids))
        classification = self.remapCollapsed(classification)
        return (classification, "@@"+info)

#------------------------------------------------------------------------------
# REGIONS

    def regionalise(self, force=False, save=True):
        """Create regions within the SOMs"""
        print "    Creating classification regions"
        if(not force):
            # first check to see that the
            ids_in_use = self.DM.getSOMDataInfo(self.PM.dbFileName)
            soms_done = []
            for b in ["mer","cov"]:
                for a in ["weights","regions"]:
                    soms_done.append(len(ids_in_use[a][b]))
            if (sum(soms_done) > 0):
                # something's been done!
                if(self.promptOnOverwrite() != 'Y'):
                    print "Operation cancelled"
                    return False
                else:
                    print "Overwriting SOM regions in db:", self.PM.dbFileName

        bids = self.BM.getBids()
        self.loadTrainingVectors()
        
        # build mer regions
        for i in self.merSoms:
            self.merSoms[i].regionalise(bids, self.kVecs)
            if(save):
                self.saveMerRegions(i)
            
        # build coverage regions
        for i in self.covSoms:
            self.covSoms[i].regionalise(bids, self.cVecs)
            if(save):
                self.saveCovRegions(i)

    def findRegionNeighbours(self, merge=False, printMergers=False):
        """Find out which regions neighbour which other regions"""
        print "    Finding region neighbours"
        mer_Ns = {}
        cov_Ns = {}
        self.collapsedMappings = {}
        # first find out who is next to whom
        for i in self.merSoms:
            for N in self.merSoms[i].findRegionNeighbours():
                if(N in mer_Ns):
                    mer_Ns[N] += 1
                else:
                    mer_Ns[N] = 1
        for i in self.covSoms:
            for N in self.covSoms[i].findRegionNeighbours():
                if(N in cov_Ns):
                    cov_Ns[N] += 1
                else:
                    cov_Ns[N] = 1
        
        # now find out who is consistently next to whom
        combined_Ns = {}
        for N in mer_Ns:
            combined_Ns[N] = mer_Ns[N]
        for N in cov_Ns:
            if(N in combined_Ns):
                combined_Ns[N] += cov_Ns[N]
                
        # now refine tis search further
        filtered_Ns = {}
        for N in combined_Ns:
            if(combined_Ns[N] >= 4):
                filtered_Ns[N] = combined_Ns[N]

        
        # now back it up with some stats
        for N in filtered_Ns:
            #print "FN", N, filtered_Ns[N],
            bin1 = self.BM.getBin(N[0])
            bin2 = self.BM.getBin(N[1])
            if(self.BM.shouldMerge(bin1, bin2)):
            #if(bin1.isSimilar(bin2)): # this test is symmetrical
                # always map down to the smaller
                self.collapsedMappings[N[1]] = N[0]
            #    print True
            #else:
            #    print False

        if(printMergers and not merge):
            ml = self.makeMergeLists(verbose=True)
            #for mml in ml:
            #    print mml
        elif(merge):
            self.merge()

    def getNeighbours(self, bids=[]):
        """Return a list of neighbours based on regions"""
        if(bids == []): return []
        ret_list = list(bids)
        for i in self.merSoms:
            ret_list.extend([i for i in self.merSoms[i].getNeighbours(bids) if i not in ret_list])
        for i in self.covSoms:
            ret_list.extend([i for i in self.covSoms[i].getNeighbours(bids) if i not in ret_list])
        return ret_list

    def validateRegions(self):
        """Basic validation of regions
        
        Classify each contig fromeach bin and see if the classification is
        correct
        """
        print "Validating regions"
        
        adds = {}
        removes = {}
        
        self.loadTrainingVectors()
        bids = self.BM.getBids()
        for bid in bids:
            adds[bid] = []
            removes[bid] = {}
            bin = self.BM.getBin(bid)
            tmp_bid = self.remapCollapsed(bid)
            print "    BID:",bid," (",tmp_bid,"):",
            total = 0
            correct = 0
            incorrect = 0
            unassigned = 0
            for row_index in bin.rowIndicies:
                total += 1
                (class_bid, info) = self.classify(row_index)
                print "{",class_bid,
                if(class_bid == 0):
                    unassigned += 1
                elif(class_bid != tmp_bid):
                    incorrect += 1
                    print info
                    removes[bid][row_index] = True
                    if(class_bid not in adds):
                        adds[class_bid] = []  
                    adds[class_bid].append(row_index)
                else:
                    correct += 1
                print "}",
            print "\nTotal %d, correct: %d, incorrect %d, unassigned %d" % (total,correct,incorrect,unassigned)
        #print adds
        #print removes
        for bid in bids:
            bin = self.BM.getBin(bid)
            bin.shuffleMembers(adds[bid], removes[bid])
        self.BM.saveBins()

    def assignmentMerge(self):
        """Determine which regions should be merged
        
        Classify each contig from each bin and see if the classification is
        correct. Keep tabs on where it's wrong and use that to choose
        which regions to merge
        """
        self.loadTrainingVectors()
        bids = self.BM.getBids()
        total_bids = len (self.BM.getBids())
        done_bids = 0
        joins = {}
        for bid in bids:
            incorrect_assignments = {}
            bin = self.BM.getBin(bid)
            done_bids += 1
            print "    classifying %d (%d) of %d" % (done_bids, len(bin.rowIndicies), total_bids)
            for row_index in bin.rowIndicies:
                (class_bid, info) = self.classify(row_index)
                if(class_bid != bid):
                    # incorrect assignment, make note of which bin it was
                    if(class_bid not in incorrect_assignments):
                        incorrect_assignments[class_bid] = 1
                    else:
                        incorrect_assignments[class_bid] += 1
            # make tuples of all the incorrect matches,
            # this way we can find the worst offenders
            for i in incorrect_assignments:
                if(incorrect_assignments[i] > 1):
                    join_key = self.makeNTuple(bid, i)
                    if(join_key not in joins):
                        joins[join_key] = incorrect_assignments[i]
                    else:
                        joins[join_key] += incorrect_assignments[i]
                    
        # sort the possible joins in descending order of number of instances
        # those with the highest number *should* be the best...
        import operator
        sorted_joins = sorted(joins.iteritems(), key=operator.itemgetter(1), reverse=True)
        for join in sorted_joins:
            if(join[1] > 2):    # don't take anything too spurious
                bid1 = self.remapCollapsed(join[0][0])
                bid2 = self.remapCollapsed(join[0][1])
                # stop circular mergers
                if(bid1 != bid2):
                    bin1 = self.BM.getBin(bid1)
                    bin2 = self.BM.getBin(bid2)
                    should_merge = self.BM.shouldMerge(bin1, bin2)
                    if(should_merge):
                        print bid1, bid2, join[1], should_merge
                        rv = self.BM.merge([bid1, bid2], saveBins=True, printInstructions=False)
                        if(rv == 2):
                            # merge happened
                            self.collapsedMappings[bid2] = bid1
                            
    def makeNTuple(self, bid1, bid2):
        """A way for making standard tuples from bids"""
        if(bid1 < bid2): return (bid1, bid2)
        return (bid2, bid1)
            
#------------------------------------------------------------------------------
# WEIGHTS 

    def buildSomWeights(self, force=False, save=True, plot=False):
        """Construct and save the weights matrix""" 
        if(not force):
            # first check to see that the
            ids_in_use = self.DM.getSOMDataInfo(self.PM.dbFileName)
            soms_done = []
            for b in ["mer","cov"]:
                for a in ["weights","regions"]:
                    soms_done.append(len(ids_in_use[a][b]))
            if (sum(soms_done) > 0):
                # something's been done!
                if(self.promptOnOverwrite() != 'Y'):
                    print "Operation cancelled"
                    return False
                else:
                    print "Overwriting SOM weights in db:", self.PM.dbFileName
        
        # now we can start
        print "Building %d SOMs of each type (coverage + kmer) with grid side %d," % (self.numSoms, self.somSide) 
        
        # get training data
        self.loadTrainingVectors()
        
        # build kmer SOMS
        for i in range(self.numSoms):
            print "\n    Training kmer SOM #%d" % i
            map = som.SOM(self.somSide,self.merDim)
            self.merSoms[i] = map
            if(i == 0 and plot):
                map.train(self.kVecs, iterations=self.somIterations, weightImgFileName="mer")
            else:
                map.train(self.kVecs, iterations=self.somIterations)
            if(save):
                self.saveMerWeights(i)

        # build coverage SOMS
        for i in range(self.numSoms):
            print "\n    Training coverage SOM #%d" % i
            map = som.SOM(self.somSide,self.covDim)
            self.covSoms[i] = map
            if(i == 0 and plot):
                map.train(self.cVecs, iterations=self.somIterations, weightImgFileName="cov")
            else:
                map.train(self.cVecs, iterations=self.somIterations)
            if(save):
                self.saveCovWeights(i)
        print "--"
        return True
    
    def whiten(self, profile):
        """Z normalize and scale profile columns"""
        v_mean = np.mean(profile, axis=0)
        v_std = np.std(profile, axis=0)
        profile = (profile-v_mean)/v_std
        v_mins = np.min(profile, axis=0)
        profile -= v_mins
        v_maxs = np.max(profile, axis=0)
        profile /= v_maxs
        return (profile, v_mean, v_std, v_mins, v_maxs)

    def whitenKVector(self, vector):
        """Z normalize and scale individual vectors"""
        vector = (vector-self.kMeans)/self.kStdevs
        vector = (vector - self.kMins)/self.kMaxs
        return np.clip(vector,0,1)

    def whitenCVector(self, vector):
        """Z normalize and scale individual vectors"""
        vector = (vector-self.cMeans)/self.cStdevs
        vector = (vector - self.cMins)/self.cMaxs
        return np.clip(vector,0,1)

#------------------------------------------------------------------------------
# MERGE BINS

    def makeMergeLists(self, verbose=False):
        """Use the collapsed mappings to build a set of merging lists"""
        working_lists = {}
        for bid in self.collapsedMappings:
            if(self.collapsedMappings[bid] not in working_lists and bid not in working_lists):
                # both new!
                tmp = [self.collapsedMappings[bid], bid]
                working_lists[self.collapsedMappings[bid]] = tmp
                working_lists[bid] = tmp
            elif(self.collapsedMappings[bid] not in working_lists):
                working_lists[self.collapsedMappings[bid]] = working_lists[bid]
                working_lists[self.collapsedMappings[bid]].append(self.collapsedMappings[bid])
            elif(bid not in working_lists):                 
                working_lists[bid] = working_lists[self.collapsedMappings[bid]]
                working_lists[bid].append(bid)
            # else both in already
        
        merge_lists = []
        used_ids = {}
        for bid in working_lists:
            if(bid not in used_ids):
                merge_lists.append(working_lists[bid])
                for inner_id in working_lists[bid]:
                    used_ids[inner_id] = True
        if(verbose):
            num_reduced = 0
            for ml in merge_lists:
                print ml
                num_reduced += (len(ml) - 1)
            print "    Merging %d into %d bins, leaving %d bins" % (num_reduced,len(merge_lists),(len(self.BM.bins.keys())-num_reduced))
        return merge_lists

    def merge(self):
        """Merge bins, keeping the soms informed of changes"""
        # self.collapsedMappings is a tree of values where key must be merged with value
        print "    Merging globally adjacent regions"
        merge_lists = self.makeMergeLists(verbose=True)
        for ml in merge_lists:
            self.BM.merge(ml, auto=True, saveBins=True, printInstructions=False)
        
        # remake the regions
        self.regionalise(force=True)
        self.findRegionNeighbours(merge=False, printMergers=False)

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def renderWeights(self, tag):
        """Render all the weights for all the SOMS"""
        for key in self.covSoms.keys():
            self.covSoms[key].renderWeights(tag+"_covWeights_"+str(key))
        for key in self.merSoms.keys():
            self.merSoms[key].renderWeights(tag+"_merWeights_"+str(key))

    def renderRegions(self, tag):
        """Render all the weights for all the SOMS"""
        # make the palette
        palette = self.BM.makeCentroidPalette()
        for key in self.covSoms.keys():
            self.covSoms[key].renderRegions(tag+"_covRegions_"+str(key), palette)
        for key in self.merSoms.keys():
            self.merSoms[key].renderRegions(tag+"_merRegions_"+str(key), palette)
    
    def renderMulti(self, tag):
        """Render large image including regions bids etc"""
        pass
    
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
        self.dataManager = mstore.GMDataManager()  # most data is saved to hdf
        self.dbFileName = dbFileName        # db containing all the data we'd like to use
        self.condition = ""                 # condition will be supplied at loading time
        # --> NOTE: ALL of the arrays in this section are in sync
        # --> each one holds information for an individual contig 
        self.indicies = np.array([])        # indicies into the data structure based on condition
        self.covProfiles = np.array([])     # coverage based coordinates
        self.transformedCP = np.array([])   # the munged data points
        self.averageCoverages = np.array([]) # average coverage across all stoits
        self.kmerSigs = np.array([])        # raw kmer signatures
        self.kmerVals = np.array([])        # PCA'd kmer sigs

        self.contigNames = np.array([])
        self.contigLengths = np.array([])
        self.contigColours = np.array([])   # calculated from kmerVals
        
        self.binIds = np.array([])          # list of bin IDs
        self.isCore = np.array([])          # True False values
        # --> end section

        # meta                
        self.validBinIds = {}               # valid bin ids -> numMembers
        self.binnedRowIndicies = {}         # dictionary of those indicies which belong to some bin
        self.restrictedRowIndicies = {}     # dictionary of those indicies which can not be binned yet
        self.numContigs = 0                 # this depends on the condition given
        self.numStoits = 0                  # this depends on the data which was parsed

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
                 loadCores=False):
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
                print "    Loading indicies (", condition,")"
            self.indicies = self.dataManager.getConditionalIndicies(self.dbFileName, condition=condition)
            self.numContigs = len(self.indicies)
            
            if(not silent):
                print "    Working with: %d contigs" % self.numContigs

            if(loadCovProfiles):
                if(verbose):
                    print "    Loading coverage profiles"
                self.covProfiles = self.dataManager.getCoverageProfiles(self.dbFileName, indicies=self.indicies)

                # work out average coverages
                self.averageCoverages = np.array([sum(i)/self.numStoits for i in self.covProfiles])

            if(loadKmerSigs):
                if(verbose):
                    print "    Loading kmer sigs"
                self.kmerSigs = self.dataManager.getKmerSigs(self.dbFileName, indicies=self.indicies)

                if(makeColours):
                    if(verbose):
                        print "    Creating colour profiles"
                    self.makeColourProfile()
                    # use HSV to RGB to generate colours
                    S = 1       # SAT and VAL remain fixed at 1. Reduce to make
                    V = 1       # Pastels if that's your preference...
                    self.contigColours = np.array([colorsys.hsv_to_rgb(val, S, V) for val in self.kmerVals])

            if(loadContigNames):
                if(verbose):
                    print "    Loading contig names"
                self.contigNames = self.dataManager.getContigNames(self.dbFileName, indicies=self.indicies)
            
            if(loadContigLengths):
                if(verbose):
                    print "    Loading contig lengths"
                self.contigLengths = self.dataManager.getContigLengths(self.dbFileName, indicies=self.indicies)
                print "    Contigs contain %d BP" % ( sum(self.contigLengths) )
            
            if(loadBins):
                if(verbose):
                    print "    Loading bins"
                self.binIds = self.dataManager.getBins(self.dbFileName, indicies=self.indicies)
                if(len(bids) != 0): # need to make sure we're not restricted in terms of bins
                    tmp_bids = self.getBinStats()
                    for bid in bids:
                        self.validBinIds[bid] = tmp_bids[bid]
                else:
                    self.validBinIds = self.getBinStats()

                # fix the binned indicies
                self.binnedRowIndicies = {}
                for i in range(len(self.indicies)):
                    if(self.binIds[i] != 0):
                        self.binnedRowIndicies[i] = True 

            if(loadCores):
                if(verbose):
                    print "    Loading core info"
                self.isCore = self.dataManager.getCores(self.dbFileName, indicies=self.indicies)
            
        except:
            print "Error loading DB:", self.dbFileName, sys.exc_info()[0]
            raise

    def reduceIndicies(self, deadRowIndicies):
        """purge indicies from the data structures
        
        Be sure that deadRowIndicies are sorted ascending
        """
        # strip out the other values        
        self.indicies = np.delete(self.indicies, deadRowIndicies, axis=0)
        self.covProfiles = np.delete(self.covProfiles, deadRowIndicies, axis=0)
        self.transformedCP = np.delete(self.transformedCP, deadRowIndicies, axis=0)
        self.contigNames = np.delete(self.contigNames, deadRowIndicies, axis=0)
        self.contigLengths = np.delete(self.contigLengths, deadRowIndicies, axis=0)
        self.contigColours = np.delete(self.contigColours, deadRowIndicies, axis=0)
        self.kmerSigs = np.delete(self.kmerSigs, deadRowIndicies, axis=0)
        self.kmerVals = np.delete(self.kmerVals, deadRowIndicies, axis=0)
        self.binIds = np.delete(self.binIds, deadRowIndicies, axis=0)
        self.isCore = np.delete(self.isCore, deadRowIndicies, axis=0)
        
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
    
    def saveBinIds(self, updates):
        """Save our bins into the DB"""
        self.dataManager.setBins(self.dbFileName, updates)
    
    def saveCores(self, updates):
        """Save our core flags into the DB"""
        self.dataManager.setCores(self.dbFileName, updates)

    def saveValidBinIds(self, updates):
        """Store the valid bin Ids and number of members
                
        updates is a dictionary which looks like:
        { tableRow : [bid , numMembers] }
        """
        self.dataManager.setBinStats(self.dbFileName, updates)
        self.setNumBins(len(updates.keys()))

    def updateValidBinIds(self, updates):
        """Store the valid bin Ids and number of members
        
        updates is a dictionary which looks like:
        { bid : numMembers }
        if numMembers == 0 then the bid is removed from the table
        if bid is not in the table yet then it is added
        otherwise it is updated
        """
        # get the current guys
        existing_bin_stats = self.dataManager.getBinStats(self.dbFileName)
        num_bins = self.getNumBins()
        # now update this dict
        for bid in updates.keys():
            if bid in existing_bin_stats:
                if updates[bid] == 0:
                    # remove this guy
                    del existing_bin_stats[bid]
                    num_bins -= 1
                else:
                    # update the count
                    existing_bin_stats[bid] = updates[bid]
            else:
                # new guy!
                existing_bin_stats[bid] = updates[bid]
                num_bins += 1
        
        # finally , save
        self.saveValidBinIds(existing_bin_stats)

#------------------------------------------------------------------------------
# DATA TRANSFORMATIONS 

    def getAverageCoverage(self, rowIndex):
        """Return the average coverage for this contig across all stoits"""
        return sum(self.transformedCP[rowIndex])/self.numStoits

    def transformCP(self, silent=False, nolog=False, min=None, max=None):
        """Do the main ransformation on the coverage profile data"""
        shrinkFn = np.log10
        if(nolog):
            shrinkFn = lambda x:x
         
        s = (self.numContigs,3)
        self.transformedCP = np.zeros(s)

        if(not silent):
            print "    Dimensionality reduction"

        # get the median distance from the origin
        unit_vectors = [(np.cos(i*2*np.pi/self.numStoits),np.sin(i*2*np.pi/self.numStoits)) for i in range(self.numStoits)]
        for i in range(len(self.indicies)):
            norm = np.linalg.norm(self.covProfiles[i])
            if(norm != 0):
                radial = shrinkFn(norm)
            else:
                radial = norm
            shifted_vector = np.array([0.0,0.0])
            flat_vector = (self.covProfiles[i] / sum(self.covProfiles[i]))
            
            for j in range(self.numStoits):
                shifted_vector[0] += unit_vectors[j][0] * flat_vector[j]
                shifted_vector[1] += unit_vectors[j][1] * flat_vector[j]

            # log scale it towards the centre
            scaling_vector = shifted_vector * self.scaleFactor
            sv_size = np.linalg.norm(scaling_vector)
            if(sv_size > 1):
                shifted_vector /= shrinkFn(sv_size)

            self.transformedCP[i,0] = shifted_vector[0]
            self.transformedCP[i,1] = shifted_vector[1]
            self.transformedCP[i,2] = radial

        if(not silent):
            print "    Reticulating splines"
            
        # finally scale the matrix to make it equal in all dimensions
        if(min is None):                
            min = np.amin(self.transformedCP, axis=0)
            max = np.amax(self.transformedCP, axis=0)
            max = max - min
            max = max / (self.scaleFactor-1)

        for i in range(0,3):
            self.transformedCP[:,i] = (self.transformedCP[:,i] -  min[i])/max[i]
            
        return(min,max)

    def makeColourProfile(self):
        """Make a colour profile based on ksig information"""
        working_data = np.array(self.kmerSigs, copy=True) 
        PCA.Center(working_data,verbose=0)
        p = PCA.PCA(working_data)
        components = p.pc()
        
        # now make the colour profile based on PC1
        self.kmerVals = np.array([float(i) for i in components[:,0]])
        
        # normalise to fit between 0 and 1
        self.kmerVals -= np.min(self.kmerVals)
        self.kmerVals /= np.max(self.kmerVals)
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
        return V_p/np.linalg.norm(V_p)
    
    def rad2deg(self, anglein):
        return 180*anglein/np.pi

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
            print "Error showing image", sys.exc_info()[0]
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
                print "Error saving image",fileName, sys.exc_info()[0]
                raise
        elif(show):
            try:
                plt.show()
            except:
                print "Error showing image", sys.exc_info()[0]
                raise
        if del_fig:
            plt.close(fig)
            del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

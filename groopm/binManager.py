#!/usr/bin/env python
###############################################################################
#                                                                             #
#    binManager.py                                                            #
#                                                                             #
#    GroopM - High level bin data management                                  #
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
from os.path import join as osp_join
from sys import exc_info, exit, stdout as sys_stdout
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
from scipy.cluster.vq import kmeans,vq

# GroopM imports
from profileManager import ProfileManager
from bin import Bin
import groopmExceptions as ge
from groopmUtils import makeSurePathExists

np_seterr(all='raise')     

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
                 loadLinks=False,
                 min=None,
                 max=None,
                 cutOff=0,
                 transform=True):
        """Load data and make bin objects"""
        # fix the condition
        condition=""
        if(cutOff != 0):
            condition="length >= %d" % cutOff
        elif(len(bids) == 0):
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
                         loadLinks=loadLinks
                        )
        
        if(makeBins):
            if transform:
                self.PM.transformCP(silent=silent, min=min, max=max)
            else:
                if self.PM.numStoits == 3:
                    self.PM.transformedCP = self.PM.covProfiles
                else:
                    print "Number of stoits != 3. You need to transform"
                    self.PM.transformCP(silent=silent, min=min, max=max)
            self.makeBins(self.getBinMembers())

    def getBinMembers(self):
        """Munge the raw data into something more usable
        
        self.PM.binIds is an array, contains 0 for unassigned rows
        By default this creates an array for the '0' bin. You need
        to ignore it later if you want to
        """
        # fill them up
        bin_members = {0:[]}
        for row_index in range(np_size(self.PM.indices)):
            try:
                bin_members[self.PM.binIds[row_index]].append(row_index)
            except KeyError:
                bin_members[self.PM.binIds[row_index]] = [row_index]

        # we need to get the largest BinId in use
        if len(bin_members) > 0:
            self.nextFreeBinId = np_max(bin_members.keys())
        
        return bin_members

    def makeBins(self, binMembers, zeroIsBin=False):
        """Make bin objects from loaded data"""
        invalid_bids = []
        for bid in binMembers:
            if bid != 0 or zeroIsBin:
                if len(binMembers[bid]) == 0:
                    invalid_bids.append(bid)
                else:
                    self.bins[bid] = Bin(np_array(binMembers[bid]), bid, self.PM.scaleFactor-1)
                    self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
        if len(invalid_bids) != 0:
            print "MT bins!"
            print invalid_bids
            exit(-1)      

    def saveBins(self, binAssignments={}):
        """Save binning results
        
        binAssignments is a hash of LOCAL row indices Vs bin ids
        { row_index : bid } 
        PM.setBinAssignments needs GLOBAL row indices
        
        We always overwrite the bins table (It is smallish)
        """

        # save the bin assignments        
        self.PM.setBinAssignments(
                                  self.getGlobalBinAssignments(binAssignments) # convert to global indices
                                  )
        # overwrite the bins table
        self.setBinStats()
        
        # we must have done something
        self.PM.setClustered()

    def getGlobalBinAssignments(self, binAssignments={}):
        """Merge the bids, raw DB indexes and core information so we can save to disk
        
        returns a hash of type:
        
        { global_index : bid }
        """
        # we need a mapping from cid (or local index) to to global index to binID
        bin_assignment_update = {}
        
        if binAssignments != {}:
            # we have been told to do only these guys
            for row_index in binAssignments:
                bin_assignment_update[self.PM.indices[row_index]] = binAssignments[row_index]

        else:
            # this are all our regularly binned guys
            for bid in self.getBids():
                for row_index in self.bins[bid].rowIndices:
                    bin_assignment_update[self.PM.indices[row_index]] = bid
            
        return bin_assignment_update

    def setBinStats(self):
        """Update / overwrite the table holding the bin stats
        
        Note that this call effectively nukes the existing table
        """
        bin_stats = {}
        for bid in self.getBids():
            # no point in saving empty bins
            if np_size(self.bins[bid].rowIndices) > 0:
                bin_stats[bid] = np_size(self.bins[bid].rowIndices)
        self.PM.setBinStats(bin_stats)

    
#------------------------------------------------------------------------------
# REMOVE ALREADY LOADED BINS

    def removeBinAndIndicies(self, bid):
        """Remove indices from the PM based on bin identity
        
        "unload" some data
        """
        # get some info
        rem_bin = self.getBin(bid)
        original_length = len(self.PM.indices)
        rem_list = np_sort(rem_bin.rowIndices)
        
        # affect the raw data in the PM
        self.PM.reduceIndicies(rem_list)
        del self.PM.validBinIds[bid]
        
        # remove the bin here
        del self.bins[bid]
        
        # now fix all the rowIndices in all the other bins
        for bid in self.getBids():
            self.bins[bid].rowIndices = self.fixRowIndexLists(original_length, np_sort(self.bins[bid].rowIndices), rem_list)


    def fixRowIndexLists(self, originalLength, oldList, remList):
        """Fix up row index lists which reference into the
        data structure after a call to reduceIndicies
        
        originalLength is the length of all possible row indices
        before the removal (ie self.indices)
        oldList is the old list of row indices
        remList is the list of indices to be removed
        
        BOTH OLD AND REM LIST MUST BE SORTED ASCENDING!
        """
        shift_down = 0;
        old_list_index = 0
        new_list = np_array([])
        for i in range(originalLength):
            if(i in remList):
                shift_down+=1
            elif(i in oldList):
                new_list = np_append(new_list, oldList[old_list_index]-shift_down)
                old_list_index += 1
        return new_list

#------------------------------------------------------------------------------
# NEIGHBOURS AND DISTANCES

    def findBinNeighbours(self, thresholdDist=50.0):
        """Construct a network of all bins and their closest neighbours"""
        num_bins = len(self.bins)
        bids = self.getBids()
        cov_centres = np_reshape([self.bins[bid].covMeans for bid in bids], (num_bins,3))
        
        # get an all vs all distance matrix
        c_dists = cdist(cov_centres, cov_centres)
        
        # reduce this to only close neighbours
        neigbour_dists = np_where(c_dists < thresholdDist, c_dists, 0.0)
        
        # now make the network
        network = {}
        outer_index = 0
        for i in range(num_bins):
            # make a structure to hold the info
            network[bids[i]] = [[bids[i]],[0.0]]
            for j in range(num_bins):
                if(neigbour_dists[i,j] != 0.0):
                    # this is a legit guy!
                    network[bids[i]][0].append(bids[j])
                    network[bids[i]][1].append(neigbour_dists[i,j])
        return network

#------------------------------------------------------------------------------
# LINKS

    def getLinkingContigs(self, bid):
        """Get all contigs and their bin IDs which link to contigs in this bin"""
        condition = ""
        bin = self.getBin(bid)
        bin2count = {}
        for row_index in bin.rowIndices:
            try: 
                #print row_index, len(self.PM.links[row_index]), self.PM.links[row_index], "===", 
                for link in self.PM.links[row_index]:
                    #print "{{{",link,"}}}",
                    try:
                        link_bid = self.PM.binIds[link[0]]
                        #print ";;", link_bid, bid ,
                        if link_bid != bid and link_bid != 0:
                            try: 
                                bin2count[link_bid] += 1.0
                            except KeyError:
                                bin2count[link_bid] = 1.0
                    except KeyError:
                        pass#print "****\n\n" 
                #print "[[[[\n\n"
            except KeyError:
                pass
        #print bin2count
        return bin2count
    
    def getConnectedBins(self, rowIndex):
        """Get a  list of bins connected to this contig"""
        ret_links = []
        for link in self.PM.links[rowIndex]:
            cid = link[0]
            try:
                bid = self.PM.binIds[cid]
            except KeyError:
                bid = 0
            ret_links.append((cid, bid, link[1])) 
        return ret_links
    
    def getAllLinks(self):
        """Return a sorted array of all links between all bins"""
        bids = self.getBids()
        # first, work out who links with whom...       
        all_links = {}
        for bid in bids:
            links = self.getLinkingContigs(bid)
            # links is a hash of type bid : num_links
            for link in links:
                key = self.makeBidKey(bid, link)
                if key not in all_links:
                    all_links[key] = links[link]
    
        # sort and return
        return sorted(all_links.iteritems(), key=itemgetter(1), reverse=True)
       
    def getWithinLinkProfiles(self):
        """Determine the average number of links between contigs for all bins"""
        bids = self.getBids()
        link_profiles = {}
        for bid in bids:
            link_profiles[bid] = self.getWithinBinLinkProfile(bid)
        return link_profiles
        
    def getWithinBinLinkProfile(self, bid):
        """Determine the average number of links between contigs in a bin"""
        bin = self.getBin(bid)
        links = []
        min_links = 1000000000
        for row_index in bin.rowIndices:
            try: 
                for link in self.PM.links[row_index]:
                    link_bid = self.PM.binIds[link[0]] 
                    if link_bid == bid:
                        links.append(link[1])
                        if link[1] < min_links:
                            min_links = link[1] 
            except KeyError:
                pass
        return (np_mean(links), np_std(links), min_links)
    
#------------------------------------------------------------------------------
# BIN REFINEMENT AND EXPANSION

    def recruitWrapper(self, inclusivity=2, step=200, saveBins=False):
        """Recuit more contigs to the bins"""
        print "Recruiting unbinned contigs"
        # make a list of all the cov and kmer vals
        num_bins = len(self.bins)
        num_expanded = 1
        total_expanded = 0
        total_binned = 0
        total_unbinned = 0
        total_contigs = len(self.PM.indices)
        shortest_binned = 1000000000          # we need to know this
        shortest_unbinned = 1000000000
        
        # we need to get a list of bin centroids
        (bin_centroid_points,
         bin_centroid_colours,
         bin_centroid_kvals,
         bids) = self.findCoreCentres(getKVals=True)
        # centroids
        tdm_centroid = np_append(bin_centroid_points,
                                 1000*np_reshape(bin_centroid_kvals,(len(bin_centroid_kvals),1)),
                                 1)
        search_tree = kdt(tdm_centroid)
        # contigs
        tdm = np_append(self.PM.transformedCP,
                        1000*np_reshape(self.PM.kmerVals,(len(self.PM.kmerVals),1)),
                        1)
        # for stats, work out number binned and unbinned and relative lengths
        unbinned = {}
        for row_index in range(len(self.PM.indices)):
            if(row_index in self.PM.binnedRowIndicies):
                if self.PM.contigLengths[row_index] < shortest_binned:
                    shortest_binned = self.PM.contigLengths[row_index] 
                total_binned += 1
            else:
                if self.PM.contigLengths[row_index] < shortest_unbinned:
                    shortest_unbinned = self.PM.contigLengths[row_index] 
                total_unbinned += 1
                unbinned[row_index] = self.PM.contigLengths[row_index] 
        
        # work out how many iterations we'll do
        if shortest_binned > shortest_unbinned:
            size_range = shortest_binned - shortest_unbinned 
            num_steps = size_range/step
            if num_steps == 0:
                steps = [shortest_unbinned]
            else:
                step_size = size_range/num_steps
                steps = [shortest_binned - i*step_size for i in range(1,num_steps)]
                steps.append(shortest_unbinned)
        else:
            steps = [shortest_unbinned]

        # talk to the user
        perc_binned = float(total_binned)/float(total_contigs)
        print "    Planned steps = ", steps
        print "    BEGIN: %0.4f" % perc_binned +"%"+" of %d requested contigs in bins" % total_contigs
        print "    %d contigs unbinned" % total_unbinned
        
        # go through the steps we decided on
        for cutoff in steps:
            print "    Recruiting contigs above: %d" % cutoff
            newly_binned = [0]
            this_step_binned = 0 
            while len(newly_binned) > 0:
                newly_binned = []
                affected_bids = []
                for row_index in unbinned:
                    if unbinned[row_index] >= cutoff:
                        # meets our criteria
                        putative_bid = int(bids[search_tree.query(tdm[row_index])[1]])
                        (covZ,merZ) = self.scoreContig(row_index, putative_bid)
                        if covZ <= inclusivity and merZ <= inclusivity:
                            # we can recruit
                            self.bins[putative_bid].rowIndices = np_append(self.bins[putative_bid].rowIndices,
                                                                            row_index
                                                                            ) 
                            affected_bids.append(putative_bid)
                            newly_binned.append(row_index)
                            this_step_binned += 1
                            total_binned += 1
                            total_expanded += 1
    
                # remove any binned contigs from the unbinned list
                for row_index in newly_binned:
                    del unbinned[row_index]

                # remake bin stats
                for bid in affected_bids:
                    self.bins[bid].makeBinDist(self.PM.transformedCP,
                                               self.PM.averageCoverages,
                                               self.PM.kmerVals,
                                               self.PM.contigLengths)      

            print "    Recruited: %d contigs" % this_step_binned

        # talk to the user
        perc_recruited = float(total_expanded)/float(total_unbinned)
        perc_binned = float(total_binned)/float(total_contigs)
        print "    Recruited %0.4f" % perc_recruited +"%"+" of %d unbinned contigs" % total_unbinned
        print "    END: %0.4f" % perc_binned +"%"+" of %d requested contigs in bins" % total_contigs
        
        # now save
        if(saveBins):
            self.saveBins()
    
    def refineWrapper(self,
                      manual=False,          # do we need to ask permission every time?
                      saveBins=False,
                      plotter=False,
                      shuffle=False,
                      links=False
                      ):
        """Iterative wrapper for the refine function"""
        if plotter:
            self.plotterRefineBins()
        if shuffle:
            print "Start automatic bin refinement"
            self.autoRefineBins(iterate=True)
            num_binned = len(self.PM.binnedRowIndicies.keys())
            print "   ",num_binned,"contigs across",len(self.bins.keys()),"cores"
            
            if saveBins:
                self.saveBins()
        if links:
            # we don't load links by default, so lets do it now
            print "    Loading links"
            self.PM.loadLinks()
            self.refineViaLinks()
            if saveBins:
                self.saveBins()

    def refineViaLinks(self):
        """Use linking information to refine bins"""
        bin_links = self.getAllLinks()
        print bin_links
                        
    def plotterRefineBins(self):
        """combine similar bins using 3d plots"""
        self.printRefinePlotterInstructions()
        self.plotBinIds()
        continue_merge = True
        while(continue_merge):
            user_option = self.promptOnPlotterRefine()
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

    def getClosestBID(self, rowIndex, searchTree, tdm, neighbourList={}, k=101, verbose=False):
        """Find the bin ID which would best describe the placement of the contig
        
        The neighbourlist is a hash of type:
        
        row_index : [neighbour, neighbour, neighbour, ...]
        
        This is built on the fly and added to or updated as need be
        """
        # find the k nearest neighbours for the query contig
        if k < 21:
            k = 21
        if k > 101:
            k = 101
        
        try:
            t_list = neighbourList[rowIndex]
            if len(t_list) > k:
                # must have made a change, we'll fix it here
                t_list = searchTree.query(tdm[rowIndex],k=k)[1]
                neighbourList[rowIndex] = t_list
        except KeyError:
            # first time, we need to do the search
            t_list = searchTree.query(tdm[rowIndex],k=k)[1]
            neighbourList[rowIndex] = t_list
            
        # calculate the distribution of neighbouring bins
        refined_t_list = {}
        for row_index in t_list:
            #print len(self.PM.binIds), "::", row_index 
            try:
                cbid = self.PM.binIds[row_index]
                if cbid != 0: # don't want these guys
                    try:
                        refined_t_list[cbid] += 1
                    except KeyError:
                        refined_t_list[cbid] = 1
            except KeyError:
                pass

        # work out the most prominent BID
        max_bid = 0
        max_count = 0
        for cbid in refined_t_list:
            if refined_t_list[cbid] > max_count:
                max_count = refined_t_list[cbid]
                max_bid = cbid
        if verbose:
            print rowIndex, "**",max_bid,"**"                        
            for cbid in refined_t_list:
                print "[", cbid, ",", refined_t_list[cbid], "]",
            print
        # we're done!
        return (max_bid, neighbourList)

    def autoRefineBins(self, iterate=False, verbose=False):
        """Automagically refine bins"""
        super_round = 1
        tdm = np_append(self.PM.transformedCP, 1000*np_reshape(self.PM.kmerVals,(len(self.PM.kmerVals),1)),1)
        neighbour_list={} # save looking things up a 1,000,000 times
        search_tree = kdt(tdm)
        while True:
            sr_contigs_reassigned = 0
            num_reassigned = -1
            round = 0
            stable_bids = {} # once a bin is stable it's stable!
            while num_reassigned != 0:
                num_reassigned = 0
                reassignment_map = {}
                moved_RIs = {}
                bids = self.getBids()
                for bid in bids:
                    bin = self.getBin(bid)
                    if bid in stable_bids:
                        try:
                            reassignment_map[bid] += list(bin.rowIndices)
                        except KeyError:
                            reassignment_map[bid] = list(bin.rowIndices)
                    else:
                        stable = True
                        for row_index in bin.rowIndices:
                            (assigned_bid, neighbour_list) = self.getClosestBID(row_index, search_tree, tdm, neighbourList=neighbour_list, verbose=verbose, k=2*bin.binSize-1)
                            if assigned_bid != bid:
                                stable = False
                                num_reassigned += 1
                                sr_contigs_reassigned += 1
                                moved_RIs[row_index] = assigned_bid 
                            
                            # keep track of where this guy lives
                            try:
                                reassignment_map[assigned_bid].append(row_index)
                            except KeyError:
                                reassignment_map[assigned_bid] = [row_index]
                                
                        if stable: # no changes this round, mark bin as stable
                            stable_bids[bid] = True
                
                # fix the lookup table
                for moved_index in moved_RIs:
                    self.PM.binIds[moved_index] = moved_RIs[moved_index]
    
                # now fix the bins
                bins_removed = 0
                for bid in bids:
                    if bid in reassignment_map:
                        self.bins[bid].rowIndices = np_array(reassignment_map[bid])
                        self.bins[bid].binSize = len(reassignment_map[bid])
                    else:
                        # empty bin
                        bins_removed += 1
                        self.deleteBins([bid], force=True, freeBinnedRowIndicies=False, saveBins=False)
                
                round += 1
                if verbose:
                    print "    Refine round %d: reassigned %d contigs, removed %d cores" % (round, num_reassigned, bins_removed)
            print "    Refine round %d complete. (%d iterations) Total contigs reassigned: %d" % (super_round, round, sr_contigs_reassigned)
            if sr_contigs_reassigned == 0 or not iterate:
                break
            super_round += 1
            
#------------------------------------------------------------------------------
# BIN UTILITIES 

    def getBids(self):
        """Return a sorted list of bin ids"""
        return sorted(self.bins.keys())

    def getCentroidProfiles(self, mode="mer"):
        """Return an array containing the centroid stats for each bin"""
        if(mode == "mer"):
            ret_vecs = np_zeros((len(self.bins)))
            outer_index = 0
            for bid in self.getBids():
                ret_vecs[outer_index] = self.bins[bid].kValMean
                outer_index += 1
            return ret_vecs
        elif(mode == "cov"):
            ret_vecs = np_zeros((len(self.bins), len(self.PM.transformedCP[0])))
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
        if(printInstructions and not auto):
            self.printSplitInstructions()

        # make some split bins
        # bids[0] holds the original bin id
        (bin_assignment_update, bids) = self.getSplitties(bid, n, mode)
        
        if(auto and saveBins):
            # charge on through
            self.deleteBins([bids[0]], force=True)  # delete the combined bin
            # save new bins
            self.saveBins(binAssignments=bin_assignment_update)
            return

        # we will need to confer with the user
        # plot some stuff
        # sort the bins by kmer val
        bid_tuples = [(tbid, self.bins[tbid].kValMean) for tbid in bids[1:]]
        bid_tuples.sort(key=itemgetter(1))
        index = 1
        for pair in bid_tuples:
            bids[index] = pair[0]
            index += 1 
        self.plotSideBySide(bids)
        
        user_option = self.promptOnSplit(n,mode)
        if(user_option == 'Y'):
            if(saveBins):
                # save the temp bins
                self.deleteBins([bids[0]], force=True)  # delete the combined bin
                # save new bins
                self.saveBins(binAssignments=bin_assignment_update)
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
        obs = np_array([])
        if(mode=='kmer'):
            obs = np_array([self.PM.kmerVals[i] for i in bin.rowIndices])
        elif(mode=='cov'):
            obs = np_array([self.PM.covProfiles[i] for i in bin.rowIndices])
        
        # do the clustering
        try:
            centroids,_ = kmeans(obs,n)
        except ValueError:
            if(verbose):
                print "Error splitting"
            return False
        idx,_ = vq(obs,centroids)

        # build some temp bins 
        # this way we can show the user what the split will look like       
        idx_sorted = np_argsort(np_array(idx))
        current_group = 0
        bids = [bid]
        bin_assignment_update = {} # row index to bin id
        holding_array = np_array([])
        split_bin = None
        for i in idx_sorted:
            if(idx[i] != current_group):
                # bin is full!
                split_bin = self.makeNewBin(holding_array)
                for row_index in holding_array:
                    bin_assignment_update[row_index] = split_bin.id
                split_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
                bids.append(split_bin.id)
                holding_array = np_array([])
                current_group = idx[i]
            holding_array = np_append(holding_array, bin.rowIndices[i])
        # do the last one
        if(np_size(holding_array) != 0):
            split_bin = self.makeNewBin(holding_array)
            for row_index in holding_array:
                bin_assignment_update[row_index] = split_bin.id  
            split_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
            bids.append(split_bin.id)

        return (bin_assignment_update, bids)

    def shouldMerge(self, bin1, bin2, ignoreCov=False, ignoreMer=False, merTol=0, confidence=0.95, verbose=False):
        """Determine whether its wise to merge two bins
        
        Perfoms a one-way anova to determine if the larger bin would be 
        significantly changed if it consumed the smaller
        
        OR does a tolerance test on kmervals. We assume that bin1 is larger than bin2
        """
        if(bin1.id != bin2.id):
            if not ignoreCov: # work out coverage distributions
                b1_c_dist = bin1.getAverageCoverageDist(self.PM.averageCoverages) 
                b2_c_dist = bin2.getAverageCoverageDist(self.PM.averageCoverages)
                c_dist_1 = b1_c_dist
                if(bin1.binSize < bin2.binSize):
                    c_dist_1 = b2_c_dist
                c_dist_2 = np_append(b2_c_dist, b1_c_dist)
                if verbose:
                    tag = "COV:"
                else:
                    tag = "" 
                cov_match = self.isSameVariance(c_dist_1, c_dist_2, confidence=confidence, tag=tag) 
            else:
                cov_match = True
                
            if not ignoreMer: # work out kmer distributions
                if not cov_match:
                    return False
                if merTol != 0:
                    # Tolerance based testing
                    upper_k_val_cut = bin1.kValMean + merTol * bin1.kValStdev
                    lower_k_val_cut = bin1.kValMean - merTol * bin1.kValStdev
                
                    if bin2.kValMean >= lower_k_val_cut and bin2.kValMean <= upper_k_val_cut:
                        mer_match = True
                    else:
                        mer_match = False
                else:                   
                    b1_k_dist = bin1.getkmerValDist(self.PM.kmerVals) 
                    b2_k_dist = bin2.getkmerValDist(self.PM.kmerVals)
                    k_dist_1 = b1_k_dist
                    if(bin1.binSize < bin2.binSize):
                        k_dist_1 = b2_k_dist
                    k_dist_2 = np_append(b2_k_dist, b1_k_dist)
                    if verbose:
                        tag = "MER: %0.4f %0.4f" % (np_mean(k_dist_2), np_std(k_dist_2))
                    else:
                        tag = "" 
                    mer_match = self.isSameVariance(k_dist_1, k_dist_2, confidence=confidence, tag=tag)
            else:
                mer_match = True
                
            return cov_match and mer_match
         
        return False

    def isSameVariance(self, dist1, dist2, confidence=0.95, tag=""):
        """Test to see if the kmerValues for two bins are the same"""
        F_cutoff =  distributions.f.ppf(confidence, 2, len(dist1)+len(dist2)-2)
        F_value = f_oneway(dist1,dist2)[0]
        if tag != "":
           print "%s [V: %f, C: %f]" % (tag, F_value, F_cutoff)
        return F_value < F_cutoff

    def merge(self, bids, auto=False, manual=False, newBid=False, saveBins=False, verbose=False, printInstructions=True):
        """Merge two or more bins
        
        It's a bit strange to have both manual and auto params
        NOTE: manual ALWAYS overrides auto. In the condensing code, auto is
        set programmaticaly, manual is always set by the user. So we listen
        to manual first
        """
        parent_bin = None

        if(printInstructions and not auto):
            self.printMergeInstructions()

        if(newBid):
            # we need to make this into a new bin
            parent_bin = makeNewBin()
            # now merge it with the first in the new list
            dead_bin = self.getBin(bids[0])
            parent_bin.consume(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths, dead_bin, verbose=verbose)
            self.deleteBins([bids[0]], force=True)
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
                tmp_bin = self.makeNewBin(np_concatenate([parent_bin.rowIndices,dead_bin.rowIndices]))
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
                some_merged = True

        if some_merged:
            # Fix up the r2b indices and bin updates
            parent_bid = parent_bin.id
            bin_assignment_update = {}
            for row_index in parent_bin.rowIndices:
                bin_assignment_update[row_index] = parent_bid
                try:
                    self.PM.binIds[row_index] = parent_bid
                except KeyError:
                    pass 
    
            if saveBins:
                self.saveBins(binAssignments=bin_assignment_update)
            
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
        bin_assignment_update = {}
        for bid in bids:
            if bid in self.bins:
                if(freeBinnedRowIndicies):
                    for row_index in self.bins[bid].rowIndices:
                        if row_index in self.PM.binnedRowIndicies:
                            del self.PM.binnedRowIndicies[row_index]
                        else:
                            print bid, row_index, "FUNG"
                        bin_assignment_update[row_index] = 0 
                del self.bins[bid]
            else:
                raise ge.BinNotFoundException("Cannot find: "+str(bid)+" in bins dicts")
            
        if(saveBins):
            self.saveBins(binAssignments=bin_assignment_update)
        return True
        
    def makeNewBin(self, rowIndices=np_array([]), bid=None):
        """Make a new bin and add to the list of existing bins"""
        if bid is None:
            self.nextFreeBinId +=1
            bid = self.nextFreeBinId
        self.bins[bid] = Bin(rowIndices, bid, self.PM.scaleFactor-1)        
        return self.bins[bid]

#------------------------------------------------------------------------------
# UI 
    
    def printRefinePlotterInstructions(self):
        raw_input( "****************************************************************\n"
                   " REFINING INSTRUCTIONS - PLEASE READ CAREFULLY\n"+
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

    def promptOnPlotterRefine(self, minimal=False):
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
            return (np_mean(np_array(Ms.values())), np_std(np_array(Ms.values())), np_median(np_array(Ss.values())), np_std(np_array(Ss.values())))
        
        else:
            cutoff = np_mean(np_array(Ms.values())) + tolerance * np_std(np_array(Ms.values()))  
            kill_list = []
            for bid in Ms:
                if(Ms[bid] > cutoff):
                    kill_list.append(bid)
            return (kill_list, cutoff)

    def findCoreCentres(self, krange=None, getKVals=False):
        """Find the point representing the centre of each core"""
        print "    Finding bin centers"
        bin_centroid_points = np_array([])
        bin_centroid_colours = np_array([])
        bin_centroid_kvals = np_array([])
        bids = np_array([])
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
                ave_kval = np_mean([self.PM.kmerVals[row_index] for row_index in self.bins[bid].rowIndices])
                if ave_kval < k_low or ave_kval > k_high:
                    add_bin = False
            if add_bin:
                bin_centroid_points = np_append(bin_centroid_points,
                                                self.bins[bid].covMeans)
                bin_centroid_colours = np_append(bin_centroid_colours, 
                                                 np_mean([
                                                          self.PM.contigColours[row_index] for row_index in 
                                                          self.bins[bid].rowIndices
                                                          ],
                                                         axis=0)
                                                 )
                if getKVals:
                    bin_centroid_kvals = np_append(bin_centroid_kvals, 
                                                   np_mean([
                                                            self.PM.kmerVals[row_index] for row_index in 
                                                            self.bins[bid].rowIndices
                                                            ],
                                                           axis=0)
                                                   )

                bids = np_append(bids, bid)
                num_added += 1
        
        if num_added != 0:
            bin_centroid_points = np_reshape(bin_centroid_points, (num_added, 3))
            bin_centroid_colours = np_reshape(bin_centroid_colours, (num_added, 3))
                
        if getKVals:
            return (bin_centroid_points, bin_centroid_colours, bin_centroid_kvals, bids)
        return (bin_centroid_points, bin_centroid_colours, bids)

    def analyseBinKVariance(self, outlierTrim=0.1, plot=False):
        """Measure within and between bin variance of kmer sigs
        
        return a list of potentially confounding kmer indices
        """
        print "    Measuring kmer type variances"        
        means = np_array([])
        stdevs = np_array([])
        bids = np_array([])
        
        # work out the mean and stdev for the kmer sigs for each bin
        for bid in self.getBids():
            bkworking = np_array([])
            for row_index in self.bins[bid].rowIndices:
                bkworking = np_append(bkworking, self.PM.kmerSigs[row_index])
            bkworking = np_reshape(bkworking, (self.bins[bid].binSize, np_size(self.PM.kmerSigs[0])))
            bids = np_append(bids, [bid])
            means = np_append(means, np_mean(bkworking, axis=0))
            stdevs = np_append(stdevs, np_std(bkworking, axis=0))
            
        means = np_reshape(means, (len(self.bins), np_size(self.PM.kmerSigs[0])))
        stdevs = np_reshape(stdevs, (len(self.bins), np_size(self.PM.kmerSigs[0])))
        
        # now work out the between and within core variances
        between = np_std(means, axis=0)
        within = np_median(stdevs, axis=0)

        B = np_arange(0, np_size(self.PM.kmerSigs[0]), 1)
        names = self.PM.getMerColNames().split(',')
        
        # we'd like to find the indices of the worst 10% for each type so we can ignore them
        # specifically, we'd like to remove the least variable between core kms and the 
        # most variable within core kms.
        sort_between_indices = np_argsort(between)
        sort_within_indices = np_argsort(within)[::-1]
        number_to_trim = int(outlierTrim* float(np_size(self.PM.kmerSigs[0])))
        
        return_indices =[]
        for i in range(0,number_to_trim):
            if(sort_between_indices[i] not in return_indices):
                return_indices.append(sort_between_indices[i])
            if(sort_within_indices[i] not in return_indices):
                return_indices.append(sort_within_indices[i]) 
        
        if(plot):
            print "BETWEEN"
            for i in range(0,number_to_trim):
                print names[sort_between_indices[i]]
            print "WITHIN" 
            for i in range(0,number_to_trim):
                print names[sort_within_indices[i]] 

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

        return return_indices

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
                stdout = open(fileName, 'w')
                self.printInner(outFormat, stdout)
            except:
                print "Error diverting stout to file:", fileName, exc_info()[0]
                raise
        else:
            self.printInner(outFormat)           

    def printInner(self, outFormat, stream=sys_stdout):
        """Print bin information to STDOUT"""
        # handle the headers first
        separator = "\t"
        if(outFormat == 'summary'):
            stream.write(separator.join(["#\"bid\"","\"totalBP\"","\"numCons\"","\"cMean\"","\"cStdev\"","\"kMean\"","\"kStdev\""])+"\n") 
        elif(outFormat == 'minimal'):
            stream.write(separator.join(["#\"bid\"","\"cid\"","\"length\""])+"\n")            
        elif(outFormat == 'full'):
            pass
        else:
            print "Error: Unrecognised format:", outFormat
            return

        for bid in self.getBids():
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
            self.bins[bid].printBin(self.PM.contigNames, self.PM.contigLengths, outFormat=outFormat, separator=separator, stream=stream)

    def plotProfileDistributions(self):
        """Plot the coverage and kmer distributions for each bin"""
        for bid in self.getBids():
            self.bins[bid].plotProfileDistributions(self.PM.transformedCP, self.PM.kmerSigs, fileName="PROFILE_"+str(bid))

    def plotBins(self, FNPrefix="BIN", sideBySide=False, folder=''):
        """Make plots of all the bins"""
        if folder != '':
            makeSurePathExists(folder)
            
        for bid in self.getBids():
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
        if(sideBySide):
            print "Plotting side by side"
            self.plotSideBySide(self.bins.keys(), tag=FNPrefix)
        else:
            print "Plotting bins"
            for bid in self.getBids():
                if folder != '':
                    self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, fileName=osp_join(folder, FNPrefix+"_"+str(bid)))
                else:
                    self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, FNPrefix+"_"+str(bid))

    def plotSideBySide(self, bids, fileName="", tag=""):
        """Plot two bins side by side in 3d"""
        fig = plt.figure()
        # we need to work out how to shape the plots
        num_plots = len(bids)
        plot_rows = float(int(np_sqrt(num_plots)))
        plot_cols = np_ceil(float(num_plots)/plot_rows)
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
                print "Error saving image:", fileName, exc_info()[0]
                raise
        elif(show):
            try:
                plt.show()
            except:
                print "Error showing image:", exc_info()[0]
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
            print "Error showing image", exc_info()[0]
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
            print "Error showing image", exc_info()[0]
            raise
        del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

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
__version__ = "0.2.2"
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

from numpy import sum as np_sum, abs as np_abs, amax as np_amax, amin as np_amin, append as np_append, arccos as np_arccos, argmin as np_argmin, argsort as np_argsort, array as np_array, ceil as np_ceil, concatenate as np_concatenate, delete as np_delete, log10 as np_log10, max as np_max, mean as np_mean, median as np_median, min as np_min, pi as np_pi, reshape as np_reshape, seterr as np_seterr, size as np_size, sort as np_sort, sqrt as np_sqrt, std as np_std, where as np_where, zeros as np_zeros, cos as np_cos, sin as np_sin
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
                 loadContigLengths=True,
                 loadLinks=False,
                 loadContigNames=True,
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
                         loadContigNames=loadContigNames,
                         loadContigLengths=loadContigLengths,
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
                        if self.bins[putative_bid].binSize > 1:
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
            if len(t_list) != k:
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
            try:
                cbid = self.PM.binIds[row_index]
                if cbid != 0: # don't want these guys
                    try:
                        refined_t_list[cbid] += 1
                    except KeyError:
                        refined_t_list[cbid] = 1
            except KeyError:
                pass

        if verbose:
            print k, refined_t_list,
            
        # work out the most prominent BID
        max_bid = 0
        max_count = 0
        for cbid in refined_t_list:
            if refined_t_list[cbid] > max_count:
                max_count = refined_t_list[cbid]
                max_bid = cbid
        if verbose:
            print self.PM.contigNames[rowIndex], "**",max_bid,max_count,"**"                        
            for cbid in refined_t_list:
                print "[", cbid, ",", refined_t_list[cbid], "]",
            print
        # we're done!
        return (max_bid, neighbourList)

    def autoRefineBins(self, iterate=False, verbose=False):
        """Automagically refine bins"""
        super_round = 1
        tdm = self.PM.transformedCP
        neighbour_list={} # save looking things up a 1,000,000 times
        search_tree = kdt(tdm)
        
        # pay attention to contig lengths when inclusing in bins
        # use grubbs test
        GT = GrubbsTester() # we need to perform grubbs test before inclusion
        
        stable_bids = {} # once a bin is stable it's stable... almost
        re_unstable = {}
        while True:
            sr_contigs_reassigned = 0
            num_reassigned = -1
            round = 0
            while num_reassigned != 0:
                num_reassigned = 0
                reassignment_map = {}
                moved_RIs = {}
                
                # make a lookup of each bins contig length distributions
                bin_c_lengths = {}
                for bid in self.getBids():
                    self.bins[bid].makeBinDist(self.PM.transformedCP, 
                                               self.PM.averageCoverages, 
                                               self.PM.kmerVals, 
                                               self.PM.contigLengths)
                    for row_index in self.bins[bid].rowIndices:
                        try:
                            bin_c_lengths[bid].append(self.PM.contigLengths[row_index])
                        except KeyError:
                            bin_c_lengths[bid] = [self.PM.contigLengths[row_index]]

                # destabilize those who've picked up new contigs
                for bid in re_unstable.keys():
                    try:
                        del stable_bids[bid]
                    except KeyError: pass
                re_unstable = {} 
                
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
                                if GT.isMaxOutlier(self.PM.contigLengths[row_index], bin_c_lengths[assigned_bid], verbose=verbose):
                                    # undo the change
                                    assigned_bid = bid
                                else:
                                    stable = False
                                    re_unstable[assigned_bid] = True
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
                if True:#verbose:
                    print "    Refine round %d: reassigned %d contigs, removed %d cores" % (round, num_reassigned, bins_removed)
            print "    Refine round %d complete. (%d iterations) Total contigs reassigned: %d" % (super_round, round, sr_contigs_reassigned)
            if sr_contigs_reassigned == 0:
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

class GrubbsTester:
    """Data and methods for performing Grubbs test
    
    cutoff values taken from qgrubs from R package outliers
    using command: qgrubbs(0.95, c(3:1002), 10)
    """
    def __init__(self):
        # cutoff values for n degress of freedom 
        # If you have 8 sample points then self.cutoffs[6]
        # is what you want!
        self.critVs = np_array([1.153118,1.462500,1.671386,1.822120,1.938135,2.031652,2.109562,2.176068,
                                2.233908,2.284953,2.330540,2.371654,2.409038,2.443272,2.474810,2.504017,
                                2.531193,2.556581,2.580388,2.602784,2.623916,2.643910,2.662873,2.680899,
                                2.698071,2.714459,2.730127,2.745132,2.759523,2.773345,2.786639,2.799440,
                                2.811782,2.823693,2.835202,2.846331,2.857105,2.867542,2.877664,2.887485,
                                2.897023,2.906293,2.915308,2.924081,2.932623,2.940946,2.949060,2.956975,
                                2.964699,2.972240,2.979608,2.986808,2.993848,3.000735,3.007474,3.014072,
                                3.020533,3.026863,3.033067,3.039150,3.045115,3.050968,3.056711,3.062349,
                                3.067885,3.073323,3.078665,3.083916,3.089077,3.094152,3.099143,3.104053,
                                3.108885,3.113640,3.118321,3.122929,3.127468,3.131939,3.136344,3.140684,
                                3.144962,3.149179,3.153337,3.157437,3.161481,3.165470,3.169405,3.173289,
                                3.177122,3.180905,3.184640,3.188327,3.191968,3.195565,3.199117,3.202627,
                                3.206094,3.209520,3.212906,3.216253,3.219562,3.222832,3.226066,3.229264,
                                3.232427,3.235555,3.238649,3.241710,3.244738,3.247735,3.250700,3.253635,
                                3.256540,3.259415,3.262261,3.265079,3.267870,3.270632,3.273368,3.276078,
                                3.278762,3.281421,3.284054,3.286663,3.289248,3.291810,3.294348,3.296863,
                                3.299356,3.301827,3.304276,3.306704,3.309110,3.311496,3.313862,3.316208,
                                3.318534,3.320840,3.323128,3.325397,3.327647,3.329879,3.332093,3.334290,
                                3.336469,3.338631,3.340776,3.342905,3.345017,3.347113,3.349193,3.351258,
                                3.353307,3.355340,3.357359,3.359363,3.361352,3.363327,3.365288,3.367235,
                                3.369167,3.371087,3.372992,3.374885,3.376764,3.378631,3.380484,3.382325,
                                3.384154,3.385970,3.387774,3.389566,3.391347,3.393116,3.394873,3.396618,
                                3.398353,3.400076,3.401789,3.403491,3.405181,3.406862,3.408532,3.410191,
                                3.411840,3.413480,3.415109,3.416728,3.418338,3.419938,3.421528,3.423109,
                                3.424681,3.426244,3.427797,3.429341,3.430877,3.432404,3.433922,3.435431,
                                3.436932,3.438424,3.439908,3.441384,3.442851,3.444311,3.445762,3.447206,
                                3.448642,3.450070,3.451490,3.452903,3.454308,3.455706,3.457096,3.458479,
                                3.459855,3.461224,3.462586,3.463940,3.465288,3.466629,3.467963,3.469290,
                                3.470611,3.471925,3.473232,3.474533,3.475828,3.477116,3.478398,3.479674,
                                3.480943,3.482206,3.483464,3.484715,3.485960,3.487199,3.488433,3.489661,
                                3.490883,3.492099,3.493309,3.494514,3.495714,3.496908,3.498096,3.499279,
                                3.500457,3.501629,3.502797,3.503958,3.505115,3.506267,3.507413,3.508555,
                                3.509691,3.510823,3.511949,3.513071,3.514188,3.515300,3.516407,3.517510,
                                3.518608,3.519701,3.520790,3.521874,3.522953,3.524028,3.525099,3.526165,
                                3.527227,3.528284,3.529337,3.530386,3.531430,3.532471,3.533507,3.534539,
                                3.535567,3.536590,3.537610,3.538626,3.539637,3.540645,3.541649,3.542648,
                                3.543644,3.544636,3.545625,3.546609,3.547590,3.548567,3.549540,3.550509,
                                3.551475,3.552437,3.553396,3.554351,3.555303,3.556251,3.557195,3.558136,
                                3.559073,3.560007,3.560938,3.561865,3.562789,3.563710,3.564627,3.565541,
                                3.566452,3.567359,3.568263,3.569164,3.570062,3.570957,3.571848,3.572737,
                                3.573622,3.574504,3.575384,3.576260,3.577133,3.578003,3.578870,3.579734,
                                3.580596,3.581454,3.582309,3.583162,3.584012,3.584859,3.585703,3.586544,
                                3.587382,3.588218,3.589051,3.589881,3.590709,3.591533,3.592355,3.593175,
                                3.593992,3.594806,3.595617,3.596426,3.597232,3.598036,3.598837,3.599636,
                                3.600432,3.601226,3.602017,3.602805,3.603592,3.604375,3.605157,3.605935,
                                3.606712,3.607486,3.608258,3.609027,3.609794,3.610559,3.611321,3.612081,
                                3.612839,3.613594,3.614347,3.615098,3.615847,3.616593,3.617338,3.618080,
                                3.618819,3.619557,3.620293,3.621026,3.621757,3.622486,3.623213,3.623938,
                                3.624661,3.625381,3.626100,3.626816,3.627531,3.628243,3.628954,3.629662,
                                3.630368,3.631073,3.631775,3.632476,3.633174,3.633871,3.634565,3.635258,
                                3.635949,3.636637,3.637324,3.638009,3.638693,3.639374,3.640053,3.640731,
                                3.641407,3.642080,3.642753,3.643423,3.644091,3.644758,3.645423,3.646086,
                                3.646747,3.647407,3.648065,3.648721,3.649375,3.650028,3.650679,3.651328,
                                3.651976,3.652621,3.653266,3.653908,3.654549,3.655188,3.655826,3.656462,
                                3.657096,3.657729,3.658360,3.658989,3.659617,3.660243,3.660868,3.661491,
                                3.662112,3.662732,3.663351,3.663968,3.664583,3.665197,3.665809,3.666420,
                                3.667029,3.667637,3.668243,3.668848,3.669451,3.670053,3.670654,3.671253,
                                3.671850,3.672446,3.673041,3.673634,3.674226,3.674816,3.675405,3.675992,
                                3.676578,3.677163,3.677746,3.678328,3.678909,3.679488,3.680066,3.680642,
                                3.681218,3.681791,3.682364,3.682935,3.683505,3.684073,3.684641,3.685206,
                                3.685771,3.686334,3.686896,3.687457,3.688016,3.688574,3.689131,3.689687,
                                3.690241,3.690794,3.691346,3.691897,3.692446,3.692995,3.693541,3.694087,
                                3.694632,3.695175,3.695717,3.696258,3.696798,3.697336,3.697874,3.698410,
                                3.698945,3.699479,3.700011,3.700543,3.701073,3.701602,3.702130,3.702657,
                                3.703183,3.703708,3.704231,3.704754,3.705275,3.705795,3.706314,3.706832,
                                3.707349,3.707865,3.708380,3.708893,3.709406,3.709917,3.710428,3.710937,
                                3.711445,3.711953,3.712459,3.712964,3.713468,3.713971,3.714473,3.714974,
                                3.715474,3.715973,3.716471,3.716967,3.717463,3.717958,3.718452,3.718945,
                                3.719437,3.719928,3.720417,3.720906,3.721394,3.721881,3.722367,3.722852,
                                3.723336,3.723819,3.724301,3.724782,3.725263,3.725742,3.726220,3.726697,
                                3.727174,3.727649,3.728124,3.728597,3.729070,3.729542,3.730013,3.730483,
                                3.730952,3.731420,3.731887,3.732354,3.732819,3.733284,3.733747,3.734210,
                                3.734672,3.735133,3.735593,3.736052,3.736511,3.736968,3.737425,3.737881,
                                3.738336,3.738790,3.739243,3.739695,3.740147,3.740598,3.741048,3.741497,
                                3.741945,3.742392,3.742839,3.743284,3.743729,3.744173,3.744617,3.745059,
                                3.745501,3.745942,3.746382,3.746821,3.747259,3.747697,3.748134,3.748570,
                                3.749005,3.749440,3.749873,3.750306,3.750739,3.751170,3.751601,3.752030,
                                3.752459,3.752888,3.753315,3.753742,3.754168,3.754594,3.755018,3.755442,
                                3.755865,3.756287,3.756709,3.757130,3.757550,3.757969,3.758388,3.758806,
                                3.759223,3.759639,3.760055,3.760470,3.760884,3.761298,3.761711,3.762123,
                                3.762535,3.762945,3.763355,3.763765,3.764174,3.764581,3.764989,3.765395,
                                3.765801,3.766207,3.766611,3.767015,3.767418,3.767821,3.768223,3.768624,
                                3.769024,3.769424,3.769823,3.770222,3.770620,3.771017,3.771413,3.771809,
                                3.772204,3.772599,3.772993,3.773386,3.773779,3.774171,3.774562,3.774953,
                                3.775343,3.775732,3.776121,3.776509,3.776897,3.777284,3.777670,3.778056,
                                3.778441,3.778825,3.779209,3.779592,3.779975,3.780357,3.780738,3.781119,
                                3.781499,3.781879,3.782258,3.782636,3.783014,3.783391,3.783768,3.784144,
                                3.784519,3.784894,3.785268,3.785642,3.786015,3.786387,3.786759,3.787130,
                                3.787501,3.787871,3.788241,3.788610,3.788978,3.789346,3.789713,3.790080,
                                3.790446,3.790812,3.791177,3.791541,3.791905,3.792269,3.792632,3.792994,
                                3.793356,3.793717,3.794077,3.794437,3.794797,3.795156,3.795514,3.795872,
                                3.796230,3.796587,3.796943,3.797299,3.797654,3.798009,3.798363,3.798717,
                                3.799070,3.799422,3.799775,3.800126,3.800477,3.800828,3.801178,3.801527,
                                3.801876,3.802225,3.802573,3.802920,3.803267,3.803614,3.803960,3.804305,
                                3.804650,3.804995,3.805338,3.805682,3.806025,3.806367,3.806709,3.807051,
                                3.807392,3.807732,3.808072,3.808412,3.808751,3.809090,3.809428,3.809765,
                                3.810102,3.810439,3.810775,3.811111,3.811446,3.811781,3.812115,3.812449,
                                3.812782,3.813115,3.813448,3.813779,3.814111,3.814442,3.814772,3.815103,
                                3.815432,3.815761,3.816090,3.816418,3.816746,3.817073,3.817400,3.817727,
                                3.818053,3.818378,3.818703,3.819028,3.819352,3.819676,3.819999,3.820322,
                                3.820645,3.820967,3.821288,3.821610,3.821930,3.822250,3.822570,3.822890,
                                3.823209,3.823527,3.823845,3.824163,3.824480,3.824797,3.825114,3.825430,
                                3.825745,3.826060,3.826375,3.826689,3.827003,3.827317,3.827630,3.827942,
                                3.828255,3.828566,3.828878,3.829189,3.829499,3.829810,3.830119,3.830429,
                                3.830738,3.831046,3.831354,3.831662,3.831969,3.832276,3.832583,3.832889,
                                3.833195,3.833500,3.833805,3.834110,3.834414,3.834718,3.835021,3.835324,
                                3.835627,3.835929,3.836231,3.836532,3.836833,3.837134,3.837434,3.837734,
                                3.838034,3.838333,3.838632,3.838930,3.839228,3.839526,3.839823,3.840120,
                                3.840417,3.840713,3.841009,3.841304,3.841599,3.841894,3.842188,3.842482,
                                3.842776,3.843069,3.843362,3.843654,3.843946,3.844238,3.844529,3.844820,
                                3.845111,3.845401,3.845691,3.845981,3.846270,3.846559,3.846847,3.847136,
                                3.847423,3.847711,3.847998,3.848285,3.848571,3.848857,3.849143,3.849428,
                                3.849713,3.849998,3.850282,3.850566,3.850850,3.851133,3.851416,3.851699,
                                3.851981,3.852263,3.852545,3.852826,3.853107,3.853388,3.853668,3.853948,
                                3.854228,3.854507,3.854786,3.855064,3.855343,3.855621,3.855898,3.856175,
                                3.856452,3.856729,3.857005,3.857281,3.857557,3.857832,3.858107,3.858382,
                                3.858656,3.858930,3.859204,3.859478,3.859751,3.860023,3.860296,3.860568,
                                3.860840,3.861111,3.861383,3.861653,3.861924,3.862194,3.862464,3.862734,
                                3.863003,3.863272,3.863541,3.863809,3.864077,3.864345,3.864613,3.864880,
                                3.865147,3.865413,3.865679,3.865945,3.866211,3.866476,3.866741,3.867006,
                                3.867271,3.867535,3.867799,3.868062,3.868325,3.868588,3.868851,3.869113,
                                3.869375,3.869637,3.869899,3.870160,3.870421,3.870681,3.870942,3.871202,
                                3.871462,3.871721,3.871980,3.872239,3.872498,3.872756,3.873014,3.873272,
                                3.873529,3.873786,3.874043,3.874300,3.874556,3.874812,3.875068,3.875324,
                                3.875579,3.875834,3.876088,3.876343,3.876597,3.876851,3.877104,3.877357]
                               )

    def isMaxOutlier(self, maxVal, compVals, verbose=False):
        """Test if the maxVal is an outlier 
        
        maxVal should NOT be included in allVals
        if len(compVals) - 1 > 1000 use the 1000 cutoff anyway
        """
        # get Z score for the maxValue 
        v = (maxVal - np_mean(compVals+[maxVal]))/np_std(compVals+[maxVal], ddof=1)
        idx = len(compVals) - 1
        if idx > 999:
            idx = 999
        if verbose:
            print np_mean(compVals+[maxVal]), np_std(compVals+[maxVal], ddof=1), maxVal, v, idx, self.critVs[idx], v > self.critVs[idx]
        return v > self.critVs[idx] 

###############################################################################
###############################################################################
###############################################################################
###############################################################################
        
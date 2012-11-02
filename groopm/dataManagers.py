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
        
        # linking row indices to bids
        self.r2b = {}

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
                         loadCores=False,
                         loadLinks=loadLinks
                        )
        
        bin_members = self.initialiseContainers()
        if(makeBins):
            if transform:
                self.PM.transformCP(silent=silent, min=min, max=max)
            else:
                if self.PM.numStoits == 3:
                    self.PM.transformedCP = self.PM.covProfiles
                else:
                    print "Number of stoits != 3. You need to transform"
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
        for row_index in range(0, np_size(self.PM.indices)):
            bin_members[self.PM.binIds[row_index]].append(row_index)
            bin_sizes[self.PM.binIds[row_index]] += self.PM.contigLengths[row_index]

        # we need to get the largest BinId in use
        bids = self.PM.getBinStats().keys()
        if(len(bids) > 0):
            self.nextFreeBinId = np_max(bids)
        
        return bin_members

    def makeBins(self, binMembers):
        """Make bin objects from loaded data"""
        for bid in self.PM.validBinIds.keys():
            self.bins[bid] = Bin(np_array(binMembers[bid]), bid, self.PM.scaleFactor-1)
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
            bin_updates[bid] = np_size(self.bins[bid].rowIndices)
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
        for row_index in range(len(self.PM.indices)):
            if row_index in self.PM.binnedRowIndicies:
                core_update[self.PM.indices[row_index]] = True

        bin_update = self.getBinUpdates()

        return (bin_update, core_update)

    def makeR2BLookup(self):
        """Get a list of row indicies to bin ids"""
        self.r2b = {}
        bids = self.getBids()
        for bid in bids:
            bin = self.getBin(bid)
            for row_index in bin.rowIndices:
                self.r2b[row_index] = bid
                
        # do the unbinned guys too
        for i in range(len(self.PM.indices)):
            if i not in self.r2b:
                self.r2b[i] = 0

    def getBinUpdates(self, c2b={}):
        """Merge the bids, raw DB indexes and core information so we can save to disk"""
        # we need a mapping from cid (or local index) to binID
        bin_update = {}
        for row_index in c2b: # load any now-unbinned critters
            bin_update[self.PM.indices[row_index]] = c2b[row_index]
            
        for bid in self.getBids():
            for row_index in self.bins[bid].rowIndices:
                bin_update[self.PM.indices[row_index]] = bid
            
        return bin_update

    def contig2bids(self):
        """Map contigIDs (row indices) to bins"""
        c2b = {}
        for bid in self.getBids():
            bin = self.bins[bid]
            for row_index in bin.rowIndices:
                c2b[row_index] = bid
            
        return c2b
        

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

    def getAllContigTrueDistances(self, binNetwork={}, c2b={}):
        """work out the "true" distance between all closely situated contigs"""
        # work out which bins lie nearby
        if binNetwork == {}:
            binNetwork = self.findBinNeighbours()
        
        # work out which contigs are in which bins
        if c2b == {}:
            c2b = self.contig2bids()
        
        all_distances = {} # hash of type cid -> [(cid, dist), (cid,dist), ... ]
        for row_index in range(len(self.PM.indices)):
            self.getContigTrueDistances(row_index, binNetwork, c2b, all_distances)
        return all_distances 

    def getContigTrueDistances(self, rowIndex, binNetwork, c2b, distances):
        """work out the "true" distance between this contig and all closely situated contigs"""
        try:
            this_bid = c2b[rowIndex]
        except KeyError:
            this_bid = 0
            pass
        if this_bid != 0:
            pass
        else:
            # we will need to generate a list of close bins 
            # using some other method
            pass

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
        sorted_Xs = np_argsort(Xs)
        sorted_Ys = np_argsort(Ys)
        sorted_Zs = np_argsort(Zs)

        x_dists = {}
        for i in range(max_index):
            for j in range(i + 1, max_index):
                dist = np_abs(Xs[sorted_Xs[i]] - Xs[sorted_Xs[j]])
                if(dist <= maxDist):
                    key = self.makeBidKey(bids[sorted_Xs[i]], bids[sorted_Xs[j]])
                    x_dists[key] = dist
                else:
                    break 
        y_dists = {}
        for i in range(max_index):
            for j in range(i + 1, max_index):
                dist = np_abs(Ys[sorted_Ys[i]] - Ys[sorted_Ys[j]])
                if(dist <= maxDist):
                    key = self.makeBidKey(bids[sorted_Ys[i]], bids[sorted_Ys[j]])
                    if(key in x_dists):
                        y_dists[key] = dist
                else:
                    break 
        actual_dists = {}
        for i in range(max_index):
            for j in range(i + 1, max_index):
                dist = np_abs(Zs[sorted_Zs[i]] - Zs[sorted_Zs[j]])
                if(dist <= maxDist):
                    key = self.makeBidKey(bids[sorted_Zs[i]], bids[sorted_Zs[j]])
                    if(key in y_dists):
                        dist = np_sqrt(y_dists[key]**2 + x_dists[key]**2 + dist**2)
                        if(dist <= maxDist):
                            actual_dists[key] = dist
                else:
                    break 

        return actual_dists


#------------------------------------------------------------------------------
# LINKS

    def getLinkingContigs(self, bid):
        """Get all contigs and their bin IDs which link to contigs in this bin"""
        condition = ""
        bin = self.getBin(bid)
        bin2count = {}
        for row_index in bin.rowIndices:
            try: 
                for link in self.PM.links[row_index]:
                    link_bid = self.r2b[link[0]] 
                    if link_bid != bid and link_bid != 0:
                        try: 
                            bin2count[link_bid] += 1.0
                        except KeyError:
                            bin2count[link_bid] = 1.0
            except KeyError:
                pass
        return bin2count
    
    def getConnectedBins(self, rowIndex, c2b):
        """Get a  list of bins connected to this contig"""
        ret_links = []
        for link in self.PM.links[rowIndex]:
            cid = link[0]
            try:
                bid = c2b[cid]
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
        self.makeR2BLookup()
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
                    link_bid = self.r2b[link[0]] 
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
            self.saveBins(doCores=False, saveBinStats=False, updateBinStats=True)
    
    def refineWrapper(self,
                      manual=False,          # do we need to ask permission every time?
                      save=False,
                      plotter=False,
                      shuffle=False,
                      ):
        """Iterative wrapper for the refine function"""
        if(plotter):
            self.plotterRefineBins()
        if(shuffle):
            self.autoRefineBins()
            self.saveBins(doCores=False, saveBinStats=True, updateBinStats=True)
                        
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

    def getClosestBID(self, rowIndex, searchTree, tdm, c2b, k=101, verbose=False):
        """Find the bin ID which would best describe the placement of the contig"""
        # find the k nearest neighbours for the query contig
        if k < 21:
            k = 21
        if k > 101:
            k = 101
        t_list = searchTree.query(tdm[rowIndex],k=k)[1]
        # calculate the distribution of neighbouring bins
        refined_t_list = {}
        for row_index in t_list:
            try:
                cbid = c2b[row_index]
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
        return max_bid

    def autoRefineBins(self, verbose=False):
        """Automagically refine bins"""
        num_reassigned = -1
        round = 0
        stable_bids = {} # once a bin is stable it's stable!
        tdm = np_append(self.PM.transformedCP, 1000*np_reshape(self.PM.kmerVals,(len(self.PM.kmerVals),1)),1)
        search_tree = kdt(tdm)
        while num_reassigned != 0:
            num_reassigned = 0
            reassignment_map = {}
            c2b = self.contig2bids()
            bids = self.getBids()
            calls = 0
            for bid in bids:
                bin = self.getBin(bid)
                if bid in stable_bids:
                    try:
                        reassignment_map[bid] += list(bin.rowIndices)
                    except KeyError:
                        reassignment_map[bid] = list(bin.rowIndices)
                else:
                    stable = True
                    #print "BID:", bid, bin.binSize 
                    for row_index in bin.rowIndices:
                        calls += 1
                        assigned_bid = self.getClosestBID(row_index, search_tree, tdm, c2b, verbose=verbose, k=2*bin.binSize-1)                        
                        if assigned_bid != bid:
                            stable = False
                            num_reassigned += 1
                        
                        # keep track of where this guy lives
                        try:
                            reassignment_map[assigned_bid].append(row_index)
                        except KeyError:
                            reassignment_map[assigned_bid] = [row_index]
                    if stable: # no changes this round, mark bin as stable
                        stable_bids[bid] = True
            
            print "Calls:", calls
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
            print "    Refine round %d: reassigned %d contigs, removed %d cores" % (round, num_reassigned, bins_removed)

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
        obs = np_array([])
        if(mode=='kmer'):
            obs = np_array([self.PM.kmerVals[i] for i in bin.rowIndices])
        elif(mode=='cov'):
            obs = np_array([self.PM.covProfiles[i] for i in bin.rowIndices])
        
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
        idx_sorted = np_argsort(np_array(idx))
        current_group = 0
        bids = [bid]
        bin_stats = {} # bin id to bin size
        bin_stats[bid]=0 # this will ensure that the old bin id will be deleted!
        bin_update = {} # row index to bin id
        holding_array = np_array([])
        split_bin = None
        for i in idx_sorted:
            if(idx[i] != current_group):
                # bin is full!
                split_bin = self.makeNewBin(holding_array)
                for row_index in holding_array:
                    bin_update[self.PM.indices[row_index]] = split_bin.id
                bin_stats[split_bin.id] = split_bin.binSize  
                split_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
                bids.append(split_bin.id)
                holding_array = np_array([])
                current_group = idx[i]
            holding_array = np_append(holding_array, bin.rowIndices[i])
        # do the last one
        if(np_size(holding_array) != 0):
            split_bin = self.makeNewBin(holding_array)
            for row_index in holding_array:
                bin_update[self.PM.indices[row_index]] = split_bin.id  
            bin_stats[split_bin.id] = split_bin.binSize  
            split_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)
            bids.append(split_bin.id)

        return (bin_stats, bin_update, bids)

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
                bin_stats[bids[i]] = 0
                bin_stats[parent_bin.id] = parent_bin.binSize
                some_merged = True

        if(saveBins and some_merged):
            self.updateBinStats(bin_stats)
            self.saveBins(doCores=False, saveBinStats=False)
            
        # now fix up the r2b indices
        if self.r2b != {}:
            parent_bid = parent_bin.id
            for row_index in self.r2b:
                try:
                    self.r2b[row_index] = parent_bid
                except KeyError:
                    pass 
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
                    for row_index in self.bins[bid].rowIndices:
                        if row_index in self.PM.binnedRowIndicies:
                            del self.PM.binnedRowIndicies[row_index]
                        else:
                            print bid, row_index, "FUNG"
                        bin_update[self.PM.indices[row_index]] = 0 
                bin_stats[bid] = 0
                del self.bins[bid]
            else:
                raise ge.BinNotFoundException("Cannot find: "+str(bid)+" in bins dicts")
            
        if(saveBins):
            self.updateBinStats(bin_stats)
            self.PM.saveBinIds(bin_update)
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
                self.printInner(outFormat)
            except:
                print "Error diverting stout to file:", fileName, exc_info()[0]
                raise
        else:
            self.printInner(outFormat)           

    def printInner(self, outFormat):
        """Print bin information to STDOUT"""
        # handle the headers first
        separator = "\t"
        if(outFormat == 'summary'):
            print separator.join(["#\"bid\"","\"totalBP\"","\"numCons\"","\"cMean\"","\"cStdev\"","\"kMean\"","\"kStdev\""]) 
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
            print "Plotting side by side"
            self.plotSideBySide(self.bins.keys(), tag=FNPrefix)
        else:
            print "Plotting bins"
            for bid in self.getBids():
                self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, fileName=FNPrefix+"_"+str(bid))

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
        self.isCore = np_array([])          # True False values
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
                 loadCores=False,
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
                if(len(bids) != 0): # need to make sure we're not restricted in terms of bins
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

            if(loadCores):
                if(verbose):
                    print "    Loading core info"
                self.isCore = self.dataManager.getCores(self.dbFileName, indices=self.indices)
                
            if(loadLinks):
                self.links = self.getLinks()
            
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
        self.isCore = np_delete(self.isCore, deadRowIndicies, axis=0)
        
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

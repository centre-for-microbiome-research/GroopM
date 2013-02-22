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

from numpy import shape as np_shape, around as np_around, argmax as np_argmax, arccos as np_cos, dot as np_dot, sum as np_sum, abs as np_abs, amax as np_amax, amin as np_amin, append as np_append, arccos as np_arccos, argmin as np_argmin, argsort as np_argsort, array as np_array, ceil as np_ceil, concatenate as np_concatenate, delete as np_delete, log10 as np_log10, max as np_max, mean as np_mean, median as np_median, min as np_min, pi as np_pi, reshape as np_reshape, seterr as np_seterr, size as np_size, sort as np_sort, sqrt as np_sqrt, std as np_std, where as np_where, zeros as np_zeros, cos as np_cos, sin as np_sin
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
from ellipsoid import EllipsoidTool

np_seterr(all='raise')     

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BinManager:
    """Class used for manipulating bins"""
    def __init__(self, dbFileName="", pm=None, minSize=10, minVol=1000000):
        # data storage
        if(dbFileName != ""):
            self.PM = ProfileManager(dbFileName)
        elif(pm is not None):
            self.PM = pm
        
        # all about bins
        self.nextFreeBinId = 0                      # increment before use!
        self.bins = {}                              # bid -> Bin
        
        # misc
        self.minSize=minSize           # Min number of contigs for a bin to be considered legit
        self.minVol=minVol             # Override on the min size, if we have this many BP
        

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

    def saveBins(self, binAssignments={}, nuke=False):
        """Save binning results
        
        binAssignments is a hash of LOCAL row indices Vs bin ids
        { row_index : bid } 
        PM.setBinAssignments needs GLOBAL row indices
        
        We always overwrite the bins table (It is smallish)
        """
        if nuke:
            self.PM.dataManager.nukeBins(self.PM.dbFileName)

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
                      links=False,
                      ignoreRanges=False
                      ):
        """Iterative wrapper for the refine function"""
        if plotter:
            self.plotterRefineBins(ignoreRanges=ignoreRanges)
        if shuffle:
            print "Start automatic bin refinement"
            self.autoRefineBins(iterate=True)
            num_binned = len(self.PM.binnedRowIndicies.keys())
            print "   ",num_binned,"contigs across",len(self.bins.keys()),"cores"
            
            if saveBins:
                self.saveBins(nuke=True)

    def plotterRefineBins(self, ignoreRanges=False):
        """combine similar bins using 3d plots"""
        ET = EllipsoidTool()
        self.printRefinePlotterInstructions()
        self.plotBinIds(ignoreRanges=ignoreRanges)
        continue_merge = True
        while(continue_merge):
            user_option = self.promptOnPlotterRefine()
            if(user_option == 'R'):
                self.plotBinIds(ignoreRanges=ignoreRanges)
            
            elif(user_option == 'P'):
                self.plotBinPoints(ignoreRanges=ignoreRanges)
            
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
                        krange = int(raw_input(" Enter kmer range (0-9): "))
                        have_range = True
                    except ValueError:
                        print "You need to enter an integer value!"
                self.plotBinIds(krange=krange, ignoreRanges=ignoreRanges)
                
            elif(user_option == 'B'):
                # print subset of bins
                have_bid = False
                bids = []
                while(not have_bid):
                    have_bid = True
                    try:
                        usr_bids = raw_input(" Enter bid(s) to plot: ")
                        bids = [int(i) for i in usr_bids.split(" ")]
                        if bids == [-1]:
                            bids = self.getBids()
                        else:
                            for bid in bids:
                                if bid not in self.bins:
                                    print "ERROR: Bin %d not found!" % bid
                                    have_bid &= False
                    except ValueError:
                        print "You need to enter an integer value!"

                self.plotSelectBins(bids, plotMers=True, ET=ET)
                
            elif(user_option == 'S'):
                # split bins
                have_bid = False
                have_parts = False
                while(not have_bid):
                    try:
                        bid = int(raw_input(" Enter bid to split: "))
                        if bid not in self.bins:
                            print "ERROR: Bin %d not found!" % bid
                        else:
                            have_bid = True
                    except ValueError:
                        print "You need to enter an integer value!"
                while(not have_parts):
                    try:
                        parts = int(raw_input(" Enter number of parts to split into: "))
                        if(parts < 2):
                            print "ERROR: Need to choose 2 or more parts"
                        else:
                            have_parts = True
                    except ValueError:
                        print "You need to enter an integer value!"
                self.split(bid, parts, mode='kmer', auto=False, saveBins=True, printInstructions=False)
            elif(user_option == 'V'):
                """plot in vicinity of a bin"""
                have_bid = False
                have_radius = False
                while(not have_bid):
                    try:
                        bid = int(raw_input(" Enter bid of interest: "))
                        if bid not in self.bins:
                            print "ERROR: Bin %d not found!" % bid
                        else:
                            have_bid = True
                    except ValueError:
                        print "You need to enter an integer value!"
                while(not have_radius):
                    try:
                        usr_radius = raw_input(" Enter radius to select from [default 100]: ")
                        if usr_radius == "":
                            radius = 100
                        else:
                            radius = int(usr_radius)
                        have_radius = True
                    except ValueError:
                        print "You need to enter an integer value!"
                
                # we need to find all points in an area about the centroid of
                # this bin
                self.bins[bid].makeBinDist(self.PM.transformedCP,
                                           self.PM.averageCoverages,
                                           self.PM.kmerVals,
                                           self.PM.contigLengths)
                # now go point shopping
                disp_vals = np_array([])
                disp_cols = np_array([])
                disp_lens = np_array([])
                num_points = 0
                seen_bids = {}
                for row_index in range(len(self.PM.indices)):
                    if np_norm(self.PM.transformedCP[row_index] - 
                              self.bins[bid].covMeans) <= 100:
                        num_points += 1
                        disp_vals = np_append(disp_vals, self.PM.transformedCP[row_index])
                        disp_lens = np_append(disp_lens, np_sqrt(self.PM.contigLengths[row_index]))
                        disp_cols = np_append(disp_cols, self.PM.contigColours[row_index])
                        try:
                            seen_bids[self.PM.binIds[row_index]].append(1)
                        except KeyError:
                            seen_bids[self.PM.binIds[row_index]] = [1]
                            
                # reshape
                disp_vals = np_reshape(disp_vals, (num_points, 3))
                disp_cols = np_reshape(disp_cols, (num_points, 3))
                
                print " Points are located in bins:"
                for seen_bid in seen_bids:
                    print "    %d - %d occurances" % (seen_bid, len(seen_bids[seen_bid]))
                
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1, projection='3d')
                ax.scatter(disp_vals[:,0],
                           disp_vals[:,1],
                           disp_vals[:,2],
                           edgecolors=disp_cols,
                           c=disp_cols,
                           s=disp_lens,
                           marker='.')
                self.bins[bid].plotOnAx(ax,
                               self.PM.transformedCP,
                               self.PM.contigColours,
                               self.PM.contigLengths,
                               ET=ET)
                try:
                    plt.show()
                except:
                    print "Error showing image:", sys.exc_info()[0]
                    raise
                plt.close(fig)
                del fig
            else:
                return
        return

    def getClosestBID(self,
                      searchIndex,
                      searchTree,
                      tdm,
                      gt2riLookup,
                      neighbourList={},
                      k=101,
                      verbose=False
                      ):
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
        if k > len(tdm):
            k = len(tdm)

        try:
            t_list = neighbourList[searchIndex]
            if len(t_list) != k:
                # must have made a change, we'll fix it here
                t_list = searchTree.query(tdm[searchIndex],k=k)[1]
                neighbourList[searchIndex] = t_list
        except KeyError:
            # first time, we need to do the search
            t_list = searchTree.query(tdm[searchIndex],k=k)[1]
            neighbourList[searchIndex] = t_list
            
        # calculate the distribution of neighbouring bins
        refined_t_list = {}
        for index in t_list:
            row_index = gt2riLookup[index]
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
            print self.PM.contigNames[searchIndex], "**",max_bid,max_count,"**"                        
            for cbid in refined_t_list:
                print "[", cbid, ",", refined_t_list[cbid], "]",
            print
        # we're done!
        return (max_bid, neighbourList, refined_t_list)

    def nukeOutliers(self, verbose=False):
        """Identify and remove small bins which contain mixed genomes
        
        Uses Grubbs testing to identify outliers
        """
        print "    Identifying possible chimeric cores"
        bids = self.getBids()
        # first we need to build a distribution!
        kval_stdev_distrb = [] 
        for bid in bids:
            self.bins[bid].makeBinDist(self.PM.transformedCP, 
                                       self.PM.averageCoverages, 
                                       self.PM.kmerVals, 
                                       self.PM.contigLengths)
            kval_stdev_distrb.append(self.bins[bid].kValStdev)
        
        # now we work out the distribution of stdevs
        stdstd = np_std(kval_stdev_distrb)
        stdmean = np_mean(kval_stdev_distrb)
        dead_bins = []
        for bid in bids:
            Z = (self.bins[bid].kValStdev - stdmean)/stdstd
            if Z > 2 and self.bins[bid].totalBP < 100000:
                 dead_bins.append(bid)
                 
        # delete the bad bins
        self.deleteBins(dead_bins,
                        force=True,
                        freeBinnedRowIndicies=True,
                        saveBins=False)
        print "    Identified %d possible chimeras leaving %d cores" % (len(dead_bins), len(self.bins))

    def makeUpperLower(self, vals, tol, decimals=0):
        """ make upper and lower cutoffs"""
        mean = np_mean(vals)
        try:
            stdev = np_std(vals)
        except FloatingPointError:
            stdev = 0
        if decimals == 0:
            return (mean + tol*stdev, mean - tol*stdev)
        else:
            return (np_around(mean + tol*stdev,decimals=decimals),
                    np_around(mean - tol*stdev,decimals=decimals))

    def doObviousMergers(self, verbose=False):
        """Merge bins which are just crying out to be merged!"""
        # use very strict params
        print "    Making obvious mergers"
        # pay attention to contig lengths when inclusing in bins
        # use grubbs test
        GT = GrubbsTester() # we need to perform grubbs test before inclusion
        
        tdm = []                # these are used in the neighbor search
        bid_2_tdm_index = {} 
        tdm_index_2_bid = {} 

        bids = self.getBids()
        K = 10                   # number of neighbours to test
        if K > len(bids):
            K = len(bids)

        # keep track of what gets merged where
        merged_bins = {}        # oldId => newId
        processed_pairs = {}    # keep track of pairs we've analysed
        bin_c_lengths = {}      # bid => [len,len,...]
        angles = {}             # store angles between contigs
        ave_covs = {}           # store average coverages

        index = 0
        for bid in bids:
            self.bins[bid].makeBinDist(self.PM.transformedCP, 
                                       self.PM.averageCoverages, 
                                       self.PM.kmerVals, 
                                       self.PM.contigLengths,
                                       merTol=2.5)  # make the merTol a little larger...
            # use a 4 dimensional vector [cov, cov, cov, mer]
            tdm.append(np_append(self.bins[bid].covMeans, [1000*self.bins[bid].kValMean]))
            
            bin_c_lengths[bid] = [self.PM.contigLengths[row_index] for row_index in self.bins[bid].rowIndices]
            
            bid_2_tdm_index[bid] = index
            tdm_index_2_bid[index] = bid
            index += 1

            # we will not process any pair twice.
            # we also wish to avoid checking if a bin will merge with itself        
            processed_pairs[self.makeBidKey(bid, bid)] = True

#-----
# ALL Vs ALL
        
        # make a search tree
        search_tree = kdt(tdm)
        tdm = np_array(tdm)
        for bid in bids:
            # get the base bid and trace the chain
            # up through mergers...
            merged_base_bid = bid
            base_bid = bid
            while merged_base_bid in merged_bins:
                merged_base_bid = merged_bins[merged_base_bid]
            base_bin = self.bins[base_bid]
            
            # get the K closest bins
            neighbor_list = [tdm_index_2_bid[i] for i in search_tree.query(tdm[bid_2_tdm_index[bid]],k=K)[1]]

            if verbose:
                print "++++++++++"
                print bid, neighbor_list

            # calculate the angles between all the contigs in
            # the base bin
            # we'd like to calculate the angles once only      
            ang_tol = 1.5  
            v_min_default = 1000000
            try:
                (base_angles, v_min_base) = angles[base_bid]
            except KeyError:
                base_angles = self.calculateAngles(base_bin.rowIndices, base_bin.rowIndices)
                try:                    
                    v_min_base = np_mean(base_angles) + ang_tol * np_std(base_angles)
                except FloatingPointError:     
                    v_min_base = v_min_default
                angles[base_bid] = (base_angles, v_min_base)
            
            # test each neighbor in turn
            for i in range(1,K):
                # get the query bid and trace the chain
                # up through mergers...
                should_merge = False     # we'll test to see if we should change this value
                query_bid = neighbor_list[i]
                merged_query_bid = query_bid
                while merged_query_bid in merged_bins:
                    merged_query_bid = merged_bins[merged_query_bid]
                    
                seen_key = self.makeBidKey(base_bid, query_bid)
                if(seen_key not in processed_pairs and
                   merged_base_bid != merged_query_bid):
                    if verbose:
                        print "-----\nbegin comparison %d (%d) %d (%d)" % (base_bid, merged_base_bid, query_bid, merged_query_bid)
                    processed_pairs[seen_key] = True
                    query_bin = self.bins[query_bid]
#-----
# CHECK ANGLES
                    # calculate the angles between all the contigs
                    # in the query bin if not already calcuated
                    try:
                        (query_angles, v_min_query) = angles[query_bid]
                    except KeyError:
                        query_angles = self.calculateAngles(query_bin.rowIndices, query_bin.rowIndices)
                        try:                    
                            v_min_query = np_mean(query_angles) + ang_tol * np_std(query_angles)
                        except FloatingPointError:     
                            v_min_query = v_min_default
                        angles[query_bid] = (query_angles, v_min_query) 
                    
                    # now calculate the angles between the base and query...
                    combined_angles = self.calculateAngles(base_bin.rowIndices, query_bin.rowIndices)

                    # work out the lowest variance in angle distribution
                    v_min = np_min([v_min_query, v_min_base])
                    if v_min == v_min_default:    # both containing single contigs!
                        should_merge = False                      

                    # work out how many of the combined angles are within the 
                    # variance
                    gooduns = len([i for i in combined_angles if i < v_min])

                    if (2*gooduns) >= len(combined_angles):
                        # these guys are in close proximity
                        should_merge = True
                    if verbose:
                        print "ANGLES vmin: %f #good: %d #comb: %d" % (v_min, gooduns, len(combined_angles)), should_merge
#-----
# CHECK AVERAGE COVERAGES
                    if should_merge:
                        should_merge = False
                        # get the average coverages for each bin
                        ave_cov_tol = 3
                        try:
                            (base_AC, base_AC_upper, base_AC_lower) = ave_covs[base_bid]
                        except KeyError:
                            base_AC = [self.PM.averageCoverages[i] for i in base_bin.rowIndices]
                            (base_AC_upper, base_AC_lower) = self.makeUpperLower(base_AC, ave_cov_tol, decimals=3)
                            ave_covs[base_bid] = (base_AC, base_AC_upper, base_AC_lower)
                        try:
                            (query_AC, query_AC_upper, query_AC_lower) = ave_covs[query_bid]
                        except KeyError:
                            query_AC = [self.PM.averageCoverages[i] for i in query_bin.rowIndices]
                            (query_AC_upper, query_AC_lower) = self.makeUpperLower(query_AC, ave_cov_tol, decimals=3)
                            ave_covs[query_bid] = (query_AC, query_AC_upper, query_AC_lower)

                        # make sure that the majority of the average coverages are kinda similar
                        gooduns1 = len([i for i in base_AC if (i <= query_AC_upper and i >= query_AC_lower)])
                        gooduns2 = -1
                        if (3*gooduns1) >= base_bin.binSize:
                             gooduns2 = len([i for i in query_AC if (i <= base_AC_upper and i >= base_AC_lower)])
                             if (3*gooduns2) >= query_bin.binSize:
                                 should_merge = True

                        if verbose:
                            print "AVE COV q: %d (%d) b: %d (%d)" % (gooduns1, base_bin.binSize, gooduns2, query_bin.binSize), should_merge
                            if should_merge == False:
                                print base_AC_lower, query_bin.cValMean, base_AC_upper, query_AC_lower, base_bin.cValMean, query_AC_upper
                                
                                if(query_bin.cValMean >= base_AC_lower and
                                   query_bin.cValMean <= base_AC_upper and
                                   base_bin.cValMean >= query_AC_lower and
                                   base_bin.cValMean <= query_AC_upper):
                                    should_merge = True
                    else:
                        continue
#-----
# CHECK LENGTHS
                    if should_merge:
                        # check to see that the contig lengths 
                        # from both bins are somewhat similar
                        #
                        # Test the smaller bin against the larger
                        if query_bin.binSize < base_bin.binSize:
                            lengths_wrong = GT.isMaxOutlier(np_median(bin_c_lengths[query_bid]),
                                                            bin_c_lengths[base_bid]
                                                            )
                        else:
                            lengths_wrong = GT.isMaxOutlier(np_median(bin_c_lengths[base_bid]),
                                                            bin_c_lengths[query_bid]
                                                            )
                        should_merge = not lengths_wrong
                        if verbose:
                            print "LENGTHS", should_merge 
                    else:
                        continue
#-----
# CHECK MER SIGS
                    if should_merge:
                        should_merge=False      
                        # we just need to check that the kmer sigs make sense
                        if(query_bin.kValMean <= base_bin.kValUpperLimit and 
                           query_bin.kValMean >= base_bin.kValLowerLimit and
                           base_bin.kValMean <= query_bin.kValUpperLimit and 
                           base_bin.kValMean >= query_bin.kValLowerLimit):
                            should_merge = True
                        if verbose:
                            print "MER1", should_merge, "[", base_bin.kValLowerLimit, query_bin.kValMean, base_bin.kValUpperLimit, "] [", query_bin.kValLowerLimit, base_bin.kValMean, query_bin.kValUpperLimit, "]" 
                    else:
                        continue

                    if should_merge:
                        should_merge=False
                        # the means make sense. Now to check the 
                        # individual contigs
                        upper_lim = query_bin.kValMean + query_bin.kValStdev
                        lower_lim = query_bin.kValMean - query_bin.kValStdev
                        num_OK = len([row_index for row_index in base_bin.rowIndices
                                      if (self.PM.kmerVals[row_index] >= lower_lim and
                                          self.PM.kmerVals[row_index] <= upper_lim)
                                      ])
                        if verbose:
                            print "MER - base %d (%d)" % (num_OK, base_bin.binSize)
                                
                        if (3*num_OK) >= base_bin.binSize:
                            upper_lim = base_bin.kValMean + base_bin.kValStdev
                            lower_lim = base_bin.kValMean - base_bin.kValStdev
                            num_OK = len([row_index for row_index in query_bin.rowIndices
                                          if (self.PM.kmerVals[row_index] >= lower_lim and
                                              self.PM.kmerVals[row_index] <= upper_lim)
                                          ])
                            if verbose:
                                print "MER - base %d (%d)" % (num_OK, query_bin.binSize)
                                    
                            if (3*num_OK) >= query_bin.binSize:
                                should_merge=True
                    else:
                        continue

                    if should_merge:
                        if merged_query_bid < merged_base_bid:
                            if verbose:
                                print "MERGE!", merged_base_bid, "=>", merged_query_bid
                            merged_bins[merged_base_bid] = merged_query_bid
                            # we just nuked the base bid
                            break
                        else:
                            if verbose:
                                print "MERGE!", merged_query_bid, "=>", merged_base_bid
                            merged_bins[merged_query_bid] = merged_base_bid

#-----
# MERGE
        # now make a bunch of mergers
        mergers = []
        processed_bids = {}     # bid => index in mergers
        for bid in merged_bins:
            
            # trace this guy until the end
            merging_bid = merged_bins[bid]
            while merging_bid in merged_bins:
                merging_bid = merged_bins[merging_bid]
            
            try:
                merge_list_id_1 = processed_bids[bid]
            except KeyError:
                merge_list_id_1 = -1
            try:
                merge_list_id_2 = processed_bids[merging_bid]
            except KeyError:
                merge_list_id_2 = -1
            
            if merge_list_id_1 == -1:
                if merge_list_id_2 == -1:
                    # all new
                    index = len(mergers)
                    processed_bids[merging_bid] = index
                    processed_bids[bid] = index
                    mergers.append([bid, merging_bid])
                else:
                    processed_bids[bid] = merge_list_id_2
                    mergers[merge_list_id_2].append(bid)
            elif merge_list_id_2 == -1:
                processed_bids[merging_bid] = merge_list_id_1
                mergers[merge_list_id_1].append(merging_bid)
            else:
                print "gonna hurt!", bid, merging_bid, merged_bins[bid], mergers[merge_list_id_1], mergers[merge_list_id_2]

        # now merge them
        num_bins_removed = 0
        for merge in mergers:
            if verbose:
                print merge
            num_bins_removed += (len(merge) - 1)
            self.merge(merge, auto=True, newBid=False, saveBins=False, verbose=False, printInstructions=False)

        print "    Removed %d cores leaving %d cores" % (num_bins_removed, len(self.bins))        
        
        return len(mergers)


    def mergeObvious(self, verbose=False):
        """Merge bins which are just crying out to be merged!"""
        
        print "    Making obvious mergers"
        # pay attention to contig lengths when inclusing in bins
        # use grubbs test
        GT = GrubbsTester() # we need to perform grubbs test before inclusion

        # Use overlapping ellipsoids to determine co-locality        
        ET = EllipsoidTool()

        tdm = []                # these are used in the neighbor search
        bid_2_tdm_index = {} 
        tdm_index_2_bid = {} 

        bids = self.getBids()
        K = 10                   # number of neighbours to test
        if K > len(bids):
            K = len(bids)

        # keep track of what gets merged where
        merged_bins = {}        # oldId => newId
        processed_pairs = {}    # keep track of pairs we've analysed
        bin_c_lengths = {}      # bid => [len,len,...]
        bin_c_ellipsoid_volumes = {}    # volume of the minimum bounding COVERAGE ellipsoid
        bin_c_ellipsoids = {}           # the matrix A representing the bins COVERAGE ellipsiod
        bin_k_ellipse_areas = {}        # area of the minimum bounding KMER ellipse
        bin_k_ellipses = {}             # the matrix A representing the bins KMER ellipse

#-----
# PREP DATA STRUCTURES

        # sort all contigs in ascending order of kSigPC
        sorted_Ks = np_argsort(self.PM.kmerVals)
        seen_bids = {}
        index = 0
        for RI in sorted_Ks:
            try:
                seen_bids[self.PM.binIds[RI]][1] = index
            except KeyError:
                seen_bids[self.PM.binIds[RI]] = [index, index]
            index += 1 

        index = 0
        for bid in bids:
            bin = self.bins[bid]
            bin.makeBinDist(self.PM.transformedCP,
                            self.PM.averageCoverages, 
                            self.PM.kmerVals, 
                            self.PM.contigLengths,
                            merTol=2.5)  # make the merTol a little larger...
            # use a 4 dimensional vector [cov, cov, cov, mer]
            tdm.append(np_append(bin.covMeans, [1000*bin.kValMean]))

            # work out the volume of the minimum bounding coverage ellipsoid and kmer ellipse
            (bin_c_ellipsoids[bid], bin_c_ellipsoid_volumes[bid]) = bin.getBoundingCEllipsoidVol(self.PM.transformedCP, ET=ET, retA=True)
            BP = np_array(zip([self.PM.kmerVals[i] for i in bin.rowIndices],
                              [self.PM.kmerVals2[i] for i in bin.rowIndices])
                          )
            (bin_k_ellipses[bid], bin_k_ellipse_areas[bid]) = bin.getBoundingKEllipseArea(BP,
                                                                                          ET=ET,
                                                                                          retA=True)
            # store the start / end location of this guy on the kmer scale
            bin.lowestK = seen_bids[bid][0]
            bin.highestK = seen_bids[bid][1]

            bin_c_lengths[bid] = [self.PM.contigLengths[row_index] for row_index in bin.rowIndices]
            
            bid_2_tdm_index[bid] = index
            tdm_index_2_bid[index] = bid
            index += 1

            # we will not process any pair twice.
            # we also wish to avoid checking if a bin will merge with itself        
            processed_pairs[self.makeBidKey(bid, bid)] = True

#-----
# ALL Vs ALL
        
        # make a search tree
        search_tree = kdt(tdm)
        tdm = np_array(tdm)
        for bid in bids:
            #verbose = bid == 23 or bid == 24 or bid == 48 or bid == 49 or bid == 50 or bid == 30 or bid == 66

            # get the base bid and trace the chain
            # up through mergers...
            merged_base_bid = bid
            base_bid = bid
            while merged_base_bid in merged_bins:
                merged_base_bid = merged_bins[merged_base_bid]
            base_bin = self.bins[base_bid]
            
            # get the K closest bins
            neighbor_list = [tdm_index_2_bid[i] for i in search_tree.query(tdm[bid_2_tdm_index[bid]],k=K)[1]]

            if verbose:
                print "++++++++++"
                print bid, neighbor_list

            # test each neighbor in turn
            for i in range(1,K):
                # get the query bid and trace the chain
                # up through mergers...
                query_bid = neighbor_list[i]
                
                #verbose = base_bid == 23 and query_bid == 24 or base_bid == 48 and query_bid == 49 or base_bid == 48 and query_bid == 50 or base_bid == 49 and query_bid == 50 or base_bid == 30 and query_bid == 66
                


                merged_query_bid = query_bid
                while merged_query_bid in merged_bins:
                    merged_query_bid = merged_bins[merged_query_bid]
                
                if verbose:
                    print "++++++++++"
                    print base_bid, query_bid, merged_base_bid, merged_query_bid 
#-----
# TIME WASTERS
                    
                # process each BID pair once only (takes care of self comparisons too!)
                seen_key = self.makeBidKey(base_bid, query_bid)
                if(seen_key in processed_pairs or
                   merged_base_bid == merged_query_bid):
                    if verbose:
                        print "TW"
                    continue
                processed_pairs[seen_key] = True
                
                query_bin = self.bins[query_bid]
#-----
# KMER ELLIPSE OVERLAP
                if (bin_k_ellipse_areas[base_bid] <= bin_k_ellipse_areas[query_bid]):
                    INTT = ET.doesIntersect2D(bin_k_ellipses[query_bid][0],
                                              bin_k_ellipses[query_bid][1],
                                              bin_k_ellipses[base_bid][0],
                                              bin_k_ellipses[base_bid][1])
                else:
                    INTT = ET.doesIntersect2D(bin_k_ellipses[base_bid][0],
                                              bin_k_ellipses[base_bid][1],
                                              bin_k_ellipses[query_bid][0],
                                              bin_k_ellipses[query_bid][1])
                if verbose:
                    
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1)
                    base_bin.plotMersOnAx(ax,
                                          self.PM.kmerVals,
                                          self.PM.kmerVals2,
                                          self.PM.contigColours,
                                          self.PM.contigLengths,
                                          ET=ET)
                    query_bin.plotMersOnAx(ax,
                                           self.PM.kmerVals,
                                           self.PM.kmerVals2,
                                           self.PM.contigColours,
                                           self.PM.contigLengths,
                                           ET=ET)
                    plt.title("MERGE: %d -> %d (%d)" % (base_bid, query_bid, INTT))                
                    plt.show()
                    plt.close(fig)
                    del fig

                if not INTT:
                    if verbose:
                        print "KINTT" 
                    continue
                    
                if False:
                    # check to see if the kmer values overlap sufficiently or not
                    olap_amounts = base_bin.overlappingKVals(self.PM.kmerVals, query_bin)
                    if not (olap_amounts[0] and (olap_amounts[1] > 0.4 or olap_amounts[2] > 0.4)):
                        if verbose:
                            print "OL", olap_amounts 
                            print base_bin.lowestK, base_bin.highestK, query_bin.lowestK, query_bin.highestK
                        continue
#-----
# CONTIG LENGTH SANITY
                # Test the smaller bin against the larger
                if query_bin.binSize < base_bin.binSize:
                    lengths_wrong = GT.isMaxOutlier(np_median(bin_c_lengths[query_bid]),
                                                    bin_c_lengths[base_bid]
                                                    )
                else:
                    lengths_wrong = GT.isMaxOutlier(np_median(bin_c_lengths[base_bid]),
                                                    bin_c_lengths[query_bid]
                                                    )
                if lengths_wrong:
                    if verbose:
                        print "LW"
                    continue
                    
#-----
# MINIMUM BOUNDING COVERAGE ELLIPSOID
                if False:
                    comb_points = np_array([self.PM.transformedCP[i] for i in 
                                            np_concatenate((base_bin.rowIndices,
                                                            query_bin.rowIndices))
                                            ])
                    (ccenter, cradii, crotation) = ET.getMinVolEllipse(comb_points)
                    comb_vol = ET.getEllipsoidVolume(cradii)
                    if np_max([bin_c_ellipsoid_volumes[base_bid], bin_c_ellipsoid_volumes[query_bid]])*4 < comb_vol:
                        if verbose:
                            print "EV", bin_c_ellipsoid_volumes[base_bid], bin_c_ellipsoid_volumes[query_bid], comb_vol
                        continue

                # determine if intersection exists
                if bin_c_ellipsoid_volumes[base_bid] <= bin_c_ellipsoid_volumes[query_bid]:
                    intersects = ET.doesIntersect3D(bin_c_ellipsoids[query_bid][0],
                                                    bin_c_ellipsoids[query_bid][1],
                                                    bin_c_ellipsoids[base_bid][0],
                                                    bin_c_ellipsoids[base_bid][1])
                else:
                    intersects = ET.doesIntersect3D(bin_c_ellipsoids[base_bid][0],
                                                    bin_c_ellipsoids[base_bid][1],
                                                    bin_c_ellipsoids[query_bid][0],
                                                    bin_c_ellipsoids[query_bid][1]) 

                if verbose:
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1, projection='3d')
                    base_bin.plotOnAx(ax, self.PM.transformedCP, self.PM.contigColours, self.PM.contigLengths, ET=ET)                
                    query_bin.plotOnAx(ax, self.PM.transformedCP, self.PM.contigColours, self.PM.contigLengths, ET=ET)
                    plt.title("MERGE: %d -> %d (%d)" % (base_bid, query_bid, INTT))                
                    plt.show()
                    plt.close(fig)
                    del fig

                if not intersects:
                    if verbose:
                        print "CINTT"
                    continue
                
                # We only get here if we're going to merge the bins
                if merged_query_bid < merged_base_bid:
                    merged_bins[merged_base_bid] = merged_query_bid
                    # we just nuked the base bid
                    break
                else:
                    merged_bins[merged_query_bid] = merged_base_bid
#-----
# MERGE
        # now make a bunch of mergers
        mergers = []
        processed_bids = {}     # bid => index in mergers
        for bid in merged_bins:
            
            # trace this guy until the end
            merging_bid = merged_bins[bid]
            while merging_bid in merged_bins:
                merging_bid = merged_bins[merging_bid]
            
            try:
                merge_list_id_1 = processed_bids[bid]
            except KeyError:
                merge_list_id_1 = -1
            try:
                merge_list_id_2 = processed_bids[merging_bid]
            except KeyError:
                merge_list_id_2 = -1
            
            if merge_list_id_1 == -1:
                if merge_list_id_2 == -1:
                    # all new
                    index = len(mergers)
                    processed_bids[merging_bid] = index
                    processed_bids[bid] = index
                    mergers.append([bid, merging_bid])
                else:
                    processed_bids[bid] = merge_list_id_2
                    mergers[merge_list_id_2].append(bid)
            elif merge_list_id_2 == -1:
                processed_bids[merging_bid] = merge_list_id_1
                mergers[merge_list_id_1].append(merging_bid)
            else:
                print "gonna hurt!", bid, merging_bid, merged_bins[bid], mergers[merge_list_id_1], mergers[merge_list_id_2]

        # now merge them
        num_bins_removed = 0
        for merge in mergers:
            self.plotSelectBins(merge, plotMers=True, ET=ET)
            if True:#verbose:
                print merge
            #num_bins_removed += (len(merge) - 1)
            #self.merge(merge, auto=True, newBid=False, saveBins=False, verbose=False, printInstructions=False)

        print "    Removed %d cores leaving %d cores" % (num_bins_removed, len(self.bins))        

    def autoRefineBins(self,
                       iterate=False,
                       mergeObvious=True,
                       removeDuds=True,
                       nukeOutliers=True,
                       verbose=False,
                       plotAfterOB=True):
        """Automagically refine bins"""

#-----
# PRELIM CULLING
        
        # identify and remove outlier bins
        if nukeOutliers:
            self.nukeOutliers()

        if mergeObvious:
            self.mergeObvious()

        if plotAfterOB:
            bids = self.getBids()
            for bid in bids:
                self.bins[bid].makeBinDist(self.PM.transformedCP, 
                                           self.PM.averageCoverages, 
                                           self.PM.kmerVals, 
                                           self.PM.contigLengths)
            self.plotBins(FNPrefix="AFTER_OB", ET=EllipsoidTool())

        return
#-----
# MAKE SEARCH TREES
        # search-based data structures
        neighbour_lists = {} # save looking things up a 1,000,000 times
        graduated_tdms = {}
        graduated_searches = {}
        gt_2_ri_lookup = {}
        ri_2_gt_lookup = {}
        
        # make 10 separate search trees based on kmer sigs
        for i in range(10):
            graduated_tdms[i] = []
            gt_2_ri_lookup[i] = {}
            ri_2_gt_lookup[i] = {}
            neighbour_lists[i] = {}
            
        for row_index in range(np_size(self.PM.indices)):
            if self.PM.binIds[row_index] != 0:
                tdm_index = int(self.PM.kmerVals[row_index]*10)
                lower = np_max([tdm_index-1, 0])
                upper = np_min([tdm_index+2, 10])
                for i in range(lower, upper):
                    gt_index = len(graduated_tdms[i])
                    ri_2_gt_lookup[i][row_index] = gt_index 
                    gt_2_ri_lookup[i][gt_index] = row_index 
                    graduated_tdms[i].append(self.PM.transformedCP[row_index])

        for i in range(10):
            #print i, len(graduated_tdms[i])
            graduated_searches[i] = kdt(graduated_tdms[i])

        # pay attention to contig lengths when including in bins
        # use grubbs test
        GT = GrubbsTester() # we need to perform grubbs test before inclusion
        
        # keep track of what gets merged where
        merged_bins = {}        # oldId => newId
        processed_pairs = {}    # keep track of pairs we've analysed
        bin_c_lengths = {}

        bids = self.getBids()


        stable_bids = {} # once a bin is stable it's stable... almost
        re_unstable = {}

        # we can get stuck in an infinite loop
        # where a contig is passed back and forth between a pair of bins. 
        all_moved_RIs = {}      # when a contig has been moved FROM a bin it can 
                                # not go BACK to THAT bin
        
        bins_changed = {}
        bin_c_lengths = {}       
        start_num_bins = len(bids)
        for bid in bids:
            self.bins[bid].covTolerance = 3
            self.bins[bid].kValTolerance = 3
            bins_changed[bid] = 1       
        
        super_round = 1
        while True:
            sr_contigs_reassigned = 0
            sr_bins_removed = 0
            num_reassigned = -1
            round = 0
            while num_reassigned != 0:
                num_reassigned = 0
                reassignment_map = {}
                moved_RIs = {}
                # make a lookup of each bins contig length distributions
                # this is used in the Grubbs testing
                bids = self.getBids()
                for bid in bids:
                    if bid in bins_changed:
                        bin_c_lengths[bid] = []
                        self.bins[bid].makeBinDist(self.PM.transformedCP, 
                                                   self.PM.averageCoverages, 
                                                   self.PM.kmerVals, 
                                                   self.PM.contigLengths)
                        for row_index in self.bins[bid].rowIndices:
                            bin_c_lengths[bid].append(self.PM.contigLengths[row_index])
                
                bins_changed = {}
                unstables = {}  # bins we destabilize by assigning contigs TO them
                num_done_now = 0
                for bid in bids:
                    bin = self.getBin(bid)
                    num_done_now = 0
                    if bid in stable_bids:
                        # We may have added something to the reassignment
                        # map for this guy already. append to be sure.
                        try:
                            reassignment_map[bid] += list(bin.rowIndices)
                        except KeyError:
                            reassignment_map[bid] = list(bin.rowIndices)
                    else:
                        stable = True   # stable till proven guilty
                        for row_index in bin.rowIndices:
                            tdm_index = np_min([int(self.PM.kmerVals[row_index]*10),9])
                            (assigned_bid,
                             neighbour_lists[tdm_index],
                             scores
                             ) = self.getClosestBID(ri_2_gt_lookup[tdm_index][row_index],
                                                    graduated_searches[tdm_index],
                                                    graduated_tdms[tdm_index],
                                                    gt_2_ri_lookup[tdm_index],
                                                    neighbourList=neighbour_lists[tdm_index],
                                                    verbose=verbose,
                                                    k=2*bin.binSize-1
                                                    )
                            assign_to_new = True    # track this var to see if we will assign or not
                            # make sure we are reassigning
                            if assigned_bid == bid:
                                assign_to_new = False
                            # make sure the contig isn't way too large for the bin
                            elif GT.isMaxOutlier(self.PM.contigLengths[row_index], bin_c_lengths[assigned_bid], verbose=verbose):
                                assign_to_new = False
                            # we need to check that we are not assigning a contig back
                            # to a bin it came from
                            elif row_index in all_moved_RIs:
                                if assigned_bid in all_moved_RIs[row_index]:
                                    assign_to_new = False
                            
                            # we have found the closest based only on coverage 
                            # profile. We need to check and see if the kmersig makes
                            # any sense
                            elif not self.bins[assigned_bid].withinLimits(self.PM.kmerVals,
                                                                          self.PM.averageCoverages,
                                                                          row_index):
                                assign_to_new = False
                                
                            if assign_to_new:
                                # we are assigning it to a new bid
                                stable = False
                                unstables[assigned_bid] = True

                                num_reassigned += 1
                                num_done_now += 1
                                sr_contigs_reassigned += 1
                                # keep this for updating the PMs bin ID hash
                                moved_RIs[row_index] = assigned_bid
                                # make note of where this guy came from
                                if row_index in all_moved_RIs:
                                    try:
                                        all_moved_RIs[row_index][bid]+=1
                                    except KeyError:
                                        all_moved_RIs[row_index][bid] = 1
                                else:
                                    all_moved_RIs[row_index] = {bid: 1}
                                
                                bins_changed[bid] = 1
                                bins_changed[assigned_bid] = 1
                                
                            else:
                                assigned_bid = bid
                            
                            # now go through and put the index where it belongs
                            try:
                                reassignment_map[assigned_bid].append(row_index)
                            except KeyError:
                                reassignment_map[assigned_bid] = [row_index]
                                
                        if stable: 
                            # no changes this round, mark bin as stable
                            stable_bids[bid] = True

                # we need to destabilise any bin who accepted contigs
                for u_bid in unstables:
                    try:
                        del stable_bids[u_bid]
                    except KeyError:
                        pass
                
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
                        #print "After AR delete:", bid, len(self.bins)
                sr_bins_removed += bins_removed
                round += 1
                if verbose:
                    print "    Refine sub-round %d: reassigned %d contigs, removed %d cores" % (round, num_reassigned, bins_removed)
                    
            print "    Refine round %d. (%d iterations, %d contigs reassigned, %d cores removed)" % (super_round, round, sr_contigs_reassigned, sr_bins_removed)
            if sr_contigs_reassigned == 0:
                break
            super_round += 1

        print "    Removed %d cores leaving %d cores" % (start_num_bins-len(self.bins), len(self.bins))        
        
        if removeDuds:    
            self.removeDuds()

    def removeDuds(self, ms=20, mv=1000000, verbose=False):
        """Run this after refining to remove scrappy leftovers"""
        print "    Removing dud cores (min %d contigs or %d bp)" % (ms, mv)
        deleters = []
        for bid in self.getBids():
            bin = self.bins[bid]
            if not self.isGoodBin(bin.totalBP, bin.binSize, ms=ms, mv=mv):
                # delete this chap!
                deleters.append(bid)
        if verbose:
            print "duds", deleters
        if len(deleters) > 0:
            self.deleteBins(deleters,
                            force=True,
                            freeBinnedRowIndicies=True,
                            saveBins=False)
        print "    Removed %d cores leaving %d cores" % (len(deleters), len(self.bins))
            
#------------------------------------------------------------------------------
# BIN UTILITIES 

    def getBids(self):
        """Return a sorted list of bin ids"""
        return sorted(self.bins.keys())

    def isGoodBin(self, totalBP, binSize, ms=0, mv=0):
        """Does this bin meet my exacting requirements?"""
        if(ms == 0):
            ms = self.minSize               # let the user choose
        if(mv == 0):
            mv = self.minVol                # let the user choose
        
        if(totalBP < mv):                   # less than the good volume
            if(binSize > ms):               # but has enough contigs
                return True
        else:                               # contains enough bp to pass regardless of number of contigs
            return True        
        return False

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
            for row_index in dead_bin.rowIndices:
                self.PM.binIds[row_index] = parent_bin.id
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
                
                for row_index in dead_bin.rowIndices:
                    self.PM.binIds[row_index] = parent_bin.id
                
                parent_bin.consume(self.PM.transformedCP,
                                   self.PM.averageCoverages,
                                   self.PM.kmerVals,
                                   self.PM.contigLengths,
                                   dead_bin,
                                   verbose=verbose)
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
                        try:
                            del self.PM.binnedRowIndicies[row_index]
                        except KeyError:
                            print bid, row_index, "FUNG"
                        self.PM.binIds[row_index] = 0
                            
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
        valid_responses = ['R','P','B','V','M','S','K','Q']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" What next? ("+vrs+") : ")
            else:
                option = raw_input(" How do you want to continue?\n" \
                                   " r = replot ids, p = replot points, b = plot one or more bins,\n" \
                                   " v = plot in vincinity of bin, m = merge, s = split,\n" \
                                   " k = set kmer range, q = quit\n" \
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

    def calculateAngles(self, rowSet1, rowSet2):
        """work out the angles between a set of contigs"""
        angles = []
        for i in rowSet1:
            for j in rowSet2:
                if i != j:
                    angles.append(self.getAngleBetween(i, j))
        return angles

    def getAngleBetween(self, rowIndex1, rowIndex2):
        """Find the angle between two contig's coverage vectors"""
        u = self.PM.covProfiles[rowIndex1]
        v = self.PM.covProfiles[rowIndex2]
        
        try:
            ac = np_arccos(np_dot(u,v)/np_norm(u)/np_norm(v))
        except FloatingPointError:
            return 0
        return ac

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

    def plotSelectBins(self, bids, plotMers=False, fileName="", plotEllipsoid=False, ET=None):
        """Plot a selection of bids in a window"""
        if plotEllipsoid and ET == None:
            ET = EllipsoidTool()
            
        fig = plt.figure()
        
        # we need to do some fancy-schmancy stuff at times!
        if plotMers:
            num_cols = 2
        else:
            num_cols = 1
            
        ax = fig.add_subplot(1,num_cols,1, projection='3d')
        for bid in bids:
            self.bins[bid].plotOnAx(ax,
                           self.PM.transformedCP,
                           self.PM.contigColours,
                           self.PM.contigLengths,
                           ET=ET,
                           printID=True
                           )
        if plotMers:
            ax.set_title('Coverage')
            ax = fig.add_subplot(1, 2, 2)
            for bid in bids:
                self.bins[bid].plotMersOnAx(ax,
                                            self.PM.kmerVals,
                                            self.PM.kmerVals2,
                                            self.PM.contigColours,
                                            self.PM.contigLengths,
                                            ET=ET,
                                            printID=True
                                            )
            ax.set_title('Kmer sig PCA')
                
        if(fileName != ""):
            try:
                fig.set_size_inches(6,6)
                plt.savefig(fileName,dpi=300)
            except:
                print "Error saving image:", fileName, exc_info()[0]
                raise
        else:
            try:
                plt.show()
            except:
                print "Error showing image:", exc_info()[0]
                raise
            
        plt.close(fig)
        del fig
        

    def plotBins(self, FNPrefix="BIN", sideBySide=False, folder='', plotEllipsoid=False, ET=None):
        """Make plots of all the bins"""
        if plotEllipsoid and ET == None:
            ET = EllipsoidTool()
        
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
                    self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, self.PM.contigLengths, fileName=osp_join(folder, FNPrefix+"_"+str(bid)), ET=ET)
                else:
                    self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, self.PM.contigLengths, FNPrefix+"_"+str(bid), ET=ET)

    def plotSideBySide(self, bids, fileName="", tag=""):
        """Plot two bins side by side in 3d"""
        ET = EllipsoidTool()
        fig = plt.figure()
        # we need to work out how to shape the plots
        num_plots = len(bids)
        plot_rows = float(int(np_sqrt(num_plots)))
        plot_cols = np_ceil(float(num_plots)/plot_rows)
        plot_num = 1
        for bid in bids:
            title = self.bins[bid].plotOnFig(fig, plot_rows, plot_cols, plot_num, self.PM.transformedCP, self.PM.contigColours, self.PM.contigLengths, ET=ET, fileName=fileName)
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

    def plotBinIds(self, krange=None, ignoreRanges=False):
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
        
        if not ignoreRanges:
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

    def plotBinPoints(self, ignoreRanges=False):
        """Render the image for validating cores"""
        (bin_centroid_points, bin_centroid_colours, bids) = self.findCoreCentres()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(bin_centroid_points[:,0], bin_centroid_points[:,1], bin_centroid_points[:,2], edgecolors=bin_centroid_colours, c=bin_centroid_colours)

        if not ignoreRanges:
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class GrubbsTester:
    """Data and methods for performing Grubbs test
    
    cutoff values taken from qgrubs from R package outliers
    using command: qgrubbs(0.99, c(3:1002), 10)
    """
    def __init__(self):
        # cutoff values for n degress of freedom 
        # If you have 8 sample points then self.cutoffs[6]
        # is what you want!
        self.critVs = np_array([1.154637,1.492500,1.748857,1.944245,2.097304,2.220833,2.323148,2.409725,
                                2.484279,2.549417,2.607020,2.658480,2.704855,2.746963,2.785445,2.820817,
                                2.853495,2.883821,2.912078,2.938503,2.963296,2.986628,3.008645,3.029473,
                                3.049223,3.067989,3.085855,3.102897,3.119180,3.134761,3.149694,3.164026,
                                3.177798,3.191049,3.203813,3.216121,3.228002,3.239482,3.250585,3.261332,
                                3.271744,3.281839,3.291634,3.301145,3.310386,3.319372,3.328114,3.336624,
                                3.344914,3.352993,3.360872,3.368558,3.376061,3.383388,3.390546,3.397543,
                                3.404385,3.411078,3.417628,3.424041,3.430321,3.436474,3.442505,3.448417,
                                3.454215,3.459902,3.465484,3.470963,3.476342,3.481626,3.486816,3.491917,
                                3.496930,3.501859,3.506706,3.511474,3.516164,3.520780,3.525324,3.529797,
                                3.534201,3.538539,3.542812,3.547022,3.551171,3.555260,3.559292,3.563266,
                                3.567186,3.571051,3.574865,3.578627,3.582340,3.586004,3.589620,3.593190,
                                3.596715,3.600196,3.603634,3.607030,3.610384,3.613698,3.616973,3.620209,
                                3.623407,3.626569,3.629695,3.632785,3.635840,3.638862,3.641851,3.644807,
                                3.647731,3.650624,3.653486,3.656319,3.659122,3.661896,3.664642,3.667360,
                                3.670050,3.672714,3.675352,3.677964,3.680551,3.683113,3.685650,3.688164,
                                3.690654,3.693121,3.695565,3.697986,3.700386,3.702764,3.705121,3.707457,
                                3.709773,3.712068,3.714344,3.716600,3.718836,3.721054,3.723253,3.725434,
                                3.727597,3.729742,3.731869,3.733979,3.736072,3.738149,3.740209,3.742253,
                                3.744281,3.746293,3.748289,3.750270,3.752237,3.754188,3.756125,3.758047,
                                3.759955,3.761849,3.763729,3.765595,3.767448,3.769287,3.771114,3.772928,
                                3.774728,3.776516,3.778292,3.780056,3.781807,3.783546,3.785274,3.786990,
                                3.788694,3.790387,3.792069,3.793740,3.795400,3.797049,3.798687,3.800315,
                                3.801932,3.803540,3.805137,3.806723,3.808300,3.809868,3.811425,3.812973,
                                3.814511,3.816040,3.817560,3.819071,3.820572,3.822065,3.823549,3.825024,
                                3.826490,3.827948,3.829397,3.830838,3.832271,3.833696,3.835112,3.836521,
                                3.837921,3.839314,3.840699,3.842076,3.843446,3.844808,3.846163,3.847510,
                                3.848850,3.850183,3.851509,3.852827,3.854139,3.855444,3.856742,3.858033,
                                3.859317,3.860595,3.861866,3.863131,3.864389,3.865640,3.866886,3.868125,
                                3.869358,3.870585,3.871805,3.873020,3.874229,3.875431,3.876628,3.877819,
                                3.879005,3.880184,3.881358,3.882526,3.883689,3.884846,3.885998,3.887144,
                                3.888285,3.889421,3.890551,3.891677,3.892797,3.893911,3.895021,3.896126,
                                3.897226,3.898321,3.899411,3.900496,3.901576,3.902651,3.903722,3.904788,
                                3.905849,3.906906,3.907958,3.909005,3.910048,3.911087,3.912121,3.913150,
                                3.914176,3.915197,3.916213,3.917226,3.918234,3.919238,3.920238,3.921233,
                                3.922225,3.923212,3.924196,3.925175,3.926151,3.927122,3.928090,3.929054,
                                3.930014,3.930970,3.931922,3.932871,3.933815,3.934756,3.935694,3.936627,
                                3.937558,3.938484,3.939407,3.940326,3.941242,3.942155,3.943064,3.943969,
                                3.944871,3.945770,3.946665,3.947557,3.948445,3.949331,3.950213,3.951091,
                                3.951967,3.952839,3.953708,3.954574,3.955437,3.956297,3.957154,3.958007,
                                3.958858,3.959705,3.960550,3.961391,3.962229,3.963065,3.963898,3.964727,
                                3.965554,3.966378,3.967199,3.968017,3.968833,3.969645,3.970455,3.971262,
                                3.972066,3.972868,3.973667,3.974463,3.975256,3.976047,3.976836,3.977621,
                                3.978404,3.979184,3.979962,3.980738,3.981510,3.982280,3.983048,3.983813,
                                3.984576,3.985336,3.986094,3.986849,3.987602,3.988353,3.989101,3.989847,
                                3.990590,3.991331,3.992070,3.992806,3.993541,3.994272,3.995002,3.995729,
                                3.996454,3.997177,3.997898,3.998616,3.999332,4.000046,4.000758,4.001468,
                                4.002175,4.002880,4.003584,4.004285,4.004984,4.005681,4.006376,4.007069,
                                4.007759,4.008448,4.009135,4.009820,4.010502,4.011183,4.011862,4.012539,
                                4.013213,4.013886,4.014557,4.015226,4.015893,4.016558,4.017222,4.017883,
                                4.018543,4.019200,4.019856,4.020510,4.021162,4.021812,4.022461,4.023108,
                                4.023752,4.024395,4.025037,4.025676,4.026314,4.026950,4.027585,4.028217,
                                4.028848,4.029477,4.030105,4.030730,4.031354,4.031977,4.032597,4.033217,
                                4.033834,4.034450,4.035064,4.035676,4.036287,4.036896,4.037504,4.038110,
                                4.038715,4.039318,4.039919,4.040519,4.041117,4.041714,4.042309,4.042903,
                                4.043495,4.044085,4.044674,4.045262,4.045848,4.046433,4.047016,4.047597,
                                4.048178,4.048756,4.049334,4.049909,4.050484,4.051057,4.051628,4.052198,
                                4.052767,4.053334,4.053900,4.054465,4.055028,4.055590,4.056150,4.056709,
                                4.057267,4.057823,4.058378,4.058932,4.059484,4.060035,4.060585,4.061133,
                                4.061680,4.062226,4.062771,4.063314,4.063856,4.064396,4.064935,4.065474,
                                4.066010,4.066546,4.067080,4.067613,4.068145,4.068676,4.069205,4.069733,
                                4.070260,4.070786,4.071310,4.071834,4.072356,4.072877,4.073396,4.073915,
                                4.074432,4.074949,4.075464,4.075977,4.076490,4.077002,4.077512,4.078022,
                                4.078530,4.079037,4.079543,4.080047,4.080551,4.081054,4.081555,4.082056,
                                4.082555,4.083053,4.083550,4.084046,4.084541,4.085035,4.085528,4.086019,
                                4.086510,4.087000,4.087488,4.087976,4.088462,4.088948,4.089432,4.089915,
                                4.090398,4.090879,4.091359,4.091839,4.092317,4.092794,4.093271,4.093746,
                                4.094220,4.094693,4.095166,4.095637,4.096107,4.096577,4.097045,4.097513,
                                4.097979,4.098445,4.098909,4.099373,4.099836,4.100297,4.100758,4.101218,
                                4.101677,4.102135,4.102592,4.103048,4.103503,4.103958,4.104411,4.104864,
                                4.105315,4.105766,4.106216,4.106665,4.107113,4.107560,4.108006,4.108452,
                                4.108896,4.109340,4.109783,4.110225,4.110666,4.111106,4.111545,4.111984,
                                4.112421,4.112858,4.113294,4.113729,4.114163,4.114597,4.115029,4.115461,
                                4.115892,4.116322,4.116751,4.117180,4.117608,4.118034,4.118460,4.118886,
                                4.119310,4.119734,4.120157,4.120579,4.121000,4.121421,4.121840,4.122259,
                                4.122677,4.123095,4.123511,4.123927,4.124342,4.124757,4.125170,4.125583,
                                4.125995,4.126406,4.126817,4.127226,4.127635,4.128044,4.128451,4.128858,
                                4.129264,4.129669,4.130074,4.130478,4.130881,4.131284,4.131685,4.132086,
                                4.132487,4.132886,4.133285,4.133683,4.134081,4.134477,4.134873,4.135269,
                                4.135663,4.136057,4.136451,4.136843,4.137235,4.137626,4.138017,4.138407,
                                4.138796,4.139184,4.139572,4.139959,4.140346,4.140732,4.141117,4.141501,
                                4.141885,4.142268,4.142651,4.143033,4.143414,4.143794,4.144174,4.144554,
                                4.144932,4.145310,4.145688,4.146064,4.146440,4.146816,4.147191,4.147565,
                                4.147938,4.148311,4.148684,4.149055,4.149427,4.149797,4.150167,4.150536,
                                4.150905,4.151273,4.151640,4.152007,4.152373,4.152739,4.153104,4.153468,
                                4.153832,4.154195,4.154558,4.154920,4.155282,4.155643,4.156003,4.156363,
                                4.156722,4.157080,4.157438,4.157796,4.158153,4.158509,4.158865,4.159220,
                                4.159574,4.159928,4.160282,4.160635,4.160987,4.161339,4.161690,4.162041,
                                4.162391,4.162740,4.163089,4.163438,4.163786,4.164133,4.164480,4.164826,
                                4.165172,4.165517,4.165862,4.166206,4.166550,4.166893,4.167235,4.167577,
                                4.167919,4.168260,4.168600,4.168940,4.169280,4.169619,4.169957,4.170295,
                                4.170632,4.170969,4.171306,4.171641,4.171977,4.172311,4.172646,4.172980,
                                4.173313,4.173646,4.173978,4.174310,4.174641,4.174972,4.175302,4.175632,
                                4.175962,4.176290,4.176619,4.176947,4.177274,4.177601,4.177927,4.178253,
                                4.178579,4.178904,4.179228,4.179552,4.179876,4.180199,4.180522,4.180844,
                                4.181166,4.181487,4.181808,4.182128,4.182448,4.182767,4.183086,4.183405,
                                4.183723,4.184040,4.184357,4.184674,4.184990,4.185306,4.185621,4.185936,
                                4.186250,4.186564,4.186878,4.187191,4.187503,4.187815,4.188127,4.188438,
                                4.188749,4.189060,4.189370,4.189679,4.189988,4.190297,4.190605,4.190913,
                                4.191220,4.191527,4.191834,4.192140,4.192446,4.192751,4.193056,4.193360,
                                4.193664,4.193968,4.194271,4.194574,4.194876,4.195178,4.195479,4.195781,
                                4.196081,4.196381,4.196681,4.196981,4.197280,4.197578,4.197877,4.198175,
                                4.198472,4.198769,4.199066,4.199362,4.199658,4.199953,4.200248,4.200543,
                                4.200837,4.201131,4.201424,4.201718,4.202010,4.202303,4.202594,4.202886,
                                4.203177,4.203468,4.203758,4.204048,4.204338,4.204627,4.204916,4.205204,
                                4.205493,4.205780,4.206068,4.206355,4.206641,4.206927,4.207213,4.207499,
                                4.207784,4.208068,4.208353,4.208637,4.208920,4.209204,4.209487,4.209769,
                                4.210051,4.210333,4.210615,4.210896,4.211176,4.211457,4.211737,4.212016,
                                4.212296,4.212575,4.212853,4.213132,4.213409,4.213687,4.213964,4.214241,
                                4.214517,4.214794,4.215069,4.215345,4.215620,4.215895,4.216169,4.216443,
                                4.216717,4.216990,4.217263,4.217536,4.217808,4.218080,4.218352,4.218623,
                                4.218894,4.219165,4.219436,4.219706,4.219975,4.220245,4.220514,4.220782,
                                4.221051,4.221319,4.221586,4.221854,4.222121,4.222388,4.222654,4.222920,
                                4.223186,4.223451,4.223716,4.223981,4.224246,4.224510,4.224774,4.225037,
                                4.225300,4.225563,4.225826,4.226088,4.226350,4.226611,4.226873,4.227134,
                                4.227394,4.227655,4.227915,4.228175,4.228434,4.228693,4.228952,4.229211,
                                4.229469,4.229727,4.229984,4.230242,4.230499,4.230755,4.231012,4.231268,
                                4.231524,4.231779,4.232034,4.232289,4.232544,4.232798,4.233052,4.233306,
                                4.233559,4.233813,4.234065,4.234318,4.234570,4.234822,4.235074,4.235325,
                                4.235576,4.235827,4.236078,4.236328,4.236578,4.236827,4.237077,4.237326,
                                4.237575,4.237823,4.238071,4.238319,4.238567,4.238815,4.239062,4.239308,
                                4.239555,4.239801,4.240047,4.240293,4.240538,4.240784,4.241028,4.241273,
                                4.241517,4.241761,4.242005,4.242249,4.242492,4.242735,4.242978,4.243220,
                                4.243462,4.243704,4.243946,4.244187,4.244428,4.244669,4.244910,4.245150,
                                4.245390,4.245630,4.245869,4.246108,4.246347,4.246586,4.246825,4.247063]
                               )

    def isMaxOutlier(self, maxVal, compVals, verbose=False):
        """Test if the maxVal is an outlier 
        
        maxVal should NOT be included in compVals
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
        
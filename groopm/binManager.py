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
__copyright__ = "Copyright 2012/2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.2.7"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################
from os.path import join as osp_join
from sys import exc_info, exit, stdout as sys_stdout
from operator import itemgetter


import matplotlib.pyplot as plt
from pylab import show
from numpy import (abs as np_abs,
                   append as np_append,
                   arange as np_arange,
                   arccos as np_arccos,
                   argmax as np_argmax,
                   argsort as np_argsort,
                   around as np_around,
                   array as np_array,
                   bincount as np_bincount,
                   ceil as np_ceil,
                   concatenate as np_concatenate,
                   dot as np_dot,
                   log10 as np_log10,
                   max as np_max,
                   mean as np_mean,
                   median as np_median,
                   min as np_min,
                   ones as np_ones,
                   reshape as np_reshape,
                   seterr as np_seterr,
                   size as np_size,
                   sort as np_sort,
                   sqrt as np_sqrt,
                   std as np_std,
                   sum as np_sum,
                   where as np_where,
                   zeros as np_zeros)

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
    def __init__(self,
                 dbFileName="",
                 pm=None,
                 minSize=10,
                 minVol=1000000):
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

    def setColorMap(self, colorMapStr):
        self.PM.setColorMap(colorMapStr)

#------------------------------------------------------------------------------
# LOADING / SAVING

    def loadBins(self,
                 timer,
                 getUnbinned=False,
                 bids=[],
                 makeBins=False,
                 silent=True,
                 loadKmerPCs=False,
                 loadRawKmers=False,
                 loadCovProfiles=True,
                 loadContigLengths=True,
                 loadLinks=False,
                 loadContigNames=True,
                 cutOff=0,
                 transform=True):
        """Load data and make bin objects"""
        # build the condition
        
        if getUnbinned:
            # get everything
            condition = "(length >= %d) " % cutOff
        else:
            # make sense of bin information
            if bids == []:
                condition = "((length >= %d) & (bid != 0))" % cutOff
            else:
                condition = "((length >= %d) & " % cutOff + " | ".join(["(bid == %d)"%bid for bid in bids])+")"

        # if we're going to make bins then we'll need kmer sigs
        if(makeBins):
            loadKmerPCs=True
            loadCovProfiles=True

        self.PM.loadData(timer,
                         condition,
                         bids=bids,
                         silent=silent,
                         loadCovProfiles=loadCovProfiles,
                         loadKmerPCs=loadKmerPCs,
                         loadRawKmers=loadRawKmers,
                         makeColors=True,
                         loadContigNames=loadContigNames,
                         loadContigLengths=loadContigLengths,
                         loadBins=True,
                         loadLinks=loadLinks
                        )

        # exit if no bins loaded
        if self.PM.numContigs == 0:
            return

        if(makeBins):
            if transform:
                self.PM.transformCP(timer, silent=silent)
            else:
                if self.PM.numStoits == 3:
                    self.PM.transformedCP = self.PM.covProfiles
                else:
                    print "Number of stoits != 3. You need to transform"
                    self.PM.transformCP(timer, silent=silent)
            if not silent:
                print "    Making bin objects"
            self.makeBins(self.getBinMembers())
            if not silent:
                print "    Loaded %d bins from database" % len(self.bins)
        if not silent:
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()

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
                    self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerNormPC1, self.PM.kmerPCs, self.PM.contigGCs, self.PM.contigLengths)
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
        # save the bin assignments
        self.PM.setBinAssignments(
                                  self.getGlobalBinAssignments(binAssignments), # convert to global indices
                                  nuke=nuke
                                  )
        # overwrite the bins table
        self.setBinStats()

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

        # create and array of tuples:
        # [(bid, size, likelyChimeric)]
        bin_stats = []
        for bid in self.getBids():
            # no point in saving empty bins
            if np_size(self.bins[bid].rowIndices) > 0:
                bin_stats.append((bid, np_size(self.bins[bid].rowIndices), self.PM.isLikelyChimeric[bid]))
        self.PM.setBinStats(bin_stats)


#------------------------------------------------------------------------------
# REMOVE ALREADY LOADED BINS

    def removeBinAndIndices(self, bid):
        """Remove indices from the PM based on bin identity

        "unload" some data
        """
        # get some info
        rem_bin = self.getBin(bid)
        original_length = len(self.PM.indices)
        rem_list = np_sort(rem_bin.rowIndices)

        # affect the raw data in the PM
        self.PM.reduceIndices(rem_list)
        del self.PM.validBinIds[bid]

        # remove the bin here
        del self.bins[bid]

        # now fix all the rowIndices in all the other bins
        for bid in self.getBids():
            self.bins[bid].rowIndices = self.fixRowIndexLists(original_length, np_sort(self.bins[bid].rowIndices), rem_list)

    def fixRowIndexLists(self, originalLength, oldList, remList):
        """Fix up row index lists which reference into the
        data structure after a call to reduceIndices

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
# LINKS

    def getLinkingContigs(self, bid):
        """Get all contigs and their bin IDs which link to contigs in this bin"""
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
# BIN UTILITIES

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

    def getBids(self):
        """Return a sorted list of bin ids"""
        return sorted(self.bins.keys())

    def getCentroidProfiles(self, mode="mer"):
        """Return an array containing the centroid stats for each bin"""
        if(mode == "mer"):
            ret_vecs = np_zeros((len(self.bins)))
            outer_index = 0
            for bid in self.getBids():
                ret_vecs[outer_index] = self.bins[bid].kValMeanNormPC1
                outer_index += 1
            return ret_vecs
        elif(mode == "cov"):
            ret_vecs = np_zeros((len(self.bins), len(self.PM.transformedCP[0])))
            outer_index = 0
            for bid in self.getBids():
                ret_vecs[outer_index] = self.bins[bid].covMedians
                outer_index += 1
            return ret_vecs
        else:
            raise ge.ModeNotAppropriateException("Mode",mode,"unknown")

    def split(self, bid, n, mode='kmer', auto=False, saveBins=False, verbose=False, printInstructions=True, use_elipses=True):
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
            print 'here!!!!'
            # charge on through
            self.deleteBins([bids[0]], force=True)  # delete the combined bin
            # save new bins
            self.saveBins(binAssignments=bin_assignment_update)
            return

        # we will need to confer with the user
        # plot some stuff
        # sort the bins by kmer val
        bid_tuples = [(tbid, self.bins[tbid].kValMeanNormPC1) for tbid in bids[1:]]
        bid_tuples.sort(key=itemgetter(1))
        index = 1
        for pair in bid_tuples:
            bids[index] = pair[0]
            index += 1

        self.plotSideBySide(bids, use_elipses=use_elipses)

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
            self.split(bid,
                       n,
                       mode='cov',
                       auto=auto,
                       saveBins=saveBins,
                       verbose=verbose,
                       printInstructions=False,
                       use_elipses=use_elipses)
        elif(user_option == 'K'):
            self.split(bid,
                       n,
                       mode='kmer',
                       auto=auto,
                       saveBins=saveBins,
                       verbose=verbose,
                       printInstructions=False,
                       use_elipses=use_elipses)
        elif(user_option == 'L'):
            self.split(bid,
                       n,
                       mode='len',
                       auto=auto,
                       saveBins=saveBins,
                       verbose=verbose,
                       printInstructions=False,
                       use_elipses=use_elipses)
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
                    self.split(bid,
                               parts,
                               mode=mode,
                               auto=auto,
                               saveBins=saveBins,
                               verbose=verbose,
                               printInstructions=False,
                               use_elipses=use_elipses)


    def getSplitties(self, bid, n, mode):
        """Return a set of split bins"""
        obs = np_array([])
        if(mode=='kmer'):
            obs = np_array([self.PM.kmerNormPC1[i] for i in self.getBin(bid).rowIndices])
        elif(mode=='cov'):
            obs = np_array([self.PM.covProfiles[i] for i in self.getBin(bid).rowIndices])
        elif(mode=='len'):
            obs = np_array([self.PM.contigLengths[i] for i in self.getBin(bid).rowIndices])

        # do the clustering
        try:
            centroids,_ = kmeans(obs,n)
        except ValueError:
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
                holding_array = holding_array.astype(int)
                split_bin = self.makeNewBin(holding_array)

                for row_index in holding_array:
                    bin_assignment_update[row_index] = split_bin.id
                split_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerNormPC1, self.PM.kmerPCs, self.PM.contigGCs, self.PM.contigLengths)
                bids.append(split_bin.id)
                holding_array = np_array([])
                current_group = idx[i]
            holding_array = np_append(holding_array, int(self.getBin(bid).rowIndices[i]))

        # do the last one
        if(np_size(holding_array) != 0):
            holding_array = holding_array.astype(int)
            split_bin = self.makeNewBin(holding_array)
            for row_index in holding_array:
                bin_assignment_update[int(row_index)] = split_bin.id
            split_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerNormPC1, self.PM.kmerPCs, self.PM.contigGCs, self.PM.contigLengths)
            bids.append(split_bin.id)

        return (bin_assignment_update, bids)

    def shouldMerge(self, bin1, bin2, ignoreCov=False, ignoreMer=False, merTol=0, confidence=0.95, verbose=False):
        """Determine whether its wise to merge two bins

        Perfoms a one-way anova to determine if the larger bin would be
        significantly changed if it consumed the smaller

        OR does a tolerance test on kmerNormPC1. We assume that bin1 is larger than bin2
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
                    upper_k_val_cut = bin1.kValMeanNormPC1 + merTol * bin1.kValStdevNormPC1
                    lower_k_val_cut = bin1.kValMeanNormPC1 - merTol * bin1.kValStdevNormPC1

                    if bin2.kValMeanNormPC1 >= lower_k_val_cut and bin2.kValMeanNormPC1 <= upper_k_val_cut:
                        mer_match = True
                    else:
                        mer_match = False
                else:
                    b1_k_dist = bin1.getkmerValDist(self.PM.kmerNormPC1)
                    b2_k_dist = bin2.getkmerValDist(self.PM.kmerNormPC1)
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

    def merge(self, bids, auto=False, manual=False, newBid=False, saveBins=False, verbose=False, printInstructions=True, use_elipses=True):
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
            parent_bin = self.makeNewBin()
            
            # now merge it with the first in the new list
            dead_bin = self.getBin(bids[0])
            for row_index in dead_bin.rowIndices:
                self.PM.binIds[row_index] = parent_bin.id
            parent_bin.consume(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerNormPC1, self.PM.kmerPCs, self.PM.contigGCs, self.PM.contigLengths, dead_bin, verbose=verbose)
            self.deleteBins([bids[0]], force=True)
        else:
            # just use the first given as the parent
            parent_bin = self.getBin(bids[0])
            
        # a merged bin specified by the user should not be considered chimeric since 
        # they are indicating both original bins were reasonable
        if manual:
            self.PM.isLikelyChimeric[parent_bin.id] = False

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
                tmp_bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerNormPC1, self.PM.kmerPCs, self.PM.contigGCs, self.PM.contigLengths)
                self.plotSideBySide([parent_bin.id,dead_bin.id,tmp_bin.id], use_elipses=use_elipses)
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
                                   self.PM.kmerNormPC1,
                                   self.PM.kmerPCs,
                                   self.PM.contigGCs,
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

    def getChimericBinIds(self):
        bids = []
        for bid in self.bins:
            if self.PM.isLikelyChimeric[bid]:
                bids.append(bid)

        return bids

    def getNonChimericBinIds(self):
        bids = []
        for bid in self.bins:
            if not self.PM.isLikelyChimeric[bid]:
                bids.append(bid)
        
        return bids

    def deleteBins(self, bids, force=False, freeBinnedRowIndices=False, saveBins=False):
        """Purge a bin from our lists"""
        if(not force):
            user_option = self.promptOnDelete(bids)
            if(user_option != 'Y'):
                return False
        bin_assignment_update = {}
        for bid in bids:
            if bid in self.bins:
                if(freeBinnedRowIndices):
                    for row_index in self.bins[bid].rowIndices:
                        try:
                            del self.PM.binnedRowIndices[row_index]
                        except KeyError:
                            print bid, row_index, "FUNG"
                        self.PM.binIds[row_index] = 0

                        bin_assignment_update[row_index] = 0
                del self.bins[bid]
                del self.PM.isLikelyChimeric[bid]
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

        self.PM.isLikelyChimeric[bid] = False
        self.bins[bid] = Bin(rowIndices, bid, self.PM.scaleFactor-1)
        return self.bins[bid]

#------------------------------------------------------------------------------
# UI

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
        valid_responses = ['Y','N','C','K','L','P']
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
                                   " k = redo but use kmer profile, l = redo but use length profile,\n" \
                                   " p = choose new number of parts\n" \
                                   " Split? ("+vrs+") : ")
            if(option.upper() in valid_responses):
                if(option.upper() == 'K' and mode.upper() == 'KMER' or option.upper() == 'C' and mode.upper() == 'COV' or option.upper() == 'L' and mode.upper() == 'LEN'):
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
    def findCoreCentres(self, gc_range=None, getKVals=False, processChimeric=False):
        """Find the point representing the centre of each core"""
        bin_centroid_points = np_array([])
        bin_centroid_colors = np_array([])
        bin_centroid_kvals = np_array([])
        bin_centroid_gc = np_array([])
        bids = np_array([])

        if gc_range is not None:
            # we only want to plot a subset of these guys
            gc_low = gc_range[0]
            gc_high = gc_range[1]
        num_added = 0
        for bid in self.getBids():
            if not processChimeric and self.PM.isLikelyChimeric[bid]:
                continue

            add_bin = True
            if gc_range is not None:
                avg_gc = np_median([self.PM.contigGCs[row_index] for row_index in self.bins[bid].rowIndices])
                if avg_gc < gc_low or avg_gc > gc_high:
                    add_bin = False
            if add_bin:
                bin_centroid_points = np_append(bin_centroid_points,
                                                self.bins[bid].covMedians)

                bin_centroid_colors = np_append(bin_centroid_colors, self.PM.colorMapGC(np_median(self.PM.contigGCs[self.bins[bid].rowIndices])))

                bin_centroid_gc = np_append(bin_centroid_gc, np_median(self.PM.contigGCs[self.bins[bid].rowIndices]))

                if getKVals:
                    bin_centroid_kvals = np_append(bin_centroid_kvals,
                                                   np_median([
                                                            self.PM.kmerNormPC1[row_index] for row_index in
                                                            self.bins[bid].rowIndices
                                                            ],
                                                           axis=0)
                                                   )

                bids = np_append(bids, bid)
                num_added += 1

        if num_added != 0:
            bin_centroid_points = np_reshape(bin_centroid_points, (num_added, 3))
            bin_centroid_colors = np_reshape(bin_centroid_colors, (num_added, 4))

        if getKVals:
            return (bin_centroid_points, bin_centroid_colors, bin_centroid_gc, bin_centroid_kvals, bids)

        return (bin_centroid_points, bin_centroid_colors, bin_centroid_gc, bids)

    def getAngleBetween(self, rowIndex1, rowIndex2, ):
        """Find the angle between two contig's coverage vectors"""
        u = self.PM.covProfiles[rowIndex1]
        v = self.PM.covProfiles[rowIndex2]
        try:
            ac = np_arccos(np_dot(u,v)/self.PM.normCoverages[rowIndex1]/self.PM.normCoverages[rowIndex1])
        except FloatingPointError:
            return 0
        return ac

    def scoreContig(self, rowIndex, bid):
        """Determine how well a particular contig fits with a bin"""
        return self.getBin(bid).scoreProfile(self.PM.kmerNormPC1[rowIndex], self.PM.transformedCP[rowIndex])

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
                (Ms[bid], Ss[bid], Rs[bid]) = self.bins[bid].getInnerVariance(self.PM.kmerNormPC1)
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
        if(outFormat == 'contigs'):
            stream.write(separator.join(["#\"bid\"","\"cid\"","\"length\"","\"GC\""])+"\n")
        elif(outFormat == 'bins'):
            header = ["\"bin id\"","\"Likely chimeric\"","\"length (bp)\"","\"# seqs\"","\"GC mean\"","\"GC std\""]
            for i in xrange(0, len(self.PM.covProfiles[0])):
                header.append("\"Coverage " + str(i+1) + " mean\"")
                header.append("\"Coverage " + str(i+1) + " std\"")
            stream.write(separator.join(header) + "\n")
        elif(outFormat == 'full'):
            pass
        else:
            print "Error: Unrecognised format:", outFormat
            return

        for bid in self.getBids():
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerNormPC1,
                                        self.PM.kmerPCs, self.PM.contigGCs, self.PM.contigLengths)
            self.bins[bid].printBin(self.PM.contigNames, self.PM.covProfiles, self.PM.contigGCs,
                                    self.PM.contigLengths, self.PM.isLikelyChimeric,
                                    outFormat=outFormat, separator=separator, stream=stream)

    def plotProfileDistributions(self):
        """Plot the coverage and kmer distributions for each bin"""
        for bid in self.getBids():
            self.bins[bid].plotProfileDistributions(self.PM.transformedCP, self.PM.kmerSigs, fileName="PROFILE_"+str(bid))

    def plotSelectBins(self,
                       bids,
                       plotMers=False,
                       fileName="",
                       plotEllipsoid=False,
                       ignoreContigLengths=False,
                       ET=None):
        """Plot a selection of bids in a window"""
        if plotEllipsoid and ET == None:
            ET = EllipsoidTool()

        # we need to do some fancy-schmancy stuff at times!
        if plotMers:
            num_cols = 2
        else:
            num_cols = 1

        fig = plt.figure(figsize=(6.5*num_cols, 6.5))
        ax = fig.add_subplot(1,num_cols,1, projection='3d')
        for i, bid in enumerate(bids):
            self.bins[bid].plotOnAx(ax,
                           self.PM.transformedCP,
                           self.PM.contigGCs,
                           self.PM.contigLengths,
                           self.PM.colorMapGC,
                           self.PM.isLikelyChimeric,
                           ET=ET,
                           printID=True,
                           ignoreContigLengths=ignoreContigLengths,
                           plotColorbar=(num_cols==1 and i==0)
                           )
        if plotMers:

            title = "Bin: %d : %d contigs : %s BP\n" % (bid, len(self.bins[bid].rowIndices), np_sum(self.PM.contigLengths[self.bins[bid].rowIndices]))
            title += "Coverage centroid: %d %d %d\n" % (np_median(self.PM.transformedCP[self.bins[bid].rowIndices,0]),
                                                                        np_median(self.PM.transformedCP[self.bins[bid].rowIndices,1]),
                                                                        np_median(self.PM.transformedCP[self.bins[bid].rowIndices,2]))
            title += "GC: median: %.4f stdev: %.4f\n" % (np_median(self.PM.contigGCs[self.bins[bid].rowIndices]), np_std(self.PM.contigGCs[self.bins[bid].rowIndices]))

            if self.PM.isLikelyChimeric[bid]:
                title += "Likely Chimeric"

            ax.set_title(title)

            ax = fig.add_subplot(1, 2, 2)
            for i, bid in enumerate(bids):
                self.bins[bid].plotMersOnAx(ax,
                                            self.PM.kmerPCs[:,0],
                                            self.PM.kmerPCs[:,1],
                                            self.PM.contigGCs,
                                            self.PM.contigLengths,
                                            self.PM.colorMapGC,
                                            ET=ET,
                                            printID=True,
                                            plotColorbar=(i==0)
                                            )
            ax.set_title('PCA of k-mer signature')

        fig.set_size_inches(6*num_cols, 6)
        if(fileName != ""):
            try:
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

    def plotMultipleBins(self,
                         bins,
                         untransformed=False,
                         semi_untransformed=False,
                         ignoreContigLengths=False,
                         squash=False):
        """plot a bunch of bins, used mainly for debugging"""
        ET = EllipsoidTool()
        fig = plt.figure()

        if untransformed or semi_untransformed:
            coords = self.PM.covProfiles
            et = None
            pc = False
        else:
            coords = self.PM.transformedCP
            et = ET
            pc = True


        if squash:
            # mix all the points together
            # we need to work out how to shape the plots
            num_plots = len(bins)
            plot_rows = float(int(np_sqrt(num_plots)))
            plot_cols = np_ceil(float(num_plots)/plot_rows)
            plot_num = 1
            for bids in bins:
                ax = fig.add_subplot(plot_rows, plot_cols, plot_num, projection='3d')
                disp_vals = np_array([])
                disp_lens = np_array([])
                gcs = np_array([])
                num_points = 0
                for bid in bids:
                    for row_index in self.bins[bid].rowIndices:
                        num_points += 1
                        disp_vals = np_append(disp_vals, coords[row_index])
                        disp_lens = np_append(disp_lens, np_sqrt(self.PM.contigLengths[row_index]))
                        gcs = np_append(gcs, self.PM.contigGCs[row_index])

                # reshape
                disp_vals = np_reshape(disp_vals, (num_points, 3))

                if ignoreContigLengths:
                    sc = ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors='none', c=gcs, cmap=self.PM.colorMapGC, vmin=0.0, vmax=1.0, s=10., marker='.')
                else:
                    sc = ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors='k', c=gcs, cmap=self.PM.colorMapGC, vmin=0.0, vmax=1.0, s=disp_lens, marker='.')
                sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect

                plot_num += 1
                
            cbar = plt.colorbar(sc, shrink=0.5)
            cbar.ax.tick_params()
            cbar.ax.set_title("% GC", size=10)
            cbar.set_ticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
            cbar.ax.set_ylim([0.15, 0.85])
            cbar.outline.set_ydata([0.15] * 2 + [0.85] * 4 + [0.15] * 3)
        else:
            # plot all separately
            # we need to work out how to shape the plots
            num_plots = len(bins) + 1
            plot_rows = float(int(np_sqrt(num_plots)))
            plot_cols = np_ceil(float(num_plots)/plot_rows)
            plot_num = 1
            ax = fig.add_subplot(plot_rows, plot_cols, plot_num, projection='3d')

            for bids in bins:
                for bid in bids:
                    self.bins[bid].plotOnAx(ax, coords, self.PM.contigGCs, self.PM.contigLengths, self.PM.colorMapGC, self.PM.isLikelyChimeric, ignoreContigLengths=ignoreContigLengths, ET=et, plotCentroid=pc)

            plot_num += 1
            if semi_untransformed:
                coords = self.PM.transformedCP
                et = ET
                pc = True

            for bids in bins:
                ax = fig.add_subplot(plot_rows, plot_cols, plot_num, projection='3d')
                plot_num += 1
                for bid in bids:
                    self.bins[bid].plotOnAx(ax, coords, self.PM.contigGCs, self.PM.contigLengths, self.PM.colorMapGC, self.PM.isLikelyChimeric, ignoreContigLengths=ignoreContigLengths, ET=et, plotCentroid=pc)

        try:
            plt.show()
        except:
            print "Error showing image:", exc_info()[0]
            raise

        plt.close(fig)
        del fig


    def plotBins(self,
                 FNPrefix="BIN",
                 sideBySide=False,
                 folder='',
                 plotEllipsoid=False,
                 ignoreContigLengths=False,
                 ET=None):
        """Make plots of all the bins"""
        if plotEllipsoid and ET == None:
            ET = EllipsoidTool()

        if folder != '':
            makeSurePathExists(folder)

        for bid in self.getBids():
            self.bins[bid].makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerNormPC1, self.PM.kmerPCs, self.PM.contigGCs, self.PM.contigLengths)

        if(sideBySide):
            print "Plotting side by side"
            self.plotSideBySide(self.bins.keys(), tag=FNPrefix, ignoreContigLengths=ignoreContigLengths)
        else:
            print "Plotting bins"
            for bid in self.getBids():
                if folder != '':
                    self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigGCs, self.PM.kmerNormPC1,
                                            self.PM.contigLengths, self.PM.colorMapGC, self.PM.isLikelyChimeric,
                                            fileName=osp_join(folder, FNPrefix+"_"+str(bid)),
                                            ignoreContigLengths=ignoreContigLengths, ET=ET)
                else:
                    self.bins[bid].plotBin(self.PM.transformedCP, self.PM.contigGCs, self.PM.kmerNormPC1,
                                              self.PM.contigLengths, self.PM.colorMapGC, self.PM.isLikelyChimeric,
                                              fileName=FNPrefix+"_"+str(bid), ignoreContigLengths=ignoreContigLengths, ET=ET)
                    
    def plotBinCoverage(self, plotEllipses=False, plotContigLengs=False, printID=False):
        """Make plots of all the bins"""

        print "Plotting first 3 stoits in untransformed coverage space"
         
        # plot contigs in coverage space
        fig = plt.figure()
  
        if plotContigLengs:
            disp_lens = np_sqrt(self.PM.contigLengths)
        else:
            disp_lens = 30
            
        # plot contigs in kmer space
        ax = fig.add_subplot(121, projection='3d')
        ax.set_xlabel('kmer PC1')
        ax.set_ylabel('kmer PC2')
        ax.set_zlabel('kmer PC3')
        ax.set_title('kmer space')
        
        sc = ax.scatter(self.PM.kmerPCs[:,0], self.PM.kmerPCs[:,1], self.PM.kmerPCs[:,2], edgecolors='k', c=self.PM.contigGCs, cmap=self.PM.colorMapGC, vmin=0.0, vmax=1.0, s=disp_lens)
        sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect
        
        if plotEllipses:
            ET = EllipsoidTool()
            for bid in self.getBids():
                row_indices = self.bins[bid].rowIndices
                (center, radii, rotation) = self.bins[bid].getBoundingEllipsoid(self.PM.kmerPCs[:, 0:3], ET=ET)
                centroid_gc = np_mean(self.PM.contigGCs[row_indices])
                centroid_color = self.PM.colorMapGC(centroid_gc)
                if printID:
                    ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=False, cageColor=centroid_color, label=self.id)
                else:
                    ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=False, cageColor=centroid_color)
        
        # plot contigs in untransformed coverage space
        ax = fig.add_subplot(122, projection='3d')
        ax.set_xlabel('coverage 1')
        ax.set_ylabel('coverage 2')
        ax.set_zlabel('coverage 3')
        ax.set_title('coverage space')
        
        sc = ax.scatter(self.PM.covProfiles[:,0], self.PM.covProfiles[:,1], self.PM.covProfiles[:,2], edgecolors='k', c=self.PM.contigGCs, cmap=self.PM.colorMapGC, vmin=0.0, vmax=1.0, s=disp_lens)
        sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect
        
        cbar = plt.colorbar(sc, shrink=0.5)
        cbar.ax.tick_params()
        cbar.ax.set_title("% GC", size=10)
        cbar.set_ticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        cbar.ax.set_ylim([0.15, 0.85])
        cbar.outline.set_ydata([0.15] * 2 + [0.85] * 4 + [0.15] * 3)

        if plotEllipses:
            ET = EllipsoidTool()
            for bid in self.getBids():
                row_indices = self.bins[bid].rowIndices
                (center, radii, rotation) = self.bins[bid].getBoundingEllipsoid(self.PM.covProfiles[:, 0:3], ET=ET)
                centroid_gc = np_mean(self.PM.contigGCs[row_indices])
                centroid_color = self.PM.colorMapGC(centroid_gc)
                if printID:
                    ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=False, cageColor=centroid_color, label=self.id)
                else:
                    ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=False, cageColor=centroid_color)

        try:
            plt.show()
            plt.close(fig)
        except:
            print "Error showing image", exc_info()[0]
            raise
        
        del fig
            
    def plotSideBySide(self, bids, fileName="", tag="", use_elipses=True, ignoreContigLengths=False):
        """Plot two bins side by side in 3d"""
        if use_elipses:
            ET = EllipsoidTool()
        else:
            ET = None
        fig = plt.figure()

        # get plot extents
        xMin = 1e6
        xMax = 0
        yMin = 1e6
        yMax = 0
        zMin = 1e6
        zMax = 0

        for bid in bids:
            x = self.PM.transformedCP[self.bins[bid].rowIndices,0]
            y = self.PM.transformedCP[self.bins[bid].rowIndices,1]
            z = self.PM.transformedCP[self.bins[bid].rowIndices,2]

            xMin = min(min(x), xMin)
            xMax = max(max(x), xMax)

            yMin = min(min(y), yMin)
            yMax = max(max(y), yMax)

            zMin = min(min(z), zMin)
            zMax = max(max(z), zMax)

        # we need to work out how to shape the plots
        num_plots = len(bids)
        plot_rows = float(int(np_sqrt(num_plots)))
        plot_cols = np_ceil(float(num_plots)/plot_rows)
        for plot_num, bid in enumerate(bids):
            title = self.bins[bid].plotOnFig(fig, plot_rows, plot_cols, plot_num+1,
                                              self.PM.transformedCP, self.PM.contigGCs, self.PM.contigLengths,
                                              self.PM.colorMapGC, self.PM.isLikelyChimeric, ET=ET, fileName=fileName,
                                              plotColorbar=(plot_num == len(bids)-1), extents=[xMin, xMax, yMin, yMax, zMin, zMax],
                                              ignoreContigLengths=ignoreContigLengths)

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

    def plotStoitNames(self, ax):
        """Plot stoit names on an existing axes"""
        self.PM.plotStoitNames(ax)

    def plotBinIds(self, gc_range=None, ignoreRanges=False, showChimeric=False):
        """Render 3d image of core ids"""
        (bin_centroid_points, bin_centroid_colors, _bin_centroid_gc, bids) = self.findCoreCentres(gc_range=gc_range,
                                                                                                  processChimeric=showChimeric)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlabel('x coverage')
        ax.set_ylabel('y coverage')
        ax.set_zlabel('z coverage')

        outer_index = 0
        for bid in bids:
            ax.text(bin_centroid_points[outer_index,0],
                    bin_centroid_points[outer_index,1],
                    bin_centroid_points[outer_index,2],
                    str(int(bid)),
                    color=bin_centroid_colors[outer_index]
                    )
            outer_index += 1

        if ignoreRanges:
            mm = np_max(bin_centroid_points, axis=0)
            ax.set_xlim3d(0, mm[0])
            ax.set_ylim3d(0, mm[1])
            ax.set_zlim3d(0, mm[2])

        else:
            self.plotStoitNames(ax)
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

    def plotBinPoints(self, ignoreRanges=False, plotColorbar=True, showChimeric=False):
        """Render the image for validating cores"""
        (bin_centroid_points, _bin_centroid_colors, bin_centroid_gc, _bids) = self.findCoreCentres(processChimeric=showChimeric)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        print bin_centroid_gc
        sc = ax.scatter(bin_centroid_points[:,0], bin_centroid_points[:,1], bin_centroid_points[:,2], edgecolors='k', c=bin_centroid_gc, cmap=self.PM.colorMapGC, vmin=0.0, vmax=1.0)
        sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect

        if plotColorbar:
            cbar = plt.colorbar(sc, shrink=0.5)
            cbar.ax.tick_params()
            cbar.ax.set_title("% GC", size=10)
            cbar.set_ticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
            cbar.ax.set_ylim([0.15, 0.85])
            cbar.outline.set_ydata([0.15] * 2 + [0.85] * 4 + [0.15] * 3)

        ax.set_xlabel('x coverage')
        ax.set_ylabel('y coverage')
        ax.set_zlabel('z coverage')

        if not ignoreRanges:
            self.plotStoitNames(ax)
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

    def plotTDist(self, scores, testPoint):
        """DEBUG: Plot the distribution of the NULL T disytribution"""
        B = np_arange(0, len(scores), 1)
        co_90 = [int(float(len(scores))*0.90)] * 100
        co_95 = [int(float(len(scores))*0.95)] * 100
        co_97 = [int(float(len(scores))*0.97)] * 100
        co_99 = [int(float(len(scores))*0.99)] * 100
        smin = np_min(scores)
        smax = np_max(scores)
        step = (smax-smin)/100
        LL = np_arange(smin, smax, step)

        fig = plt.figure()

        plt.plot(co_90, LL, 'g-')
        plt.plot(co_95, LL, 'g-')
        plt.plot(co_97, LL, 'g-')
        plt.plot(co_99, LL, 'g-')

        plt.plot(B, scores, 'r-')
        plt.plot(B, [testPoint]*len(scores), 'b-')
        plt.show()
        plt.close(fig)
        del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

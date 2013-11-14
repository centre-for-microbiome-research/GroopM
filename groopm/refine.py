#!/usr/bin/env python
###############################################################################
#                                                                             #
#    refine.py                                                                #
#                                                                             #
#    A collection of classes / methods used when bin refinment and expansion  #
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
__version__ = "0.2.8"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Alpha"

###############################################################################

import sys
from sys import stdout as sys_stdout

import math

from colorsys import hsv_to_rgb as htr
import matplotlib.pyplot as plt
from numpy import (abs as np_abs,
                   append as np_append,
                   arange as np_arange,
                   arccos as np_arccos,
                   argmin as np_argmin,
                   argsort as np_argsort,
                   around as np_around,
                   array as np_array,
                   concatenate as np_concatenate,
                   copy as np_copy,
                   dot as np_dot,
                   max as np_max,
                   mean as np_mean,
                   median as np_median,
                   min as np_min,
                   ones as np_ones,
                   reshape as np_reshape,
                   seterr as np_seterr,
                   sqrt as np_sqrt,
                   std as np_std,
                   sum as np_sum,
                   where as np_where,
                   zeros as np_zeros)
from numpy.linalg import norm as np_norm
from numpy.random import (randint as randint,
                          shuffle as shuffle)

from scipy.spatial import KDTree as kdt
from scipy.spatial.distance import cdist, squareform, pdist

# GroopM imports
from binManager import BinManager
from ellipsoid import EllipsoidTool
from PCA import PCA, Center
import groopmExceptions as ge
from som import SOM
np_seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# dimension of our SOM
SOMDIM = 4 # trans cov  + 1st pca kmers

class RefineEngine:
    """Workhorse wrapper for bin refinement"""
    def __init__(self,
                 timer,
                 BM=None,
                 dbFileName="",
                 transform=True,
                 getUnbinned=False,
                 loadContigNames=False,
                 cutOff=0,
                 bids=[]):

        # worker classes
        if BM is None:
            # make our own ones from scratch
            self.BM = BinManager(dbFileName=dbFileName)
            self.BM.loadBins(timer,
                             bids=bids,
                             makeBins=True,
                             silent=False,
                             cutOff=cutOff,
                             loadContigNames=loadContigNames,
                             getUnbinned=getUnbinned,
                             transform=transform)
        else:
            self.BM = BM

        self.PM = self.BM.PM

        # make pretty ellipses for fun and profit
        self.ET = EllipsoidTool()

        # pay attention to contig lengths when including in bins use grubbs test
        self.GT = GrubbsTester()

        self.transform = transform        # are we going to transform the data

#------------------------------------------------------------------------------
# REFINING

    def refineBins(self, timer, auto=False, saveBins=False, plotFinal="", gf=""):
        """Iterative wrapper for the refine function"""
        if self.transform:      # do we want to plot to 1000*1000*1000
            ignoreRanges=False
        else:
            ignoreRanges=True

        if auto:
            print "    Start automatic bin refinement"
            num_binned = len(self.PM.binnedRowIndices.keys())
            perc = "%.2f" % round((float(num_binned)/float(self.PM.numContigs))*100,2)
            print "   ",num_binned,"contigs across",len(self.BM.bins.keys()),"cores (",perc,"% )"

            graph = self.autoRefineBins(timer, makeGraph=gf!="")
            if graph is not None:
                print "    Writing graph to:", gf
                try:
                    with open(gf, "w") as gv_fh:
                        gv_fh.write(graph)
                except:
                    print "Error writing graph to:", gf

            num_binned = len(self.PM.binnedRowIndices.keys())
            perc = "%.2f" % round((float(num_binned)/float(self.PM.numContigs))*100,2)
            print "   ",num_binned,"contigs across",len(self.BM.bins.keys()),"cores (",perc,"% )"

            if plotFinal != "":
                bids = self.BM.getBids()
                for bid in bids:
                    self.BM.bins[bid].makeBinDist(self.PM.transformedCP,
                                               self.PM.averageCoverages,
                                               self.PM.kmerNormPC1,
                                               self.PM.kmerPCs,
                                               self.PM.contigGCs,
                                               self.PM.contigLengths)
                self.BM.plotBins(FNPrefix=plotFinal, ET=self.ET)

            if saveBins:
                self.BM.saveBins(nuke=True)
        else:
            self.plotterRefineBins(ignoreRanges=ignoreRanges)

    def plotterRefineBins(self, ignoreRanges=False):
        """combine similar bins using 3d plots"""
        use_elipses = True
        show_chimeric_bins = False
        ET = self.ET
        self.printRefinePlotterInstructions()
        #self.BM.plotBinIds(ignoreRanges=ignoreRanges)
        continue_merge = True
        while(continue_merge):
            user_option = self.promptOnPlotterRefine()

            if(user_option == 'Q'):
                print '\nBye!'
                return

            elif(user_option == 'C'):
                print "Select colormap:"
                print "  1. HSV"
                print "  2. Accent"
                print "  3. Blues"
                print "  4. Spectral"
                print "  5. Grayscale"
                print "  6. Discrete (14 colors)"
                print "  7. Discrete paired (14 colors)"

                bValid = False
                while(not bValid):
                    try:
                        colormap_id = int(raw_input(" Enter colormap number (e.g., 1): "))
                        if colormap_id < 1 or colormap_id > 7:
                            raise ValueError('Invalid colormap id.')
                        bValid = True
                    except ValueError:
                        print "Colormap must be specified as a number between 1 and 7."

                if colormap_id == 1:
                    self.PM.setColorMap('HSV')
                elif colormap_id == 2:
                    self.PM.setColorMap('Accent')
                elif colormap_id == 3:
                    self.PM.setColorMap('Blues')
                elif colormap_id == 4:
                    self.PM.setColorMap('Spectral')
                elif colormap_id == 5:
                    self.PM.setColorMap('Grayscale')
                elif colormap_id == 6:
                    self.PM.setColorMap('Discrete')
                elif colormap_id == 7:
                    self.PM.setColorMap('DiscretePaired')

            elif(user_option == 'E'):
                if use_elipses:
                    ET = None
                    use_elipses = False
                    print "\nEllipses off"
                else:
                    ET = self.ET
                    use_elipses = True
                    print "\nEllipses on"

            elif(user_option == 'X'):
                if show_chimeric_bins:
                    show_chimeric_bins = False
                    print "\nHiding likely chimeric bins."
                else:
                    show_chimeric_bins = True
                    print "\nShowing likely chimeric bins."

            elif(user_option == 'R'):
                self.BM.plotBinIds(ignoreRanges=ignoreRanges, showChimeric=show_chimeric_bins)

            elif(user_option == 'P'):
                self.BM.plotBinPoints(ignoreRanges=ignoreRanges, showChimeric=show_chimeric_bins)
            elif(user_option == 'U'):
                self.BM.plotBinCoverage(plotEllipses = use_elipses)
            elif(user_option == 'M'):
                # merge bins
                merge_bids = self.BM.getPlotterMergeIds()
                if(not len(merge_bids) == 0):
                    self.BM.merge(merge_bids,
                                  auto=False,
                                  manual=True,
                                  newBid=False,
                                  saveBins=True,
                                  verbose=False,
                                  printInstructions=False,
                                  use_elipses=use_elipses)

            elif(user_option == 'G'):
                # display a subset only!
                have_range = False
                while(not have_range):
                    try:
                        gc_range_str = raw_input(" Enter GC range to examine (e.g., 0.5-0.6): ")

                        if '-' not in gc_range_str:
                            raise ValueError('Incorrectly formatted GC range.')
                        else:
                            values = gc_range_str.split('-')
                            start = float(values[0])
                            end = float(values[1])
                            gc_range=[start,end]

                        have_range = True
                    except ValueError:
                        print "GC ranges must be entered as 'a-b' (e.g., 0.5-0.6)."
                self.BM.plotBinIds(gc_range=gc_range, ignoreRanges=ignoreRanges)

            elif(user_option == 'B'):
                # print subset of bins
                have_bid = False
                bids = []
                while(not have_bid):
                    have_bid = True
                    try:
                        usr_bids = raw_input(" Enter 'space' seperated bin id(s) to plot: ")
                        bids = [int(i) for i in usr_bids.split(" ")]
                        if bids == [-1]:
                            bids = self.BM.getBids()
                        else:
                            for bid in bids:
                                if bid not in self.BM.bins:
                                    print "ERROR: Bin %d not found!" % bid
                                    have_bid &= False
                    except ValueError:
                        print "You need to enter an integer value!"

                if len(bids) > 0:
                    self.BM.plotSelectBins(bids, plotMers=True, ET=ET)

            elif(user_option == 'S'):
                # split bins
                have_bid = False
                have_parts = False
                while(not have_bid):
                    try:
                        bid = int(raw_input(" Enter bid to split: "))
                        if bid not in self.BM.bins:
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
                self.BM.split(bid,
                              parts,
                              mode='kmer',
                              auto=False,
                              saveBins=True,
                              printInstructions=False,
                              use_elipses=use_elipses)

            elif(user_option == 'V'):
                """plot in vicinity of a bin"""
                have_bid = False
                have_radius = False
                while(not have_bid):
                    try:
                        bid = int(raw_input(" Enter bid of interest: "))
                        if bid not in self.BM.bins:
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
                self.BM.bins[bid].makeBinDist(self.PM.transformedCP,
                                              self.PM.averageCoverages,
                                              self.PM.kmerNormPC1,
                                              self.PM.kmerPCs,
                                              self.PM.contigGCs,
                                              self.PM.contigLengths)
                # now go point shopping
                disp_vals = np_array([])
                disp_gc = np_array([])
                disp_lens = np_array([])
                num_points = 0
                seen_bids = {}
                for row_index in range(len(self.PM.indices)):
                    if np_norm(self.PM.transformedCP[row_index] -
                               self.BM.bins[bid].covMedians) <= radius:
                        num_points += 1
                        disp_vals = np_append(disp_vals, self.PM.transformedCP[row_index])
                        disp_lens = np_append(disp_lens, np_sqrt(self.PM.contigLengths[row_index]))
                        disp_gc = np_append(disp_gc, self.PM.contigGCs[row_index])
                        try:
                            seen_bids[self.PM.binIds[row_index]].append(1)
                        except KeyError:
                            seen_bids[self.PM.binIds[row_index]] = [1]

                # reshape
                disp_vals = np_reshape(disp_vals, (num_points, 3))

                print " Points are located in bins:"
                for seen_bid in seen_bids:
                    print "    %d - %d occurances" % (seen_bid, len(seen_bids[seen_bid]))

                fig = plt.figure()
                ax = fig.add_subplot(1,1,1, projection='3d')
                sc = ax.scatter(disp_vals[:,0],
                           disp_vals[:,1],
                           disp_vals[:,2],
                           edgecolors='k',
                           c=disp_gc,
                           cmap=self.PM.colorMapGC,
                           s=disp_lens,
                           marker='.')
                sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect

                self.BM.bins[bid].plotOnAx(ax,
                                           self.PM.transformedCP,
                                           self.PM.contigGCs,
                                           self.PM.contigLengths,
                                           self.PM.colorMapGC,
                                           self.PM.isLikelyChimeric,
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

    def autoRefineBins(self,
                       timer,
                       mergeSimilarBins=True,
                       removeDuds=False,
                       markLikelyChimeric=True,
                       shuffleRefine=True,
                       verbose=False,
                       plotAfterOB=False,
                       makeGraph=False):
        """Automagically refine bins"""

        sys_stdout.flush()
        graph = None
        if makeGraph:
            # lets make a graph
            graph = [{}, []]

            for bid in self.BM.getBids():
                bin = self.BM.bins[bid]
                centroid_gc = np_mean(self.PM.contigGCs[bin.rowIndices])
                centroid_color = self.PM.colorMapGC(centroid_gc)
                ncc = [int(i) for i in centroid_color * 255]
                hex_color = '#%02x%02x%02x' % (ncc[0], ncc[1], ncc[2])
                graph[0][bid] = '\t%d [fontcolor="%s" color="%s"];\n' % (bid, hex_color, hex_color)

        # identify and remove outlier bins
        if markLikelyChimeric:
            nuked = self.markLikelyChimericBins()
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()

            if makeGraph:
                # make all these uys point at the deleted bin
                if len(nuked) > 0:
                    graph[0][-1] = '\t"OUTLIERS" [fontcolor="#FF0000" color="#000000"];\n'
                for bid in nuked:
                    graph[1].append('\t%d -> "OUTLIERS";' % (bid))

        # merge bins together
        if mergeSimilarBins:
            self.mergeSimilarBins(graph=graph, verbose=False)
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()

        if plotAfterOB:
            bids = self.BM.getBids()
            for bid in bids:
                self.BM.bins[bid].makeBinDist(self.PM.transformedCP,
                                              self.PM.averageCoverages,
                                              self.PM.kmerNormPC1,
                                              self.PM.kmerPCs,
                                              self.PM.contigGCs,
                                              self.PM.contigLengths)
            self.BM.plotBins(FNPrefix="AFTER_OB", ET=self.ET)
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()

        if shuffleRefine:
            nuked = self.shuffleRefineContigs(timer)
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()
            if makeGraph:
                # Make sure we know these guys were deleted
                if len(nuked) > 0:
                    graph[0][-2] = '\t"SHUFFLED" [fontcolor="#FF0000" color="#000000"];\n'
                for bid in nuked:
                    graph[1].append('\t%d -> "SHUFFLED";' % (bid))

        if removeDuds:
            nuked = self.removeDuds()
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()
            if makeGraph:
                # Make sure we know these guys were deleted
                if len(nuked) > 0:
                    graph[0][-3] = '\t"DUDS" [fontcolor="#FF0000" color="#000000"];\n'
                for bid in nuked:
                    graph[1].append('\t%d -> "DUDS";' % (bid))

        if makeGraph:
            return self.writeGV(graph)

    def markLikelyChimericBins(self, verbose=False):
        """ Identify bins which contain mixed genomes based on GC.
              Small bins are nuked, large bins are flagged as chimeric. """
        print "    Identifying possible chimeric bins"
        sys_stdout.flush()

        # first we need to build a distribution!
        gc_stdev_distrb = []
        for bid in self.BM.getBids():
            self.BM.bins[bid].makeBinDist(self.PM.transformedCP,
                                       self.PM.averageCoverages,
                                       self.PM.kmerNormPC1,
                                       self.PM.kmerPCs,
                                       self.PM.contigGCs,
                                       self.PM.contigLengths)
            gc_stdev_distrb.append(self.BM.bins[bid].gcStdev)

        # now we work out the distribution of stdevs
        stdstd = np_std(gc_stdev_distrb)
        stdmean = np_mean(gc_stdev_distrb)

        dead_bins = []
        num_chimeric_bins = 0
        for bid in self.BM.getBids():
            Z = (self.BM.bins[bid].gcStdev - stdmean)/stdstd
            if Z > 2:
                if self.BM.bins[bid].totalBP < 100000:
                    dead_bins.append(bid)
                else:
                    self.PM.isLikelyChimeric[bid] = True
                num_chimeric_bins += 1
            else:
                self.PM.isLikelyChimeric[bid] = False

        # delete the bad bins
        self.BM.deleteBins(dead_bins,
                           force=True,
                           freeBinnedRowIndices=True,
                           saveBins=False)

        print "    Identified %d likely chimeric bin(s), removed %d small chimeric bin(s)" % (num_chimeric_bins, len(dead_bins))
        print "    %s" % ",".join(str(u) for u in dead_bins)
        return dead_bins

    def mergeSimilarBins(self, verbose=False, graph=None, silent=False):
        """Merge bins which are just crying out to be merged!"""
        kCutMedian, kCutStd = self.getKCut()
        cCutMedian, cCutStd = self.getCCut()

        orig_num_bins = len(self.BM.getNonChimericBinIds())

        if not silent:
            print "    Merging similar bins (%d) with kCut %0.2f (+/-%0.3f) cCut %0.2f (+/-%0.3f)" % (orig_num_bins, kCutMedian, kCutStd, cCutMedian, cCutStd)

        # identify merging groups and then merge them
        mergers = self.findMergeGroups(kCutMedian, kCutStd, cCutMedian, cCutStd, verbose=verbose)

        num_bins_removed = 0
        for merge in mergers:
            bins_removed = self.combineMergers(merge, kCutMedian, kCutStd, cCutMedian, cCutStd, graph=graph)
            num_bins_removed += len(bins_removed)
        if not silent:
            print "    Merged %d of %d cores leaving %d cores total" % (num_bins_removed, orig_num_bins, len(self.BM.getNonChimericBinIds()))

        return num_bins_removed

    def findMergeGroups(self, kCutMedian, kCutStd, cCutMedian, cCutStd, verbose=False):
        """Identify groups of contigs which could be merged"""
        cov_tdm = []                # these are used in the neighbor search
        kmer_tdm = []

        bid_2_tdm_index = {}
        tdm_index_2_bid = {}

        K = 6                   # number of neighbours to test,
                                # 6 seems reasonable since this allows
                                # for 1 neighbour per side in a 3D space
        if K > len(self.BM.getNonChimericBinIds()):
            K = len(self.BM.getNonChimericBinIds())

        # keep track of what gets merged where
        merged_bins = {}                # oldId => newId
        processed_pairs = {}            # keep track of pairs we've analysed
        bin_c_lengths = {}              # bid => [len,len,...]
        bin_c_ellipsoid_volumes = {}    # volume of the minimum bounding COVERAGE ellipsoid
        bin_c_ellipsoids = {}           # the matrix A representing the bins COVERAGE ellipsiod
        bin_k_ellipse_areas = {}        # area of the minimum bounding KMER ellipse
        bin_k_ellipses = {}             # the matrix A representing the bins KMER ellipse

#-----
# PREP DATA STRUCTURES

        index = 0
        for bid in self.BM.getNonChimericBinIds():
            bin = self.BM.bins[bid]

            bin.makeBinDist(self.PM.transformedCP,
                            self.PM.averageCoverages,
                            self.PM.kmerNormPC1,
                            self.PM.kmerPCs,
                            self.PM.contigGCs,
                            self.PM.contigLengths)

            # build coverage and kmer vectors for each bin (normalized as required)
            cov_tdm.append(bin.covMedians)
            kmer_tdm.append(bin.kMedian)

            # work out the volume of the minimum bounding coverage ellipsoid and kmer ellipse
            (bin_c_ellipsoids[bid], bin_c_ellipsoid_volumes[bid]) = bin.getBoundingCEllipsoidVol(self.PM.transformedCP,
                                                                                                 ET=self.ET,
                                                                                                 retA=True)

            BP = self.PM.kmerPCs[bin.rowIndices,0:3]
            (bin_k_ellipses[bid], bin_k_ellipse_areas[bid]) = bin.getBoundingKEllipseArea(BP,
                                                                                          ET=self.ET,
                                                                                          retA=True)

            bin_c_lengths[bid] = self.PM.contigLengths[bin.rowIndices]

            bid_2_tdm_index[bid] = index
            tdm_index_2_bid[index] = bid
            index += 1

            # we will not process any pair twice.
            # we also wish to avoid checking if a bin will merge with itself
            processed_pairs[self.BM.makeBidKey(bid, bid)] = True

#-----
# ALL Vs ALL

        # make a search tree from whitened coverage medians and kmer medians
        cp_cov_tdm = np_copy(cov_tdm)
        c_mean_tdm = np_mean(cp_cov_tdm, axis=0)
        c_std_tdm = np_std(cp_cov_tdm, axis=0)
        c_std_tdm += np_where(c_std_tdm == 0, 1, 0) # make sure std dev is never zero
        c_whiten_tdm = (cp_cov_tdm - c_mean_tdm) / c_std_tdm

        cov_search_tree = kdt(c_whiten_tdm)
        kmer_search_tree = kdt(kmer_tdm)
        for bid in self.BM.getNonChimericBinIds():
            # get the base bid and trace the chain up through mergers...
            merged_base_bid = bid
            base_bid = bid
            while merged_base_bid in merged_bins:
                merged_base_bid = merged_bins[merged_base_bid]
            base_bin = self.BM.bins[base_bid]

            # get the K closest bins in coverage and kmer space
            cov_neighbor_list = [tdm_index_2_bid[i] for i in cov_search_tree.query(c_whiten_tdm[bid_2_tdm_index[bid]], k=K)[1]]
            kmer_neighbor_list = [tdm_index_2_bid[i] for i in kmer_search_tree.query(kmer_tdm[bid_2_tdm_index[bid]], k=K)[1]]
            common_neighbors = set(cov_neighbor_list).intersection(set(kmer_neighbor_list))

            if verbose:
                print "++++++++++"
                print bid, cov_neighbor_list
                print bid, kmer_neighbor_list
                print bid, common_neighbors

            # test each neighbor in turn
            for i, neighbor_index in enumerate(common_neighbors):
                # get the query bid and trace the chain up through mergers...
                query_bid = neighbor_index
                merged_query_bid = query_bid
                while merged_query_bid in merged_bins:
                    merged_query_bid = merged_bins[merged_query_bid]

                if verbose:
                    print "++++++++++"
                    print base_bid, query_bid, merged_base_bid, merged_query_bid
#-----
# TIME WASTERS

                # process each BID pair once only (takes care of self comparisons too!)
                seen_key = self.BM.makeBidKey(base_bid, query_bid)
                if(seen_key in processed_pairs or merged_base_bid == merged_query_bid):
                    if verbose:
                        print "TW"
                    continue
                processed_pairs[seen_key] = True

                query_bin = self.BM.bins[query_bid]

#-----
# CONTIG LENGTH SANITY
                # Test the smaller bin against the larger
                if query_bin.binSize < base_bin.binSize:
                    lengths_wrong = self.GT.isMaxOutlier(np_median(bin_c_lengths[query_bid]),
                                                         bin_c_lengths[base_bid]
                                                         )
                else:
                    lengths_wrong = self.GT.isMaxOutlier(np_median(bin_c_lengths[base_bid]),
                                                         bin_c_lengths[query_bid]
                                                         )
                if lengths_wrong:
                    if verbose:
                        print "LW"
                    continue

#-----
# K and C SPACE SIMILARITY CHECK
                # If the bins are highly similar in their coverage and kmer distances
                # compared to other core bins than just merge them now
                k_dist_bw = self.kDistBetweenBins(base_bin, query_bin)
                c_dist_bw = self.cDistBetweenBins(base_bin, query_bin)

                if verbose:
                    print 'k_dist_bw, c_dist_bw'
                    print k_dist_bw, c_dist_bw
                    print '---------------------'


                if k_dist_bw < kCutMedian and c_dist_bw < cCutMedian:
                    if verbose:
                        print 'MERGED'
                        print '---------------------'

                    if merged_query_bid < merged_base_bid:
                        merged_bins[merged_base_bid] = merged_query_bid
                        break # we just nuked the base bid
                    else:
                        merged_bins[merged_query_bid] = merged_base_bid
                        continue

#-----
# KMER ELLIPSE OVERLAP
                if (bin_k_ellipse_areas[base_bid] <= bin_k_ellipse_areas[query_bid]):
                    INTT = self.ET.doesIntersect3D(bin_k_ellipses[query_bid][0],
                                                   bin_k_ellipses[query_bid][1],
                                                   bin_k_ellipses[base_bid][0],
                                                   bin_k_ellipses[base_bid][1])
                else:
                    INTT = self.ET.doesIntersect3D(bin_k_ellipses[base_bid][0],
                                                   bin_k_ellipses[base_bid][1],
                                                   bin_k_ellipses[query_bid][0],
                                                   bin_k_ellipses[query_bid][1])

                if verbose:
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1)
                    base_bin.plotMersOnAx(ax,
                                          self.PM.kmerPCs[:,0],
                                          self.PM.kmerPCs[:,1],
                                          self.PM.contigGCs,
                                          self.PM.contigLengths,
                                          self.PM.colorMapGC,
                                          ET=self.ET)
                    query_bin.plotMersOnAx(ax,
                                           self.PM.kmerPCs[:,0],
                                           self.PM.kmerPCs[:,1],
                                           self.PM.contigGCs,
                                           self.PM.contigLengths,
                                           self.PM.colorMapGC,
                                           ET=self.ET)
                    plt.title("MERGE: %d -> %d (%d)" % (base_bid, query_bid, INTT))
                    plt.show()
                    plt.close(fig)
                    del fig

                if not INTT:
                    if verbose:
                        print "KINTT"
                    continue
#-----
# MINIMUM BOUNDING COVERAGE ELLIPSOID

                # determine if intersection exists
                if bin_c_ellipsoid_volumes[base_bid] <= bin_c_ellipsoid_volumes[query_bid]:
                    intersects = self.ET.doesIntersect3D(bin_c_ellipsoids[query_bid][0],
                                                         bin_c_ellipsoids[query_bid][1],
                                                         bin_c_ellipsoids[base_bid][0],
                                                         bin_c_ellipsoids[base_bid][1])
                else:
                    intersects = self.ET.doesIntersect3D(bin_c_ellipsoids[base_bid][0],
                                                         bin_c_ellipsoids[base_bid][1],
                                                         bin_c_ellipsoids[query_bid][0],
                                                         bin_c_ellipsoids[query_bid][1])

                if verbose:
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1, projection='3d')
                    base_bin.plotOnAx(ax, self.PM.transformedCP, self.PM.contigGCs, self.PM.contigLengths, self.PM.colorMapGC, self.PM.isLikelyChimeric, ET=self.ET)
                    query_bin.plotOnAx(ax, self.PM.transformedCP, self.PM.contigGCs, self.PM.contigLengths, self.PM.colorMapGC, self.PM.isLikelyChimeric, ET=self.ET)
                    plt.title("MERGE: %d -> %d (%d)" % (base_bid, query_bid, intersects))
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
# CREATE FINAL MERGE GROUPS

        # now make a bunch of possible mergers
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

        return mergers

    def combineMergers(self, bidList, kCutMedian, kCutStd, cCutMedian, cCutStd, graph=None):
        """Merge similar bins in the given list"""
        merged_bids = []

        while len(bidList) > 1:
            # sort bins by length in bp
            length_and_bid = []
            for bid in bidList:
                length_and_bid.append([self.BM.getBin(bid).totalBP,bid])

            length_and_bid.sort(reverse=True)
            sorted_bid = [x[1] for x in length_and_bid]

            # find distance from largest bin to each putative bin fragment
            cur_bid = sorted_bid[0]
            cur_bin = self.BM.getBin(cur_bid)

            dists = []
            for i in xrange(1, len(sorted_bid)):
                frag_bid = sorted_bid[i]
                frag_bin = self.BM.getBin(frag_bid)

                k_dist_bw = self.kDistBetweenBins(cur_bin, frag_bin)
                c_dist_bw = self.cDistBetweenBins(cur_bin, frag_bin)

                z_dist = (k_dist_bw - kCutMedian) / kCutStd
                z_dist += (c_dist_bw - cCutMedian) / cCutStd

                dists.append([z_dist, k_dist_bw, c_dist_bw, frag_bid])

            # check if closest fragment should be merged
            dists.sort()

            closest_frag_kdist = dists[0][1]
            closest_frag_cdist = dists[0][2]
            closest_frag_bid = dists[0][3]

            if closest_frag_kdist < (kCutMedian + 2*kCutStd) and closest_frag_cdist < (cCutMedian + 2*cCutStd):
                # merge bins
                merged_bids.append(closest_frag_bid)
                bidList.remove(closest_frag_bid)
                self.BM.merge([cur_bid, closest_frag_bid],
                              auto=True,
                              manual=False,
                              newBid=False,
                              saveBins=False,
                              verbose=False,
                              printInstructions=False,
                              use_elipses=False)
            else:
                bidList.remove(cur_bid)

        return merged_bids

    def small2indices(self, index, side):
        """Return the indices of the comparative items
        when given an index into a condensed distance matrix
        """
        step = 0
        while index >= (side-step):
            index = index - side + step
            step += 1
        return (step, step + index + 1)

    def combineMergersMike(self, bidList, kCut, cCut, graph=None):
        """Try to merge similar bins in the given list"""

        merged_bids = []
        too_big = 10000

        # PCA kmers to find out who is most similar to whom
        (bin_mer_PCAs, mer_con_PCAs) = self.rePCA(bidList, doBoth=True)
        side = len(bidList)
        sq_dists = cdist(bin_mer_PCAs, bin_mer_PCAs, 'cityblock')
        dists = squareform(sq_dists)

        # raw coverage averages for each bin
        raw_coverage_centroids = {}

        while True:
            # find the closest pair
            closest = np_argmin(dists)
            if dists[closest] == too_big:
                break
            (i,j) = self.PM.small2indices(closest, side-1)
            bid1 = bidList[i]
            bid2 = bidList[j]
            should_merge = False

            # test if the mer dist is teensy tiny.
            # this is a time saver...
            k_diff = np_median(cdist(self.PM.kmerPCs[self.BM.bins[bid1].rowIndices], self.PM.kmerPCs[self.BM.bins[bid2].rowIndices], 'cityblock'))

            #if VVB:
            #    print bid1, bid2, k_diff,
            if k_diff <= kCut:
                should_merge = True
            else:
                # not teensy tiny, but is it still OK?
                (test, null) = self.testMergeMer(bid1, bid2, mer_con_PCAs)
                should_merge = test < null

            if should_merge:
                # now we know we should merge based on mers,
                # test if the bins lie in the same coverage region
                # get the angle between the two bins
                try:
                    c1 = raw_coverage_centroids[bid1]
                except KeyError:
                    c1 = np_mean([self.PM.covProfiles[row_index] for row_index in self.BM.bins[bid1].rowIndices], axis=0)
                    raw_coverage_centroids[bid1] = c1

                try:
                    c2 = raw_coverage_centroids[bid2]
                except KeyError:
                    c2 = np_mean([self.PM.covProfiles[row_index] for row_index in self.BM.bins[bid2].rowIndices], axis=0)
                    raw_coverage_centroids[bid2] = c2
                try:
                    ang = np_arccos(np_dot(c1,c2) / np_norm(c1) / np_norm(c2))
                except FloatingPointError:
                    ang = 0.0

                if ang > cCut:
                    # if the angle between is not teensy timy
                    (test, null) = self.testMergeCoverage(bid1, bid2)
                    if test >= null:
                        should_merge = False

            if should_merge:
                # get these before we merge!
                if graph is not None:
                    graph[1].append("\t%d -> %d;" % (bid2, bid1))
                b1_size = self.BM.bins[bid1].binSize
                b2_size = self.BM.bins[bid2].binSize

                # merge!
                merged_bids.append(bid2)
                self.BM.merge([bid1, bid2],
                              auto=True,
                              manual=False,
                              newBid=False,
                              saveBins=False,
                              verbose=False,
                              printInstructions=False,
                              use_elipses=False)

                # we use the weighted average of the two previous pca positions
                # to determine where the newly merged bin should reside
                bin_mer_PCAs[i] = (bin_mer_PCAs[i] * b1_size + bin_mer_PCAs[j] * b2_size) / (b1_size + b2_size)
                raw_coverage_centroids[bid1] = (raw_coverage_centroids[bid1] * b1_size + raw_coverage_centroids[bid2] * b2_size) / (b1_size + b2_size)

                # re-calc the distances
                new_dists = cdist([bin_mer_PCAs[i]], bin_mer_PCAs, 'cityblock')

                # we need to fix the distance matrix
                sq_dists[j,:] = too_big
                sq_dists[:,j] = too_big
                sq_dists[j,j] = 0.0
                sq_dists[i,:] = np_where(sq_dists[i,:] == too_big, sq_dists[i,:], new_dists)
                sq_dists[:,i] = np_where(sq_dists[:,i] == too_big, sq_dists[:,i], new_dists)
                dists = squareform(sq_dists)
            else:
                # we won't check this again
                sq_dists[i,j] = too_big
                sq_dists[j,i] = too_big
                dists = squareform(sq_dists)

        return merged_bids

    def buildSOM(self,
                 timer,
                 maskBoundaries=False,
                 defineBins=False,
                 retrain=False,
                 render=False,
                 silent=False,
                 animateFilePrefix=""):
        """Build, train and return a SOM for the given bids"""
        # produce the actual training data
        # this is the bin centroid values
        bids = self.BM.getBids()
        #bids = self.BM.getNonChimericBinIds()
        training_data = np_zeros((len(bids), SOMDIM))
        i = 0
        for bid in bids:
            training_data[i,:-1] = np_mean(self.PM.transformedCP[self.BM.bins[bid].rowIndices], axis=0)
            training_data[i,-1] = np_mean(self.PM.kmerNormPC1[self.BM.bins[bid].rowIndices], axis=0)
            i += 1

        som_side = np_max([100, len(bids)*5])

        # normalise the data so it fits between 0 and 1
        # but make sure that the max global CP and mer values are
        # used to scale
        minz = np_zeros((SOMDIM))
        minz[:-1] = np_min(self.PM.transformedCP, axis=0)
        minz[-1] = np_min(self.PM.kmerNormPC1, axis=0)

        maxz = np_zeros((SOMDIM))
        maxz[:-1] = np_max(self.PM.transformedCP, axis=0)
        maxz[-1] = np_max(self.PM.kmerNormPC1, axis=0)

        maxz -= minz
        training_data -= minz
        training_data /= maxz

        # during training, use the max and min vals of
        # the actual training data
        tminz = np_min(training_data, axis=0)
        tmaxz = np_max(training_data, axis=0)

        # set training in motion
        SS = SOM(som_side, SOMDIM, lc=tminz, uc=tmaxz)
        SS.train(training_data,
                 influenceRate=0.15,
                 iterations=800,
                 silent=silent,
                 weightImgFileNamePrefix=animateFilePrefix)
        print "    --"
        print "    %s" % timer.getTimeStamp()
        if render:
            SS.renderWeights("S1")

        if maskBoundaries:
            if not silent:
                print "    Creating boundary mask"
            # make a boundary mask
            if render:
                SS.makeBoundaryMask(plotMaskFile="S2.png")
            else:
                SS.makeBoundaryMask()
            SS.maskBoundaries()
        if defineBins:
            # assign regions on som surface to specific bins
            if not silent:
                print "    Defining bin regions"
            SS.defineBinRegions(bids, training_data, render=render)
            if render:
                SS.renderBoundaryMask("S5.png")
        if maskBoundaries:
            # mask out regions where we don't like it
            if not silent:
                print "    Masking SOM classifier"
            SS.maskBoundaries(addNoise=False, doFlat=True)
        if render:
            SS.renderWeights("S6")

        print "    %s" % timer.getTimeStamp()
        if retrain:
            # retrain bin regions using contigs from the bin
            if not silent:
                print "    Retraining SOM classifier"
            for i in range(len(bids)):
                bid = bids[i]
                sys_stdout.write("\r    Retraining on bin: %d (%d of %d)" % (bid, i+1, len(bids)))
                sys_stdout.flush()
                self.retrainSOM(SS,
                                bid,
                                SS.makeBinMask(training_data[i]),
                                som_side,
                                minz,
                                maxz,
                                silent=silent,
                                render=render)
            if render:
                SS.renderWeights("gg")
            print "    --"

        if render:
            SS.renderWeights("S7")
            pass

        return (SS, minz, maxz, som_side)

    def retrainSOM(self,
                   SS,
                   bid,
                   maskPoints,
                   som_side,
                   minz,
                   maxz,
                   silent=False,
                   render=False):
        """Further training of a SOM built using bin means"""
        bin = self.BM.bins[bid]

        # make a training set of just this node's contigs
        block = np_zeros((bin.binSize,SOMDIM))
        block[:,:-1] = self.PM.transformedCP[bin.rowIndices]
        block[:,-1] = self.PM.kmerNormPC1[bin.rowIndices]

        # global normalisation
        block -= minz
        block /= maxz

        # now we'd like to centre the weights and mask within an
        # appropriately sized square
        min_p = np_min(maskPoints.keys(), axis=0)
        max_p = np_max(maskPoints.keys(), axis=0)
        diffs = max_p - min_p
        small_side = np_min(diffs)
        sweights = np_copy(SS.weights.nodes[min_p[0]:min_p[0]+diffs[0]+1,min_p[1]:min_p[1]+diffs[1]+1])
        #SS.weights.renderSurface("C_%d.png"%bid, nodes=sweights)

        # shift and mask out all other bins
        shifted_mask_points = {}
        shifted_bin_mask = np_ones((diffs[0]+1,diffs[1]+1))
        for (r,c) in maskPoints.keys():
            shift = maskPoints[(r,c)] - min_p
            shifted_bin_mask[shift[0],shift[1]] = 0
            shifted_mask_points[(shift[0], shift[1])] = shift
        SS.maskBoundaries(weights=sweights, mask=shifted_bin_mask)

        # train on the set of dummy weights
        sweights = SS.train(block,
                            weights=sweights,
                            iterations=50,
                            mask=shifted_mask_points,
                            radius=small_side/3,
                            influenceRate=0.1)

        #SS.weights.renderSurface("D_%d.png"%bid, nodes=sweights)
        # update the torusMesh values appropriately
        for (r,c) in maskPoints.keys():
            shift = maskPoints[(r,c)] - min_p
            SS.weights.nodes[r,c] = sweights[shift[0], shift[1]]
        SS.weights.fixFlatNodes()

        if render:
            SS.renderWeights("S_%d"%bid)


    def shuffleRefineContigs(self, timer, inclusivity=2):
        """refine bins by shuffling contigs around"""
        print "    Start shuffle refinement"

        # first, build a SOM
        bids = self.BM.getBids()
        bin_c_lengths = {}      # bid => [len,len,...]
        for bid in bids:
            bin_c_lengths[bid] = [self.PM.contigLengths[row_index] for row_index in self.BM.bins[bid].rowIndices]

        (SS, minz, maxz, side) = self.buildSOM(timer,
                                               maskBoundaries=True,
                                               defineBins=True,
                                               retrain=True)

        print "    %s" % timer.getTimeStamp()

        # now do the shuffle refinement, keep an eye out for
        new_assignments = {}
        wrongs = {}
        news = {}
        rights = {}
        nones = {}

        # we load all contigs into the block
        block = np_zeros((len(self.PM.transformedCP),SOMDIM))
        block[:,:-1] = self.PM.transformedCP
        block[:,-1] = self.PM.kmerNormPC1

        # apply sane normalisation
        block -= minz
        block /= maxz

        for i in range(len(self.PM.indices)):
            assigned = False
            old_bid = self.PM.binIds[i]
            putative_bid = SS.classifyContig(block[i])
            if putative_bid == old_bid:
                # assigned to the old bin
                # nothing much to do here...
                assigned = True
                try:
                    new_assignments[old_bid].append(i)
                except KeyError:
                    new_assignments[old_bid] = [i]
                try:
                    rights[old_bid] += 1
                except KeyError:
                    rights[old_bid] = 1

            elif self.BM.bins[putative_bid].binSize > 1:
                # stats f**k up on single contig bins, soz...
                length_wrong = self.GT.isMaxOutlier(self.PM.contigLengths[i],
                                                    bin_c_lengths[putative_bid]
                                                    )
                if not length_wrong:
                    # fits length cutoff
                    (covZ,merZ) = self.BM.scoreContig(i, putative_bid)
                    if covZ <= inclusivity and merZ <= inclusivity:
                        # we can recruit
                        try:
                            new_assignments[putative_bid].append(i)
                        except KeyError:
                            new_assignments[putative_bid] = [i]
                        assigned = True

                        #----------------------
                        if old_bid != 0:
                            if putative_bid != old_bid:
                                try:
                                    wrongs[old_bid] += 1
                                except KeyError:
                                    wrongs[old_bid] = 1
                            else:
                                try:
                                    rights[putative_bid] += 1
                                except KeyError:
                                    rights[putative_bid] = 1
                        else:
                            try:
                                news[putative_bid] += 1
                            except KeyError:
                                news[putative_bid] = 1
                        #----------------------
            if not assigned:
                # could not put it anwhere
                # assign to the old bin
                try:
                    new_assignments[old_bid].append(i)
                except KeyError:
                    new_assignments[old_bid] = [i]
                try:
                    nones[old_bid] += 1
                except KeyError:
                    nones[old_bid] = 1


	if False:
            print "    ------------------------------------------------------"
            print "     BID    ORIG    CHGE    SAME    NEWS    NONE    TOTAL"
            print "    ------------------------------------------------------"
            for bid in bids:
                print "   %4d    %5d   " % (bid, self.BM.bins[bid].binSize),
                if bid in wrongs:
                    print "%04d   " % wrongs[bid],
                else:
                    print "0000   ",
                if bid in rights:
                    print "%04d   " % rights[bid],
                else:
                    print "0000   ",
                if bid in news:
                    print "%04d   " % news[bid],
                else:
                    print "0000   ",
                if bid in nones:
                    print "%04d   " % nones[bid],
                else:
                    print "0000   ",
                print "%04d   " % len(new_assignments[bid])
            print "\n    ---------------------------------------------"

        # now get ready for saving.
        # first, we nuke all non-chimeric bins
        self.BM.deleteBins(bids, force=True, freeBinnedRowIndices=False, saveBins=False)
        self.BM.bins = {}

        # these are profile manager variables. We will overwrite
        # these here so that everything stays in sync..
        self.PM.binIds = np_zeros((len(self.PM.indices))) # list of bin IDs
        self.PM.validBinIds = {}              # { bid : numMembers }
        self.PM.binnedRowIndices = {}         # dictionary of those indices which belong to some bin
        self.PM.restrictedRowIndices = {}     # dictionary of those indices which can not be binned yet
        self.PM.isLikelyChimeric = {}

        # now we rebuild all the bins but with the new assignments
        for bid in new_assignments:
            if bid != 0:
                row_indices = np_array(new_assignments[bid])
                self.BM.makeNewBin(rowIndices=row_indices, bid=bid)
                self.PM.validBinIds[bid] = len(row_indices)
                for row_index in row_indices:
                    self.PM.binIds[row_index] = bid
                    self.PM.binnedRowIndices[row_index] = True

        # recheck bins for likely chimeric bins
        self.markLikelyChimericBins()

        return []

    def removeDuds(self, ms=20, mv=1000000, verbose=False):
        """Run this after refining to remove scrappy leftovers"""
        print "    Removing dud cores (min %d contigs or %d bp)" % (ms, mv)
        deleters = []
        for bid in self.BM.getBids():
            self.BM.bins[bid]
            if not self.BM.isGoodBin(bin.totalBP, bin.binSize, ms=ms, mv=mv):
                # delete this chap!
                deleters.append(bid)
        if verbose:
            print "duds", deleters
        if len(deleters) > 0:
            self.BM.deleteBins(deleters,
                               force=True,
                               freeBinnedRowIndices=True,
                               saveBins=False)
        print "    Removed %d cores leaving %d cores" % (len(deleters), len(self.BM.bins))
        return deleters

#------------------------------------------------------------------------------
# UTILITIES

    def rePCA(self,
              bidList,
              mode='mer',
              doBoth=False,
              doContigs=False,
              addZeros=False,
              addOnes=False):
        """Re-calculate PCA coords for a collection of bins

        mode may be 'mer', 'cov' or 'trans'

        If do contigs is set then it returns a n X 2 array
        of PC1, PC2 for all contigs in all bins in bidList
        The ordering is bidList -> rowIndices.
        IE. n = sum(rowIndices) for all bins in bidList

        If doContigs is not set then we average across all the contigs
        in each bin and return a n X 2 array where n is the number of
        bins in bidList

        If doBoth is set then we do both. However, we return a dict of type:
        {row_index : np_array([PC1, PC2])}

        addZeros adds a zero vector as the final entry in the PCA output
        addOnes adds a ones vector as the second last entry.
        addOnes implies addZeros
        """
        signal = []
        both_tmp = {}
        both_ret = {}
        num_ss = 0
        if mode == 'mer':
            data = self.PM.kmerSigs
        elif mode == 'cov':
            data = self.PM.covProfiles
        elif mode == 'trans':
            data = self.PM.transformedCP
        else:
            raise ge.ModeNotAppropriateException("Invlaid mode " + type)

        pc_len = len(data[self.BM.bins[bidList[0]].rowIndices[0]])
        if doBoth:
            # the eventual goal is to produce a dict of RI -> kPCA
            # we just need to shuffle things about for a second here...
            for bid in bidList:
                for row_index in self.BM.bins[bid].rowIndices:
                    signal.append(data[row_index])
                    both_tmp[num_ss] = row_index
                    num_ss += 1
        else:
            for bid in bidList:
                for row_index in self.BM.bins[bid].rowIndices:
                    signal.append(data[row_index])
                    num_ss += 1

        if addOnes:
            addZeros = True
            signal.append(np_ones(pc_len))
            num_ss += 1

        if addZeros:
            signal.append(np_zeros(pc_len))
            num_ss += 1

        signal = np_reshape(signal, (num_ss, pc_len))


        # do the PCA analysis
        Center(signal,verbose=0)
        p = PCA(signal)
        components = p.pc()

        # now make the color profile based on PC1
        PC1 = np_array([float(i) for i in components[:,0]])
        PC2 = np_array([float(i) for i in components[:,1]])

        # normalise to fit between 0 and 1
        PC1 -= np_min(PC1)
        PC1 /= np_max(PC1)
        PC2 -= np_min(PC2)
        PC2 /= np_max(PC2)

        if doContigs:
            return np_reshape([[PC1[i], PC2[i]] for i in range(len(PC1))],
                              (num_ss,2))
        if doBoth:
            # make the actual return dict
            for i in range(num_ss):
                both_ret[both_tmp[i]] = np_array([PC1[i], PC2[i]])

        # else work out the average for each bin
        ml_2d = np_array([])
        index_start = 0
        for bid in bidList:
            nri = len(self.BM.bins[bid].rowIndices)
            mPCA = np_mean(np_reshape([[PC1[i], PC2[i]] for i in range(index_start, index_start+nri)],
                                      (nri,2)),
                           axis=0)

            index_start += nri
            ml_2d = np_append(ml_2d, mPCA)

        if addOnes:
            # second last index is the 1
            # last index is the 0 - PCA stylez!
            ml_2d = np_append(ml_2d, [PC1[-2], PC2[-2]])
            ml_2d = np_append(ml_2d, [PC1[-1], PC2[-1]])
            if doBoth:
                return (np_reshape(ml_2d, (len(bidList)+2,2)), both_ret)
            return np_reshape(ml_2d, (len(bidList)+2,2))
        elif addZeros:
            # last index is the 0 - PCA stylez!
            ml_2d = np_append(ml_2d, [PC1[-1], PC2[-1]])
            if doBoth:
                return (np_reshape(ml_2d, (len(bidList)+1,2)), both_ret)
            return np_reshape(ml_2d, (len(bidList)+1,2))
        else:
            if doBoth:
                return (np_reshape(ml_2d, (len(bidList),2)), both_ret)
            return np_reshape(ml_2d, (len(bidList),2))

    def getKCut(self):
        """Work out the easy cutoff for kmerVal distance"""
        median_k_vals = []
        for bid in self.BM.getNonChimericBinIds():
            if len(self.BM.getBin(bid).rowIndices) > 1:
                kdist = self.kDist(self.BM.getBin(bid).rowIndices)
                if kdist != None:
                    median_k_vals.append(kdist)

        return np_median(median_k_vals), np_std(median_k_vals)

    def kDist(self, row_indices):
        bin_k_vals = self.PM.kmerPCs[row_indices]
        k_dist = pdist(bin_k_vals, 'cityblock')
        if len(k_dist) > 0:
            return np_median(k_dist)

        return None

    def kDistMergedBins(self, bin1, bin2):
        merged_indices = np_concatenate((bin1.rowIndices, bin2.rowIndices))
        return self.kDist(merged_indices)

    def kDistBetweenBins(self, bin1, bin2):
        return np_median(cdist(self.PM.kmerPCs[bin1.rowIndices], self.PM.kmerPCs[bin2.rowIndices], 'cityblock'))

    def getEvenlySpacedPtsZ(self, row_indices, sample_size):
        # select samples evenly along Z-axis of coverage space
        # (this makes the program deterministic while getting a good 'random' spread of points)
        sorted_indices = np_argsort(self.PM.transformedCP[row_indices, -1])
        step_size = float(len(row_indices)) / sample_size
        si = []
        index = 0.0
        for _i in xrange(0, sample_size):
            si.append(row_indices[sorted_indices[int(index)]])
            index += step_size

        return si

    def getCCut(self):
        """Work out the easy cutoff for coverage angle difference"""
        median_angles = []
        for bid in self.BM.getNonChimericBinIds():
            if len(self.BM.getBin(bid).rowIndices) > 1:
                cdistance = self.cDist(self.BM.getBin(bid).rowIndices)
                median_angles.append(cdistance)

        return np_median(median_angles), np_std(median_angles)

    def cDist(self, row_indices):
        max_in_bin = 100

        if len(row_indices) > max_in_bin:
            sample_size = max_in_bin
            si = self.getEvenlySpacedPtsZ(row_indices, max_in_bin)
        else:
            sample_size = len(row_indices)
            si = row_indices

        # select a few at random
        angles = []
        for i in range(sample_size):
            for j in range(i+1, sample_size):
                r1 = si[i]
                r2 = si[j]
                try:
                    ang = np_arccos(np_dot(self.PM.covProfiles[r1],self.PM.covProfiles[r2]) /
                                             (self.PM.normCoverages[r1]*self.PM.normCoverages[r2]))
                    angles.append(ang)
                except FloatingPointError:
                    pass

        return np_median(angles)

    def cDistMergedBins(self, bin1, bin2):
        merged_indices = np_concatenate((bin1.rowIndices, bin2.rowIndices))
        return self.cDist(merged_indices)

    def cDistBetweenBins(self, bin1, bin2):
        max_in_bin = 100

        if len(bin1.rowIndices) > max_in_bin:
            indices1 = self.getEvenlySpacedPtsZ(bin1.rowIndices, max_in_bin)
        else:
            indices1 = bin1.rowIndices

        if len(bin2.rowIndices) > max_in_bin:
            indices2 = self.getEvenlySpacedPtsZ(bin2.rowIndices, max_in_bin)
        else:
            indices2 = bin2.rowIndices

        angles = []
        for i in xrange(0, min(len(bin1.rowIndices), max_in_bin)):
            r1 = indices1[i]

            for j in xrange(0, min(len(bin2.rowIndices), max_in_bin)):
                r2 = indices2[j]
                try:
                    ang = np_arccos(np_dot(self.PM.covProfiles[r1], self.PM.covProfiles[r2]) /
                                             (self.PM.normCoverages[r1]*self.PM.normCoverages[r2]))
                    angles.append(ang)
                except FloatingPointError:
                    pass

        return np_median(angles)

#-----------------------------
# MERGE TESTING BASED ON KMERS

    def testMergeMer(self,
                     bid1,
                     bid2,
                     merPCAs,
                     maxSample = 200,
                     confidence=0.97,
                     verbose=False):
        """Determine if a merge makes sense in mer land"""
        null_loops = 99
        front_RIs = self.BM.bins[bid1].rowIndices
        rear_RIs = self.BM.bins[bid2].rowIndices

        front_RIs_size = len(front_RIs)
        rear_RIs_size = len(rear_RIs)

        # generate a NULL distribution of alpha splits
        x = np_concatenate([front_RIs, rear_RIs])
        shuffle(x)
        x_size = len(x)

        # sub sample the total space
        if maxSample != 0:
            front_sample_size = np_min([maxSample, front_RIs_size])
            rear_sample_size = np_min([maxSample, rear_RIs_size])
        else:
            front_sample_size = front_RIs_size
            rear_sample_size = rear_RIs_size

        # we do different things here depending on the size of the samples
        if front_sample_size < 64:
            # We should brute force this guy
            F_funct = self.calculateMerAlphaPart
        else:
            # take a sampling approach
            F_funct = self.calculateMerAlphaPartSampleNoSame

        if rear_sample_size < 64:
            R_funct = self.calculateMerAlphaPart
        else:
            R_funct = self.calculateMerAlphaPartSampleNoSame

        if front_sample_size * rear_sample_size < 2000:
            C_funct = self.calculateMerAlphaPart
        else:
            C_funct = self.calculateMerAlphaPartSample

        merAlphas={}
        # test the given split against it
        if maxSample != 0:
            shuffle(front_RIs)
            shuffle(rear_RIs)
            FRI = front_RIs[:front_sample_size]
            RRI = rear_RIs[:rear_sample_size]
        else:
            FRI = front_RIs
            RRI = rear_RIs
        test_T_score = C_funct(FRI,RRI,merPCAs,merAlphas) / ( F_funct(FRI,FRI,merPCAs,merAlphas) + R_funct(RRI,RRI,merPCAs,merAlphas))

        T_scores = []
        index_array = np_arange(x_size)
        fsi = np_arange(front_sample_size)
        rsi = np_arange(front_sample_size, front_sample_size+rear_sample_size)
        for i in range(null_loops):
            # select two sets of bins at random
            shuffle(index_array)
            FRI = x[index_array[fsi]]
            RRI = x[index_array[rsi]]
            T_scores.append(C_funct(FRI,RRI,merPCAs,merAlphas) / ( F_funct(FRI,FRI,merPCAs,merAlphas) + R_funct(RRI,RRI,merPCAs,merAlphas)))

        T_scores = sorted(T_scores)
        index = int(np_around(float(null_loops+1)*confidence))

        return (test_T_score, T_scores[index])

    def calculateMerAlphaPartSampleNoSame(self, RI1, RI2, profile, alphas, sampleSize=2000):
        """Calculate one part of an alpha score

        Subsample of all vs all, assumes RI1 == RI2"""
        lr1 = len(RI1)
        lr2 = len(RI2)
        score = 0.0
        samples_selected = {}
        num_samples_complete = 0
        while num_samples_complete < sampleSize:
            combo = (randint(lr1), randint(lr2))
            while (combo in samples_selected) or (combo[0] == combo[1]):
                combo = (randint(lr1), randint(lr2))
            samples_selected[combo] = True
            r1 = RI1[combo[0]]
            r2 = RI2[combo[1]]
            num_samples_complete += 1
            try:
                # try to pull the value out
                score += alphas[(r1,r2)]
            except KeyError:
                try:
                    # try to pull the value out
                    score += alphas[(r1,r2)]
                except KeyError:
                    # calculate it if it's not there...
                    # inline this because we call it alot!
                    dist = np_sum(np_abs(profile[r1] - profile[r2]))
                    alphas[(r1,r2)] = dist
                    score += dist
        return score / sampleSize

    def calculateMerAlphaPartSample(self, RI1, RI2, profile, alphas, sampleSize=2000):
        """Calculate one part of an alpha score

        Subsample of all vs all, assumes RI1 != RI2
        """
        lr1 = len(RI1)
        lr2 = len(RI2)
        score = 0.0
        samples_selected = {}
        num_samples_complete = 0
        while num_samples_complete < sampleSize:
            combo = (randint(lr1), randint(lr2))
            while combo in samples_selected:
                combo = (randint(lr1), randint(lr2))
            samples_selected[combo] = True
            r1 = RI1[combo[0]]
            r2 = RI2[combo[1]]
            num_samples_complete += 1
            try:
                # try to pull the value out
                score += alphas[(r1,r2)]
            except KeyError:
                try:
                    # try to pull the value out
                    score += alphas[(r1,r2)]
                except KeyError:
                    # calculate it if it's not there...
                    # inline this because we call it alot!
                    dist = np_sum(np_abs(profile[r1] - profile[r2]))
                    alphas[(r1,r2)] = dist
                    score += dist
        return score / sampleSize

    def calculateMerAlphaPart(self, RI1, RI2, profile, alphas):
        """Calculate one part of an alpha score

        All vs all complete
        """
        lr1 = len(RI1)
        lr2 = len(RI2)
        score = 0.0
        for r1 in RI1:
            for r2 in RI2:
                try:
                    score += alphas[(r1,r2)]
                except KeyError:
                    try:
                        score += alphas[(r2,r1)]
                    except KeyError:
                        dist = np_sum(np_abs(profile[r1] - profile[r2]))
                        alphas[(r1,r2)] = dist
                        score += dist
        return score / (lr1*lr2)


    def calculateMerAlphaTScore(self, RI1, RI2, profile, alphas):
        """Measure the goodness of the separation into lists

        specifically, calculate: t = INTER_DIST / (INTRA_LIST1_DISTANCE + INTRA_LIST2_DISTANCE)
        """
        lr1 = len(RI1)
        lr2 = len(RI2)
        R1_intras = 0.0
        R2_intras = 0.0
        inters = 0.0
        for i in range(lr1):
            for j in range(i+1, lr1):
                try:
                    # try to pull the value out
                    R1_intras += alphas[(RI1[i], RI1[j])]
                except KeyError:
                    try:
                        # try to pull the value out
                        R1_intras += alphas[(RI1[j], RI1[i])]
                    except KeyError:
                        # calculate it if it's not there...
                        # inline this because we call it alot!
                        dist = np_sum(np_abs(profile[RI1[i]] - profile[RI1[j]]))
                        alphas[(RI1[i],RI1[j])] = dist
                        R1_intras += dist
        R1_intras /= ((lr1 - 1)*lr1 / 2)

        for i in range(lr2):
            for j in range(i+1, lr2):
                try:
                    R2_intras += alphas[(RI2[i], RI2[j])]
                except KeyError:
                    try:
                        R2_intras += alphas[(RI2[j], RI2[i])]
                    except KeyError:
                        dist = np_sum(np_abs(profile[RI2[i]] - profile[RI2[j]]))
                        alphas[(RI2[i],RI2[j])] = dist
                        R2_intras += dist
        R2_intras /= ((lr2 - 1)*lr2 / 2)

        for r1 in RI1:
            for r2 in RI2:
                try:
                    inters += alphas[(r1,r2)]
                except KeyError:
                    try:
                        inters += alphas[(r2,r1)]
                    except KeyError:
                        dist = np_sum(np_abs(profile[r1] - profile[r2]))
                        alphas[(r1,r2)] = dist
                        inters += dist
        inters /= (lr1*lr2)
        return inters/(R1_intras + R2_intras)


#-----------------------------
# MERGE TESTING BASED ON COVERAGE

    def testMergeCoverage(self,
                          bid1,
                          bid2,
                          maxSample = 200,
                          confidence=0.97,
                          verbose=False):
        """Determine if a merge based on kmer PCAs makes sense in coverage land

        Calculates a t-score which measures the coverage separation of the two groups
        Higher t-scores indicate greater separation between the groups.
        """
        null_loops = 99
        front_RIs = self.BM.bins[bid1].rowIndices
        rear_RIs = self.BM.bins[bid2].rowIndices

        front_RIs_size = len(front_RIs)
        rear_RIs_size = len(rear_RIs)

        # generate a NULL distribution of alpha splits
        x = np_concatenate([front_RIs, rear_RIs])
        shuffle(x)
        x_size = len(x)

        # sub sample the total space
        if maxSample != 0:
            front_sample_size = np_min([maxSample, front_RIs_size])
            rear_sample_size = np_min([maxSample, rear_RIs_size])
        else:
            front_sample_size = front_RIs_size
            rear_sample_size = rear_RIs_size

        # we do different things here depending on the size of the samples
        if front_sample_size < 64:
            # We should brute force this guy
            F_funct = self.calculateCovAlphaPart
        else:
            # take a sampling approach
            F_funct = self.calculateCovAlphaPartSampleNoSame

        if rear_sample_size < 64:
            R_funct = self.calculateCovAlphaPart
        else:
            R_funct = self.calculateCovAlphaPartSampleNoSame

        if front_sample_size * rear_sample_size < 2000:
            C_funct = self.calculateCovAlphaPart
        else:
            C_funct = self.calculateCovAlphaPartSample
        covAlphas={}

        # test the given split against it
        if maxSample != 0:
            shuffle(front_RIs)
            shuffle(rear_RIs)
            FRI = front_RIs[:front_sample_size]
            RRI = rear_RIs[:rear_sample_size]
        else:
            FRI = front_RIs
            RRI = rear_RIs
        test_T_score = C_funct(FRI,RRI,covAlphas) / ( F_funct(FRI,FRI,covAlphas) + R_funct(RRI,RRI,covAlphas))

        T_scores = []
        index_array = np_arange(x_size)
        fsi = np_arange(front_sample_size)
        rsi = np_arange(front_sample_size, front_sample_size+rear_sample_size)
        for i in range(null_loops):
            # select two sets of bins at random
            shuffle(index_array)
            FRI = x[index_array[fsi]]
            RRI = x[index_array[rsi]]
            T_scores.append(C_funct(FRI,RRI,covAlphas) / ( F_funct(FRI,FRI,covAlphas) + R_funct(RRI,RRI,covAlphas)))

        T_scores = sorted(T_scores)
        index = int(np_around(float(null_loops+1)*confidence))

        return (test_T_score, T_scores[index])

    def calculateCovAlphaPartSampleNoSame(self, RI1, RI2, alphas, sampleSize=2000):
        """Calculate one part of an alpha score

        Subsample of all vs all, assumes RI1 == RI2
        """
        lr1 = len(RI1)
        lr2 = len(RI2)
        score = 0.0
        samples_selected = {}
        num_samples_complete = 0
        while num_samples_complete < sampleSize:
            combo = (randint(lr1), randint(lr2))
            while (combo in samples_selected) or (combo[0] == combo[1]):
                combo = (randint(lr1), randint(lr2))
            samples_selected[combo] = True
            r1 = RI1[combo[0]]
            r2 = RI2[combo[1]]
            num_samples_complete += 1
            try:
                # try to pull the value out
                score += alphas[(r1,r2)]
            except KeyError:
                try:
                    # try to pull the value out
                    score += alphas[(r1,r2)]
                except KeyError:
                    # calculate it if it's not there...
                    # inline this becuase we call it alot!
                    try:
                        ang = np_arccos(np_dot(self.PM.covProfiles[r1],self.PM.covProfiles[r2]) /
                                        (self.PM.normCoverages[r1]*self.PM.normCoverages[r2]))
                    except FloatingPointError:
                        ang = 0.0
                    alphas[(r1,r2)] = ang
                    score += ang
        return score / sampleSize

    def calculateCovAlphaPartSample(self, RI1, RI2, alphas, sampleSize=2000):
        """Calculate one part of an alpha score

        Subsample of all vs all, assumes RI1 != RI2
        """
        lr1 = len(RI1)
        lr2 = len(RI2)
        score = 0.0
        samples_selected = {}
        num_samples_complete = 0
        while num_samples_complete < sampleSize:
            combo = (randint(lr1), randint(lr2))
            while combo in samples_selected:
                combo = (randint(lr1), randint(lr2))
            samples_selected[combo] = True
            r1 = RI1[combo[0]]
            r2 = RI2[combo[1]]
            num_samples_complete += 1
            try:
                # try to pull the value out
                score += alphas[(r1,r2)]
            except KeyError:
                try:
                    # try to pull the value out
                    score += alphas[(r1,r2)]
                except KeyError:
                    # calculate it if it's not there...
                    # inline this becuase we call it alot!
                    try:
                        ang = np_arccos(np_dot(self.PM.covProfiles[r1],self.PM.covProfiles[r2]) /
                                        (self.PM.normCoverages[r1]*self.PM.normCoverages[r2]))
                    except FloatingPointError:
                        ang = 0.0
                    alphas[(r1,r2)] = ang
                    score += ang
        return score / sampleSize

    def calculateCovAlphaPart(self, RI1, RI2, alphas):
        """Calculate one part of an alpha score

        All vs all
        """
        lr1 = len(RI1)
        lr2 = len(RI2)
        score = 0.0
        for r1 in RI1:
            for r2 in RI2:
                try:
                    score += alphas[(r1,r2)]
                except KeyError:
                    try:
                        score += alphas[(r2,r1)]
                    except KeyError:
                        try:
                            ang = np_arccos(np_dot(self.PM.covProfiles[r1],self.PM.covProfiles[r2]) /
                                            (self.PM.normCoverages[r1]*self.PM.normCoverages[r2]))
                        except FloatingPointError:
                            ang = 0.0
                        alphas[(r1,r2)] = ang
                        score += ang
        return score / (lr1*lr2)


    def calculateCovAlphaTScore(self, RI1, RI2, alphas):
        """Measure the goodness of the separation into lists

        specifically, calculate: t = INTER_DIST / (INTRA_LIST1_DISTANCE + INTRA_LIST2_DISTANCE)
        """
        lr1 = len(RI1)
        lr2 = len(RI2)
        R1_intras = 0.0
        R2_intras = 0.0
        inters = 0.0
        for i in range(lr1):
            for j in range(i+1, lr1):
                try:
                    # try to pull the value out
                    R1_intras += alphas[(RI1[i], RI1[j])]
                except KeyError:
                    try:
                        # try to pull the value out
                        R1_intras += alphas[(RI1[j], RI1[i])]
                    except KeyError:
                        # calculate it if it's not there...
                        # inline this becuase we call it alot!
                        try:
                            ang = np_arccos(np_dot(self.PM.covProfiles[RI1[i]],self.PM.covProfiles[RI1[j]]) /
                                            (self.PM.normCoverages[RI1[i]]*self.PM.normCoverages[RI1[j]]))
                        except FloatingPointError:
                            ang = 0.0
                        alphas[(RI1[i],RI1[j])] = ang
                        R1_intras += ang
        R1_intras /= ((lr1 - 1)*lr1 / 2)

        for i in range(lr2):
            for j in range(i+1, lr2):
                try:
                    R2_intras += alphas[(RI2[i], RI2[j])]
                except KeyError:
                    try:
                        R2_intras += alphas[(RI2[j], RI2[i])]
                    except KeyError:
                        try:
                            ang = np_arccos(np_dot(self.PM.covProfiles[RI2[i]],self.PM.covProfiles[RI2[j]]) /
                                            (self.PM.normCoverages[RI2[i]]*self.PM.normCoverages[RI2[j]]))
                        except FloatingPointError:
                            ang = 0.0
                        alphas[(RI2[i],RI2[j])] = ang
                        R2_intras += ang
        R2_intras /= ((lr2 - 1)*lr2 / 2)
        for r1 in RI1:
            for r2 in RI2:
                try:
                    inters += alphas[(r1,r2)]
                except KeyError:
                    try:
                        inters += alphas[(r2,r1)]
                    except KeyError:
                        try:
                            ang = np_arccos(np_dot(self.PM.covProfiles[r1],self.PM.covProfiles[r2]) /
                                            (self.PM.normCoverages[r1]*self.PM.normCoverages[r2]))
                        except FloatingPointError:
                            ang = 0.0
                        alphas[(r1,r2)] = ang
                        inters += ang
        inters /= (lr1*lr2)
        return inters/(R1_intras + R2_intras)

#------------------------------------------------------------------------------
# RECRUITMENT

    def recruitWrapper(self, timer, inclusivity=2, step=200, nukeAll=False, saveBins=False):
        """Recuit more contigs to the bins"""
        print "Recruiting unbinned contigs"

        # make a list of all the cov and kmer vals
        total_expanded = 0
        total_binned = 0
        total_unbinned = 0
        bin_c_lengths = {}
        total_contigs = len(self.PM.indices)
        shortest_binned = 1000000000          # we need to know this
        shortest_unbinned = 1000000000

        # for stats, work out number binned and unbinned and relative lengths
        unbinned = {}
        for row_index in range(len(self.PM.indices)):
            if(row_index in self.PM.binnedRowIndices):
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

        # build the classifier on all the existing bins
        (SS, minz, maxz, side) = self.buildSOM(timer,
                                               maskBoundaries=True,
                                               defineBins=True,
                                               retrain=True)

        print "    %s" % timer.getTimeStamp()

        # go through the steps we decided on
        affected_bids = list(np_copy(self.BM.getBids()))
        for cutoff in steps:
            # work out the bin length, mer, etc stats
            for bid in affected_bids:
                bin_c_lengths[bid] = [self.PM.contigLengths[row_index] for row_index in self.BM.bins[bid].rowIndices]
                self.BM.bins[bid].makeBinDist(self.PM.transformedCP,
                                           self.PM.averageCoverages,
                                           self.PM.kmerNormPC1,
                                           self.PM.kmerPCs,
                                           self.PM.contigGCs,
                                           self.PM.contigLengths)
            affected_bids = []
            this_step_binned = 0
            new_binned = []

            # load the unbinned guys into a block
            unbinned_rows = []
            unbinned_lens = []
            for row_index in unbinned:
                if unbinned[row_index] >= cutoff:
                    unbinned_rows.append(row_index)
                    unbinned_lens.append(unbinned[row_index])
            block = np_zeros((len(unbinned_rows),SOMDIM))
            block[:,:-1] = self.PM.transformedCP[unbinned_rows]
            block[:,-1] = self.PM.kmerNormPC1[unbinned_rows]
            # apply sane normalisation
            block -= minz
            block /= maxz

            print "    Recruiting contigs above: %d (%d contigs)" % (cutoff, len(unbinned_rows))

            for i in range(len(unbinned_rows)):
                putative_bid = SS.classifyContig(block[i])
                if self.BM.bins[putative_bid].binSize > 1:
                    # stats f**k up on single contig bins, soz...
                    length_wrong = self.GT.isMaxOutlier(unbinned_lens[i],
                                                        bin_c_lengths[putative_bid]
                                                        )
                    if not length_wrong:
                        # fits length cutoff
                        (covZ,merZ) = self.BM.scoreContig(unbinned_rows[i], putative_bid)
                        if covZ <= inclusivity and merZ <= inclusivity:
                            # we can recruit
                            self.BM.bins[putative_bid].rowIndices = np_append(self.BM.bins[putative_bid].rowIndices,
                                                                              unbinned_rows[i]
                                                                              )
                            affected_bids.append(putative_bid)
                            this_step_binned += 1
                            total_binned += 1
                            total_expanded += 1
                            new_binned.append(unbinned_rows[i])

            # need only check this guy once
            for row_index in new_binned:
                del unbinned[row_index]

            print "    Recruited: %d contigs" % this_step_binned
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()

        # talk to the user
        perc_recruited = float(total_expanded)/float(total_unbinned)
        perc_binned = float(total_binned)/float(total_contigs)
        print "    Recruited %0.4f" % perc_recruited +"%"+" of %d unbinned contigs" % total_unbinned
        print "    END: %0.4f" % perc_binned +"%"+" of %d requested contigs in bins" % total_contigs
        print "    %s" % timer.getTimeStamp()
        sys_stdout.flush()

        # now save
        if(saveBins):
            print "Saving bins"
            self.BM.saveBins()

#------------------------------------------------------------------------------
# UI and IMAGE RENDERING

    def writeGV(self, graph):
        """Output a valid graphviz dot file"""
        op = "digraph refine {\n"
        # render nodes
        for bid in graph[0].keys():
            op += graph[0][bid]
        # render edges
        op += "\n".join(graph[1])
        op += "\n};\n"
        return op

    def printRefinePlotterInstructions(self):
        raw_input( "****************************************************************\n"
                   " REFINING INSTRUCTIONS - PLEASE READ CAREFULLY\n"+
                   "****************************************************************\n"
                   " You have chosen to refine in plotter mode. Congratulations!\n"
                   " You will be shown a 3d plot of all the bins, colored by kmer\n"
                   " profile. Bin Ids in close proximity and similar color may need\n"
                   " to be merged. Conversely, you can split bins which appear chimeric\n"
                   " Follow the instructions to merge or split these bins\n\n"
                   " Good Luck!\n\n"
                   " Press return to continue...")
        print "****************************************************************"

    def promptOnPlotterRefine(self, minimal=False):
        """Find out what the user wishes to do next when refining bins"""
        input_not_ok = True
        valid_responses = ['R','P','G','U','B','V','M','S', 'C','E','X','Q']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" What next? ("+vrs+") : ")
            else:
                option = raw_input("\n Please choose from the following options:\n" \
                                   "------------------------------------------------------------\n" \
                                   " r = plot entire space using bin ids\n" \
                                   " p = plot entire space with bins as points\n" \
                                   " g = plot entire space for bins within a specific GC range\n" \
                                   " u = plot all contigs in untransformed coverage space (first 3 stoits only)\n"
                                   " b = plot one or more bins\n" \
                                   " v = plot all contigs in vincinity of bin\n" \
                                   " m = merge two or more bins\n" \
                                   " s = split a bin into multiple pieces\n" \
                                   " c = change colormap\n" \
                                   " e = toggle elipses (default = on)\n" \
                                   " x = toggle chimeric bins (default = hidden)\n" \
                                   " q = quit\n" \
                                   "------------------------------------------------------------\n" \
                                   " What next? ("+vrs+") : ")
            if(option.upper() in valid_responses):
                return option.upper()
            else:
                print "Error, unrecognised choice '"+option+"'"
                minimal=True

    def PCA2Col(self, PCAs):
        """Convert a set of PCA coords into a color"""
        # use HSV to RGB to generate colors
        S = 1       # SAT and VAL remain fixed at 1. Reduce to make
        V = 1       # Pastels if that's your preference...
        return np_reshape(np_array([htr(val, S, V) for val in PCAs[:,0]]),
                          (len(PCAs),3))

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

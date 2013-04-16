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
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.2.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Alpha"

###############################################################################

from sys import exc_info, exit, stdout as sys_stdout
from Queue import Queue
from operator import itemgetter
import readline

from colorsys import hsv_to_rgb as htr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

from numpy import arange as numpy_arange, copy as np_copy, arange as np_arange, ravel as np_ravel, ones as np_ones, eye as np_eye, shape as np_shape, around as np_around, argmax as np_argmax, arccos as np_cos, dot as np_dot, sum as np_sum, abs as np_abs, amax as np_amax, amin as np_amin, append as np_append, arccos as np_arccos, argmin as np_argmin, argsort as np_argsort, array as np_array, ceil as np_ceil, concatenate as np_concatenate, delete as np_delete, log10 as np_log10, max as np_max, mean as np_mean, median as np_median, min as np_min, pi as np_pi, reshape as np_reshape, seterr as np_seterr, size as np_size, sort as np_sort, sqrt as np_sqrt, std as np_std, where as np_where, zeros as np_zeros, cos as np_cos, sin as np_sin
from numpy.linalg import norm as np_norm 
from numpy.random import shuffle as shuffle, randint as randint

from scipy.spatial import KDTree as kdt
from scipy.cluster.vq import kmeans,vq,whiten,kmeans2
from scipy.spatial.distance import cdist, pdist, squareform

# GroopM imports
from binManager import BinManager
from bin import Bin
from ellipsoid import EllipsoidTool
from PCA import PCA, Center
import groopmTimekeeper as gtime
import groopmExceptions as ge

np_seterr(all='raise')     

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class RefineEngine:
    """Workhorse wrapper for bin refinement"""
    def __init__(self,
                 timer,
                 BM=None,
                 dbFileName="",
                 transform=True,
                 getUnbinned=False,
                 loadContigNames=False,
                 bids=[]
                 ):
        # worker classes
        if BM is None:
            # make our own ones from scratch
            self.BM = BinManager(dbFileName=dbFileName)
            self.BM.loadBins(timer,
                             bids=bids,
                             makeBins=True,
                             silent=False,
                             loadContigNames=loadContigNames,
                             getUnbinned=getUnbinned,
                             transform=transform,
                             loadRawKmers=True)
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
                                               self.PM.kmerVals, 
                                               self.PM.contigLengths)
                self.BM.plotBins(FNPrefix=plotFinal, ET=self.ET)
            
            if saveBins:
                self.BM.saveBins(nuke=True)
        else:
            self.plotterRefineBins(ignoreRanges=ignoreRanges)

    def plotterRefineBins(self, ignoreRanges=False):
        """combine similar bins using 3d plots"""
        use_elipses = True
        ET = self.ET
        self.printRefinePlotterInstructions()
        #self.BM.plotBinIds(ignoreRanges=ignoreRanges)
        continue_merge = True
        while(continue_merge):
            user_option = self.promptOnPlotterRefine()
            
            if(user_option == 'Q'):
                print 'Hasta luego mi amigo'
                return

            elif(user_option == 'E'):
                if use_elipses:
                    ET = None
                    use_elipses = False
                    print "Ellipses off"
                else:
                    ET = self.ET
                    use_elipses = True
                    print "Ellipses on"
            
            elif(user_option == 'R'):
                self.BM.plotBinIds(ignoreRanges=ignoreRanges)
            
            elif(user_option == 'P'):
                self.BM.plotBinPoints(ignoreRanges=ignoreRanges)
            
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
                self.BM.plotBinIds(krange=krange, ignoreRanges=ignoreRanges)
                
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
                            bids = self.BM.getBids()
                        else:
                            for bid in bids:
                                if bid not in self.BM.bins:
                                    print "ERROR: Bin %d not found!" % bid
                                    have_bid &= False
                    except ValueError:
                        print "You need to enter an integer value!"

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
                               self.BM.bins[bid].covMeans) <= radius:
                        num_points += 1
                        disp_vals = np_append(disp_vals, self.PM.transformedCP[row_index])
                        disp_lens = np_append(disp_lens, np_sqrt(self.PM.contigLengths[row_index]))
                        disp_cols = np_append(disp_cols, self.PM.contigColors[row_index])
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
                self.BM.bins[bid].plotOnAx(ax,
                                           self.PM.transformedCP,
                                           self.PM.contigColors,
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

    def autoRefineBins(self,
                       timer,
                       mergeSimilarBins=True,
                       removeDuds=True,
                       nukeOutliers=True,
                       shuffleRefine=False,
                       verbose=False,
                       plotAfterOB=False,
                       makeGraph=False):
        """Automagically refine bins"""

        sys_stdout.flush()
        kCut = self.getKCut()
        cCut = self.getCCut()
        graph = None
        if makeGraph:
            # lets make a graph
            graph = [{}, []]
            
            for bid in self.BM.getBids():
                bin = self.BM.bins[bid]
                centroid_color = np_mean([self.PM.contigColors[row_index] for row_index in bin.rowIndices],
                              axis=0)
                ncc = [int(i) for i in centroid_color * 255]
                hex_color = '#%02x%02x%02x' % (ncc[0], ncc[1], ncc[2])
                graph[0][bid] = '\t%d [fontcolor="%s" color="%s"];\n' % (bid, hex_color, hex_color)
        
        # identify and remove outlier bins
        if nukeOutliers:
            nuked = self.nukeOutliers()
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
            self.mergeSimilarBins(kCut, cCut, graph=graph)
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()

        if plotAfterOB:
            bids = self.BM.getBids()
            for bid in bids:
                self.BM.bins[bid].makeBinDist(self.PM.transformedCP, 
                                              self.PM.averageCoverages, 
                                              self.PM.kmerVals, 
                                              self.PM.contigLengths)
            self.BM.plotBins(FNPrefix="AFTER_OB", ET=self.ET)
            print "    %s" % timer.getTimeStamp()
            sys_stdout.flush()

        if shuffleRefine:
            nuked = self.shuffleRefineConitgs()
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

    def nukeOutliers(self, verbose=False):
        """Identify and remove small bins which contain mixed genomes
        
        Uses Grubbs testing to identify outliers
        """
        print "    Identifying possible chimeric cores"
        sys_stdout.flush()
        
        bids = self.BM.getBids()
        # first we need to build a distribution!
        kval_stdev_distrb = [] 
        for bid in bids:
            self.BM.bins[bid].makeBinDist(self.PM.transformedCP, 
                                       self.PM.averageCoverages, 
                                       self.PM.kmerVals, 
                                       self.PM.contigLengths)
            kval_stdev_distrb.append(self.BM.bins[bid].kValStdev)
        
        # now we work out the distribution of stdevs
        stdstd = np_std(kval_stdev_distrb)
        stdmean = np_mean(kval_stdev_distrb)
        dead_bins = []
        for bid in bids:
            Z = (self.BM.bins[bid].kValStdev - stdmean)/stdstd
            if Z > 2 and self.BM.bins[bid].totalBP < 100000:
                 dead_bins.append(bid)
                 
        # delete the bad bins
        self.BM.deleteBins(dead_bins,
                           force=True,
                           freeBinnedRowIndices=True,
                           saveBins=False)
        print "    Identified %d possible chimeras leaving %d cores" % (len(dead_bins), len(self.BM.bins))
        return dead_bins

    def mergeSimilarBins(self, kCut, cCut, bids=[], verbose=False, graph=None):
        """Merge bins which are just crying out to be merged!"""
        print "    Merging similar bins with kCut %0.4f cCut %0.4f" % (kCut,cCut)
        if bids == []:
            bids = self.BM.getBids()
        # identify merging groups
        mergers = self.findMergeGroups(verbose=verbose)
        num_bins_removed = 0
        # and then merge them
        for merge in mergers:
            bins_removed = self.combineMergers(merge, kCut, cCut, graph=graph)
            num_bins_removed += len(bins_removed)
        print "    Merged %d cores leaving %d cores" % (num_bins_removed, len(self.BM.bins))        
        return num_bins_removed                
        
    def findMergeGroups(self, bids=[], verbose=False):
        """Identify groups of contigs which could be merged"""
        tdm = []                # these are used in the neighbor search
        bid_2_tdm_index = {} 
        tdm_index_2_bid = {} 

        if bids == []:
            bids = self.BM.getBids()
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
            bin = self.BM.bins[bid]
            bin.makeBinDist(self.PM.transformedCP,
                            self.PM.averageCoverages, 
                            self.PM.kmerVals, 
                            self.PM.contigLengths,
                            merTol=2.5)  # make the merTol a little larger...
            
            # use a 4 dimensional vector [cov, cov, cov, mer]
            tdm.append(np_append(bin.covMeans, [1000*bin.kValMean]))

            # work out the volume of the minimum bounding coverage ellipsoid and kmer ellipse
            (bin_c_ellipsoids[bid], bin_c_ellipsoid_volumes[bid]) = bin.getBoundingCEllipsoidVol(self.PM.transformedCP, ET=self.ET, retA=True)
            BP = np_array(zip([self.PM.kmerVals[i] for i in bin.rowIndices],
                              [self.PM.kmerVals2[i] for i in bin.rowIndices])
                          )
            (bin_k_ellipses[bid], bin_k_ellipse_areas[bid]) = bin.getBoundingKEllipseArea(BP,
                                                                                          ET=self.ET,
                                                                                          retA=True)
            bin_c_lengths[bid] = [self.PM.contigLengths[row_index] for row_index in bin.rowIndices]
            
            bid_2_tdm_index[bid] = index
            tdm_index_2_bid[index] = bid
            index += 1

            # we will not process any pair twice.
            # we also wish to avoid checking if a bin will merge with itself        
            processed_pairs[self.BM.makeBidKey(bid, bid)] = True
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
            base_bin = self.BM.bins[base_bid]
            
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
                if(seen_key in processed_pairs or
                   merged_base_bid == merged_query_bid):
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
# KMER ELLIPSE OVERLAP
                if (bin_k_ellipse_areas[base_bid] <= bin_k_ellipse_areas[query_bid]):
                    INTT = self.ET.doesIntersect2D(bin_k_ellipses[query_bid][0],
                                                   bin_k_ellipses[query_bid][1],
                                                   bin_k_ellipses[base_bid][0],
                                                   bin_k_ellipses[base_bid][1])
                else:
                    INTT = self.ET.doesIntersect2D(bin_k_ellipses[base_bid][0],
                                                   bin_k_ellipses[base_bid][1],
                                                   bin_k_ellipses[query_bid][0],
                                                   bin_k_ellipses[query_bid][1])
                if verbose:
                    
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 1, 1)
                    base_bin.plotMersOnAx(ax,
                                          self.PM.kmerVals,
                                          self.PM.kmerVals2,
                                          self.PM.contigColors,
                                          self.PM.contigLengths,
                                          ET=self.ET)
                    query_bin.plotMersOnAx(ax,
                                           self.PM.kmerVals,
                                           self.PM.kmerVals2,
                                           self.PM.contigColors,
                                           self.PM.contigLengths,
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
                    base_bin.plotOnAx(ax, self.PM.transformedCP, self.PM.contigColors, self.PM.contigLengths, ET=self.ET)                
                    query_bin.plotOnAx(ax, self.PM.transformedCP, self.PM.contigColors, self.PM.contigLengths, ET=self.ET)
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

    def combineMergers(self, bidList, kCut, cCut, graph=None):
        """Try to merge similar bins in the given list"""

        merged_bids = []
        too_big = 10000

        # PCA kmers to find out who is most similar to whom
        (bin_mer_PCAs, mer_con_PCAs) = self.rePCA(bidList, doBoth=True)
        side = len(bidList)
        sq_dists = cdist(bin_mer_PCAs, bin_mer_PCAs)
        dists = squareform(sq_dists)
        
        # raw coverage averages for each bin
        raw_coverage_centroids = {}
        
        while True:
            # find the closest pair 
            closest = np_argmin(dists)
            if dists[closest] == too_big:
                break
            (i,j) = self.small2indices(closest, side-1)
            bid1 = bidList[i]
            bid2 = bidList[j]
            should_merge = False                
            # test if the mer dist is teensy tiny.
            # this is a time saver...
            k_diff = np_abs(self.BM.bins[bid1].kValMean - self.BM.bins[bid2].kValMean)
            #if VVB:
            #    print bid1, bid2, k_diff, 
            mers_OK = False
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
                new_dists = cdist([bin_mer_PCAs[i]], bin_mer_PCAs)
                
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
        
        #self.BM.plotMultipleBins([[h] for h in bidList])
        return merged_bids

    def shuffleRefineConitgs(self):
        """refine bins by shuffling contigs around"""
        print "    Start shuffle refinement"
        start_num_bins = len(self.BM.bins)
        bids = self.BM.getBids()
        
        # assume bin distributions have been made!
        # get an array of bin centroids
        bin_centroids = np_array([])
        for bid in bids:
           bin_centroids = np_append(bin_centroids,
                                     np_mean([self.PM.covProfiles[row_index] for row_index in self.BM.bins[bid].rowIndices], axis=0))
        print np_shape(bin_centroids)
        print (len(self.BM.bins), self.PM.numStoits)
        bin_centroids = np_reshape(bin_centroids, (len(self.BM.bins), self.PM.numStoits))
        
        for bid in bids:
            print "===================="
            print bid
            bin = self.BM.bins[bid]
            for row_index in bin.rowIndices:
                i = 0
                for cent in bin_centroids:
                    #print bids[i], row_index, self.BM.getAngleBetweenVectors(self.PM.covProfiles[row_index], cent)
                    i += 1
                print "---"
        
    
        print "    Removed %d cores leaving %d cores" % (start_num_bins-len(self.BM.bins), len(self.BM.bins))    
        return []    

    def removeDuds(self, ms=20, mv=1000000, verbose=False):
        """Run this after refining to remove scrappy leftovers"""
        print "    Removing dud cores (min %d contigs or %d bp)" % (ms, mv)
        deleters = []
        for bid in self.BM.getBids():
            bin = self.BM.bins[bid]
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

    def small2indices(self, index, side):
        """Return the indices of the comparative items
        when given an index into a condensed distance matrix
        """
        step = 0
        while index >= (side-step):
            index = index - side + step 
            step += 1
        return (step, step + index + 1)
    
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
        mean_k_vals = []
        for bid in self.BM.getBids():
            bin = self.BM.bins[bid]
            bin_k_vals = [[self.PM.kmerVals[row_index]] for row_index in bin.rowIndices]
            k_dist = squareform(cdist(bin_k_vals, bin_k_vals))
            if len(k_dist) > 0:
                mean_k_vals.append(np_mean(k_dist))
        
        return np_mean(mean_k_vals)# + np_std(mean_k_vals)  

    def getCCut(self):
        """Work out the easy cutoff for coverage angle difference"""
        mean_angles = []
        max_in_bin = 100 # at most XXX contigs per bin
        for bid in self.BM.getBids():
            bin = self.BM.bins[bid]
            if bin.binSize < max_in_bin:
                sample_size = bin.binSize
            else:  
                sample_size = max_in_bin
            # select a few at random
            si = np_arange(bin.binSize)
            shuffle(si)
            for i in range(sample_size):
                for j in range(i+1, sample_size):
                    r1 = bin.rowIndices[si[i]]
                    r2 = bin.rowIndices[si[j]]
                    try:
                        ang = np_arccos(np_dot(self.PM.covProfiles[r1],self.PM.covProfiles[r2]) /
                                                 self.PM.normCoverages[r1]/self.PM.normCoverages[r2])
                    except FloatingPointError:
                        ang = 0.0
                    mean_angles.append(ang)
                    
        return np_mean(mean_angles)

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
        index_array = numpy_arange(x_size)
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
        index_array = numpy_arange(x_size)
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
                                        self.PM.normCoverages[r1]/self.PM.normCoverages[r2])
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
                                        self.PM.normCoverages[r1]/self.PM.normCoverages[r2])
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
                                            self.PM.normCoverages[r1]/self.PM.normCoverages[r2])
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
                                            self.PM.normCoverages[RI1[i]]/self.PM.normCoverages[RI1[j]])
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
                                            self.PM.normCoverages[RI2[i]]/self.PM.normCoverages[RI2[j]])
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
                                            self.PM.normCoverages[r1]/self.PM.normCoverages[r2])
                        except FloatingPointError:
                            ang = 0.0
                        alphas[(r1,r2)] = ang
                        inters += ang
        inters /= (lr1*lr2)
        return inters/(R1_intras + R2_intras)

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


#------------------------------------------------------------------------------
# RECRUITMENT

    def recruitWrapper(self, timer, inclusivity=2, step=200, saveBins=False):
        """Recuit more contigs to the bins"""
        print "Recruiting unbinned contigs"
        # make a list of all the cov and kmer vals
        num_bins = len(self.BM.bins)
        num_expanded = 1
        total_expanded = 0
        total_binned = 0
        total_unbinned = 0
        total_contigs = len(self.PM.indices)
        shortest_binned = 1000000000          # we need to know this
        shortest_unbinned = 1000000000
        
        # we need to get a list of bin centroids
        (bin_centroid_points,
         bin_centroid_colors,
         bin_centroid_kvals,
         bids) = self.BM.findCoreCentres(getKVals=True)
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
                        if self.BM.bins[putative_bid].binSize > 1:
                            (covZ,merZ) = self.BM.scoreContig(row_index, putative_bid)
                            if covZ <= inclusivity and merZ <= inclusivity:
                                # we can recruit
                                self.BM.bins[putative_bid].rowIndices = np_append(self.BM.bins[putative_bid].rowIndices,
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
                    self.BM.bins[bid].makeBinDist(self.PM.transformedCP,
                                               self.PM.averageCoverages,
                                               self.PM.kmerVals,
                                               self.PM.contigLengths)      

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
        valid_responses = ['R','P','B','V','M','S','E', 'K','Q']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" What next? ("+vrs+") : ")
            else:
                option = raw_input("\n Please choose from the following options:\n" \
                                   "------------------------------------------------------------\n" \
                                   " r = replot entire space using bin ids\n" \
                                   " p = replot entire space with bins as points\n" \
                                   " k = replot entire space for bins within kmer range\n" \
                                   " b = plot one or more bins\n" \
                                   " v = plot all contigs in vincinity of bin\n" \
                                   " m = merge two or more bins\n" \
                                   " s = split a bin into multiple pieces\n" \
                                   " e = toggle elipses\n" \
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

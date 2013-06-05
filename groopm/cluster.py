#!/usr/bin/env python
###############################################################################
#                                                                             #
#    cluster.py                                                               #
#                                                                             #
#    A collection of classes / methods used when clustering contigs           #
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

import time

from sys import exc_info, exit, stdout

from colorsys import hsv_to_rgb as htr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure
from numpy import (abs as np_abs,
                   allclose as np_allclose,
                   append as np_append,
                   arange as np_arange,
                   argmax as np_argmax,
                   argsort as np_argsort,
                   around as np_around,
                   array as np_array,
                   bincount as np_bincount,
                   column_stack as np_col_stack,
                   concatenate as np_concatenate,
                   copy as np_copy,
                   cos as np_cos,
                   delete as np_delete,
                   fill_diagonal as np_fill_diagonal,
                   finfo as np_finfo,
                   hypot as np_hypot,
                   log10 as np_log10,
                   max as np_max,
                   mean as np_mean,
                   median as np_median,
                   min as np_min,
                   newaxis as np_newaxis,
                   ones as np_ones,
                   pi as np_pi,
                   reshape as np_reshape,
                   seterr as np_seterr,
                   seterr as np_seterr,
                   shape as np_shape,
                   sin as np_sin,
                   size as np_size,
                   sort as np_sort,
                   square as np_square,
                   std as np_std,
                   sum as np_sum,
                   sqrt as np_sqrt,
                   tile as np_tile,
                   unravel_index as np_unravel_index,
                   where as np_where,
                   zeros as np_zeros)
from numpy.random import randn as np_randn
import scipy.ndimage as ndi
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.cluster.vq import kmeans,vq
from scipy.misc import imsave

# GroopM imports
from profileManager import ProfileManager
from binManager import BinManager, CenterFinder
import groopmTimekeeper as gtime
import refine
from PCA import PCA, Center
from groopmExceptions import *

np_seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ClusterEngine:
    """Top level interface for clustering contigs"""
    def __init__(self,
                 dbFileName,
                 plot=False,
                 finalPlot=False,
                 force=False,
                 numImgMaps=1,
                 minSize=5,
                 minVol=1000000,
                 squish=False):
        # worker classes
        self.PM = ProfileManager(dbFileName, squish=squish) # store our data
        self.BM = BinManager(pm=self.PM, minSize=minSize, minVol=minVol)

        # heat maps
        self.numImgMaps = numImgMaps
        self.imageMaps = np_zeros((self.numImgMaps,self.PM.scaleFactor,self.PM.scaleFactor))
        self.blurredMaps = np_zeros((self.numImgMaps,self.PM.scaleFactor,self.PM.scaleFactor))

        # we need a way to reference from the imageMaps back onto the transformed data
        self.im2RowIndices = {}

        # When blurring the raw image maps I chose a radius to suit my data, you can vary this as you like
        self.blurRadius = 2
        self.span = 45                  # amount we can travel about when determining "hot spots"


        self.HP = HoughPartitioner()

        # misc
        self.forceWriting = force
        self.debugPlots = plot
        self.finalPlot = finalPlot
        self.imageCounter = 1           # when we print many images
        self.roundNumber = 0            # how many times have we tried to make a bin?

    def promptOnOverwrite(self, minimal=False):
        """Check that the user is ok with possibly overwriting the DB"""
        if(self.PM.isClustered()):
            if(not self.forceWriting):
                input_not_ok = True
                valid_responses = ['Y','N']
                vrs = ",".join([str.lower(str(x)) for x in valid_responses])
                while(input_not_ok):
                    if(minimal):
                        option = raw_input(" Overwrite? ("+vrs+") : ")
                    else:
                        option = raw_input(" ****WARNING**** Database: '"+self.PM.dbFileName+"' has already been clustered.\n" \
                                           " If you continue you *MAY* overwrite existing bins!\n" \
                                           " Overwrite? ("+vrs+") : ")
                    if(option.upper() in valid_responses):
                        print "****************************************************************"
                        if(option.upper() == "N"):
                            print "Operation cancelled"
                            return False
                        else:
                            break
                    else:
                        print "Error, unrecognised choice '"+option.upper()+"'"
                        minimal = True
            print "Will Overwrite database",self.PM.dbFileName
        return True

#------------------------------------------------------------------------------
# CORE CONSTRUCTION AND MANAGEMENT

    def makeCores(self, timer, coreCut, gf=""):
        """Cluster the contigs to make bin cores"""
        # check that the user is OK with nuking stuff...
        if(not self.promptOnOverwrite()):
            return False

        RE = refine.RefineEngine(timer, BM=self.BM)

        # get some data
        self.PM.loadData(timer, loadRawKmers=True, condition="length >= "+str(coreCut))
        print "    %s" % timer.getTimeStamp()

        # transform the data
        print "Apply data transformations"
        self.PM.transformCP(timer)
        # plot the transformed space (if we've been asked to...)
        if(self.debugPlots):
            self.PM.renderTransCPData()
        print "    %s" % timer.getTimeStamp()

        # cluster and bin!
        print "Create cores"
        cum_contigs_used_good = self.initialiseCores(RE)
        print "    %s" % timer.getTimeStamp()

        # condense cores
        print "Refine cores [begin: %d]" % len(self.BM.bins)
        if self.finalPlot:
            prfx = "CORE"
        else:
            prfx = ""
        RE.refineBins(timer, auto=True, saveBins=False, plotFinal=prfx, gf=gf)

        # Now save all the stuff to disk!
        print "Saving bins"
        self.BM.saveBins(nuke=True)
        print "    %s" % timer.getTimeStamp()

    def initialiseCores(self, RE):
        """Process contigs and form CORE bins"""
        num_below_cutoff = 0            # how many consecutive attempts have produced small bins
        breakout_point = 100            # how many will we allow before we stop this loop

        # First we need to find the centers of each blob.
        # We can make a heat map and look for hot spots
        self.populateImageMaps()
        sub_counter = 0
        print "     .... .... .... .... .... .... .... .... .... ...."
        print "%4d" % sub_counter,
        new_line_counter = 0
        num_bins = 0

        prev_putative_clusters = []
        while(num_below_cutoff < breakout_point):
            stdout.flush()

            # apply a gaussian blur to each image map to make hot spots
            # stand out more from the background
            self.blurMaps()

            # now search for the "hottest" spots on the blurred map
            # and check for possible bin centroids
            bids_made = []
            putative_clusters = self.findNewClusterCenters()

            if(putative_clusters is None):
                break
            else:
                partitions = putative_clusters[0]
                [max_blur_value, max_x, max_y] = putative_clusters[1]
                self.roundNumber += 1
                sub_round_number = 1
                for center_row_indices in partitions:
                    # some of these row indices may have been eaten in a call to
                    # bin.recruit. We need to fix this now!
                    cri = np_array([ri for ri in center_row_indices if (ri not in self.PM.binnedRowIndices and ri not in self.PM.restrictedRowIndices)])
                    if len(cri) > 0:
                        center_row_indices = cri
                        total_BP = np_sum(self.PM.contigLengths[center_row_indices])
                        bin_size = len(center_row_indices)
                    else:
                        total_BP = 0
                        bin_size = 0

                    if self.BM.isGoodBin(total_BP, bin_size):   # Can we trust very small bins?.
                        # time to make a bin
                        bin = self.BM.makeNewBin(rowIndices=center_row_indices)

                        # work out the distribution in points in this bin
                        bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)

                        # append this bins list of mapped rowIndices to the main list
                        bids_made.append(bin.id)
                        num_bins += 1
                        self.updatePostBin(bin)
                        num_below_cutoff = 0

                        if(self.debugPlots):
                            bin.plotBin(self.PM.transformedCP, self.PM.contigColors, self.PM.kmerVals, self.PM.contigLengths, fileName="FRESH_"+str(self.imageCounter))
                            self.imageCounter += 1
                            self.plotHeat("HM_%d.%d.png" % (self.roundNumber, sub_round_number), max=max_blur_value, x=max_x, y=max_y)
                            sub_round_number += 1
                    else:
                        # this partition was too small, restrict these guys we don't run across them again
                        self.restrictRowIndices(center_row_indices)
                        num_below_cutoff += 1

                # did we do anything?
                num_bids_made = len(bids_made)
                if(num_bids_made == 0):
                    num_below_cutoff += 1
                    # nuke the lot!
                    for row_indices in partitions:
                        self.restrictRowIndices(row_indices)

                # try to merge these guys here...
                if num_bids_made > 1:
                    RE.mergeSimilarBins(None,
                                        None,
                                        bids=bids_made,
                                        loose=2.,
                                        silent=True)

                # do some post processing
                for bid in bids_made:
                    try:
                        bin = self.BM.getBin(bid)

                        # recruit more contigs
                        bin.recruit(self.PM.transformedCP,
                                    self.PM.averageCoverages,
                                    self.PM.kmerVals,
                                    self.PM.contigLengths,
                                    self.im2RowIndices,
                                    self.PM.binnedRowIndices,
                                    self.PM.restrictedRowIndices
                                    )
                        self.updatePostBin(bin)

                        new_line_counter += 1
                        print "% 4d" % bin.binSize,

                        # make the printing prettier
                        if(new_line_counter > 9):
                            new_line_counter = 0
                            sub_counter += 10
                            print "\n%4d" % sub_counter,

                        bin.plotBin(self.PM.transformedCP, self.PM.contigColors, self.PM.kmerVals, self.PM.contigLengths, fileName="P_BIN_%d"%(bin.id)) #***slow plot!

                    except BinNotFoundException: pass

        print "\n     .... .... .... .... .... .... .... .... .... ...."

    def findNewClusterCenters(self):
        """Find a putative cluster"""
        inRange = lambda x,l,u : x >= l and x < u

        # we work from the top view as this has the base clustering
        max_index = np_argmax(self.blurredMaps[0])
        max_value = self.blurredMaps[0].ravel()[max_index]

        max_x = int(max_index/self.PM.scaleFactor)
        max_y = max_index - self.PM.scaleFactor*max_x
        max_z = -1

        ret_values = [max_value, max_x, max_y]

        start_span = int(1.5 * self.span)
        span_len = 2*start_span+1

        if(self.debugPlots):
            self.plotRegion(max_x,max_y,max_z, fileName="Image_"+str(self.imageCounter), tag="column", column=True)
            self.imageCounter += 1

        # make a 3d grid to hold the values
        working_block = np_zeros((span_len, span_len, self.PM.scaleFactor))

        # go through the entire column
        (x_lower, x_upper) = self.makeCoordRanges(max_x, start_span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, start_span)
        super_putative_row_indices = []
        for p in self.im2RowIndices:
            if inRange(p[0],x_lower,x_upper) and inRange(p[1],y_lower,y_upper):
                for row_index in self.im2RowIndices[p]:
                    # check that the point is real and that it has not yet been binned
                    if row_index not in self.PM.binnedRowIndices and row_index not in self.PM.restrictedRowIndices:
                        # this is an unassigned point.
                        multiplier = np_log10(self.PM.contigLengths[row_index])
                        self.incrementAboutPoint3D(working_block, p[0]-x_lower, p[1]-y_lower, p[2],multiplier=multiplier)
                        super_putative_row_indices.append(row_index)

        # blur and find the highest value
        bwb = ndi.gaussian_filter(working_block, 8)#self.blurRadius)
        densest_index = np_unravel_index(np_argmax(bwb), (np_shape(bwb)))
        max_x = densest_index[0] + x_lower
        max_y = densest_index[1] + y_lower
        max_z = densest_index[2]

        # now get the basic color of this dense point
        putative_center_row_indices = []

        (x_lower, x_upper) = self.makeCoordRanges(max_x, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, self.span)
        (z_lower, z_upper) = self.makeCoordRanges(max_z, 2*self.span)

        for row_index in super_putative_row_indices:
            p = np_around(self.PM.transformedCP[row_index])
            if inRange(p[0],x_lower,x_upper) and inRange(p[1],y_lower,y_upper) and inRange(p[2],z_lower,z_upper):
                # we are within the range!
                putative_center_row_indices.append(row_index)

        putative_center_row_indices = np_array(putative_center_row_indices)

        # make sure we have something to go on here
        if(np_size(putative_center_row_indices) == 0):
            # it's all over!
            return None

        if(np_size(putative_center_row_indices) == 1):
            # get out of here but keep trying
            # the calling function may restrict these indices
            return [[np_array(putative_center_row_indices)], ret_values]
        else:
            total_BP = np_sum(self.PM.contigLengths[putative_center_row_indices])
            if not self.BM.isGoodBin(total_BP, len(putative_center_row_indices), ms=5): # Can we trust very small bins?.
                # get out of here but keep trying
                # the calling function should restrict these indices
                return [[np_array(putative_center_row_indices)], ret_values]
            else:
                putative_clusters = self.smartTwoWayContraction(putative_center_row_indices)
                if putative_clusters == None:
                  return None

                return [putative_clusters, ret_values]

    def smartPartDevo(self, rowIndices):
      """Partition a collection of contigs into 'core' groups"""
      from scipy.spatial import KDTree as kdt
      from scipy.cluster.vq import kmeans, whiten, vq

      # sanity check that there is enough data here to try a determine 'core' groups
      total_BP = np_sum(self.PM.contigLengths[rowIndices])
      if not self.BM.isGoodBin(total_BP, len(rowIndices), ms=5): # Can we trust very small bins?.
          # get out of here but keep trying
          # the calling function should restrict these indices
          return [np_array(rowIndices)]

      # try clustering dat using simple k-means
      k_max = 100

      # whiten data so each dimension is given equal weight
      c_whiten_dat = whiten(self.PM.transformedCP[rowIndices])

      # identify contigs that are 'noise'
      k_noise = 5
      search_tree = kdt(c_whiten_dat)
      dist_kth_contig = []
      for index in xrange(len(c_whiten_dat)):
          # get the K closest contigs
          k = np_min([k_noise, len(c_whiten_dat)-1])
          dist = search_tree.query(c_whiten_dat[index], k=k_noise)[0][k_noise-1]
          dist_kth_contig.append(dist)

      noise_threshold = np_mean(dist_kth_contig) + 2*np_std(dist_kth_contig)
      core_contigs = []
      core_contig_indices = []
      num_noisy_contig = 0
      for index in xrange(len(c_whiten_dat)):
          # get the K closest contigs
          k = np_min([k_noise, len(c_whiten_dat)-1])
          dist_kth_contig = search_tree.query(c_whiten_dat[index], k=k_noise)[0][k_noise-1]
          if dist_kth_contig < noise_threshold:
            core_contigs.append(c_whiten_dat[index])
            core_contig_indices.append(index)
          else:
            num_noisy_contig += 1

      core_contigs = np_array(core_contigs)

      # determine number of clusters in data
      avg_silhouette_max = -1
      k_best = 0
      for k in xrange(2, k_max):
        codebook, distortion = kmeans(core_contigs, k, iter=10, thresh=1e-05)
        cluster_index, dist = vq(core_contigs, codebook)

        # calculate CH index
        within = 0
        between = 0
        for i in xrange(len(core_contigs)):
          cluster_idI = cluster_index[i]
          for j in xrange(i+1, len(core_contigs)):
            cluster_idJ = cluster_index[j]
            if cluster_idI == cluster_idJ:
              within += np_sum((core_contigs[i] - core_contigs[j])**2)
            else:
              between += np_sum((core_contigs[i] - core_contigs[j])**2)

        ch_index = (between / (k - 1)) / (within / (len(core_contigs)-k))

        # calculate silhouette index
        avg_silhouette = 0
        for i in xrange(len(core_contigs)):
          cluster_idI = cluster_index[i]

          dist = [[] for x in xrange(k)]
          for j in xrange(0, len(core_contigs)):
            cluster_idJ = cluster_index[j]
            dist[cluster_idJ].append(np_sum(np_abs(core_contigs[i] - core_contigs[j])))

          a_i = np_mean(dist[cluster_idI])
          b_i = 1e10
          for j in xrange(0, k):
            if j != cluster_idI:
              mean_dist = np_mean(dist[j])
              if mean_dist < b_i:
                b_i = mean_dist

          avg_silhouette += (b_i - a_i) / max(a_i, b_i)

        avg_silhouette /= len(core_contigs)

        if avg_silhouette > avg_silhouette_max:
          avg_silhouette_max = avg_silhouette
          k_best = k

        # check if clustering results are becoming worse and we can bail on
        # trying larger k values
        if abs(k-k_best) > 3 and k > 10:
          break

      # determine final clustering with optimal K value
      codebook, distortion = kmeans(core_contigs, k_best, iter=100, thresh=1e-05)
      cluster_index, dist = vq(core_contigs, codebook)

      # create paritioning of data into 'core' clusters
      ret_parts = [[] for x in xrange(k_best)]
      for i in xrange(len(core_contigs)):
        ret_parts[cluster_index[i]].append(rowIndices[core_contig_indices[i]])
      ret_parts = [np_array(x) for x in ret_parts]

      return np_array(ret_parts)

    def smartTwoWayContraction(self, rowIndices):
      """Partition a collection of contigs into 'core' groups"""
      num_iterations = 10

      k_eps = np_max([0.05 * len(rowIndices), np_min([10, len(rowIndices)-1])])

      k_move_perc = 0.3
      c_move_perc = 0.2

      # sanity check that there is enough data here to try a determine 'core' groups
      total_BP = np_sum(self.PM.contigLengths[rowIndices])
      if not self.BM.isGoodBin(total_BP, len(rowIndices), ms=5): # Can we trust very small bins?.
          # get out of here but keep trying
          # the calling function should restrict these indices
          return [np_array(rowIndices)]

      # make a copy of the data we'll be munging
      k_dat = np_copy(np_col_stack((self.PM.kmerVals[rowIndices], self.PM.kmerVals2[rowIndices])))
      c_dat = np_copy(self.PM.transformedCP[rowIndices])
      l_dat = np_copy(self.PM.contigLengths[rowIndices])
      col_dat = np_copy(self.PM.contigColors[rowIndices])
      row_indices = np_copy(rowIndices)

      # calculate radius threshold in whitened transformed coverage space
      c_std = np_std(c_dat, axis=0)
      try:
        c_whiten_dat = (c_dat-np_mean(c_dat, axis=0)) / c_std
      except:
        # data has zero standard deviation ?
        return [np_array(row_indices)]

      dist_matrix = squareform(pdist(c_whiten_dat))
      c_radius = np_median(np_sort(dist_matrix)[:,k_eps-1])

      # calculate radius threshold in whitened kmer space
      k_std = np_std(k_dat, axis=0)
      try:
        k_whiten_dat = k_dat / k_std
      except:
        # data has zero standard deviation ?
        return [np_array(row_indices)]

      dist_matrix = squareform(pdist(k_whiten_dat))
      k_radius = np_median(np_sort(dist_matrix)[:,k_eps-1])

      # perform two-way contraction magic
      iter = 0
      while iter < num_iterations:
        iter += 1

        if True:
          if iter == 1:
            try:
              self.cluster_num
            except:
              self.cluster_num = 0

            self.cluster_num += 1

          # use HSV to RGB to generate colors
          S = 1       # SAT and VAL remain fixed at 1. Reduce to make
          V = 1       # Pastels if that's your preference...
          num_points = len(row_indices)
          disp_cols = np_array([htr(val, S, V) for val in k_dat[:,0]]).reshape(((num_points, 3)))

          fig = plt.figure()
          ax = fig.add_subplot(111, projection='3d')

          ax.scatter(c_dat[:,0],
                     c_dat[:,1],
                     c_dat[:,2],
                     edgecolors=disp_cols,
                     c=disp_cols,
                     marker='.')

          title = "Points: " + str(len(c_dat[:,0]))
          plt.title(title)

          if iter == 1:
            xlim = [np_min(c_dat[:,0]) - 0.05*np_min(c_dat[:,0]), np_max(c_dat[:,0]) + 0.05*np_max(c_dat[:,0])]
            ylim = [np_min(c_dat[:,1]) - 0.05*np_min(c_dat[:,1]), np_max(c_dat[:,1]) + 0.05*np_max(c_dat[:,1])]
            zlim = [np_min(c_dat[:,2]) - 0.05*np_min(c_dat[:,2]), np_max(c_dat[:,2]) + 0.05*np_max(c_dat[:,2])]

          ax.set_xlim(xlim)
          ax.set_ylim(ylim)
          ax.set_zlim(zlim)

          fig.set_size_inches(6,6)

          #fileName = "../../images/gh_%d_%d" % (self.cluster_num, iter)
          fileName = "gh_%d_%d" % (self.cluster_num, iter)
          plt.savefig(fileName + '.png',dpi=96)

          plt.close(fig)
          del fig

        # calculate distance matrices
        c_dist_matrix = squareform(pdist(c_whiten_dat))
        k_dist_matrix = squareform(pdist(k_whiten_dat))

        # find nearest neighbours to each point in normalized transformed coverage space,
        # and use this to converage a point's kmer profile
        new_k_dat = np_zeros(k_dat.shape)
        k_putative_noise = set()
        for index, row in enumerate(c_dist_matrix):
          neigbhours = np_where(row <= c_radius)[0]
          if len(neigbhours) > np_max([1, 0.1*k_eps]):
            neigbhours = neigbhours[1:] # ignore self match
          else:
            # no neighbours, so continue on
            k_putative_noise.add(index)
            new_k_dat[index] = k_dat[index]
            continue

          # use distance between kmer profiles as weights for moving similar
          # points towards each other; a minimum distance based of the kmer
          # radius is used to avoid zeros and ensure all neighbours provide
          # some weight
          neighbour_dist = k_dist_matrix[index][neigbhours] + 0.1 * k_radius

          # move point towards neighbours using inverse distance weighting
          inv_dist = 1.0 / neighbour_dist
          sum_inv_dist = np_sum(inv_dist)
          neighbour_weights = inv_dist / sum_inv_dist
          new_k_dat[index] = (1-k_move_perc) * k_dat[index] + k_move_perc * np_sum( (k_dat[neigbhours].T * neighbour_weights).T, axis = 0 )

        k_dat = new_k_dat

        # find nearest neighbours to each point in normalized kmer space,
        # and use this to converage a point's coverage profile
        new_c_dat = np_zeros(c_dat.shape)
        c_putative_noise = set()
        for index, row in enumerate(k_dist_matrix):
          neigbhours = np_where(row <= k_radius)[0]
          if len(neigbhours) > np_max([1, 0.1*k_eps]):
            neigbhours = neigbhours[1:] # ignore self match
          else:
            # no neighbours, so continue on
            c_putative_noise.add(index)
            new_c_dat[index] = c_dat[index]
            continue

          # use distance between coverage profiles as weights for moving similar
          # points towards each other; a minimum distance based of the kmer
          # radius is used to avoid zeros and ensure all neighbours provide
          # some weight
          neighbour_dist = c_dist_matrix[index][neigbhours] + 0.1 * c_radius

          # move point towards neighbours using inverse distance weighting
          inv_dist = 1.0 / neighbour_dist
          sum_inv_dist = np_sum(inv_dist)
          neighbour_weights = inv_dist / sum_inv_dist
          new_c_dat[index] = (1-c_move_perc) * c_dat[index] + c_move_perc * np_sum( (c_dat[neigbhours].T * neighbour_weights).T, axis = 0 )

        c_dat = new_c_dat

        # remove points that have no neighbours, unless they are long enough to be interesting
        noise = []
        putative_noise = c_putative_noise.intersection(k_putative_noise)
        #putative_noise = c_putative_noise.union(k_putative_noise)
        for index in putative_noise:
          if not self.BM.isGoodBin(l_dat[index], 1):
            noise.append(index)

        if len(noise) > 0:
          c_dat = np_delete(c_dat, noise, axis = 0)
          k_dat = np_delete(k_dat, noise, axis = 0)
          c_whiten_dat = np_delete(c_whiten_dat, noise, axis = 0)
          k_whiten_dat = np_delete(k_whiten_dat, noise, axis = 0)
          l_dat = np_delete(l_dat, noise, axis = 0)
          col_dat = np_delete(col_dat, noise, axis = 0)
          row_indices = np_delete(row_indices, noise, axis = 0)

        if len(row_indices) == 0:
          return None

      # perform hough transform clustering
      if True:
        self.HP.hc += 1
        data = k_dat[:,0]
        data -= np_min(data)
        try:
            data /= np_max(data)
        except FloatingPointError:
            pass

        #k_partitions = self.HP.houghPartition(data)
        k_partitions = self.HP.houghPartition(data, l_dat, imgTag="MER")

        if(len(k_partitions) == 0):
          return None

        partitions = []

        #-----------------------
        # GRID
        fig = plt.figure()

        orig_k_dat = self.PM.kmerVals[rowIndices]
        orig_k2_dat = self.PM.kmerVals2[rowIndices]
        orig_c_dat = self.PM.transformedCP[rowIndices][:,2]/10
        orig_l_dat = np_sqrt(self.PM.contigLengths[rowIndices])
        orig_col_dat = self.PM.contigColors[rowIndices]

        ax = plt.subplot(221)
        plt.xlabel("PCA1")
        plt.ylabel("PCA2")
        
        ax.scatter(orig_k_dat, orig_k2_dat, edgecolors=orig_col_dat, c=orig_col_dat, s=orig_l_dat)
        
        ax = plt.subplot(223)
        plt.title("%s contigs" % len(rowIndices))
        plt.xlabel("MER PARTS")
        plt.ylabel("COV PARTS")
        
        ax.scatter(orig_k_dat, orig_c_dat, edgecolors=orig_col_dat, c=orig_col_dat, s=orig_l_dat)

        
        c_plot_data = np_copy(c_dat[:,2])/10
        k_plot_data = np_copy(k_dat[:,0])
        k2_plot_data = np_copy(k_dat[:,1])
        k_sorted_indices = np_argsort(k_plot_data)
        c_sorted_indices = np_argsort(c_plot_data)
        k_plot_data = k_plot_data[k_sorted_indices]
        c_plot_data = c_plot_data[c_sorted_indices]
        c_max = np_max(c_plot_data) * 1.1
        k_max = np_max(k_plot_data) * 1.1
        c_min = np_min(c_plot_data) * 0.9
        k_min = np_min(k_plot_data) * 0.9
        k_eps = (k_max - k_min) / len(row_indices)
        c_eps = (c_max - c_min) / len(row_indices)
        cols=self.PM.contigColors[row_indices]
        lens = np_sqrt(self.PM.contigLengths[row_indices])

        start = 0
        k_lines = []
        k_sizes = [len(p) for p in k_partitions]
        for k in range(len(k_sizes)-1):
            k_lines.append(k_plot_data[k_sizes[k]+start]+k_eps)
            start += k_sizes[k]

        k_temp = {}
        c_temp = {}
        for ii in range(len(row_indices)):
            k_temp[row_indices[k_sorted_indices[ii]]] = ii
            c_temp[row_indices[c_sorted_indices[ii]]] = ii
        schooched_c = []
        schooched_k = []
        for ri in row_indices:
            schooched_k.append(k_temp[ri])
            schooched_c.append(c_temp[ri])
        schooched_k = np_array(schooched_k)
        schooched_c = np_array(schooched_c)

        ax = plt.subplot(222)
        plt.xlabel("PCA1")
        plt.ylabel("PCA2")

        ax.scatter(k_plot_data[schooched_k], k2_plot_data[schooched_k], edgecolors=cols, c=disp_cols, s=lens)

        ax = plt.subplot(224)
        plt.title("%s contigs" % len(row_indices))
        plt.xlabel("MER PARTS")
        plt.ylabel("COV PARTS")

        ax.scatter(k_plot_data[schooched_k], c_plot_data[schooched_c], edgecolors=cols, c=disp_cols, s=lens)

        for k in k_lines:
            plt.plot([k,k], [c_min, c_max], 'b-')

        pc = 0
        for k_part in k_partitions:
            pc += 1
            k_sep_indices = row_indices[k_part]
            part_bp = np_sum(l_dat[k_part])
            if self.BM.isGoodBin(part_bp, len(k_part), ms=5):

                data = np_copy(c_dat[k_part,2]/10)
                l_data = np_copy(l_dat[k_part])
                data -= np_min(data)
                try:
                    data /= np_max(data)
                except FloatingPointError:
                    pass

                c_partitions = self.HP.houghPartition(data, l_data, imgTag="COV")

                #-----
                # GRID
                c_plot_data = self.PM.transformedCP[k_sep_indices][:,2]/10
                c_plot_data = c_plot_data[np_argsort(c_plot_data)]

                start = 0
                c_lines = []
                c_sizes = [len(p) for p in c_partitions]
                for c in range(len(c_sizes)-1):
                    c_lines.append(c_plot_data[c_sizes[c]+start]+c_eps)
                    start += c_sizes[c]

                if pc == 1:
                    k_line_min = k_min
                else:
                    k_line_min = k_lines[pc-2]

                if pc == len(k_partitions):
                    k_line_max = k_max
                else:
                    k_line_max = k_lines[pc-1]


                for c in c_lines:
                    plt.plot([k_line_min,k_line_max], [c, c], 'g-')

                for c_part in c_partitions:
                    partitions.append(np_array(k_part[c_part]))

        ax.set_xlim(k_min, k_max)
        ax.set_ylim(c_min, c_max)
        fig.set_size_inches(12,12)
        plt.savefig("%d_GRID" % self.HP.hc,dpi=300)
        plt.close()
        del fig

        if len(partitions) == 0:
          return None

        ret_parts = []
        for p in partitions:
            ret_parts.append(np_array(row_indices[p]))

        return np_array(ret_parts)

      else:
        return [np_array(row_indices)]

    def smartPart(self, rowIndices):
        """Partition a collection of contigs into 'core' groups"""
        partitions = []
        total_BP = np_sum(self.PM.contigLengths[rowIndices])
        if not self.BM.isGoodBin(total_BP, len(rowIndices), ms=5): # Can we trust very small bins?.
            # get out of here but keep trying
            # the calling function should restrict these indices
            return [np_array(rowIndices)]
        else:
            # we've got a few good guys here, partition them up!
            self.HP.hc += 1

            data = np_copy(self.PM.kmerVals[rowIndices])
            data -= np_min(data)
            try:
                data /= np_max(data)
            except FloatingPointError:
                pass

            k_partitions = self.HP.houghPartition(data, self.PM.contigLengths[rowIndices], imgTag="MER")
            print "\n==========KKKKK==========", len(k_partitions)

            if(len(k_partitions) == 0):
                return None

            # for new plot!
            self.plotIndices(rowIndices, fileName="%d_CUT" % self.HP.hc)

            c_plot_data = self.PM.transformedCP[rowIndices][:,2]/10
            k_plot_data = np_copy(self.PM.kmerVals[rowIndices])
            k_sorted_indices = np_argsort(k_plot_data)
            c_sorted_indices = np_argsort(c_plot_data)
            k_plot_data = k_plot_data[k_sorted_indices]
            c_plot_data = c_plot_data[c_sorted_indices]
            c_max = np_max(c_plot_data) * 1.1
            k_max = np_max(k_plot_data) * 1.1
            c_min = np_min(c_plot_data) * 0.9
            k_min = np_min(k_plot_data) * 0.9
            k_eps = (k_max - k_min) / len(rowIndices)
            c_eps = (c_max - c_min) / len(rowIndices)

            start = 0
            k_lines = []
            k_sizes = [len(p) for p in k_partitions]
            for k in range(len(k_sizes)-1):
                k_lines.append(k_plot_data[k_sizes[k]+start]+k_eps)
                start += k_sizes[k]

            k_temp = {}
            c_temp = {}
            for ii in range(len(rowIndices)):
                k_temp[rowIndices[k_sorted_indices[ii]]] = ii
                c_temp[rowIndices[c_sorted_indices[ii]]] = ii
            schooched_c = []
            schooched_k = []
            for ri in rowIndices:
                schooched_k.append(k_temp[ri])
                schooched_c.append(c_temp[ri])
            schooched_k = np_array(schooched_k)
            schooched_c = np_array(schooched_c)

            cols=self.PM.contigColors[rowIndices]
            lens = np_sqrt(self.PM.contigLengths[rowIndices])

            fig = plt.figure()
            ax = plt.subplot(111)

            ax.scatter(k_plot_data[schooched_k], c_plot_data[schooched_c], edgecolors=cols, c=cols, s=lens)

            for k in k_lines:
                plt.plot([k,k], [c_min, c_max], 'b-')

            pc = 0
            for p in k_partitions:
                pc += 1
                k_sep_indices = rowIndices[p]
                part_bp = np_sum(self.PM.contigLengths[rowIndices])
                if self.BM.isGoodBin(part_bp, len(k_sep_indices), ms=5):
                    if pc == 1:
                        k_line_min = k_min
                    else:
                        k_line_min = k_lines[pc-2]

                    if pc == len(k_partitions):
                        k_line_max = k_max
                    else:
                        k_line_max = k_lines[pc-1]

                    data = np_copy(self.PM.transformedCP[k_sep_indices][:,2])
                    data -= np_min(data)
                    try:
                        data /= np_max(data)
                    except FloatingPointError:
                        pass

                    c_partitions = self.HP.houghPartition(data, self.PM.contigLengths[rowIndices], imgTag="COV(%0.2f-%0.2f)" % (k_line_min,k_line_max))
                    c_plot_data = self.PM.transformedCP[k_sep_indices][:,2]/10
                    c_plot_data = c_plot_data[np_argsort(c_plot_data)]

                    start = 0
                    c_lines = []
                    c_sizes = [len(p) for p in c_partitions]
                    for c in range(len(c_sizes)-1):
                        c_lines.append(c_plot_data[c_sizes[c]+start]+c_eps)
                        start += c_sizes[c]

                    for c in c_lines:
                        plt.plot([k_line_min,k_line_max], [c, c], 'g-')

                    print "\n==========CCCCCC==========", pc, len(c_partitions)
                    for p in c_partitions:
                        partitions.append(np_array(k_sep_indices[p]))


            ax.set_xlim(k_min, k_max)
            ax.set_ylim(c_min, c_max)
            fig.set_size_inches(6,6)
            plt.savefig("%d_GRID" % self.HP.hc,dpi=300)
            plt.close()
            del fig

            return partitions

    def partitionVals(self, vals, tag=None, stdevCutoff=0.04, maxSpread=0.15):
        """Work out where shifts in kmer/coverage vals happen"""
        cf = CenterFinder()
        partitions = []
        working_list = np_copy(vals)
        fix_dict = dict(zip(range(len(working_list)),range(len(working_list))))
        while(len(working_list) > 2):
            c_index = cf.findArrayCenter(working_list, tag=tag)
            expanded_indices = cf.expandSelection(c_index, working_list, stdevCutoff=stdevCutoff, maxSpread=maxSpread, tag=tag)
            # fix any munges from previous deletes
            morphed_indices = [fix_dict[i] for i in expanded_indices]
            partitions.append(morphed_indices)
            # shunt the indices to remove down!
            shunted_indices = []
            for offset, index in enumerate(expanded_indices):
                shunted_indices.append(index - offset)

            # make an updated working list and fix the fix dict
            nwl = []
            nfd = {}
            shifter = 0
            for i in range(len(working_list) - len(shunted_indices)):
                if(len(shunted_indices) > 0):
                    if(i >= shunted_indices[0]):
                        tmp = shunted_indices.pop(0)
                        shifter += 1
                        # consume any and all conseqs
                        while(len(shunted_indices) > 0):
                            if(shunted_indices[0] == tmp):
                                shunted_indices.pop(0)
                                shifter += 1
                            else:
                                break

                nfd[i] = fix_dict[i + shifter]
                nwl.append(working_list[i + shifter])

            fix_dict = nfd
            working_list = np_array(nwl)

        if(len(working_list) > 0):
            partitions.append(fix_dict.values())
        return partitions

#------------------------------------------------------------------------------
# DATA MAP MANAGEMENT

    def populateImageMaps(self):
        """Load the transformed data into the main image maps"""
        # reset these guys... JIC
        self.imageMaps = np_zeros((self.numImgMaps,self.PM.scaleFactor,self.PM.scaleFactor))
        self.im2RowIndices = {}

        # add to the grid wherever we find a contig
        row_index = -1
        for point in np_around(self.PM.transformedCP):
            row_index += 1
            # can only bin things once!
            if row_index not in self.PM.binnedRowIndices and row_index not in self.PM.restrictedRowIndices:
                # add to the row_index dict so we can relate the
                # map back to individual points later
                p = tuple(point)
                try:
                    self.im2RowIndices[p].append(row_index)
                except KeyError:
                    self.im2RowIndices[p] = [row_index]

                # now increment in the grid
                # for each point we encounter we incrmement
                # it's position + the positions to each side
                # and touching each corner
                self.incrementViaRowIndex(row_index, p)

    def incrementViaRowIndex(self, rowIndex, point=None):
        """Wrapper to increment about point"""
        if(point is None):
            point = tuple(np_around(self.PM.transformedCP[rowIndex]))
        multiplier = np_log10(self.PM.contigLengths[rowIndex])
        self.incrementAboutPoint(0, point[0], point[1], multiplier=multiplier)
        if(self.numImgMaps > 1):
            self.incrementAboutPoint(1, self.PM.scaleFactor - point[2] - 1, point[1], multiplier=multiplier)
            self.incrementAboutPoint(2, self.PM.scaleFactor - point[2] - 1, self.PM.scaleFactor - point[0] - 1, multiplier=multiplier)

    def decrementViaRowIndex(self, rowIndex, point=None):
        """Wrapper to decrement about point"""
        if(point is None):
            point = tuple(np_around(self.PM.transformedCP[rowIndex]))
        #px = point[0]
        #py = point[1]
        #pz = point[2]
        multiplier = np_log10(self.PM.contigLengths[rowIndex])
        self.decrementAboutPoint(0, point[0], point[1], multiplier=multiplier)
        if(self.numImgMaps > 1):
            self.decrementAboutPoint(1, self.PM.scaleFactor - point[2] - 1, point[1], multiplier=multiplier)
            self.decrementAboutPoint(2, self.PM.scaleFactor - point[2] - 1, self.PM.scaleFactor - point[0] - 1, multiplier=multiplier)

    def incrementAboutPoint(self, view_index, px, py, valP=1, valS=0.6, valC=0.2, multiplier=1):
        """Increment value at a point in the 2D image maps

        Increment point by valP, increment neighbouring points at the
        sides and corners of the target point by valS and valC

        multiplier is proportional to the contigs length
        """
        valP *= multiplier
        valS *= multiplier
        valC *= multiplier
        if px > 0:
            if py > 0:
                self.imageMaps[view_index,px-1,py-1] += valC      # Top left corner
            self.imageMaps[view_index,px-1,py] += valS            # Top
            if py < self.PM.scaleFactor-1:
                self.imageMaps[view_index,px-1,py+1] += valC      # Top right corner

        if py > 0:
            self.imageMaps[view_index,px,py-1] += valS            # Left side
        self.imageMaps[view_index,px,py] += valP                  # Point
        if py < self.PM.scaleFactor-1:
            self.imageMaps[view_index,px,py+1] += valS            # Right side

        if px < self.PM.scaleFactor-1:
            if py > 0:
                self.imageMaps[view_index,px+1,py-1] += valC      # Bottom left corner
            self.imageMaps[view_index,px+1,py] += valS            # Bottom
            if py < self.PM.scaleFactor-1:
                self.imageMaps[view_index,px+1,py+1] += valC      # Bottom right corner

    def decrementAboutPoint(self, view_index, px, py, valP=1, valS=0.6, valC=0.2, multiplier=1):
        """Decrement value at a point in the 2D image maps

        multiplier is proportional to the contigs length
        """
        valP *= multiplier
        valS *= multiplier
        valC *= multiplier
        if px > 0:
            if py > 0:
                self.safeDecrement(self.imageMaps[view_index], px-1, py-1, valC) # Top left corner
            self.safeDecrement(self.imageMaps[view_index], px-1, py, valS)       # Top
            if py < self.PM.scaleFactor-1:
                self.safeDecrement(self.imageMaps[view_index], px-1, py+1, valC) # Top right corner

        if py > 0:
            self.safeDecrement(self.imageMaps[view_index], px, py-1, valS)       # Left side
        self.safeDecrement(self.imageMaps[view_index], px, py, valP)             # Point
        if py < self.PM.scaleFactor-1:
            self.safeDecrement(self.imageMaps[view_index], px, py+1, valS)       # Right side

        if px < self.PM.scaleFactor-1:
            if py > 0:
                self.safeDecrement(self.imageMaps[view_index], px+1, py-1, valC) # Bottom left corner
            self.safeDecrement(self.imageMaps[view_index], px+1, py, valS)       # Bottom
            if py < self.PM.scaleFactor-1:
                self.safeDecrement(self.imageMaps[view_index], px+1, py+1, valC) # Bottom right corner

    def safeDecrement(self, map, px, py, value):
        """Decrement a value and make sure it's not negative or something shitty"""
        map[px][py] -= value
        if map[px][py] < np_finfo(float).eps:
            map[px][py] = 0

    def incrementAboutPoint3D(self, workingBlock, px, py, pz, vals=(6.4,4.9,2.5,1.6), multiplier=1):
        """Increment a point found in a 3D column

        used when finding the centroid of a hot area
        update the 26 points which surround the centre point
        z spans the height of the entire column, x and y have been offset to
        match the column subspace

        multiplier is proportional to the contigs length
        """
        valsM = [x*multiplier for x in vals]
        # top slice
        if pz < self.PM.scaleFactor-1:
            self.subIncrement3D(workingBlock, px, py, pz+1, valsM, 1)

        # center slice
        self.subIncrement3D(workingBlock, px, py, pz, valsM, 0)

        # bottom slice
        if pz > 0:
            self.subIncrement3D(workingBlock, px, py, pz-1, valsM, 1)

    def subIncrement3D(self, workingBlock, px, py, pz, vals, offset):
        """AUX: Called from incrementAboutPoint3D does but one slice

        multiplier is proportional to the contigs length
        """
        # get the size of the working block
        shape = np_shape(workingBlock)
        if px > 0:
            if py > 0:
                workingBlock[px-1,py-1,pz] += vals[offset + 2]      # Top left corner
            workingBlock[px-1,py,pz] += vals[offset + 1]            # Top
            if py < shape[1]-1:
                workingBlock[px-1,py+1,pz] += vals[offset + 2]      # Top right corner

        if py > 0:
            workingBlock[px,py-1,pz] += vals[offset + 1]            # Left side
        workingBlock[px,py,pz] += vals[offset]                      # Point
        if py < shape[1]-1:
            workingBlock[px,py+1,pz] += vals[offset + 1]            # Right side

        if px < shape[0]-1:
            if py > 0:
                workingBlock[px+1,py-1,pz] += vals[offset + 2]      # Bottom left corner
            workingBlock[px+1,py,pz] += vals[offset + 1]            # Bottom
            if py < shape[1]-1:
                workingBlock[px+1,py+1,pz] += vals[offset + 2]      # Bottom right corner

    def blurMaps(self):
        """Blur the 2D image maps"""
        self.blurredMaps = np_zeros((self.numImgMaps,self.PM.scaleFactor,self.PM.scaleFactor))
        for i in range(self.numImgMaps): # top, front and side
            self.blurredMaps[i,:,:] = ndi.gaussian_filter(self.imageMaps[i,:,:], 8)#self.blurRadius)

    def makeCoordRanges(self, pos, span):
        """Make search ranges which won't go out of bounds"""
        lower = pos-span
        upper = pos+span+1
        if(lower < 0):
            lower = 0
        if(upper > self.PM.scaleFactor):
            upper = self.PM.scaleFactor
        return (lower, upper)

    def updatePostBin(self, bin):
        """Update data structures after assigning contigs to a new bin"""
        bid = bin.id
        for row_index in bin.rowIndices:
            self.PM.binIds[row_index] = bid
            self.setRowIndexAssigned(row_index)

    def setRowIndexAssigned(self, rowIndex):
        """fix the data structures to indicate that rowIndex belongs to a bin

        Use only during initial core creation
        """
        if(rowIndex not in self.PM.restrictedRowIndices and rowIndex not in self.PM.binnedRowIndices):
            self.PM.binnedRowIndices[rowIndex] = True
            # now update the image map, decrement
            self.decrementViaRowIndex(rowIndex)

    def setRowIndexUnassigned(self, rowIndex):
        """fix the data structures to indicate that rowIndex no longer belongs to a bin

        Use only during initial core creation
        """
        if(rowIndex in self.PM.restrictedRowIndices and rowIndex not in self.PM.binnedRowIndices):
            del self.PM.binnedRowIndices[rowIndex]
            # now update the image map, increment
            self.incrementViaRowIndex(rowIndex)

    def restrictRowIndices(self, indices):
        """Add these indices to the restricted list"""
        for row_index in indices:
            # check that it's not binned or already restricted
            if(row_index not in self.PM.restrictedRowIndices and row_index not in self.PM.binnedRowIndices):
                self.PM.restrictedRowIndices[row_index] = True
                # now update the image map, decrement
                self.decrementViaRowIndex(row_index)

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING

    def plotIndices(self, rowIndices, fileName=""):
        """Plot these contigs in transformed space"""
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        disp_vals = np_array([])
        disp_cols = np_array([])
        disp_lens = np_array([])
        num_points = 0
        for row_index in rowIndices:
            num_points += 1
            disp_vals = np_append(disp_vals, self.PM.transformedCP[row_index])
            disp_cols = np_append(disp_cols, self.PM.contigColors[row_index])
            disp_lens = np_append(disp_lens, np_sqrt(self.PM.contigLengths[row_index]))

        # reshape
        disp_vals = np_reshape(disp_vals, (num_points, 3))
        disp_cols = np_reshape(disp_cols, (num_points, 3))

        ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors=disp_cols, c=disp_cols, s=disp_lens, marker='.')

        if(fileName != ""):
            fig.set_size_inches(6,6)
            plt.savefig(fileName,dpi=300)
        elif(show):
            plt.show()

        plt.close(fig)
        del fig

    def plotRegion(self, px, py, pz, fileName="", tag="", column=False):
        """Plot the region surrounding a point """
        import matplotlib as mpl
        disp_vals = np_array([])
        disp_cols = np_array([])
        num_points = 0
        # plot all points within span
        (z_lower, z_upper) = self.makeCoordRanges(pz, self.span)
        if(column):
            z_lower = 0
            z_upper = self.PM.scaleFactor - 1

        (x_lower, x_upper) = self.makeCoordRanges(px, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(py, self.span)
        for z in range(z_lower, z_upper):
            realz = self.PM.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.im2RowIndices):
                        for row_index in self.im2RowIndices[(x,y,realz)]:
                            if row_index not in self.PM.binnedRowIndices and row_index not in self.PM.restrictedRowIndices:
                                num_points += 1
                                disp_vals = np_append(disp_vals, self.PM.transformedCP[row_index])
                                disp_cols = np_append(disp_cols, self.PM.contigColors[row_index])

        # make a black mark at the max values
        small_span = self.span/2
        (x_lower, x_upper) = self.makeCoordRanges(px, small_span)
        (y_lower, y_upper) = self.makeCoordRanges(py, small_span)
        (z_lower, z_upper) = self.makeCoordRanges(pz, small_span)
        for z in range(z_lower, z_upper):
            realz = self.PM.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.im2RowIndices):
                        for row_index in self.im2RowIndices[(x,y,realz)]:
                            if row_index not in self.PM.binnedRowIndices and row_index not in self.PM.restrictedRowIndices:
                                num_points += 1
                                disp_vals = np_append(disp_vals, self.PM.transformedCP[row_index])
                                disp_cols = np_append(disp_cols, htr(0,0,0))
        # reshape
        disp_vals = np_reshape(disp_vals, (num_points, 3))
        disp_cols = np_reshape(disp_cols, (num_points, 3))

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        cm = mpl.colors.LinearSegmentedColormap('my_colormap', disp_cols, 1024)
        result = ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors=disp_cols, c=disp_cols, cmap=cm, marker='.')
        title = str.join(" ", ["Focus at: (",str(px), str(py), str(self.PM.scaleFactor - pz - 1),")\n",tag])
        plt.title(title)

        if(fileName != ""):
            fig.set_size_inches(6,6)
            plt.savefig(fileName,dpi=300)
        elif(show):
            plt.show()

        plt.close(fig)
        del fig

    def plotHeat(self, fileName = "", max=-1, x=-1, y=-1):
        """Print the main heat maps

        Useful for debugging
        """
        fig = plt.figure()
        images = []
        ax = None
        if(self.numImgMaps == 1):
            ax = fig.add_subplot(121)
            images.append(ax.imshow(self.blurredMaps[0,:,:]**0.5))
            if(max > 0):
                title = "Max value: %f (%f, %f)" % (max, x, y)
                plt.title(title)
        else:
            ax = fig.add_subplot(231)
            images.append(ax.imshow(self.blurredMaps[0,:,:]**0.5))
            if(max > 0):
                title = str.join(" ", ["Max value:",str(max)])
                plt.title(title)
            ax = fig.add_subplot(232)
            images.append(ax.imshow(self.blurredMaps[1,:,:]**0.5))
            ax = fig.add_subplot(233)
            images.append(ax.imshow(self.blurredMaps[2,:,:]**0.5))

        if(self.numImgMaps == 1):
            ax = fig.add_subplot(122)
            images.append(ax.imshow(self.imageMaps[0,:,:]**0.5))
        else:
            ax = fig.add_subplot(234)
            images.append(ax.imshow(self.imageMaps[0,:,:]**0.5))
            ax = fig.add_subplot(235)
            images.append(ax.imshow(self.imageMaps[1,:,:]**0.5))
            ax = fig.add_subplot(236)
            images.append(ax.imshow(self.imageMaps[2,:,:]**0.5))

        if(fileName != ""):
            if(self.numImgMaps == 1):
                fig.set_size_inches(12,6)
            else:
                fig.set_size_inches(18,18)

            plt.savefig(fileName,dpi=300)
        elif(show):
            plt.show()

        plt.close(fig)
        del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class HoughPartitioner:
    def __init__(self):
        self.hc = 0

    def houghPartition(self, dAta, lengths, level=0, side="C", imgTag=None):
        d_len_raw = int(len(dAta))
        if d_len_raw < 2:
            return np_array([[0]])

        if d_len_raw < 3:
            return np_array([[0,1]])

        # fudge the data to make longer contigs have more say in the
        # diff line we'll be making. This way we may be able to avoid lumping
        # super long contigs in with the short riff raff by accident.
        wobble = 0.05
        
        # we give you at least one point, but we give you 
        # more points if you're longer. Let's say 1 point per 10,000bp
        scales_per = {}
        for i in range(len(dAta)):
            scales_per[i] = int((lengths[i] - 1)/10000) + 1
            
        spread_data = []       
        spread2real = {}
        #real2spread = {}
        
        j = 0
        for i in range(len(dAta)):
            spread_data.append(dAta[i])
            spread2real[j] = i
#            try:
#                real2spread[i].append(j)
#            except KeyError:
#                real2spread[i] = [j]
#            rep = scales[int(np_log10(lengths[i]))]
            rep = scales_per[i] 
            j += 1 
            for k in range(1, rep):
                spread_data.append(wobble * np_randn() + dAta[i])
                spread2real[j] = i
#                try:
#                    real2spread[i].append(j)
#                except KeyError:
#                    real2spread[i] = [j]
                j += 1 
            
        # Force the data to fill the space
        spread_data = np_array(spread_data)
        d_len = len(spread_data)
        sorted_indices_raw = np_argsort(dAta)
        sorted_indices_spread = np_argsort(spread_data)
        data = spread_data[sorted_indices_spread]
        data -= np_min(data)    # shift to 0 but DO NOT scale to 1

        # make an array of "real" indices sorted by their spread data values.
        # NOTE: this array may / will contain multiple copies of each real index
        spread_sorted_reals = np_array([spread2real[ii] for ii in sorted_indices_spread])
        seens = {}  # keep track of how many times we've seen this guy
        # now we'd like to find out where the center of this guy lies.
        # Just take the mean
        real_center_in_spread = {}  # the center of the blob of "real" vales in the spread data
        for ii in range(len(spread_sorted_reals)):
            ri = spread_sorted_reals[ii]
            try:
                seens[ri].append(ii)
            except KeyError:
                seens[ri] = [ii]
        for ri in spread_sorted_reals:
            real_center_in_spread[ri] = np_mean(seens[ri])

        o_cutz = {}
        for ri in spread_sorted_reals:
            o_cutz[real_center_in_spread[ri]] = ri

        # work out weightings
        # we want to know how much each value differs from it's neighbours
        back_diffs = [float(data[i] - data[i-1]) for i in range(1,d_len)]
        diffs = [back_diffs[0]]
        for i in range(len(back_diffs)-1):
            diffs.append((back_diffs[i] + back_diffs[i+1])/2)
        diffs.append(back_diffs[-1])
        diffs = np_array(diffs)**2

        # replace the data array by the sum of it's diffs
        for i in range(1, d_len):
            diffs[i] += diffs[i-1]

        # scale to fit between 0 and len(diffs)
        # HT works better on a square
        diffs -= np_min(diffs)
        try:
            diffs /= np_max(diffs)
        except FloatingPointError:
            pass
        diffs *= len(diffs)

        t_data = np_array(zip(diffs, np_arange(d_len)))
        im_shape = (int(np_max(t_data, axis=0)[0]+1), d_len)

        (m, c, accumulator) = self.hough(t_data.astype(float), im_shape)
        # find the most prominent line

        # create the line we found and see which of the original points lie on
        # the line
        found_line = self.points2Line(np_array([[c,0],[m*im_shape[1]+c,im_shape[1]]]), im_shape[1], im_shape[0], 3)

        # we need to protect against the data line crossing
        # in and out of the "found line"
        in_block = False
        block_starts = []
        block_lens = []
        ii = -1
        for p in t_data.astype('int'):
            ii += 1
            if tuple(p) in found_line:
                if not in_block:
                    in_block = True
                    block_starts.append(ii)
            else:
                if in_block:
                    in_block = False
                    block_lens.append(ii - block_starts[-1] + 1)

        if in_block:
            # finishing block
            block_lens.append(ii - block_starts[-1] + 1)

        if imgTag is not None:
            #print "%d_%s_%s_%d DL: %d BL: %d" % (self.hc, imgTag, side, level, len(data), len(block_lens))
            # make a pretty picture
            fff = np_ones(im_shape) * 255
            for p in found_line.keys():
                fff[p[0],p[1]] = 220
            for p in t_data:
                fff[p[0],p[1]] = 0
            # scale so colors look sharper
            accumulator -= np_min(accumulator)
            accumulator /= np_max(accumulator)
            accumulator *= 255

            #imsave("%d_%s_%s_%d.png" % (self.hc, imgTag, side, level), np_concatenate([accumulator,fff]))

        if len(block_lens) == 0:
            return np_array([np_arange(len(dAta))])

        # work out what we'll keep and what we'll refine more
        longest_block = np_argmax(block_lens)
        spread_start = block_starts[longest_block]
        spread_end =  block_lens[longest_block] + spread_start

        left = []
        selected = []
        right = []

        for ii in range(spread_start):
            try:
                left.append(o_cutz[ii])
            except KeyError:
                pass
        for ii in range(spread_start, spread_end):
            try:
                selected.append(o_cutz[ii])
            except KeyError:
                pass
        for ii in range(spread_end, d_len):
            try:
                right.append(o_cutz[ii])
            except KeyError:
                pass

        left = np_array(left)
        selected = np_array(selected)
        right = np_array(right)

        rets = []
        if len(left) > 0:
            left_p = self.houghPartition(dAta[left], lengths[left], side="%sL" %side, level=level+1, imgTag=imgTag)
            for A in left_p:
                rets.append(np_array([left[i] for i in A]))

        if len(selected) > 0:
            rets.append(selected)

        if len(right) > 0:
            right_p = self.houghPartition(dAta[right], lengths[right], side="%sR" %side, level=level+1, imgTag=imgTag)
            for A in right_p:
                rets.append(np_array([right[i] for i in A]))

        if False:
            # this was kinds hard to write and
            # is still an idea worth pursuing.
            # can it for now...
            flat_data = []
            for i in range(len(dAta)):
                for k in range(scales[int(np_log10(lengths[i]))]):
                    flat_data.append(dAta[i])

            sorted_indices_flat = np_argsort(flat_data)
            data = np_array(flat_data)[sorted_indices_flat]
            data -= np_min(data)

            # work out weightings
            # we want to know how much each value differs from it's neighbours
            back_diffs = [float(data[i] - data[i-1]) for i in range(1,d_len)]
            diffs = [back_diffs[0]]
            for i in range(len(back_diffs)-1):
                diffs.append((back_diffs[i] + back_diffs[i+1])/2)
            diffs.append(back_diffs[-1])
            diffs = np_array(diffs)**2

            # replace the data array by the sum of it's diffs
            for i in range(1, d_len):
                diffs[i] += diffs[i-1]

            diffs -= np_min(diffs)
            try:
                diffs /= np_max(diffs)
            except FloatingPointError:
                pass
            diffs *= len(diffs)

            for ret in rets:
                spread_ret = []
                # these are real indices
                for index in ret:
                    # now get spread indices
                    for spread_index in real2spread[index]:
                        spread_ret.append(spread_index)
                # examine diffs[spread_ret],

        return np_array(rets)

    def points2Line(self, points, xIndexLim, yIndexLim, thickness):
        """Draw a thick line between a series of points"""
        line_points = []
        num_points = len(points)
        for i in range(1, num_points):
            # draw a line between this point and the last point
            x_gap = float(np_abs(points[i-1,1] - points[i,1]))
            y_gap = float(np_abs(points[i-1,0] - points[i,0]))
            largest_gap = np_max([x_gap, y_gap])

            if points[i-1,0] >= points[i,0]:
                Ys = [int(j) for j in np_around((np_arange(largest_gap)*y_gap/largest_gap) + points[i,0])]
            elif points[i,0] > points[i-1,0]:
                Ys = [int(j) for j in np_around((np_arange(largest_gap)*y_gap/largest_gap) + points[i-1,0])][::-1]
            if points[i-1,1] >= points[i,1]:
                Xs = [int(j) for j in np_around((np_arange(largest_gap)*x_gap/largest_gap) + points[i,1])]
            elif points[i,1] > points[i-1,1]:
                Xs = [int(j) for j in np_around((np_arange(largest_gap)*x_gap/largest_gap) + points[i-1,1])][::-1]

            for p in zip(Ys, Xs):
                line_points.append(p)

            # now make the line thicker
            thick_points = {}
            for point in line_points:
                for y in range(np_max([point[0]-thickness, 0]),np_min([point[0]+thickness+1,yIndexLim])):
                    for x in range(np_max([point[1]-thickness, 0]),np_min([point[1]+thickness+1,xIndexLim])):
                        thick_points[(y,x)] = 1

        return thick_points

    def hough(self, data, imShape):
        """Calculate Hough transform

        Data is a 2D numpy array"""
        (rows, cols) = imShape
        d_len = len(data)
        half_rows = rows/2
        if half_rows == 0:
            half_rows = 1
            rows = 2
        rmax = np_hypot(rows, cols)
        dr = rmax / (half_rows)
        dth = np_pi / cols
        accumulator = np_ones((rows * cols))*255

        """
        For speed we numpify this loop. I just keep this here
        so that I can remember what it is I am actually doing...

        Note that this needs the accumulator to be a 2D array
        ie.
        accumulator = np_ones((rows, cols))*255

        for p in data:
            for theta_index in range(cols):
                th = dth * theta_index
                r = p[1]*cos(th) + p[0]*sin(th)
                iry = half_rows + int(r/dr)
                accumulator[iry, theta_index] -= 1
        """
        cos_sin_array = np_array(zip([np_sin(dth * theta_index) for theta_index in range(cols)],
                                     [np_cos(dth * theta_index) for theta_index in range(cols)]))
        Rs = np_array(np_sum(np_reshape([p * cos_sin_array for p in data], (d_len*cols,2)),
                             axis=1)/dr).astype('int') + half_rows
        Cs = np_array(range(cols)*d_len)
        flat_indices = Rs * cols + Cs

        # update the accumulator with integer decrements
        decrements = np_bincount(flat_indices.astype('int'))
        index = 0
        for d in decrements:
            accumulator[index] -= d
            index += 1

        minindex = accumulator.argmin()

        # find the liniest line
        min_row = int(minindex/cols)
        min_col = minindex - (min_row*cols)
        theta = float(min_col) * dth
        rad = float(min_row - half_rows)*dr

        # now de hough!
        try:
            m = -1. * np_cos(theta) / np_sin(theta)
            c = rad / np_sin(theta)
        except FloatingPointError:
            # when we are trying to do a perfectly
            # straight line
            m = 0.
            c = rad

        # rounding errors suck
        if np_allclose([m], [0.]):
            m = 0.
        if np_allclose([c], [0.]):
            c = 0.

        if m < 0:
            m = 0.

        return (m, c, accumulator.reshape((rows, cols)))

###############################################################################
###############################################################################
###############################################################################
###############################################################################


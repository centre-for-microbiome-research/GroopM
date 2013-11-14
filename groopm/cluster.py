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
__copyright__ = "Copyright 2012/2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.2.10.12"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Alpha"

###############################################################################

from sys import stdout, exit

from colorsys import hsv_to_rgb as htr
import matplotlib.pyplot as plt
from pylab import show
from numpy import (abs as np_abs,
                   allclose as np_allclose,
                   append as np_append,
                   arange as np_arange,
                   argmax as np_argmax,
                   argsort as np_argsort,
                   around as np_around,
                   array as np_array,
                   bincount as np_bincount,
                   concatenate as np_concatenate,
                   copy as np_copy,
                   cos as np_cos,
                   delete as np_delete,
                   finfo as np_finfo,
                   hypot as np_hypot,
                   inf as np_inf,
                   log as np_log,
                   log10 as np_log10,
                   max as np_max,
                   mean as np_mean,
                   median as np_median,
                   min as np_min,
                   ones as np_ones,
                   pi as np_pi,
                   prod as np_prod,
                   reshape as np_reshape,
                   seterr as np_seterr,
                   shape as np_shape,
                   sin as np_sin,
                   size as np_size,
                   sort as np_sort,
                   std as np_std,
                   sum as np_sum,
                   sqrt as np_sqrt,
                   unravel_index as np_unravel_index,
                   where as np_where,
                   zeros as np_zeros)
from numpy.linalg import norm as np_norm
import scipy.ndimage as ndi
from scipy.spatial.distance import pdist, squareform, cityblock, euclidean
from scipy.misc import imsave

# GroopM imports
from profileManager import ProfileManager
from binManager import BinManager
from refine import GrubbsTester, RefineEngine
from PCA import PCA, Center
from groopmExceptions import BinNotFoundException

np_seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ClusterEngine:
    """Top level interface for clustering contigs"""
    def __init__(self,
                 dbFileName,
                 timer,
                 plot=0,
                 finalPlot=False,
                 force=False,
                 numImgMaps=1,
                 minSize=5,
                 minVol=1000000):

        # worker classes
        self.PM = ProfileManager(dbFileName) # store our data
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

        # Misc tools we'll need
        self.timer = timer
        self.RE = RefineEngine(self.timer, BM=self.BM)   # Basic refinement techniques
        self.HP = HoughPartitioner()                     # Finding cluster centers for unknown K
        self.GT = GrubbsTester()                         # Test length conformity

        # misc
        self.forceWriting = force
        self.debugPlots = plot
        self.finalPlot = finalPlot
        self.imageCounter = 1           # when we print many images
        self.roundNumber = 0            # how many times have we tried to make a bin?
        self.subRoundNumber = 0         # measure sub rounds too!
        self.TSpan = 0.                 # dist from centre to the corner

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

    def makeCores(self, coreCut, gf="", kmerThreshold=0.2, coverageThreshold=0.05):
        """Cluster the contigs to make bin cores"""
        # check that the user is OK with nuking stuff...
        if(not self.promptOnOverwrite()):
            return False

        # get some data
        self.PM.loadData(self.timer, "length >= "+str(coreCut))
        print "    %s" % self.timer.getTimeStamp()

        # transform the data
        print "    Loading transformed data"
        self.PM.transformCP(self.timer)
        # plot the transformed space (if we've been asked to...)
        #if(self.debugPlots >= 3):
        #    self.PM.renderTransCPData()

        # now we can make this guy
        self.TSpan = np_mean([np_norm(self.PM.corners[i] - self.PM.TCentre) for i in range(self.PM.numStoits)])

        print "    %s" % self.timer.getTimeStamp()

        # cluster and bin!
        print "Create cores"
        self.initialiseCores(kmerThreshold, coverageThreshold)
        print "    %s" % self.timer.getTimeStamp()

        # condense cores
        print "Refine cores [begin: %d]" % len(self.BM.bins)
        if self.finalPlot:
            prfx = "CORE"
        else:
            prfx = ""
        self.RE.refineBins(self.timer, auto=True, saveBins=False, plotFinal=prfx, gf=gf)

        # Now save all the stuff to disk!
        print "Saving bins"
        self.BM.saveBins(nuke=True)
        print "    %s" % self.timer.getTimeStamp()

    def initialiseCores(self, kmerThreshold, coverageThreshold):
        """Process contigs and form CORE bins"""
        num_below_cutoff = 0            # how many consecutive attempts have produced small bins
        breakout_point = 30            # how many will we allow before we stop this loop

        # First we need to find the centers of each blob.
        # We can make a heat map and look for hot spots
        self.populateImageMaps()
        sub_counter = 0
        print "     .... .... .... .... .... .... .... .... .... ...."
        print "%4d" % sub_counter,
        new_line_counter = 0
        num_bins = 0

        while(num_below_cutoff < breakout_point):
            stdout.flush()

            # apply a gaussian blur to each image map to make hot spots
            # stand out more from the background
            self.blurMaps()

            # now search for the "hottest" spots on the blurred map
            # and check for possible bin centroids
            putative_clusters = self.findNewClusterCenters(kmerThreshold, coverageThreshold)

            if(putative_clusters is None):
                break
            else:
                bids_made = []
                partitions = putative_clusters[0]
                [max_x, max_y] = putative_clusters[1]
                self.roundNumber += 1
                self.subRoundNumber = 1

                for center_row_indices in partitions:
                    total_BP = np_sum(self.PM.contigLengths[center_row_indices])
                    bin_size = len(center_row_indices)

                    if self.BM.isGoodBin(total_BP, bin_size, ms=3, mv=10000):   # Can we trust very small bins?.
                        # time to make a bin
                        bin = self.BM.makeNewBin(rowIndices=center_row_indices)

                        # work out the distribution in points in this bin
                        bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerNormPC1, self.PM.kmerPCs, self.PM.contigGCs, self.PM.contigLengths)

                        # append this bins list of mapped rowIndices to the main list
                        bids_made.append(bin.id)
                        num_bins += 1
                        self.updatePostBin(bin)

                        if(self.debugPlots >= 2):
                            bin.plotBin(self.PM.transformedCP, self.PM.contigGCs, self.PM.kmerNormPC1,
                                        self.PM.contigLengths, self.PM.colorMapGC, self.PM.isLikelyChimeric,
                                        fileName="FRESH_"+str(self.imageCounter))

                            self.imageCounter += 1
                            self.subRoundNumber += 1
                    else:
                        # this partition was too small, restrict these guys we don't run across them again
                        self.restrictRowIndices(center_row_indices)

                # did we do anything?
                num_bids_made = len(bids_made)
                if(num_bids_made == 0):
                    num_below_cutoff += 1
                    # nuke the lot!
                    for row_indices in partitions:
                        self.restrictRowIndices(row_indices)

                # do some post processing
                for bid in bids_made:
                    try:
                        bin = self.BM.getBin(bid)

                        # recruit more contigs
                        bin.recruit(self.PM,
                                    self.GT,
                                    self.im2RowIndices
                                    )
                        self.updatePostBin(bin)

                        new_line_counter += 1
                        print "% 4d" % bin.binSize,

                        # make the printing prettier
                        if(new_line_counter > 9):
                            new_line_counter = 0
                            sub_counter += 10
                            print "\n%4d" % sub_counter,

                        if(self.debugPlots >= 1):
                            #***slow plot!
                            bin.plotBin(self.PM.transformedCP, self.PM.contigGCs, self.PM.kmerNormPC1, self.PM.contigLengths, self.PM.colorMapGC, self.PM.isLikelyChimeric, fileName="CORE_BIN_%d"%(bin.id))

                    except BinNotFoundException: pass

        print "\n     .... .... .... .... .... .... .... .... .... ...."

    def findNewClusterCenters(self, kmerThreshold, coverageThreshold):
        """Find a putative cluster"""
        inRange = lambda x,l,u : x >= l and x < u

        # we work from the top view as this has the base clustering
        max_index = np_argmax(self.blurredMaps[0])
        max_value = self.blurredMaps[0].ravel()[max_index]

        max_x = int(max_index/self.PM.scaleFactor)
        max_y = max_index - self.PM.scaleFactor*max_x

        ret_values = [max_x, max_y]

        if(self.debugPlots >= 2):
            self.plotHeat("HM_%d.%d.png" % (self.roundNumber+1, self.subRoundNumber), x=max_x, y=max_y)

        start_span = int(1.5 * self.span)
        span_len = 2*start_span+1

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

        # z span should be greater at the corners and shallower at the center
        z_span_min = 100
        z_span_max = 400
        dist_from_centre = euclidean(self.PM.TCentre[:2], np_array([max_x,max_y]))
        z_span = (z_span_max - z_span_min)/self.PM.transRadius*dist_from_centre + z_span_min

        (x_lower, x_upper) = self.makeCoordRanges(max_x, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, self.span)
        (z_lower, z_upper) = self.makeCoordRanges(max_z, z_span)

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
                putative_clusters = self.twoWayContraction(putative_center_row_indices,
                                                           [max_x, max_y],
                                                           kmerThreshold,
                                                           coverageThreshold)
                if putative_clusters == None:
                    return None

                return [putative_clusters, ret_values]

    def twoWayContraction(self, rowIndices, positionInPlane, kmerThreshold, coverageThreshold):
        """Partition a collection of contigs into 'core' groups"""

        # sanity check that there is enough data here to try a determine 'core' groups
        total_BP = np_sum(self.PM.contigLengths[rowIndices])
        if not self.BM.isGoodBin(total_BP, len(rowIndices), ms=5): # Can we trust very small bins?.
            # get out of here but keep trying
            # the calling function should restrict these indices
            return [np_array(rowIndices)]

        # make a copy of the data we'll be munging
        k_dat = np_copy(self.PM.kmerPCs[rowIndices])
        #c_dat = np_copy(self.PM.covProfiles[rowIndices])
        c_dat = np_copy(self.PM.transformedCP[rowIndices])
        l_dat = np_copy(self.PM.contigLengths[rowIndices])
        if self.debugPlots >= 2:
            n_dat = np_copy(self.PM.contigNames[rowIndices])

        row_indices = np_copy(rowIndices)

        # calculate shortest distance to a corner (this isn't currently used)
        #min_dist_to_corner = 1e9
        #for corner in self.PM.corners:
        #    dist = cityblock(corner[:,[0,1]], positionInPlane)
        #   min_dist_to_corner = min(dist, min_dist_to_corner)

        # calculate radius threshold in whitened transformed coverage space
        #eps_neighbours = np_max([0.05 * len(rowIndices), np_min([10, len(rowIndices)-1])])
        eps_neighbours = np_min([10, int(len(rowIndices)/2)])

        # calculate mean and std in coverage space for whitening data
        c_mean = np_mean(c_dat, axis=0)
        c_std = np_std(c_dat, axis=0)
        c_std += np_where(c_std == 0, 1, 0) # make sure std dev is never zero
        c_whiten_dat = (c_dat-c_mean) / c_std
        c_whiten_dist = pdist(c_whiten_dat, 'cityblock')

        try:
            c_dist_matrix = squareform(c_whiten_dist)
        except MemoryError:
            print "\n"
            print '*******************************************************************************'
            print '*********************************    ERROR    *********************************'
            print '*******************************************************************************'
            print 'GroopM is attempting to do some maths on a putative bin which contains:'
            print
            print '\t\t%d contigs'  % (len(rowIndices))
            print
            print 'This has caused your machine to run out of memory.'
            print 'The most likely cause is that your samples are very different from each other.'
            print 'You can confirm this by running:'
            print
            print '\t\tgroopm explore -m allcontigs %s' % self.PM.dbFileName
            print
            print 'If you notice only vertical "spears" of contigs at the corners of the plot then'
            print 'this means that your samples are very different and you are not getting a good'
            print 'mapping from all samples to all contigs. You may get more mileage by assembling'
            print 'and binning your samples separately.'
            print
            print 'If you notice "clouds" of contigs then congratulations! You have found a bug.'
            print 'Please let me know at mike@mikeimelfort.com or via github.com/minillinim/GroopM'
            print
            print 'GroopM is aborting... sorry'
            print
            print '*******************************************************************************'
            print "\n"
            exit(-1)

        c_radius = np_median(np_sort(c_dist_matrix)[:,eps_neighbours-1])

        # calculate radius threshold in kmer space
        k_dist = pdist(k_dat, 'cityblock')
        k_dist_matrix = squareform(k_dist)
        k_radius = np_median(np_sort(k_dist_matrix)[:,eps_neighbours-1])

        # calculate convergence criteria
        k_converged = kmerThreshold * 30.0 #5e-2 * np_mean(k_dist)
        c_converged = coverageThreshold * 3.4  # 5e-2 * np_mean(c_whiten_dist)
        k_delt = 0.
        c_delt = 0.
        max_iterations = 1  # don't worry about params. Just do once

        k_move_perc = 0.1
        c_move_perc = 0.1

        # perform two-way contraction of kmer and coverage space
        iter = 0
        while iter < max_iterations:
            iter += 1

            if self.debugPlots >= 2:
                if iter == 1:
                    try:
                        self.cluster_num
                    except:
                        self.cluster_num = 0

                    self.cluster_num += 1

                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(c_dat[:,0],
                           c_dat[:,1],
                           c_dat[:,2],
                           edgecolors='k',
                           c=self.PM.contigGCs[row_indices],
                           cmap=self.PM.colorMapGC,
                           vmin=0.0, vmax=1.0,
                           s=np_sqrt(l_dat),
                           marker='.')

                title = "Points: %s Cd: %f Kd: %f" % (str(len(c_dat[:,0])), c_delt, k_delt)
                plt.title(title)

                if iter == 1:
                    xlim = [np_min(c_dat[:,0]) - 0.05*np_min(c_dat[:,0]), np_max(c_dat[:,0]) + 0.05*np_max(c_dat[:,0])]
                    ylim = [np_min(c_dat[:,1]) - 0.05*np_min(c_dat[:,1]), np_max(c_dat[:,1]) + 0.05*np_max(c_dat[:,1])]
                    zlim = [np_min(c_dat[:,2]) - 0.05*np_min(c_dat[:,2]), np_max(c_dat[:,2]) + 0.05*np_max(c_dat[:,2])]

                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                ax.set_zlim(zlim)

                fig.set_size_inches(6,6)

                fileName = "gh_%d_%d" % (self.cluster_num, iter)
                plt.savefig(fileName + '.png',dpi=96)

                plt.close(fig)
                del fig

            # calculate distance matrices
            c_dist_matrix = squareform(pdist(c_whiten_dat, 'cityblock'))
            k_dist_matrix = squareform(pdist(k_dat, 'cityblock'))

            # find nearest neighbours to each point in whitened coverage space,
            # and use this to converage a point's kmer profile
            new_k_dat = np_zeros(k_dat.shape)
            k_putative_noise = set()
            k_deltas = []
            for index, row in enumerate(c_dist_matrix):
                neigbhours = np_where(row <= c_radius)[0]
                if len(neigbhours) > np_max([1, 0.1*eps_neighbours]):
                    neigbhours = neigbhours[1:] # ignore self match
                else:
                    # extremely few neighbours so mark this as noise
                    k_putative_noise.add(index)
                    new_k_dat[index] = k_dat[index]
                    continue

                # use distance between kmer profiles as weights for moving similar
                # points towards each other; a minimum distance based on the kmer
                # radius is used to avoid zeros and ensure all neighbours provide
                # some weight
                neighbour_dist = k_dist_matrix[index][neigbhours] + 0.1 * k_radius

                # move point towards neighbours using inverse distance weighting
                try:
                    inv_dist = 1.0 / neighbour_dist
                except FloatingPointError:
                    inv_dist = 0.
                sum_inv_dist = np_sum(inv_dist)
                try:
                    neighbour_weights = inv_dist / sum_inv_dist
                except FloatingPointError:
                    neighbour_weights = 0.
                new_k_dat[index] = (1-k_move_perc) * k_dat[index] + k_move_perc * np_sum( (k_dat[neigbhours].T * neighbour_weights).T, axis = 0 )

                k_deltas.append(cityblock(k_dat[index], new_k_dat[index]))

            k_dat = new_k_dat

            # find nearest neighbours to each point in kmer space,
            # and use this to converage a point's coverage profile
            new_c_dat = np_zeros(c_dat.shape)
            c_putative_noise = set()
            c_deltas = []
            for index, row in enumerate(k_dist_matrix):
                neigbhours = np_where(row <= k_radius)[0]
                if len(neigbhours) > np_max([1, 0.1*eps_neighbours]):
                    neigbhours = neigbhours[1:] # ignore self match
                else:
                    # extremely few neighbours so mark this as nois
                    c_putative_noise.add(index)
                    new_c_dat[index] = c_dat[index]
                    continue

                # use distance between whitened coverage profiles as weights for moving similar
                # points towards each other; a minimum distance based of the coverage
                # radius is used to avoid zeros and ensure all neighbours provide some weight
                neighbour_dist = c_dist_matrix[index][neigbhours] + 0.1 * c_radius

                # move point towards neighbours using inverse distance weighting
                inv_dist = 1.0 / neighbour_dist
                sum_inv_dist = np_sum(inv_dist)
                neighbour_weights = inv_dist / sum_inv_dist
                new_c_dat[index] = (1-c_move_perc) * c_dat[index] + c_move_perc * np_sum( (c_dat[neigbhours].T * neighbour_weights).T, axis = 0 )
                new_c_whiten_dat = (1-c_move_perc) * c_whiten_dat[index] + c_move_perc * np_sum( (c_whiten_dat[neigbhours].T * neighbour_weights).T, axis = 0 )

                c_deltas.append(cityblock(c_whiten_dat[index], new_c_whiten_dat))

            c_dat = new_c_dat

            # remove points that have no or few neighbours in both spaces,
            # unless they are long enough to be of interest
            noise = []
            putative_noise = c_putative_noise.intersection(k_putative_noise)
            for index in putative_noise:
                if not self.BM.isGoodBin(l_dat[index], 1):
                    noise.append(index)

            if len(noise) > 0:
                c_dat = np_delete(c_dat, noise, axis = 0)
                k_dat = np_delete(k_dat, noise, axis = 0)
                l_dat = np_delete(l_dat, noise, axis = 0)
                row_indices = np_delete(row_indices, noise, axis = 0)
                if self.debugPlots >= 2:
                    #print "Noise deleting_%d_%d:\n" % (self.cluster_num, iter), n_dat[noise]
                    n_dat = np_delete(n_dat, noise, axis = 0)

            # get whitened version of modified coverage data (using original transformation)
            c_whiten_dat = (c_dat-c_mean) / c_std

            if len(row_indices) == 0:
                return None

            # check for convergence
            k_delt = np_mean(k_deltas)
            c_delt = np_mean(c_deltas)
            if k_delt < k_converged or c_delt < c_converged:
                break

        if self.debugPlots >= 2:

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(c_dat[:,0],
                       c_dat[:,1],
                       c_dat[:,2],
                       edgecolors='k',
                       c=self.PM.contigGCs[row_indices],
                       cmap=self.PM.colorMapGC,
                       vmin=0.0, vmax=1.0,
                       s=np_sqrt(l_dat),
                       marker='.')

            title = "Points: %s Cd: %f Kd: %f" % (str(len(c_dat[:,0])), c_delt, k_delt)
            plt.title(title)

            if iter == 1:
                xlim = [np_min(c_dat[:,0]) - 0.05*np_min(c_dat[:,0]), np_max(c_dat[:,0]) + 0.05*np_max(c_dat[:,0])]
                ylim = [np_min(c_dat[:,1]) - 0.05*np_min(c_dat[:,1]), np_max(c_dat[:,1]) + 0.05*np_max(c_dat[:,1])]
                zlim = [np_min(c_dat[:,2]) - 0.05*np_min(c_dat[:,2]), np_max(c_dat[:,2]) + 0.05*np_max(c_dat[:,2])]

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_zlim(zlim)

            fig.set_size_inches(6,6)

            fileName = "gh_%d_final" % (self.cluster_num)
            plt.savefig(fileName + '.png',dpi=96)

            plt.close(fig)
            del fig

        # test total GC spread
        # if it's very small we won't bother doing hough stuff
        gcs = self.PM.contigGCs[row_indices]
        gc_spread = np_max(gcs) - np_min(gcs)
        self.HP.hc += 1
        if gc_spread <= 0.15:
            # no hough
            (k_partitions, k_keeps) = (np_array([np_arange(len(row_indices))]),np_array([True]))
        else:
            # perform hough transform clustering
            if self.debugPlots >= 3:
                (k_partitions, k_keeps) = self.HP.houghPartition(k_dat[:,0], l_dat, imgTag="MER")
            else:
                (k_partitions, k_keeps) = self.HP.houghPartition(k_dat[:,0], l_dat)

        if(len(k_partitions) == 0):
            return None

        partitions = []
        k_sizes = [len(p) for p in k_partitions]

        #-----------------------
        # GRID
        if self.debugPlots >= 2:
            c_max = np_max(c_dat[:,2]/10) * 1.1
            k_max = np_max(k_dat[:,0]) * 1.1
            c_min = np_min(c_dat[:,2]/10) * 0.9
            k_min = np_min(k_dat[:,0]) * 0.9

            k_index_sort = np_argsort(k_dat[:,0])
            start = 0
            k_lines = []
            for k in range(len(k_sizes)-1):
                k_lines.append(k_dat[k_index_sort,0][k_sizes[k]+start])
                start += k_sizes[k]

            fig = plt.figure()

            orig_k_dat = self.PM.kmerPCs[rowIndices,0]
            orig_k2_dat = self.PM.kmerPCs[rowIndices,1]
            orig_c_dat = self.PM.transformedCP[rowIndices][:,2]/10
            orig_l_dat = np_sqrt(self.PM.contigLengths[rowIndices])

            ax = plt.subplot(221)
            plt.xlabel("PCA1")
            plt.ylabel("PCA2")

            from matplotlib.patches import Rectangle
            alpha = 0.35
            ax.scatter(orig_k_dat, orig_k2_dat, edgecolors='none', c=self.PM.contigGCs[rowIndices], cmap=self.PM.colorMapGC, vmin=0.0, vmax=1.0, s=orig_l_dat, zorder=10, alpha=alpha)
            XX = ax.get_xlim()
            YY = ax.get_ylim()
            ax.add_patch(Rectangle((XX[0], YY[0]),XX[1]-XX[0],YY[1]-YY[0],facecolor='#000000'))

            ax = plt.subplot(223)
            plt.title("%s contigs" % len(rowIndices))
            plt.xlabel("MER PARTS")
            plt.ylabel("COV PARTS")
            ax.scatter(orig_k_dat, orig_c_dat, edgecolors='none', c=self.PM.contigGCs[rowIndices], cmap=self.PM.colorMapGC, s=orig_l_dat, zorder=10, alpha=alpha)
            XX = ax.get_xlim()
            YY = ax.get_ylim()
            ax.add_patch(Rectangle((XX[0], YY[0]),XX[1]-XX[0],YY[1]-YY[0],facecolor='#000000'))

            lens = np_sqrt(self.PM.contigLengths[row_indices])

            ax = plt.subplot(222)
            plt.xlabel("PCA1")
            plt.ylabel("PCA2")
            ax.scatter(k_dat[:,0],
                       k_dat[:,1],
                       edgecolors='none',
                       c=self.PM.contigGCs[row_indices],
                       cmap=self.PM.colorMapGC,
                       vmin=0.0,
                       vmax=1.0,
                       s=lens,
                       zorder=10,
                       alpha=alpha)
            XX = ax.get_xlim()
            YY = ax.get_ylim()
            ax.add_patch(Rectangle((XX[0], YY[0]),XX[1]-XX[0],YY[1]-YY[0],facecolor='#000000'))

            ax = plt.subplot(224)
            plt.title("%s contigs" % len(row_indices))
            plt.xlabel("MER PARTS")
            plt.ylabel("COV PARTS")
            ax.scatter(k_dat[:,0],
                       c_dat[:,2]/10,
                       edgecolors='none',
                       c=self.PM.contigGCs[row_indices],
                       cmap=self.PM.colorMapGC,
                       vmin=0.0,
                       vmax=1.0,
                       s=lens,
                       zorder=10,
                       alpha=alpha)

            ax.set_xlim(k_min, k_max)
            ax.set_ylim(c_min, c_max)
            XX = ax.get_xlim()
            YY = ax.get_ylim()
            ax.add_patch(Rectangle((XX[0], YY[0]),XX[1]-XX[0],YY[1]-YY[0],facecolor='#000000'))

            k_line_cols = []
            for k in range(len(k_sizes)):
                if k == 0:
                    if k_keeps[k]:
                        k_line_cols = ['r-']
                    else:
                        k_line_cols = ['r--']
                elif k == len(k_sizes) - 1:
                    if k_keeps[k]:
                        k_line_cols[-1] = 'r-'
                    elif not k_keeps[k-1]:
                        k_line_cols[-1] = 'r--'
                else:
                    if k_keeps[k]:
                        k_line_cols[k-1] = 'r-'
                        k_line_cols.append('r-')
                    else:
                        k_line_cols.append('r--')

            for k in range(len(k_lines)):
                plt.plot([k_lines[k],k_lines[k]], [c_min, c_max], k_line_cols[k], zorder=11)

        pc = 0
        for k in range(len(k_sizes)):
            pc += 1
            if k_keeps[k]:
                k_part = k_partitions[k]
                part_bp = np_sum(l_dat[k_part])
                if self.BM.isGoodBin(part_bp, len(k_part), ms=5):

                    # select just the subset of covs for this kmer range
                    data = np_copy(c_dat[k_part])

                    # PCA the subset and cluster on the 1st component
                    Center(data,verbose=0)
                    p = PCA(data)
                    components = p.pc()
                    data = np_array([float(i) for i in components[:,0]])

                    l_data = np_copy(l_dat[k_part])

                    # The PCA may reverse the ordering. So we just check here quickly
                    if self.debugPlots >= 3:
                        (c_partitions, c_keeps) = self.HP.houghPartition(data,
                                                                         l_data,
                                                                         imgTag="COV_%s" % (pc-1),
                                                                         gData=self.PM.transformedCP[row_indices[k_part],2],
                                                                         gCut=10.)
                    else:
                        (c_partitions, c_keeps) = self.HP.houghPartition(data,
                                                                         l_data,
                                                                         gData=self.PM.transformedCP[row_indices[k_part],2],
                                                                         gCut=10.)


                    if self.debugPlots >= 2:
                        #-----
                        # GRID
                        c_sorted_data = np_copy(c_dat[k_part,2])/10.
                        c_sorted_data = c_sorted_data[np_argsort(c_sorted_data)]

                        start = 0
                        c_lines = []
                        if k_part[c_partitions[0]][0] > k_part[c_partitions[-1]][0]:
                            c_sizes = [len(p) for p in c_partitions][::-1]
                            cc_keeps = c_keeps[::-1]
                        else:
                            c_sizes = [len(p) for p in c_partitions]
                            cc_keeps = c_keeps

                        for c in range(len(c_sizes)-1):
                            c_lines.append(c_sorted_data[c_sizes[c]+start])
                            start += c_sizes[c]

                        c_line_cols = []
                        for c in range(len(c_sizes)):
                            if c == 0:
                                if cc_keeps[c]:
                                    c_line_cols = ['r-']
                                else:
                                    c_line_cols = ['r--']
                            elif c == len(c_sizes) - 1:
                                if cc_keeps[c]:
                                    c_line_cols[-1] = 'r-'
                                elif not cc_keeps[c-1]:
                                    c_line_cols[-1] = 'r--'
                            else:
                                if cc_keeps[c]:
                                    c_line_cols[c-1] = 'r-'
                                    c_line_cols.append('r-')
                                else:
                                    c_line_cols.append('r--')

                        if pc == 1:
                            k_line_min = k_min
                        else:
                            k_line_min = k_lines[pc-2]

                        if pc == len(k_partitions):
                            k_line_max = k_max
                        else:
                            k_line_max = k_lines[pc-1]

                        for c in range(len(c_lines)):
                            plt.plot([k_line_min,k_line_max], [c_lines[c], c_lines[c]], c_line_cols[c], zorder=11)

                    for c in range(len(c_partitions)):
                        if c_keeps[c]:
                            c_part = c_partitions[c]
                            partitions.append(np_array(k_part[c_part]))

        if self.debugPlots >= 2:
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
        disp_lens = np_array([])
        num_points = 0
        for row_index in rowIndices:
            num_points += 1
            disp_vals = np_append(disp_vals, self.PM.transformedCP[row_index])
            disp_lens = np_append(disp_lens, np_sqrt(self.PM.contigLengths[row_index]))

        # reshape
        disp_vals = np_reshape(disp_vals, (num_points, 3))

        ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors='k', c=self.PM.contigGCs[rowIndices], cmap=self.PM.colorMapGC, vmin=0.0, vmax=1.0, s=disp_lens, marker='.')

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
                                disp_cols = np_append(disp_cols, self.PM.colorMapGC(self.PM.contigGCs[row_index]))

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
        sc = ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors=disp_cols, c=disp_cols, cmap=cm, marker='.')
        sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect
        title = str.join(" ", ["Focus at: (",str(px), str(py), str(self.PM.scaleFactor - pz - 1),")\n",tag])
        plt.title(title)

        if(fileName != ""):
            fig.set_size_inches(6,6)
            plt.savefig(fileName,dpi=300)
        elif(show):
            plt.show()

        plt.close(fig)
        del fig

    def plotHeat(self, fileName = "", x=-1, y=-1):
        """Print the main heat maps

        Useful for debugging
        """
        fig = plt.figure()
        images = []
        ax = None
        if(self.numImgMaps == 1):
            ax = fig.add_subplot(121)
            ax.imshow(self.blurredMaps[0,:,:]**0.5)
            if x != -1:
                ax.scatter([y],
                           [x],
                           edgecolors='#ff69b4',
                           c='none',
                           s=5000.,
                           marker='.')
        else:
            ax = fig.add_subplot(231)
            images.append(ax.imshow(self.blurredMaps[0,:,:]**0.5))
            ax = fig.add_subplot(232)
            images.append(ax.imshow(self.blurredMaps[1,:,:]**0.5))
            ax = fig.add_subplot(233)
            images.append(ax.imshow(self.blurredMaps[2,:,:]**0.5))

        if(self.numImgMaps == 1):
            ax = fig.add_subplot(122)
            ax.imshow(self.imageMaps[0,:,:]**0.5)
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

    def houghPartition(self,
                       dAta,            # data to cluster with
                       lData,           # lengths for filling in about long contigs
                       imgTag=None,     # make debug plots
                       gData=None,      # gap data to check partitions with
                       gCut=0           # gap cutoff to stop merging of partitions
                       ):
        """Use the hough transform to find k clusters for some unknown value of K"""
        d_len_raw = int(len(dAta))
        if d_len_raw < 2:
            return np_array([[0]])

        if d_len_raw < 3:
            return np_array([[0,1]])

        #----------------------------------------------------------------------
        # prep the data
        #
        sorted_indices_raw = np_argsort(dAta)
        nUm_points = len(dAta)

        # fudge the data to make longer contigs have more say in the
        # diff line we'll be making. This way we may be able to avoid lumping
        # super long contigs in with the short riff raff by accident.
        scales_per = {}
        spread_data = []
        spread2real = {}
        j = 0

        # all points get at least one point, but long ones get more
        # let's say 1 point per 2000bp
        for i in range(len(dAta)):
            real_index = sorted_indices_raw[i]
            rep = int((lData[real_index] - 1.)/5000.) + 1
            scales_per[real_index] = rep
            if rep == 1:
                spread_data.append(dAta[real_index])
                spread2real[j] = real_index
                j += 1
            else:
                # rep >= 2
                if i == 0:
                    left_stop = dAta[real_index]
                else:
                    left_stop = (dAta[real_index] + dAta[sorted_indices_raw[i-1]])/2.

                if i == nUm_points-1:
                    right_stop = dAta[real_index]
                else:
                    right_stop = (dAta[real_index] + dAta[sorted_indices_raw[i+1]])/2.

                spread_jump = (right_stop - left_stop) / (rep + 1.)


                for ii in range(rep):
                    spread_data.append(left_stop + ((ii + 1) * spread_jump))
                    spread2real[j] = real_index
                    j += 1

        # Force the data to fill the space
        data = np_array(spread_data)
        data -= np_min(data)    # shift to 0 but DO NOT scale to 1
        d_len = len(data)
        scale = np_max(dAta) - np_min(dAta)

        # we want to know how much each value differs from it's neighbours
        back_diffs = [float(data[i] - data[i-1]) for i in range(1,d_len)]
        diffs = [back_diffs[0]]
        for i in range(len(back_diffs)-1):
            diffs.append((back_diffs[i] + back_diffs[i+1])/2)
        diffs.append(back_diffs[-1])
        diffs = np_array(diffs)**2  # square it! Makes things more betterrer

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

        # make it 2D
        t_data = np_array(zip(diffs, np_arange(d_len)))
        im_shape = (int(np_max(t_data, axis=0)[0]+1), d_len)

        #----------------------------------------------------------------------
        # Apply hough transformation and find the most prominent line
        #
        # rets should be an array of arrays of real indices
        # IE. indices into the array dAta
        rets = self.recursiveSelect(t_data,
                                    im_shape,
                                    spread2real,
                                    0,
                                    d_len,
                                    {},
                                    imgTag=imgTag)

        #----------------------------------------------------------------------
        # Squish things up
        #
        # build a flat data set similar to the gradiated data set
        real2spread = {}
        j = 0
        flat_data = []
        if len(rets) > 1:
            for i in range(len(dAta)):
                real_index = sorted_indices_raw[i]
                rep = scales_per[real_index]
                for k in range(rep):
                    flat_data.append(dAta[real_index])
                    try:
                        real2spread[real_index].append(j)
                    except KeyError:
                        real2spread[real_index] = [j]
                    j += 1
            data -= np_min(flat_data)
            back_diffs = [float(data[i] - data[i-1]) for i in range(1,d_len)]
            diffs = [back_diffs[0]]
            for i in range(len(back_diffs)-1):
                diffs.append((back_diffs[i] + back_diffs[i+1])/2)
            diffs.append(back_diffs[-1])
            diffs = np_array(diffs)**2
            for i in range(1, d_len):
                diffs[i] += diffs[i-1]
            diffs -= np_min(diffs)
            try:
                diffs /= np_max(diffs)
            except FloatingPointError:
                pass
            diffs *= len(diffs)

            # diffs is now the same size as the gradiated data sent through
            # to hough in level 0. We wish to find the gradients of the lines
            # returned by recursive partitioning
            gradients = []
            for ret in rets:
                sis = []
                for ii in ret:
                    for si in real2spread[ii]:
                        sis.append(diffs[si])
                l_sis = len(sis)
                if l_sis == 1:
                    gradients.append(-1)
                else:
                    sis = sorted(sis)
                    gradients.append((sis[-1] - sis[0])/l_sis)

            gradients = np_array(gradients)

            # get all the -1 gradients and make them equal to the larger
            # of their neighbours
            last = 0.
            for g in range(len(gradients)):
                if gradients[g] == -1:
                    h = g + 1
                    next = 0.
                    # find the next not -1 gradient
                    while(h < len(gradients)):
                        if gradients[h] != -1:
                            # use this one
                            next = gradients[h]
                            break
                        h += 1
                    gradients[g] = np_max([last, next])
                else:
                    last = gradients[g]

        else:
            gradients=np_array([0.])

        keeps = np_where(gradients >= 1, False, True)

        squished_rets = []
        squished_keeps = []
        last_squished = []

        if gData is not None:
            # measure the gaps between potential squishables
            last_max = -1
            gaps = []
            for i in range(len(rets)):
                tmp_gs = []
                for ii in rets[i]:
                    tmp_gs.append(gData[ii])    # collate the gData for this partition

                A = np_min(tmp_gs)              # find it's boundaries
                B = np_max(tmp_gs)
                if last_max > 0:
                    gaps.append(A - last_max)   # measure the gap
                last_max = B
            gap_info = [True] + [ggg < gCut for ggg in gaps]
        else:
            gap_info = [True] * len(keeps)

        for i in range(len(rets)):
            if keeps[i]:
                # check and see if we are allowed to squish
                if gap_info[i]:
                    # squish OK
                    for ii in rets[i]:
                        last_squished.append(ii)
                else:
                    # shouldn't squish
                    # push on the last fella
                    if len(last_squished) > 0:
                        squished_rets.append(np_array(last_squished))
                        squished_keeps.append(True)
                    # remake the last_squished holder
                    last_squished = []
                    for ii in rets[i]:
                        last_squished.append(ii)
            else:
                if len(last_squished) > 0:
                    squished_rets.append(np_array(last_squished))
                    last_squished = []
                    squished_keeps.append(True)

                # add the dud
                squished_rets.append(rets[i])
                squished_keeps.append(False)
        if len(last_squished) > 0:
            squished_rets.append(np_array(last_squished))
            squished_keeps.append(True)

        return (np_array(squished_rets), np_array(squished_keeps))

    def recursiveSelect(self,
                        tData,
                        imShape,
                        spread2real,
                        startRange,
                        endRange,
                        assigned,
                        level=0,
                        side="C",
                        imgTag=None):
        """Recursively select clusters from the data"""
        d_len = len(tData)
        (m, ret_point, accumulator) = self.houghTransform(tData.astype(float)[startRange:endRange,:], imShape)

        if m == np_inf:
            # this is a vertical line through ret_point
            start_p = [0, ret_point[1]]
            end_p = [imShape[0], ret_point[1]]
        else:
            # find out where the line crosses axes
            y_int = ret_point[0] - m * ret_point[1]
            if m == 0:
                # line is flat
                x_int = np_inf
            else:
                x_int = ret_point[1] - ret_point[0] / m

            if np_abs(y_int) < np_abs(x_int):
                # m line passes above (0,0)
                start_p = [y_int, 0]
                end_p = [imShape[1]*m + y_int, imShape[1]]
            else:
                # m line passes below (0,0)
                start_p = [0, x_int]
                end_p = [imShape[0], imShape[0]/m + x_int]

        # draw a nice thick line over the top of the data
        # found_line is a set of points
        found_line = self.points2Line(np_array([start_p,end_p]), imShape[1], imShape[0], 5)

        # make an image if we're that way inclined
        if imgTag is not None:
            # make a pretty picture
            fff = np_ones(imShape) * 255
            for p in found_line.keys():
                fff[p[0],p[1]] = 220
            for p in tData:
                fff[p[0],p[1]] = 0

            # scale so colors look sharper
            accumulator -= np_min(accumulator)
            accumulator /= np_max(accumulator)
            accumulator *= 255

            imsave("%d_%s_%s_%d.png" % (self.hc, imgTag, side, level), np_concatenate([accumulator,fff]))
            print "%d_%s_%s_%d.png" % (self.hc, imgTag, side, level)

        # see which points lie on the line
        # we need to protect against the data line crossing
        # in and out of the "found line"
        in_block = False
        block_starts = []
        block_lens = []
        for ii in np_arange(startRange, endRange):
            p = tData.astype('int')[ii]
            if tuple(p) in found_line:
                if not in_block:
                    in_block = True
                    block_starts.append(ii)
            else:
                if in_block:
                    in_block = False
                    block_lens.append(ii - block_starts[-1])

        if in_block:
            # finishing block
            block_lens.append(endRange - block_starts[-1])

        # check to see the line hit something
        if len(block_lens) == 0:
            tmp = {}
            for ii in np_arange(startRange, endRange):
                real_index = spread2real[ii]
                if real_index not in assigned:
                    tmp[real_index] = None
                    assigned[real_index] = None
            centre = np_array(tmp.keys())
            if len(centre) > 0:
                return np_array([centre])
            # nuffin
            return np_array([])

        # get the start and end indices in the longest block found
        longest_block = np_argmax(block_lens)
        spread_start = block_starts[longest_block]
        spread_end =  block_lens[longest_block] + spread_start  # 1 after the end of the block

        # select all the guys with their centres between the start and end
        # this is the end of the line for these guys so we fill centre with
        # "real" indices.
        tmp = {}
        for ii in np_arange(spread_start, spread_end):
            real_index = spread2real[ii]
            if real_index not in assigned:
                tmp[real_index] = None
                assigned[real_index] = None
        centre = np_array(tmp.keys())

        rets = []

        # recursive call for leftmost indices
        if (spread_start - startRange) > 0:
            if (spread_start - startRange) < 3:
                # end of the line for left expansion, give up "real" indices
                # select all the guys with their centres to the left of the start
                tmp = {}
                for ii in np_arange(startRange, spread_start):
                    real_index = spread2real[ii]
                    if real_index not in assigned:
                        tmp[real_index] = None
                        assigned[real_index] = None

                if len(tmp.keys()) > 0:
                    rets.append(np_array(tmp.keys()))

            else:
                # otherwise we keep working with ranges
                left_p = self.recursiveSelect(tData,
                                              imShape,
                                              spread2real,
                                              startRange,
                                              spread_start,
                                              assigned,
                                              level=level+1,
                                              side="%sL" %side,
                                              imgTag=imgTag)
                for L in left_p:
                    if len(L) > 0:
                        rets.append(L)

        # add the centre in
        if len(centre) > 0:
            rets.append(centre)

        # recursive call for rightmost indices
        if (endRange - spread_end) > 0:
            if (endRange - spread_end) < 3:
                # end of the line for left expansion, give up "real" indices
                # select all the guys with their centres right of the end
                tmp = {}
                for ii in np_arange(spread_end, endRange):
                    real_index = spread2real[ii]
                    if real_index not in assigned:
                        tmp[real_index] = None
                        assigned[real_index] = None

                if len(tmp.keys()) > 0:
                    rets.append(np_array(tmp.keys()))
            else:
                right_p = self.recursiveSelect(tData,
                                               imShape,
                                               spread2real,
                                               spread_end,
                                               endRange,
                                               assigned,
                                               level=level+1,
                                               side="%sR" %side,
                                               imgTag=imgTag)
                for R in right_p:
                    if len(R) > 0:
                        rets.append(R)
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

    def houghTransform(self, data, imShape):
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
        # first get a point our found line passes through
        Ps = []
        for p in range(len(data)):
            point = data[p]
            th = dth * min_col
            r = point[1]*np_cos(th) + point[0]*np_sin(th)
            iry = half_rows + int(r/dr)
            # this point does
            if iry == min_row:
                Ps.append(data[p])
        # take the average of all of em'
        ret_point = np_mean(Ps, axis=0)

        # get the gradient
        if theta != 0:
            m = -1. * np_cos(theta) / np_sin(theta)
        else:
            m = np_inf

        # rounding errors suck
        if np_allclose([m], [0.]):
            m = 0.

        return (m, ret_point, accumulator.reshape((rows, cols)))

###############################################################################
###############################################################################
###############################################################################
###############################################################################

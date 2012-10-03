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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import sys
import time

import colorsys
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

import numpy as np
import scipy.ndimage as ndi
from scipy.spatial.distance import cdist

# GroopM imports
import PCA
import dataManagers
import bin
import som
import torusMesh

np.seterr(all='raise')      

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ClusterEngine:
    """Top level interface for clustering contigs"""
    def __init__(self, dbFileName, plot=False, force=False, numImgMaps=1):
        # worker classes
        self.PM = dataManagers.ProfileManager(dbFileName) # store our data
        self.BM = dataManagers.BinManager(pm=self.PM)   # store our bins
    
        # heat maps
        self.numImgMaps = numImgMaps
        self.imageMaps = np.zeros((self.numImgMaps,self.PM.scaleFactor,self.PM.scaleFactor))
        self.blurredMaps = np.zeros((self.numImgMaps,self.PM.scaleFactor,self.PM.scaleFactor))
        
        # we need a way to reference from the imageMaps back onto the transformed data
        self.im2RowIndicies = {}  
        
        # When blurring the raw image maps I chose a radius to suit my data, you can vary this as you like
        self.blurRadius = 2
        self.span = 30                  # amount we can travel about when determining "hot spots"
        
        # misc
        self.minSize=10                 # Min number of contigs for a bin to be considered legit
        self.minVol=1000000             # Override on the min size, if we have this many BP
        self.forceWriting = force
        self.debugPlots = plot
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
            print "Overwriting database",self.PM.dbFileName
            self.PM.dataManager.nukeBins(self.PM.dbFileName)
        return True

#------------------------------------------------------------------------------
# BIN EXPANSION USING SOMS

    def expandBins(self, force=False, plot=False):
        """Expand the bins"""
        pass
    
#------------------------------------------------------------------------------
# CORE CONSTRUCTION AND MANAGEMENT
        
    def makeCores(self, coreCut, minSize, minVol):
        """Cluster the contigs to make bin cores"""
        # check that the user is OK with nuking stuff...
        if(not self.promptOnOverwrite()):
            return False

        self.minVol = minVol
        self.minSize = minSize

        # get some data
        t0 = time.time()
        self.PM.loadData(condition="length >= "+str(coreCut))
        t1 = time.time()
        print "    THIS: [",self.secondsToStr(t1-t0),"]\tTOTAL: [",self.secondsToStr(t1-t0),"]"
        
        # transform the data
        print "Apply data transformations"
        self.PM.transformCP()
        # plot the transformed space (if we've been asked to...)
        if(self.debugPlots):
            self.PM.renderTransCPData()
        t2 = time.time()
        print "    THIS: [",self.secondsToStr(t2-t1),"]\tTOTAL: [",self.secondsToStr(t2-t0),"]"
        
        # cluster and bin!
        print "Create cores"
        cum_contigs_used_good = self.initialiseCores()
        t3 = time.time()
        print "    THIS: [",self.secondsToStr(t3-t2),"]\tTOTAL: [",self.secondsToStr(t3-t0),"]"

        # Now save all the stuff to disk!
        print "Saving bins"
        self.BM.saveBins(doCores=True, saveBinStats=True)
        t4 = time.time()
        print "    THIS: [",self.secondsToStr(t4-t3),"]\tTOTAL: [",self.secondsToStr(t4-t0),"]"

    def initialiseCores(self):
        """Process contigs and form CORE bins"""
        num_below_cutoff = 0            # how many consecutive attempts have produced small bins
        breakout_point = 50             # how many will we allow before we stop this loop
        
        # First we need to find the centers of each blob.
        # We can make a heat map and look for hot spots
        self.populateImageMaps()
        sub_counter = 0
        print "    .... .... .... .... .... .... .... .... .... ...."
        print "%03d" % sub_counter,
        new_line_counter = 0
        num_bins = 0
        ss=0
        while(num_below_cutoff < breakout_point):
            #if(num_bins > 70):
            #    break
            sys.stdout.flush()
            # apply a gaussian blur to each image map to make hot spots
            # stand out more from the background 
            self.blurMaps()
    
            # now search for the "hottest" spots on the blurred map
            # and check for possible bin centroids
            ss += 200
            putative_clusters = self.findNewClusterCenters(ss=ss)
            if(putative_clusters is None):
                break
            else:
                bids_made = []
                partitions = putative_clusters[0]
                [max_blur_value, max_x, max_y] = putative_clusters[1]
                self.roundNumber += 1
                sub_round_number = 1
                for center_row_indicies in partitions:
                    total_BP = sum([self.PM.contigLengths[i] for i in center_row_indicies])
                    num_contigs = len(center_row_indicies)
                    #MM__print "Round: %d tBP: %d tC: %d" % (sub_round_number, total_BP, num_contigs)
                    if self.isGoodBin(total_BP, num_contigs, ms=5):   # Can we trust very small bins?.

                        # time to make a bin
                        bin = self.BM.makeNewBin(rowIndicies=center_row_indicies)
                        #MM__print "NEW:", total_BP, len(center_row_indicies)
                        # work out the distribution in points in this bin
                        bin.makeBinDist(self.PM.transformedCP, self.PM.averageCoverages, self.PM.kmerVals, self.PM.contigLengths)     

                        # Plot?
                        if(self.debugPlots):          
                            bin.plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, fileName="Image_"+str(self.imageCounter))
                            self.imageCounter += 1

                        # recruit more contigs
                        bin_size = bin.recruit(self.PM.transformedCP,
                                               self.PM.averageCoverages,
                                               self.PM.kmerVals,
                                               self.PM.contigLengths, 
                                               self.im2RowIndicies, 
                                               self.PM.binnedRowIndicies, 
                                               self.PM.restrictedRowIndicies
                                               )

                        if(self.debugPlots):
                            self.plotHeat("HM_%d.%d.png" % (self.roundNumber, sub_round_number), max=max_blur_value, x=max_x, y=max_y)
                            sub_round_number += 1

                        if(self.isGoodBin(self, bin.totalBP, bin_size)):
                            # Plot?
                            bids_made.append(bin.id)
                            num_bins += 1
                            if(True):#self.debugPlots):          
                                bin.plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, fileName="P_BIN_%d"%(bin.id))

                            # append this bins list of mapped rowIndicies to the main list
                            self.updatePostBin(bin)
                            num_below_cutoff = 0
                            print "%04d"%bin_size,
                        else:
                            # we just throw these indicies away for now
                            self.restrictRowIndicies(bin.rowIndicies)
                            self.BM.deleteBins([bin.id], force=True)
                            num_below_cutoff += 1
                            print str(bin_size).rjust(4,'X'),
        
                        # make the printing prettier
                        new_line_counter += 1
                        if(new_line_counter > 9):
                            new_line_counter = 0
                            sub_counter += 10
                            print "\n%03d" % sub_counter,
                    else:
                        # this partition was too small, restrict these guys we don't run across them again
                        self.restrictRowIndicies(center_row_indicies)
                
                # did we do anything?
                num_bids_made = len(bids_made)
                if(num_bids_made == 0):
                    # nuke the lot!
                    for row_indicies in partitions:
                        self.restrictRowIndicies(row_indicies)
                else:#if(False):
                    # now is as good a time as any to see if we can merge these guys
                    #MM__print "\n"
                    bids_consumed = {}
                    for i in range(num_bids_made):
                        base_bid = bids_made[i]
                        if(base_bid not in bids_consumed):
                            base_bin = self.BM.getBin(base_bid)
                            for j in range(i+1, num_bids_made):
                                query_bid = bids_made[j]
                                if(query_bid not in bids_consumed):
                                    #MM__print base_bid,query_bid,
                                    if(self.BM.shouldMerge(base_bin, self.BM.getBin(query_bid))):
                                        # merge!
                                        #MM__print "OK"
                                        self.BM.merge([base_bid,query_bid], auto=True, manual=False, newBid=False, saveBins=False, verbose=False, printInstructions=False)
                                        bids_consumed[query_bid] = True
                                        if(self.debugPlots):          
                                            base_bin.plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, fileName="C_BIN_%d"%(bin.id))
#MM__                                    else:
#MM__                                        print "NOPE"
        print "\n    .... .... .... .... .... .... .... .... .... ...."

        # condense bins
        self.BM.autoCondenseBins(2*self.span)

        # neaten up the bins
        #self.removeOutliersWrapper()
        
        num_binned = len(self.PM.binnedRowIndicies.keys())
        perc = "%.2f" % round((float(num_binned)/float(self.PM.numContigs))*100,2)
        print "\n   ",num_binned,"contigs are distributed across",len(self.BM.bins.keys()),"cores (",perc,"% )"

    def isGoodBin(self, totalBP, binSize, ms=0):
        """Does this bin meet my exacting requirements?"""
        if(ms == 0):
            ms = self.minSize               # let the user choose
        if(totalBP < self.minVol):          # less than the good volume
            if(binSize > ms):               # but has enough contigs
                return True
        else:                               # contains enough bp to pass regardless of number of contigs
            return True        
        return False
    
    def findNewClusterCenters(self, ss=0):
        """Find a putative cluster"""
        
        inRange = lambda x,l,u : x >= l and x < u

        # we work from the top view as this has the base clustering
        max_index = np.argmax(self.blurredMaps[0])
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
        working_block = np.zeros((span_len, span_len, self.PM.scaleFactor))
        
        # go through the entire column
        (x_lower, x_upper) = self.makeCoordRanges(max_x, start_span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, start_span)
        super_putative_row_indicies = []
        for p in self.im2RowIndicies:
            if inRange(p[0],x_lower,x_upper) and inRange(p[1],y_lower,y_upper):
                for row_index in self.im2RowIndicies[p]: 
                    # check that the point is real and that it has not yet been binned
                    if row_index not in self.PM.binnedRowIndicies and row_index not in self.PM.restrictedRowIndicies:
                        # this is an unassigned point. 
                        multiplier = np.log10(self.PM.contigLengths[row_index])
                        self.incrementAboutPoint3D(working_block, p[0]-x_lower, p[1]-y_lower, p[2],multiplier=multiplier)
                        super_putative_row_indicies.append(row_index)
    
        # blur and find the highest value
        bwb = ndi.gaussian_filter(working_block, 8)#self.blurRadius)
        densest_index = np.unravel_index(np.argmax(bwb), (np.shape(bwb)))
        max_x = densest_index[0] + x_lower
        max_y = densest_index[1] + y_lower
        max_z = densest_index[2]

                    
        # now get the basic color of this dense point
        putative_center_row_indicies = []

        (x_lower, x_upper) = self.makeCoordRanges(max_x, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, self.span)
        (z_lower, z_upper) = self.makeCoordRanges(max_z, 2*self.span)

        for row_index in super_putative_row_indicies:
            p = np.around(self.PM.transformedCP[row_index])
            if inRange(p[0],x_lower,x_upper) and inRange(p[1],y_lower,y_upper) and inRange(p[2],z_lower,z_upper):  
                # we are within the range!
                putative_center_row_indicies.append(row_index)
         
        # make sure we have something to go on here
        if(np.size(putative_center_row_indicies) == 0):
            # it's all over!
            return None
        elif(np.size(putative_center_row_indicies) == 1):
            # get out of here but keep trying
            # the callig function should restrict these indicies
            return [[np.array(putative_center_row_indicies)], ret_values]
        else:
            total_BP = sum([self.PM.contigLengths[i] for i in putative_center_row_indicies])
            if not self.isGoodBin(total_BP, len(putative_center_row_indicies), ms=5):   # Can we trust very small bins?.
                # get out of here but keep trying
                # the calling function should restrict these indicies
                return [[np.array(putative_center_row_indicies)], ret_values]
            else:
                # we've got a few good guys here, partition them up!
                # shift these guys around a bit
                center_k_vals = np.array([self.PM.kmerVals[i] for i in putative_center_row_indicies])
                centre_transformed_CP = np.reshape(np.array([self.PM.transformedCP[i] for i in putative_center_row_indicies]),(len(putative_center_row_indicies),3))
                c_cols = np.array([self.PM.contigColours[i] for i in putative_center_row_indicies])
                #MM__print "PRE SHIFT: %d" % len(center_k_vals)
                shifted_k_vals = self.shiftVals(center_k_vals, centre_transformed_CP, c_cols, ss=ss)
                if(len(shifted_k_vals) == 0):
                    return None
                
                sk_partitions = self.partitionVals(shifted_k_vals)
                #MM__print "# PARTS: %d" % len(sk_partitions)
                
                if(len(sk_partitions) == 0):
                    return None
                else:
                    ret_ps = [np.array([putative_center_row_indicies[i] for i in p]) for p in sk_partitions]
                    return [ret_ps, ret_values]

    def shiftVals(self, kVals, tCP, cols, ss=0):
        """Use principles from self organisation to auto cluster these points
        
        assumes len(kVals) == len(cVals)
        """
        v_size = len(kVals)
        fig = plt.figure()
        if(False):
            plt.clf()
            ax = plt.subplot(111, projection='3d')
            ax.scatter(tCP[:,0], tCP[:,1], tCP[:,2], edgecolors=cols, c=cols, marker='.')
            ax.azim = 45
            ax.elev = 45
            filename="original"
            try:
                fig.set_size_inches(4,4)
                plt.savefig(filename,dpi=100)
            except:
                print "Error showing image", sys.exc_info()[0]
                raise
            
###############################################################################
###############################################################################
###############################################################################
  #
  # BEWARE! HERE BE DRAGONS!
  #
  # MAGIC NUMBERS ABOUND! WHY ARE THEY HERE? HOW WERE THEY CHOSEN?
  #  
  # ONE DAY I WILL EXPLAIN, UNTIL THEN REST ASSURED THAT:
  #
  # A. THEY (KINDA) WORK WELL ENOUGH
  # B. IT TAKES ABOUT THREE DAYS OF LONG HARD TRIAL AND ERROR TO DERIVE THEM
  # C. THEY BASICALLY ENCODE MY INTUITION 
  #
###############################################################################
###############################################################################
###############################################################################
  #
  # first we normalise the entire space. We would like a cube with a diagonal of 128px
  #
        tCP -= np.min(tCP, axis=0)
        try:
            tCP /= (np.max(tCP, axis=0)/70) # 70 is the side of said cube
        except FloatingPointError:
            # must be a zero in there somwhere...
            return []
  #
  # want to add a long range repulsive force which tries to push all the points
  # apart, should be quadratically decreasing, think gravity!
  # This force, if left to its owen devices will push the points
  # apart until they fill the spherical space.
  # [0 3], [128, 0]
  #
        global_repulsive_force = lambda x : 3 - 0.000183105*np.square(x)
  #
  # Add to this another repulsive force which is proportional to the
  # kmer distance. In combination with the global_repulsive_force it
  # will force similar kmers to group together, by forcing dissimilar kmers
  # to move apart
  # [0 -1], [0.05 0], [1 1]
  #
        kmer_force = lambda x : -18.9473684211*np.square(x) + 20.9473684211*x -1.0
  #
  # as kmers won;t change, we can just calculate this force once.
  #
        k_dm = np.reshape(np.tile(kVals, v_size),(v_size,v_size))
        k_dm_full = k_dm - k_dm.T
        k_dm = np.abs(k_dm_full)
        kf = kmer_force(k_dm)
  #
  # [0 -2], [40 0], [128 1]
  # cov_force = lambda x : -0.0003018466*np.square(x) + 0.0620738636*x + -2.0
  # [0 -2], [25 0], [128 2]
  #
        cov_force = lambda x : -0.0004733010*np.square(x) + 0.0918325243*x -2.0
        tCP_dm = cdist(tCP,tCP,'euclidean')
        cf = cov_force(tCP_dm)
   #
###############################################################################
###############################################################################
###############################################################################
        
        counter = ss
        true_counter = 0
        total_movement = 0.0
        movement_grad = 0.0
        sg = 0
        while(True):
            if(False):
                if(true_counter % 4 == 0):
                    plt.clf()
                    ax = plt.subplot(111, projection='3d')
                    ax.scatter(tCP[:,0], tCP[:,1], tCP[:,2], edgecolors=cols, c=cols, marker='.')
                    ax.azim = 45
                    ax.elev = 45
                    plt.title("Total movement %0.4f \nGradient: %0.4f" % (total_movement, movement_grad))
                    filename="squeeze_%d.png" % counter
                    try:
                        fig.set_size_inches(4,4)
                        plt.savefig(filename,dpi=100)
                    except:
                        print "Error showing image", sys.exc_info()[0]
                        raise
                               
            true_counter += 1
            counter += 1
            
            # get the two lists of values as points and work out how far apart the points are
            tCP_dm = cdist(tCP,tCP,'euclidean')
            
            #k_dm = np.reshape(np.tile(kVals, v_size),(v_size,v_size))
            #k_dm = np.abs(k_dm - k_dm.T)

            # get the inverse of the distance matrix so we can normalise
            #tCP_idm = np.array(tCP_dm, copy=True)
            # no dividing by 0
            #np.fill_diagonal(tCP_idm, 1.0)
            tCP_idm = np.where(tCP_dm > 0, tCP_dm, 1.0)
            # now work out the residual forces
            # work out the force between points
            # we will need to normalise for distance below but we'll
            # just do it now...
            forces = (global_repulsive_force(tCP_dm) + kf + cf)/tCP_idm
            
            #forces = np.where(forces > 0, forces, 0)
            #forces *= (1-k_dm)
            
            # get a list of vectors pointing from one vector to the
            # other, we should normalise this but we'll be multiplying 
            # it  by the forces vector above so...
            # THIS IS UNREADABLE, BUT IT'S THE FASTEST I COULD MAKE IT GO
            pointers = np.reshape(np.reshape(np.repeat(tCP.ravel(),v_size), (3*v_size,v_size)).T - np.tile(tCP, v_size), (v_size,v_size,3))
            
            # now we can calculate the effect on each point being made by each other point
            movement = np.mean(pointers * forces[:,:,np.newaxis], axis=0)     
            movement_grad = total_movement 
            total_movement = np.sum(np.sqrt(np.sum(np.square(movement),axis=1)))
            movement_grad = np.abs(movement_grad - total_movement)
            tCP += movement 
            tCP -= np.min(tCP, axis=0)
            tCP /= (np.max(tCP, axis=0)/70)

            # break out if we stop condensing
            if(sg == 0):
                sg = movement_grad
            if(movement_grad < 0.1 and true_counter > 30):
                break
        #MM__print "Start: %0.4f, End: %0.4F, Rounds: %d" % (sg, movement_grad, true_counter)
        
        # Now do a quick few rounds with just kmers, to pull things a little tighter...
        # [0 -3], [0.15 0], [1 0.5] make this guy stronger
        kmer_force = lambda x : -19.4117647059*np.square(x) + 22.9117647059*x -3.0
        kf = kmer_force(k_dm)   
        # [0 2], [128, 0] # make this guy weaker
        global_repulsive_force = lambda x : 2 - 0.00012207*np.square(x)
        for l in range(15):
            if(False):
                plt.clf()
                ax = plt.subplot(111, projection='3d')
                ax.scatter(tCP[:,0], tCP[:,1], tCP[:,2], edgecolors=cols, c=cols, marker='.')
                ax.azim = 45
                ax.elev = 45
                plt.title("Total movement %0.4f \nGradient: %0.4f" % (total_movement, movement_grad))
                filename="squeeze_%d_FINAL.png" % counter
                try:
                    fig.set_size_inches(4,4)
                    plt.savefig(filename,dpi=100)
                except:
                    print "Error showing image", sys.exc_info()[0]
                    raise
                               
            counter += 1
            tCP_dm = cdist(tCP,tCP,'euclidean')
            tCP_idm = np.where(tCP_dm > 0, tCP_dm, 1.0)
            # no coverage force!
            forces = (global_repulsive_force(tCP_dm) + kf)/tCP_idm
            pointers = np.reshape(np.reshape(np.repeat(tCP.ravel(),v_size), (3*v_size,v_size)).T - np.tile(tCP, v_size), (v_size,v_size,3))
            movement = np.mean(pointers * forces[:,:, np.newaxis], axis=0)     
            movement_grad = total_movement 
            total_movement = np.sum(np.sqrt(np.sum(np.square(movement),axis=1)))
            movement_grad = np.abs(movement_grad - total_movement)
            tCP += movement 
        #MM__print "R2: %0.4F" % movement_grad
        
        # Now fade the values together, this will make kmer selection easier
        k_blur_factor = lambda k,c : 1 - np.square(k)/0.0009 - np.square(c)/25
        tCP_dm = cdist(tCP,tCP,'euclidean')
        
        blurs = k_blur_factor(k_dm, tCP_dm)
        # no negative numbers and no self interference!
        blurs = np.where(blurs > 0, blurs, 0.0)
        np.fill_diagonal(blurs,0.0)

        min_kVals = np.min(kVals)
        max_kVals = np.max(kVals)
        k_range = max_kVals - min_kVals
        min_kVals += min_kVals*k_range/5 
        max_kVals -= max_kVals*k_range/5
        
        # get the direction the colours should move in
        effect = np.sum(blurs*k_dm_full, axis=1)
        kVals += effect
        
        # renorm within the same zone
        kVals -= np.min(kVals)                  # shift to 0
        kVals *= ((max_kVals - min_kVals) / np.max(kVals))
        kVals += min_kVals

        if(False):
            # remake the colours
            S = 1
            V = 1
            cols = np.array([colorsys.hsv_to_rgb(val, S, V) for val in kVals])
            
            plt.clf()
            ax = plt.subplot(111, projection='3d')
            ax.scatter(tCP[:,0], tCP[:,1], tCP[:,2], edgecolors=cols, c=cols, marker='.')
            ax.azim = 45
            ax.elev = 45
            plt.title("Total movement %0.4f \nGradient: %0.4f" % (total_movement, movement_grad))
            filename="squeeze_%d_BLUR.png" % counter
            try:
                fig.set_size_inches(4,4)
                plt.savefig(filename,dpi=100)
            except:
                print "Error showing image", sys.exc_info()[0]
                raise

        del fig 
        return kVals
            
    def expandSelection(self, startIndex, vals, stdevCutoff=0.05, maxSpread=0.1):
        """Expand a selection left and right from a staring index in a list of values
        
        Keep expanding unless the stdev of the values goes above the cutoff
        Return a list of indices into the original list
        """
        ret_list = [startIndex]   # this is what we will give back
        start_val = vals[startIndex]
        value_store = [start_val]
        
        sorted_indicies = np.argsort(vals)
        max_index = len(vals)
        
        # set the upper and lower to point to the position
        # where the start resides 
        lower_index = 0
        upper_index = 0
        for i in range(max_index):
            if(sorted_indicies[i] == startIndex):
                break
            lower_index += 1
            upper_index += 1
        do_lower = True
        do_upper = True
        max_index -= 1
        
        while(do_lower or do_upper):
            if(do_lower):
                do_lower = False
                if(lower_index > 0):
                    try_val = vals[sorted_indicies[lower_index - 1]]
                    if(np.abs(try_val - start_val) < maxSpread):
                        try_array = value_store + [try_val]
                        if(np.std(try_array) < stdevCutoff):
                            value_store = try_array
                            lower_index -= 1
                            ret_list.append(sorted_indicies[lower_index])
                            do_lower = True
            if(do_upper):
                do_upper = False
                if(upper_index < max_index):
                    try_val = vals[sorted_indicies[upper_index + 1]]
                    if(np.abs(try_val - start_val) < maxSpread):
                        try_array = value_store + [try_val]
                        if(np.std(try_array) < stdevCutoff):
                            value_store = try_array
                            upper_index += 1
                            ret_list.append(sorted_indicies[upper_index])
                            do_upper = True
        return sorted(ret_list)

    def partitionVals(self, vals, stdevCutoff=0.04, maxSpread=0.15):
        """Work out where shifts in kmer/coverage vals happen"""
        partitions = []
        working_list = list(vals)
        fix_dict = dict(zip(range(len(working_list)),range(len(working_list))))
        while(len(working_list) > 2):
            cf = CenterFinder()
            c_index = cf.findArrayCenter(working_list)
            expanded_indicies = self.expandSelection(c_index, working_list, stdevCutoff=stdevCutoff, maxSpread=maxSpread)
            # fix any munges from previous deletes
            morphed_indicies = [fix_dict[i] for i in expanded_indicies]
            partitions.append(morphed_indicies)
            # shunt the indicies to remove down!
            shunted_indicies = []
            for offset, index in enumerate(expanded_indicies):
                shunted_indicies.append(index - offset)

            #print "FD:", fix_dict 
            #print "EI:", expanded_indicies
            #print "MI:", morphed_indicies
            #print "SI:", shunted_indicies
            
            # make an updated working list and fix the fix dict
            nwl = []
            nfd = {}
            shifter = 0
            for i in range(len(working_list) - len(shunted_indicies)):
                #print "================="
                if(len(shunted_indicies) > 0):
                    #print i, shunted_indicies[0], shifter
                    if(i >= shunted_indicies[0]):
                        tmp = shunted_indicies.pop(0)
                        shifter += 1
                        # consume any and all conseqs
                        while(len(shunted_indicies) > 0):
                            if(shunted_indicies[0] == tmp):
                                shunted_indicies.pop(0)
                                shifter += 1
                            else:
                                break
                #else:
                #    print i, "_", shifter

                nfd[i] = fix_dict[i + shifter]
                nwl.append(working_list[i + shifter])

                #print nfd
                #print nwl
                
            fix_dict = nfd
            working_list = nwl
            
        if(len(working_list) > 0):
            partitions.append(fix_dict.values())       
        return partitions

#------------------------------------------------------------------------------
# CORE MANAGEMENT 

    def condenseCores(self, auto=False):
        """Itterative wrapper for the BinManager method"""
        condensing_round = 0
        num_cores_condensed = 0
        while True: # do while loop anyone?
            condensing_round += 1
            (num_cores_condensed,continue_merge) = self.BM.condenseBins(verbose=True,
                                                                        auto=auto      
                                                                       )
            if(num_cores_condensed == 0):
                break
            else:
                print "    Core condensing round:", condensing_round, "Incorporated", num_cores_condensed, "cores into larger cores"
        
        num_binned = len(self.PM.binnedRowIndicies.keys())
        perc = "%.2f" % round((float(num_binned)/float(self.PM.numContigs))*100,2)
        print "   ",num_binned,"contigs are distributed across",len(self.BM.bins.keys()),"cores (",perc,"% )"
            
        return 

    def removeOutliersWrapper(self, mode="kmer"):
        """remove the outliers for all bins"""
        print "    Removing outliers"
        for bid in self.BM.bins:
            self.removeOutliers(bid, mode=mode)

    def removeOutliers(self, bid, fixBinnedRI=True, mode="kmer"):
        """remove outliers for a single bin"""
        dead_row_indicies = self.BM.bins[bid].findOutliers(self.PM.transformedCP, self.PM.kmerVals, mode=mode)
        if(len(dead_row_indicies)>0):
            if(fixBinnedRI):
                for row_index in dead_row_indicies:
                    self.setRowIndexUnassigned(row_index)
            self.BM.bins[bid].purge(dead_row_indicies,
                                    self.PM.transformedCP,
                                    self.PM.averageCoverages,
                                    self.PM.kmerVals,
                                    self.PM.contigLengths,
                                    self.PM.kmerVals)
        
#------------------------------------------------------------------------------
# DATA MAP MANAGEMENT 

    def populateImageMaps(self):
        """Load the transformed data into the main image maps"""
        # reset these guys... JIC
        self.imageMaps = np.zeros((self.numImgMaps,self.PM.scaleFactor,self.PM.scaleFactor))
        self.im2RowIndicies = {}
        
        # add to the grid wherever we find a contig
        row_index = -1
        for point in np.around(self.PM.transformedCP):
            row_index += 1
            # can only bin things once!
            if row_index not in self.PM.binnedRowIndicies and row_index not in self.PM.restrictedRowIndicies:
                # add to the row_index dict so we can relate the 
                # map back to individual points later
                p = tuple(point)
                if p in self.im2RowIndicies:
                    self.im2RowIndicies[p].append(row_index)
                else:
                    self.im2RowIndicies[p] = [row_index]
                
                # now increment in the grid
                # for each point we encounter we incrmement
                # it's position + the positions to each side
                # and touching each corner
                self.incrementViaRowIndex(row_index, p)

    def incrementViaRowIndex(self, rowIndex, point=None):
        """Wrapper to increment about point"""
        if(point is None):
            point = tuple(np.around(self.PM.transformedCP[rowIndex]))
        #px = point[0]
        #py = point[1]
        #pz = point[2]
        multiplier = np.log10(self.PM.contigLengths[rowIndex])
        self.incrementAboutPoint(0, point[0], point[1], multiplier=multiplier)
        if(self.numImgMaps > 1):
            self.incrementAboutPoint(1, self.PM.scaleFactor - point[2] - 1, point[1], multiplier=multiplier)
            self.incrementAboutPoint(2, self.PM.scaleFactor - point[2] - 1, self.PM.scaleFactor - point[0] - 1, multiplier=multiplier)

    def decrementViaRowIndex(self, rowIndex, point=None):
        """Wrapper to decrement about point"""
        if(point is None):
            point = tuple(np.around(self.PM.transformedCP[rowIndex]))
        #px = point[0]
        #py = point[1]
        #pz = point[2]
        multiplier = np.log10(self.PM.contigLengths[rowIndex])
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
        if map[px][py] < np.finfo(float).eps:
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
        shape = np.shape(workingBlock)
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
        self.blurredMaps = np.zeros((self.numImgMaps,self.PM.scaleFactor,self.PM.scaleFactor))
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
        for row_index in bin.rowIndicies:
            self.setRowIndexAssigned(row_index)
            
    def setRowIndexAssigned(self, rowIndex):
        """fix the data structures to indicate that rowIndex belongs to a bin
        
        Use only during initial core creation
        """        
        if(rowIndex not in self.PM.restrictedRowIndicies and rowIndex not in self.PM.binnedRowIndicies):
            self.PM.binnedRowIndicies[rowIndex] = True
            # now update the image map, decrement
            self.decrementViaRowIndex(rowIndex)

    def setRowIndexUnassigned(self, rowIndex):
        """fix the data structures to indicate that rowIndex no longer belongs to a bin
        
        Use only during initial core creation
        """
        if(rowIndex in self.PM.restrictedRowIndicies and rowIndex not in self.PM.binnedRowIndicies):
            del self.PM.binnedRowIndicies[rowIndex]
            # now update the image map, increment
            self.incrementViaRowIndex(rowIndex)

    def restrictRowIndicies(self, indicies):
        """Add these indicies to the restricted list"""
        for row_index in indicies:
            # check that it's not binned or already restricted
            if(row_index not in self.PM.restrictedRowIndicies and row_index not in self.PM.binnedRowIndicies):
                self.PM.restrictedRowIndicies[row_index] = True
                # now update the image map, decrement
                self.decrementViaRowIndex(row_index)
        
#------------------------------------------------------------------------------
# MISC 

    def secondsToStr(self, t):
        rediv = lambda ll,b : list(divmod(ll[0],b)) + ll[1:]
        return "%d:%02d:%02d.%03d" % tuple(reduce(rediv,[[t*1000,],1000,60,60]))
    
#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def plotRegion(self, px, py, pz, fileName="", tag="", column=False):
        """Plot the region surrounding a point """
        disp_vals = np.array([])
        disp_cols = np.array([])
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
                    if((x,y,realz) in self.im2RowIndicies):
                        for row_index in self.im2RowIndicies[(x,y,realz)]:
                            if row_index not in self.PM.binnedRowIndicies and row_index not in self.PM.restrictedRowIndicies:
                                num_points += 1
                                disp_vals = np.append(disp_vals, self.PM.transformedCP[row_index])
                                disp_cols = np.append(disp_cols, self.PM.contigColours[row_index])
        
        # make a black mark at the max values
        small_span = self.span/2
        (x_lower, x_upper) = self.makeCoordRanges(px, small_span)
        (y_lower, y_upper) = self.makeCoordRanges(py, small_span)
        (z_lower, z_upper) = self.makeCoordRanges(pz, small_span)
        for z in range(z_lower, z_upper):
            realz = self.PM.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.im2RowIndicies):
                        for row_index in self.im2RowIndicies[(x,y,realz)]:
                            if row_index not in self.PM.binnedRowIndicies and row_index not in self.PM.restrictedRowIndicies:
                                num_points += 1
                                disp_vals = np.append(disp_vals, self.PM.transformedCP[row_index])
                                disp_cols = np.append(disp_cols, colorsys.hsv_to_rgb(0,0,0))
        # reshape
        disp_vals = np.reshape(disp_vals, (num_points, 3))
        disp_cols = np.reshape(disp_cols, (num_points, 3))
        
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

class CenterFinder:
    """When a plain old mean won't cut it

    Uses a bouncing ball algorithm. Imagine walking along a "path",
    (through the array) hitting a ball into the air each time you
    come across a value. Gravity is bringing the ball down. If we plot
    the height of the ball vs array index then the highest the ball
    reaches is the index in the center of the densest part of the array 
    """
    def __init__(self): pass
    
    def findArrayCenter(self, vals):
        """Find the center of the numpy array vals, return the index of the center"""
        # parameters
        current_val_max = -1
        delta = 0
        bounce_amount = 0.1
        height = 0
        last_val= 0

        working = np.array([])
        final_index = -1
        
        # sort and normalise between 0 -> 1
        sorted_indicies = np.argsort(vals)
        vals_sorted = [vals[i] for i in sorted_indicies]
        vals_sorted -= vals_sorted[0]
        if(vals_sorted[-1] != 0):
            vals_sorted /= vals_sorted[-1]        

        # run through in one direction
        for val in vals_sorted:
            # calculate delta
            delta = val - last_val
            # reduce the current value according to the delta value
            height = self.reduceViaDelta(height, bounce_amount, delta)
            # bounce the ball up
            height += bounce_amount
            
            # store the height
            working = np.append(working, height)
            final_index += 1

            # save the last val            
            last_val = val

        current_val_max = -1
        height = 0
        last_val = 0
        
        # run through in the reverse direction
        vals_sorted = vals_sorted[::-1]
        for val in vals_sorted:
            if last_val == 0:
                delta = 0
            else:
                delta = last_val - val
            height = self.reduceViaDelta(height, bounce_amount, delta)
            height += bounce_amount
            # add to the old heights
            working[final_index] += height
            final_index -= 1
            last_val = val

        # find the original index!
        return sorted_indicies[np.argmax(working)]
    
    def reduceViaDelta(self, height, bounce_amount, delta):
        """Reduce the height of the 'ball'"""
        perc = (delta / bounce_amount)**0.5
        if(perc > 1):
            #print height, delta, 1, " H: ", 0
            return 0
        #print height, delta, (1-perc), " H: ", (height * (1-perc)) 
        return height * (1-perc)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

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
        """Expand bins using SOMs"""
        SM = dataManagers.SOMManager(self.BM, load=True)
        SM.regionalise(force=True, save=True)
        SM.findRegionNeighbours(merge=False, printMergers=True)
        SM.validateRegions()
        return
        SM = dataManagers.SOMManager(self.BM, load=True)
        SM.regionalise(force=True, save=True)
        SM.findRegionNeighbours(merge=True, printMergers=True)
        SM.assignmentMerge()
        SM.validateRegions()
        return
    
#------------------------------------------------------------------------------
# CORE CONSTRUCTION AND MANAGEMENT
        
    def makeCores(self, coreCut, minSize, minVol):
        """Cluster the contigs to make bin cores"""
        # check that the user is OK with nuking stuff...
        if(not self.promptOnOverwrite()):
            return False

        # get some data
        t0 = time.time()
        self.PM.loadData(condition="length >= "+str(coreCut))
        t1 = time.time()
        print "    THIS: [",self.secondsToStr(t1-t0),"]\tTOTAL: [",self.secondsToStr(t1-t0),"]"
        
        # transform the data
        print "Apply data transformations"
        self.PM.transformCP()
        # plot the transformed space (if we've been asked to...)
#        if(self.debugPlots):
#            self.PM.renderTransCPData()
        t2 = time.time()
        print "    THIS: [",self.secondsToStr(t2-t1),"]\tTOTAL: [",self.secondsToStr(t2-t0),"]"
        
        # cluster and bin!
        print "Create cores"
        cum_contigs_used_good = self.initialiseCores(minVol)
        t3 = time.time()
        print "    THIS: [",self.secondsToStr(t3-t2),"]\tTOTAL: [",self.secondsToStr(t3-t0),"]"

        # Now save all the stuff to disk!
        print "Saving bins"
        self.BM.saveBins(doCores=True, saveBinStats=True)
        t4 = time.time()
        print "    THIS: [",self.secondsToStr(t4-t3),"]\tTOTAL: [",self.secondsToStr(t4-t0),"]"

    def initialiseCores(self, minVol):
        """Process contigs and form CORE bins"""
        small_bin_cutoff = 10           # anything with less than this many contigs is not a real bin
        num_below_cutoff = 0            # how many consecutive attempts have produced small bins
        breakout_point = 50             # how many will we allow before we stop this loop
        
        # First we need to find the centers of each blob.
        # We can make a heat map and look for hot spots
        self.populateImageMaps()
        sub_counter = 0
        print "    .... .... .... .... .... .... .... .... .... ...."
        print "%03d" % sub_counter,
        new_line_counter = 0
        while(num_below_cutoff < breakout_point):
            #if(self.numBins > 3):
            #    break
            sys.stdout.flush()
            # apply a gaussian blur to each image map to make hot spots
            # stand out more from the background 
            self.blurMaps()
    
            # now search for the "hottest" spots on the blurred map
            # this is a new bin centroid
            (center_row_indicies, (max_blur_value, max_x, max_y)) = self.findNewClusterCenter()
            if(np.size(center_row_indicies) == 0):
                break
            else:
                # make sure we got something
                self.roundNumber += 1
                
                # time to make a bin
                bin = self.BM.makeNewBin(rowIndicies=center_row_indicies)
                
                # work out the distribution in points in this bin
                bin.makeBinDist(self.PM.transformedCP, self.PM.kmerSigs)     
                
                # Plot?
                if(self.debugPlots):          
                    bin.plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, fileName="Image_"+str(self.imageCounter))
                    self.imageCounter += 1
        
                # neaten up the bin
                #self.removeOutliers(bin.id, fixBinnedRI=False)
                
                # recruit more contigs
                bin_size = bin.recruit(self.PM.transformedCP, self.PM.kmerSigs, self.im2RowIndicies, self.PM.binnedRowIndicies)

                is_good_bin = False
                if(bin.calcTotalSize(self.PM.contigLengths) < minVol):    # less than the good volume
                    if(bin_size > small_bin_cutoff):                      # but has enough contigs
                        is_good_bin = True
                else:                                                     # contains enough bp to pass regardless of number of contigs
                    is_good_bin = True

                if(self.debugPlots):
                    self.plotHeat("HM_"+str(self.roundNumber)+".png", max=max_blur_value, x=max_x, y=max_y)
                    
                if(is_good_bin):
                    # Plot?
                    if(self.debugPlots):          
                        bin.plotBin(self.PM.transformedCP, self.PM.contigColours, self.PM.kmerVals, fileName="P_BIN_"+str(self.imageCounter))
                        self.imageCounter += 1

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
        
        print "\n    .... .... .... .... .... .... .... .... .... ...."

        # neaten up the bins
        self.removeOutliersWrapper()
        
        # remove possible chmimeras        
        self.BM.removeChimeras()
        
        num_binned = len(self.PM.binnedRowIndicies.keys())
        perc = "%.2f" % round((float(num_binned)/float(self.PM.numContigs))*100,2)
        print "\n   ",num_binned,"contigs are distributed across",len(self.BM.bins.keys()),"cores (",perc,"% )"

    def removeOutliersWrapper(self, mode="kmer"):
        """remove the outliers for all bins"""
        print "    Removing outliers"
        for bid in self.BM.bins:
            self.removeOutliers(bid, mode=mode)

    def removeOutliers(self, bid, fixBinnedRI=True, mode="kmer"):
        """remove outliers for a single bin"""
        dead_row_indicies = self.BM.bins[bid].findOutliers(self.PM.transformedCP, self.PM.kmerSigs, mode=mode)
        if(len(dead_row_indicies)>0):
            if(fixBinnedRI):
                for row_index in dead_row_indicies:
                    self.setRowIndexUnassigned(row_index)
            self.BM.bins[bid].purge(dead_row_indicies,
                                    self.PM.transformedCP,
                                    self.PM.kmerSigs,
                                    self.PM.contigLengths,
                                    self.PM.kmerVals)

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

        return
        # there's still a lot of background signal to remove
        # we wish to remove 90% of the data, this will leave just the really hot spots
        # Make a histogram of the data (use the top face)
        [vals,points] = np.histogram(np.reshape(self.blurredMaps[0,:,:], (self.PM.scaleFactor, self.PM.scaleFactor,1)), 50)
        total = np.sum(vals)*0.80
        lop_index = 1       # where we lop off the low values
        for val in vals:
            total -= val
            if total <= 0:
                break
            lop_index += 1
        lop_val = points[lop_index]

        # remove these low values and down normalise so that the largest value is equal to exactly 1
        for i in range (self.numImgMaps): # top, front and side
            self.blurredMaps[i,:,:] = np.where(self.blurredMaps[i,:,:] >= lop_val, self.blurredMaps[i,:,:], 0)/lop_val

    def makeCoordRanges(self, pos, span):
        """Make search ranges which won't go out of bounds"""
        lower = pos-span
        upper = pos+span+1
        if(lower < 0):
            lower = 0
        if(upper > self.PM.scaleFactor):
            upper = self.PM.scaleFactor
        return (lower, upper)

    def findNewClusterCenter(self):
        """Find a putative cluster"""
        # we work from the top view as this has the base clustering
        max_index = np.argmax(self.blurredMaps[0])
        max_value = self.blurredMaps[0].ravel()[max_index]

        max_x = int(max_index/self.PM.scaleFactor)
        max_y = max_index - self.PM.scaleFactor*max_x
        max_z = -1

        ret_values = (max_value, max_x, max_y)

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
        for z in range(self.PM.scaleFactor):
            realz = self.PM.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    # check that the point is real and that it has not yet been binned
                    if((x,y,realz) in self.im2RowIndicies):
                        for row_index in self.im2RowIndicies[(x,y,realz)]:
                            if row_index not in self.PM.binnedRowIndicies and row_index not in self.PM.restrictedRowIndicies:
                                # this is an unassigned point. 
                                multiplier = np.log10(self.PM.contigLengths[row_index])
                                self.incrementAboutPoint3D(working_block, x-x_lower, y-y_lower, z,multiplier=multiplier)
                                super_putative_row_indicies.append(row_index)

        # blur and find the highest value
        bwb = ndi.gaussian_filter(working_block, self.blurRadius)
        densest_index = np.unravel_index(np.argmax(bwb), (np.shape(bwb)))
        max_x = densest_index[0] + x_lower
        max_y = densest_index[1] + y_lower
        max_z = densest_index[2]
                    
        # now get the basic color of this dense point
        center_k_vals = np.array([])
        center_contig_lengths = np.array([]) 
        putative_center_row_indicies = []

        (x_lower, x_upper) = self.makeCoordRanges(max_x, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, self.span)
        (z_lower, z_upper) = self.makeCoordRanges(max_z, 2*self.span)
        
        # fix this
        tmp = z_lower
        z_lower = self.PM.scaleFactor - z_upper - 1
        z_upper = self.PM.scaleFactor - tmp - 1
        for row_index in super_putative_row_indicies:
            point = np.around(self.PM.transformedCP[row_index])
            if(point[0] >= x_lower and point[0] < x_upper):
                if(point[1] >= y_lower and point[1] < y_upper):
                    if(point[2] >= z_lower and point[2] < z_upper):
                        # we are within the range!
                        center_k_vals = np.append(center_k_vals, self.PM.kmerVals[row_index])
                        center_contig_lengths = np.append(center_contig_lengths, self.PM.contigLengths[row_index])
                        putative_center_row_indicies.append(row_index)

        # make sure we have something to go on here
        if(np.size(center_k_vals) == 0):
            print "GLUG"
            return (np.array([]), (-1,-1,-1))
        elif(np.size(center_k_vals) == 1):
            return (np.array([putative_center_row_indicies[0]]), ret_values)
        else:
            # find the kmer val which lies at the center
            # more specifically, we need to find the kmer vall which
            # is giving the strongest signal
            # partition of the colours in this column
            partitions = self.partitionKvals(center_k_vals)
            # now go and score these guys
            chosen_partiton_index = -1
            max_p_score = 0
            print len(partitions), [len(p) for p in partitions]
            for i in range(len(partitions)):
                # find the greatest number of base pairs
                score = sum([center_contig_lengths[j] for j in partitions[i]])
                if(score > max_p_score):
                    max_p_score = score
                    chosen_partiton_index = i
    
            if(chosen_partiton_index == -1):
                print "flug", partitions
                return (np.array([]), (-1,-1,-1))
            
            # restrict the fluff...
            for i in range(len(partitions)):
                if(len(partitions[i]) < 5 and i != chosen_partiton_index):
                    self.restrictRowIndicies([putative_center_row_indicies[j] for j in partitions[i]])
            
            if(self.debugPlots):
                median_kVal = np.median([center_k_vals[i] for i in partitions[chosen_partiton_index]])
                self.plotRegion(max_x,max_y,max_z, fileName="Image_"+str(self.imageCounter), tag="first approx (%f)" % median_kVal)
                self.imageCounter += 1

            center_row_indicies = np.array([putative_center_row_indicies[i] for i in partitions[chosen_partiton_index]])
            return (center_row_indicies, ret_values)

    def expandSelection(self, startIndex, vals, stdevCutoff=0.05, maxSpread=0.05):
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

    def partitionKvals(self, kVals, stdevCutoff=0.04):
        """Work out where shifts in kmer vals happen"""
        partitions = []
        working_list = list(kVals)
        fix_dict = dict(zip(range(len(working_list)),range(len(working_list))))
        while(len(working_list) > 2):
            cf = CenterFinder()
            c_index = cf.findArrayCenter(working_list)
            expanded_indicies = self.expandSelection(c_index, working_list, stdevCutoff=stdevCutoff)
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
    
    def Ablur(self, blur, density, incAtPoint, index, offset, size):
        """AUX: Used when finding the densest point in a small block"""
        point = index + offset;
        if(point >= 0 and point < size):
            blur[point] += incAtPoint[abs(offset)] * density[index]
            
    def updatePostBin(self, bin):
        """Update data structures after assigning contigs to a new bin"""
        for row_index in bin.rowIndicies:
            self.setRowIndexAssigned(row_index)
            
    def setRowIndexAssigned(self, rowIndex):
        """fix the data structures to indicate that rowIndex belongs to a bin
        
        Use only during initial core creation
        """        
        self.PM.binnedRowIndicies[rowIndex] = True
        
        # now update the image map, decrement
        self.decrementViaRowIndex(rowIndex)

    def setRowIndexUnassigned(self, rowIndex):
        """fix the data structures to indicate that rowIndex no longer belongs to a bin
        
        Use only during initial core creation
        """
        del self.PM.binnedRowIndicies[rowIndex]
        
        # now update the image map, increment
        self.incrementViaRowIndex(rowIndex)

    def restrictRowIndicies(self, indicies):
        """Add these indicies to the restricted list"""
        for row_index in indicies:
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

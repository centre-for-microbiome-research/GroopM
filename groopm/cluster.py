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
    def __init__(self, dbFileName, plot=False, force=False):
        # worker classes
        self.PM = dataManagers.ProfileManager(dbFileName) # store our data
        self.BM = dataManagers.BinManager(pm=self.PM)   # store our bins
    
        # heat maps
        self.imageMaps = np.zeros((3,self.PM.scaleFactor,self.PM.scaleFactor))
        self.blurredMaps = np.zeros((3,self.PM.scaleFactor,self.PM.scaleFactor))
        
        # we need a way to reference from the imageMaps back onto the transformed data
        self.im2RowIndicies = {}  
        
        # When blurring the raw image maps I chose a radius to suit my data, you can vary this as you like
        self.blurRadius = 12
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

    def expandBins(self, force=False):
        """Expand bins using SOMs"""
        SM = dataManagers.SOMManager(self.PM, self.BM, load=True)
        SM.renderWeights("test")
        return
        SM = dataManagers.SOMManager(self.PM, self.BM)
        SM.buildSomWeights(force)
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
        if(self.debugPlots):
            self.PM.renderTransCPData()
        t2 = time.time()
        print "    THIS: [",self.secondsToStr(t2-t1),"]\tTOTAL: [",self.secondsToStr(t2-t0),"]"
        
        # cluster and bin!
        print "Create cores"
        cum_contigs_used_good = self.initialiseCores(minVol)
        t3 = time.time()
        print "    THIS: [",self.secondsToStr(t3-t2),"]\tTOTAL: [",self.secondsToStr(t3-t0),"]"
        
        # now we assume that some true bins may be separated across two cores
        # try to condense things a little
        #self.BM.plotBins(FNPrefix="PRE_")
        print "Condense cores"
        #self.condenseCores()
        t4 = time.time()
        print "    THIS: [",self.secondsToStr(t4-t3),"]\tTOTAL: [",self.secondsToStr(t4-t0),"]"
        #self.BM.plotBins()

        # Now save all the stuff to disk!
        print "Saving bins"
        self.BM.saveBins(doCores=True, saveBinStats=True)
        t5 = time.time()
        print "    THIS: [",self.secondsToStr(t5-t4),"]\tTOTAL: [",self.secondsToStr(t5-t0),"]"

    def initialiseCores(self, minVol):
        """Process contigs and form CORE bins"""
        small_bin_cutoff = 10           # anything with less than this many contigs is not a real bin
        num_below_cutoff = 0            # how many consecutive attempts have produced small bins
        breakout_point = 50             # how many will we allow before we stop this loop
        
        # First we need to find the centers of each blob.
        # We can make a heat map and look for hot spots
        self.populateImageMaps()
        print "    .... .... .... .... .... .... .... .... .... ...."
        print "   ",
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
            (center_row_indicies, max_blur_value) = self.findNewClusterCenter()
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
                    bin.plotBin(self.PM.transformedCP, self.PM.contigColours, fileName="Image_"+str(self.imageCounter))
                    self.imageCounter += 1
        
                # neaten up the bin
                self.removeOutliers(bin.id, fixBinnedRI=False)
                
                # recruit more contigs
                bin_size = bin.recruit(self.PM.transformedCP, self.PM.kmerSigs, self.im2RowIndicies, self.PM.binnedRowIndicies)

                is_good_bin = False
                if(bin.calcTotalSize(self.PM.contigLengths) < minVol):    # less than the good volume
                    if(bin_size > small_bin_cutoff):                      # but has enough contigs
                        is_good_bin = True
                else:                                                     # contains enough bp to pass regardless of number of contigs
                    is_good_bin = True
                if(is_good_bin):
                    # Plot?
                    if(self.debugPlots):          
                        bin.plotBin(self.PM.transformedCP, self.PM.contigColours, fileName="Image_"+str(self.imageCounter))
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
                    print "\n   ",
                                        
                if(self.debugPlots):
                    self.plotHeat("HM_"+str(self.roundNumber)+".png", max=max_blur_value)
        
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
                                    self.PM.contigColours)

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
        self.imageMaps = np.zeros((3,self.PM.scaleFactor,self.PM.scaleFactor))
        self.im2RowIndicies = {}
        
        # add to the grid wherever we find a contig
        row_index = -1
        for point in np.around(self.PM.transformedCP):
            row_index += 1

            # can only bin things once!
            if row_index not in self.PM.binnedRowIndicies and row_index not in self.PM.restrictedRowIndicies:
                # readability
                px = point[0]
                py = point[1]
                pz = point[2]
                
                # add to the row_index dict so we can relate the 
                # map back to individual points later
                if (px,py,pz) in self.im2RowIndicies:
                    self.im2RowIndicies[(px,py,pz)].append(row_index)
                else:
                    self.im2RowIndicies[(px,py,pz)] = [row_index]
                
                # now increment in the grid
                # for each point we encounter we incrmement
                # it's position + the positions to each side
                # and touching each corner
                multiplier = np.log10(self.PM.contigLengths[row_index])
                self.incrementAboutPoint(0, px, py, multiplier=multiplier)
                self.incrementAboutPoint(1, self.PM.scaleFactor - pz - 1, py, multiplier=multiplier)
                self.incrementAboutPoint(2, self.PM.scaleFactor - pz - 1, self.PM.scaleFactor - px - 1, multiplier=multiplier)

    def incrementViaRowIndex(self, rowIndex):
        """Wrapper to increment about point"""
        point = np.around(self.PM.transformedCP[rowIndex])
        # readability
        px = point[0]
        py = point[1]
        pz = point[2]
        multiplier = np.log10(self.PM.contigLengths[rowIndex])
        self.incrementAboutPoint(0, px, py, multiplier=multiplier)
        self.incrementAboutPoint(1, self.PM.scaleFactor - pz - 1, py, multiplier=multiplier)
        self.incrementAboutPoint(2, self.PM.scaleFactor - pz - 1, self.PM.scaleFactor - px - 1, multiplier=multiplier)

    def decrementViaRowIndex(self, rowIndex):
        """Wrapper to decrement about point"""
        point = np.around(self.PM.transformedCP[rowIndex])
        # readability
        px = point[0]
        py = point[1]
        pz = point[2]
        multiplier = np.log10(self.PM.contigLengths[rowIndex])
        self.decrementAboutPoint(0, px, py, multiplier=multiplier)
        self.decrementAboutPoint(1, self.PM.scaleFactor - pz - 1, py, multiplier=multiplier)
        self.decrementAboutPoint(2, self.PM.scaleFactor - pz - 1, self.PM.scaleFactor - px - 1, multiplier=multiplier)

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
                self.imageMaps[view_index,px-1,py-1] -= valC      # Top left corner
                if self.imageMaps[view_index,px-1,py-1] < np.finfo(float).eps:
                    self.imageMaps[view_index,px-1,py-1] = 0
                
            self.imageMaps[view_index,px-1,py] -= valS            # Top
            if self.imageMaps[view_index,px-1,py] < np.finfo(float).eps:
                self.imageMaps[view_index,px-1,py] = 0
            if py < self.PM.scaleFactor-1:             
                self.imageMaps[view_index,px-1,py+1] -= valC      # Top right corner
                if self.imageMaps[view_index,px-1,py+1] < np.finfo(float).eps:
                    self.imageMaps[view_index,px-1,py+1] = 0

        if py > 0:
            self.imageMaps[view_index,px,py-1] -= valS            # Left side
            if self.imageMaps[view_index,px,py-1] < np.finfo(float).eps:
                self.imageMaps[view_index,px,py-1] = 0
            
        self.imageMaps[view_index,px,py] -= valP                  # Point
        if self.imageMaps[view_index,px,py] < np.finfo(float).eps:
            self.imageMaps[view_index,px,py] = 0
        if py < self.PM.scaleFactor-1:             
            self.imageMaps[view_index,px,py+1] -= valS            # Right side
            if self.imageMaps[view_index,px,py+1] < np.finfo(float).eps:
                self.imageMaps[view_index,px,py+1] = 0

        if px < self.PM.scaleFactor-1:
            if py > 0:
                self.imageMaps[view_index,px+1,py-1] -= valC      # Bottom left corner
                if self.imageMaps[view_index,px+1,py-1] < np.finfo(float).eps:
                    self.imageMaps[view_index,px+1,py-1] = 0
            self.imageMaps[view_index,px+1,py] -= valS            # Bottom
            if self.imageMaps[view_index,px+1,py] < np.finfo(float).eps:
                self.imageMaps[view_index,px+1,py] = 0
            if py < self.PM.scaleFactor-1:             
                self.imageMaps[view_index,px+1,py+1] -= valC      # Bottom right corner
                if self.imageMaps[view_index,px+1,py+1] < np.finfo(float).eps:
                    self.imageMaps[view_index,px+1,py+1] = 0

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
        self.blurredMaps = np.zeros((3,self.PM.scaleFactor,self.PM.scaleFactor))
        #self.maxMaps = np.zeros((3,self.PM.scaleFactor,self.PM.scaleFactor))
        
        for i in range (0,3): # top, front and side
            self.blurredMaps[i,:,:] = ndi.gaussian_filter(self.imageMaps[i,:,:]**0.5, (self.blurRadius,self.blurRadius)) 

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
        for i in range (0,3): # top, front and side
            self.blurredMaps[i,:,:] = np.where(self.blurredMaps[i,:,:] >= lop_val, self.blurredMaps[i,:,:], 0)/lop_val

    def makeCoordRanges(self, pos, span):
        """Make search ranges which won't go out of bounds"""
        lower = pos-span
        upper = pos+span+1
        if(lower < 0):
            lower = 0
        if(upper >= self.PM.scaleFactor):
            upper = self.PM.scaleFactor - 1
        return (lower, upper)

    def findNewClusterCenter(self):
        """Find a putative cluster"""
        # we work from the top view as this has the base clustering
        max_index = np.argmax(self.blurredMaps[0])
        max_value = self.blurredMaps[0].ravel()[max_index]

        max_x = int(max_index/self.PM.scaleFactor)
        max_y = max_index - self.PM.scaleFactor*max_x
        max_z = -1

        this_span = int(1.5 * self.span)
        span_len = 2*this_span+1
        
        # work out the region this max value lives in
        x_density = np.zeros(span_len)
        x_offset = max_x - this_span
        
        y_density = np.zeros(span_len)
        y_offset = max_y - this_span
        
        if(self.debugPlots):
            self.plotRegion(max_x,max_y,max_z, fileName="Image_"+str(self.imageCounter), tag="column", column=True)
            self.imageCounter += 1

        # make a 3d grid to hold the values
        working_block = np.zeros((span_len, span_len, self.PM.scaleFactor))
        
        # go through the entire column
        (x_lower, x_upper) = self.makeCoordRanges(max_x, this_span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, this_span)
        for z in range(0, self.PM.scaleFactor):
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

        # blur and find the highest value
        bwb = ndi.gaussian_filter(working_block, self.blurRadius)
        
        densest_index = np.unravel_index(np.argmax(bwb), (np.shape(bwb)))
        max_x = densest_index[0] + x_lower
        max_y = densest_index[1] + y_lower
        max_z = densest_index[2]
       
        if(self.debugPlots):
            self.plotRegion(max_x,max_y,max_z, fileName="Image_"+str(self.imageCounter), tag="first approx")
            self.imageCounter += 1

        # now get the basic color of this dense point
        (x_lower, x_upper) = self.makeCoordRanges(max_x, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, self.span)
        (z_lower, z_upper) = self.makeCoordRanges(max_z, self.span)
        center_values = np.array([])
        cv_colours = np.array([])
        c_inc = 0
        for z in range(z_lower, z_upper):
            realz = self.PM.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.im2RowIndicies):
                        for row_index in self.im2RowIndicies[(x,y,realz)]:
                            if row_index not in self.PM.binnedRowIndicies and row_index not in self.PM.restrictedRowIndicies:
                                center_values = np.append(center_values, self.PM.kmerSigs[row_index])
                                cv_colours = np.append(cv_colours, self.PM.contigColours[row_index])
                                c_inc += 1

        # make sure we have something to go on here
        if(np.size(center_values) == 0 or c_inc < 2):
            return (np.array([]), -1)

        # reshape these guys!
        center_values = np.reshape(center_values, (c_inc, np.size(self.PM.kmerSigs[0])))
        cv_colours = np.reshape(cv_colours, (c_inc, 3))
        
        # transform them into one dimensional points
        oneD_center_values = np.zeros(c_inc)
        working_centres = np.array(center_values, copy=True) 
        PCA.Center(working_centres,verbose=0)
        p = PCA.PCA(working_centres)
        components = p.pc()
        index = 0
        for point in components:
            oneD_center_values[index] = components[index,0]
            index += 1

        if(False):
            plt.figure(1)
            plt.subplot(111)
            cm = mpl.colors.LinearSegmentedColormap('my_colormap', cv_colours, 1024)
            plt.scatter(components[:,0], components[:,1], edgecolors=cv_colours, c=cv_colours, cmap=cm, marker='.')
            plt.show()

        # find the sig which lies at the center
        cf = CenterFinder()
        oneD_center_values -= np.min(oneD_center_values)
        oneD_center_values /= np.max(oneD_center_values)
        centroid_sig = center_values[cf.findArrayCenter(oneD_center_values)]
        
        # now we need to work out how close to this sig we need to be...
        dists = np.array([])
        for sig in center_values:
            dists = np.append(dists, np.linalg.norm(sig - centroid_sig))
            
        # get a first approximation on the upper cutoff
        tol = 2
        dists = np.sort(dists)
        sub_dists = dists[0:int(np.size(dists)/4):1]
        if(not np.size(sub_dists) > 0):
            return (np.array([]), -1)
        upper_dist = np.mean(sub_dists) + tol * np.std(sub_dists)

        # now refine (expand) based on the first approximation
        sub_dists = np.array([x for x in dists if x < upper_dist])
        if(np.size(sub_dists) == 0):
            return (np.array([]), -1)                
            
        upper_dist = np.mean(sub_dists) + tol * np.std(sub_dists)

        # now scoot out around this point and soak up similar points
        # get all the real rowIndicies of these points so we can use them in
        # the primary data map
        center_row_indicies = np.array([])
        for z in range(z_lower, z_upper):
            realz = self.PM.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    # check that the point is real and that it has not yet been binned
                    if((x,y,realz) in self.im2RowIndicies):
                        for row_index in self.im2RowIndicies[(x,y,realz)]:
                            if(row_index not in center_row_indicies) and (row_index not in self.PM.binnedRowIndicies and row_index not in self.PM.restrictedRowIndicies):
                                # make sure the kmer sig is close enough
                                dist = np.linalg.norm(self.PM.kmerSigs[row_index] - centroid_sig)
                                if(dist < upper_dist):
                                    center_row_indicies = np.append(center_row_indicies, row_index)
        if(np.size(center_row_indicies) > 0):
            return (center_row_indicies, max_value)
        
        return (np.array([]), -1)

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
        
        # now update the image map, decrement
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
    
    def plotHeat(self, fileName = "", max=-1):
        """Print the main heat maps
        
        Useful for debugging
        """
        fig = plt.figure()
        images = []

        ax = fig.add_subplot(231)
        images.append(ax.imshow(self.blurredMaps[0,:,:]**0.5))
        if(max > 0):
            title = str.join(" ", ["Max value:",str(max)])
            plt.title(title)
        ax = fig.add_subplot(232)
        images.append(ax.imshow(self.blurredMaps[1,:,:]**0.5))
        ax = fig.add_subplot(233)
        images.append(ax.imshow(self.blurredMaps[2,:,:]**0.5))

        ax = fig.add_subplot(234)
        images.append(ax.imshow(self.imageMaps[0,:,:]**0.5))
        ax = fig.add_subplot(235)
        images.append(ax.imshow(self.imageMaps[1,:,:]**0.5))
        ax = fig.add_subplot(236)
        images.append(ax.imshow(self.imageMaps[2,:,:]**0.5))
        
        if(fileName != ""):
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
    
    NOTE: Assumes the data is normalised between 0 and 1!
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
        # run through in one direction
        vals_sorted = np.sort(vals)
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
        
        max_index = np.argmax(working)
        vals_sorted = np.sort(vals_sorted)
        max_value = vals_sorted[max_index]
        
        # find the original index!
        index = 0
        for val in vals:
            if(val == max_value):
                return index
            index += 1
        return -1
    
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

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
import mstore
import binUtils

np.seterr(all='raise')      
###############################################################################
###############################################################################
###############################################################################
###############################################################################
class ClusterEngine:
    """Top level interface for clustering contigs"""
    def __init__(self, dbFileName, plot=False, outFile="", force=False):
        # worker classes
        self.DB = mstore.DataBlob(dbFileName)
        self.CB = ClusterBlob(self.DB, debugPlots=plot)
    
        # misc
        self.plot = plot
        self.outFile = outFile
        self.forceWriting = force
        
    def makeCores(self, coreCut, minSize, minVol):
        """Cluster the contigs to make bin cores"""
        # check that the user is OK with nuking stuff...
        if(not self.promptForOverwrite()):
            return False

        # get some data
        t0 = time.time()
        print "Load data"
        self.DB.loadData(condition="length >= "+str(coreCut))
        t1 = time.time()
        print "\tTHIS: [",self.secondsToStr(t1-t0),"]\tTOTAL: [",self.secondsToStr(t1-t0),"]"
        
        # transform the data
        print "Apply data transformations"
        self.DB.transformCP()
        # plot the transformed space (if we've been asked to...)
        if(self.plot):
            self.DB.renderTransCPData()
        t2 = time.time()
        print "\tTHIS: [",self.secondsToStr(t2-t1),"]\tTOTAL: [",self.secondsToStr(t2-t0),"]"
        
        # cluster and bin!
        print "Create cores"
        cum_contigs_used_good = self.CB.initialiseCores(minVol)
        t3 = time.time()
        print "\tTHIS: [",self.secondsToStr(t3-t2),"]\tTOTAL: [",self.secondsToStr(t3-t0),"]"
        
        # now we assume that some true bins may be separated across two cores
        # try to condense things a little
        self.CB.plotBins(FNPrefix="PRE_")
        if(True):
            print "Condense cores"
            changed = True
            con_round = 1
            while(changed):
                changed = self.CB.condenseCores(cum_contigs_used_good, con_round, minSize)
                con_round+=1
        t4 = time.time()
        print "\tTHIS: [",self.secondsToStr(t4-t3),"]\tTOTAL: [",self.secondsToStr(t4-t0),"]"

        self.CB.plotBins()        
        # Now save all the stuff to disk!
        print "Saving bins"
        self.saveBins()
        t5 = time.time()
        print "\tTHIS: [",self.secondsToStr(t5-t4),"]\tTOTAL: [",self.secondsToStr(t5-t0),"]"

    def saveBins(self):
        """Save binning results"""
        (c2b_update, core_update) = self.CB.getCoreBinUpdates()
        self.DB.saveBins(c2b_update)
        self.DB.saveCores(core_update)
        self.DB.setClustered()
        # Merge bids and number of members so we can save to disk
        bin_updates = {}
        for bid in self.CB.bins:
            bin_updates[bid] = np.size(self.CB.bins[bid].indicies)
        self.DB.saveBinIds(bin_updates)
        
    def expandBins(self):
        """Load cores and expand bins"""
        # check that the user is OK with nuking stuff...
        if(not self.promptForOverwrite()):
            return False
        
        # get some data
        t0 = time.time()
        print "Load data"
        self.DB.loadData(condition="length >= 10000")#(length >= 4000 ) & (length <= 4300)")
        t1 = time.time()
        print "\tTHIS: [",self.secondsToStr(t1-t0),"]\tTOTAL: [",self.secondsToStr(t1-t0),"]"

        # Now we use SOMs to classify the remaininfg contigs
        print "Start SOM classification"
        t2 = time.time()
        print "\tTHIS: [",self.secondsToStr(t2-t1),"]\tTOTAL: [",self.secondsToStr(t2-t0),"]"

    def promptForOverwrite(self):
        """Check that the user is ok with possibly overwriting the DB"""
        if(not self.forceWriting):
            if(self.DB.isClustered()):
                option = raw_input(" ****WARNING**** Database: '"+self.DB.dbFileName+"' has already been clustered.\n" \
                                   " If you continue you *MAY* overwrite existing bins!\n" \
                                   " Overwrite? (y,n) : ")
                print "****************************************************************"
                if(option.upper() != "Y"):
                    print "Operation cancelled"
                    return False
                else:
                    print "Overwriting database",self.DB.dbFileName
                    self.DB.dataManager.nukeBins(self.DB.dbFileName)
        elif(self.DB.isClustered()):
            print "Overwriting database",self.DB.dbFileName
            self.DB.dataManager.nukeBins(self.DB.dbFileName)
        return True
    
    def secondsToStr(self, t):
        rediv = lambda ll,b : list(divmod(ll[0],b)) + ll[1:]
        return "%d:%02d:%02d.%03d" % tuple(reduce(rediv,[[t*1000,],1000,60,60]))
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################
class ClusterBlob:
    """Main class for performing contig clustering
    
    All the bits and bobs you'll need to cluster and bin out 
    pre-transformed primary data
    """    
    def __init__(self, dataBlob, debugPlots=False):
        # See DataTransformer for details about these variables
        self.DB = dataBlob
        self.scaleFactor = self.DB.scaleFactor
        
        # get enough memory for three heat maps
        self.imageMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        self.blurredMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        #self.maxMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        
        # we need a way to reference from the imageMaps back onto the transformed data
        self.mappedIndicies = {}
        self.binnedIndicies = {}

        # store our bins
        self.nextFreeBinId = 0          # increment before use!
        self.bins = {}                  # bins we kinda trust
        self.badBins = {}               # bins we don't kinda trust
        
        # housekeeping / parameters
        self.roundNumber = 0            # how many times have we done this?
        self.span = 30                  # amount we can travel about when determining "hot spots"

        # When blurring the raw image maps I chose a radius to suit my data, you can vary this as you like
        self.blurRadius = 12
        
        self.debugPlots = debugPlots
        self.imageCounter = 1           # when we print many images

    def initialiseCores(self, minVol):
        """Process contigs and form CORE bins"""
        small_bin_cutoff = 10           # anything with less than this many contigs is not a real bin
        num_below_cutoff = 0            # how many consecutive attempts have produced small bins
        breakout_point = 50             # how many will we allow before we stop this loop
        tmp_num_bins = 10000            # gotta keep count, high numbered bins ar baaad!
        
        cum_contigs_used_good = 0
        cum_contigs_used_bad = 0
        # First we need to find the centers of each blob.
        # We can make a heat map and look for hot spots
        self.populateImageMaps()
        print "\t",
        nl_counter = 0
        while(num_below_cutoff < breakout_point):
            #if(self.numBins > 3):
            #    break
            sys.stdout.flush()
            # apply a gaussian blur to each image map to make hot spots
            # stand out more from the background 
            self.blurMaps()
    
            # now search for the "hottest" spots on the blurred map
            # this is a new bin centroid
            (center_indicies, max_blur_value) = self.findNewClusterCenter()
            if(np.size(center_indicies) == 0):
                break
            else:
                # make sure we got something
                self.roundNumber += 1
                
                # time to make a bin
                tmp_num_bins += 1
                bin = binUtils.Bin(center_indicies, self.DB.kmerSigs, tmp_num_bins)
                
                # work out the distribution in points in this bin
                bin.makeBinDist(self.DB.transformedCP, self.DB.kmerSigs)     
                
                # Plot?
                if(self.debugPlots):          
                    bin.plotBin(self.DB.transformedCP, self.DB.contigColours, fileName="Image_"+str(self.imageCounter), tag="Initial")
                    self.imageCounter += 1
        
                # make the bin more gooder
                is_good_bin = True
                bin_size = bin.recruit(self.DB.transformedCP, self.DB.kmerSigs, self.mappedIndicies, self.binnedIndicies)
                if(bin.calcTotalSize(self.DB.contigLengths) < minVol):    # less than the good volume
                    if(bin_size < small_bin_cutoff):
                        is_good_bin = False
                        cum_contigs_used_bad += bin_size
                        self.badBins[tmp_num_bins] = bin
                        num_below_cutoff += 1
                        print "-",
                    # else bin is large enough!
                if(is_good_bin):
                    # make this bin legit!
                    cum_contigs_used_good += bin_size
                    self.nextFreeBinId += 1
                    bin.id = self.nextFreeBinId 
                    self.bins[self.nextFreeBinId] = bin
                    # Plot?
                    if(self.debugPlots):          
                        bin.plotBin(self.DB.transformedCP, self.DB.contigColours, fileName="Image_"+str(self.imageCounter), tag="CORE")
                        self.imageCounter += 1
                    num_below_cutoff = 0
                    print "+",

                # make the printing prettier
                nl_counter += 1
                if(nl_counter > 9):
                    nl_counter = 0
                    print "\n\t",
                                        
                if(self.debugPlots):
                    self.plotHeat("3X3_"+str(self.roundNumber)+".png", max=max_blur_value)

                # append this bins list of mapped indicies to the main list
                self.updatePostBin(bin)
        print ""
        perc = "%.2f" % round((float(cum_contigs_used_good)/float(self.DB.numContigs))*100,2)
        print "\t",cum_contigs_used_good,"contigs are distributed across",len(self.bins),"cores (",perc,"% )"
        print "\t",cum_contigs_used_bad,"contigs are distributed across",len(self.badBins),"pseudo cores"
        
        return cum_contigs_used_good

    def condenseCores(self, cumContigsUsedGood, con_round, minSize):
        """combine similar CORE bins"""
        print "\tCondense, round:", con_round
        changed = False
        mer_stdevs = 5
        cov_stdevs = 4
        consumed_indicies = {} # who is getting consumed by who?

        # go through all the bins, sorted according to kmer profile
        num_cores_consumed = 0
        num_cores_upgraded = 0
        num_cores_deleted = 0
        contigs_upgraded = 0
        contigs_freed = 0
        all_bins_sorted = sorted(self.bins.values() + self.badBins.values())

        for subject_index in range(0, len(all_bins_sorted)):
            if(subject_index not in consumed_indicies):
                subject_bin = all_bins_sorted[subject_index]
                subject_upper = subject_bin.kDistMean + mer_stdevs * subject_bin.kDistStdev  
                query_index = subject_index+1
                while(query_index < len(all_bins_sorted)):
                    if(query_index not in consumed_indicies):
                        query_bin = all_bins_sorted[query_index]
                        query_upper = query_bin.kDistMean + mer_stdevs * query_bin.kDistStdev 
                        # only worth comparing if their kmersigs are similar
                        k_dists = np.array([])
                        continue_check = False
                        # pick the bin with the highest kDistStdev
                        if(query_upper > subject_upper):
                            for index in subject_bin.indicies:
                                k_dists = np.append(k_dists, query_bin.getKDist(self.DB.kmerSigs[index]))
                            if np.median(k_dists) <= query_upper:
                                continue_check = True
                        else:
                            for index in query_bin.indicies:
                                k_dists = np.append(k_dists, subject_bin.getKDist(self.DB.kmerSigs[index]))
                            if np.median(k_dists) <= subject_upper:
                                continue_check = True

                        if(continue_check):
                            if(subject_bin.isSimilar(query_bin, stdevs=cov_stdevs)):
                                consumed_indicies[query_index] = subject_index

                    query_index += 1
            subject_index += 1
        
        # first we worry about making all the bins ok, then we delete the consumed guys
        dead_bids = []
        for dead_index in consumed_indicies.keys():
            num_cores_consumed += 1
            if dead_index in self.badBins:
                contigs_upgraded += all_bins_sorted[dead_index].binSize
            print all_bins_sorted[dead_index].id,"consumed by:",all_bins_sorted[consumed_indicies[dead_index]].id
            all_bins_sorted[consumed_indicies[dead_index]].consume(self.DB.transformedCP,
                                                                       self.DB.kmerSigs,
                                                                       self.DB.contigLengths,
                                                                       all_bins_sorted[dead_index]
                                                                       )
            # save this guy for the killing!
            dead_bids.append(all_bins_sorted[dead_index].id)
            changed = True
            
        # remove the consumed bins
        for bid in dead_bids:
            if(bid in self.badBins):
                del self.badBins[bid]
            elif(bid in self.bins):
                del self.bins[bid]

        print "\tIncorporated:",num_cores_consumed,"smaller cores into larger ones"

        if(not changed):
            # consuming done, now upgrade or delete "bad" bins 
            # move all bad_bins with greater than minSize contigs into self.bins.
            dead_cores = []
            for bid in self.badBins.keys():
                if(bid not in dead_bids):
                    if(self.badBins[bid].binSize > minSize):
                        # move it to the good bins pile
                        num_cores_upgraded += 1
                        contigs_upgraded += self.badBins[bid].binSize
                        self.nextFreeBinId += 1 
                        self.badBins[bid].id = self.nextFreeBinId
                        self.bins[self.nextFreeBinId] = self.badBins[bid]
                        dead_cores.append(bid)
                    else:
                        # we need to free these indicies!
                        for index in self.badBins[bid].indicies:
                             del self.binnedIndicies[index]
                             contigs_freed += 1
                        num_cores_deleted += 1
                        dead_cores.append(bid)
            
            # delete the dead cores
            for dead_core in dead_cores:
                del self.badBins[dead_core]
                
            # now we prettify the bids!
            bid_upgrades = []
            for bid in self.bins.keys():
                if bid > self.nextFreeBinId:
                    bid_upgrades.append(bid)
            for bid in bid_upgrades:
                self.nextFreeBinId += 1 
                self.bins[bid].id = self.nextFreeBinId
                self.bins[self.nextFreeBinId] = self.bins[bid]
                del self.bins[bid] 
            
            self.numBins = len(self.bins)
            
            print "\tUpgraded:",num_cores_upgraded,"psuedo cores"
            print "\tDeleted:",num_cores_deleted,"psuedo cores"
            cumContigsUsedGood += contigs_upgraded
            perc = "%.2f" % round((float(cumContigsUsedGood)/float(self.DB.numContigs))*100,2)
            print "\t",(cumContigsUsedGood),"contigs are distributed across",self.numBins,"cores (",perc,"% )"
            print "\t",contigs_freed,"contigs free'd"
        return changed

    def getCoreBinUpdates(self):
        """Merge the bids, raw DB indexes and core information so we can save to disk"""
        core_update = dict(zip(self.DB.indicies, [False]*np.size(self.DB.indicies)))
        bin_update = dict(zip(self.DB.indicies, [0]*np.size(self.DB.indicies)))

        # we need a mapping from cid (or local index) to binID
        c2b = dict(zip(range(0,np.size(self.DB.indicies)), [0]*np.size(self.DB.indicies)))
        for bid in self.bins:
            for index in self.bins[bid].indicies:
                c2b[index] = bid
        
        # at this stage, all bins are cores
        for index in range(0, self.DB.numContigs):
            if index in self.binnedIndicies:
                bin_update[self.DB.indicies[index]] = c2b[index]
                core_update[self.DB.indicies[index]] = True

        return (bin_update, core_update)

    def populateImageMaps(self):
        """Load the transformed data into the main image maps"""
        # reset these guys... JIC
        self.imageMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        self.mappedIndicies = {}
        
        # add to the grid wherever we find a contig
        index = -1
        for point in np.around(self.DB.transformedCP):
            index += 1

            # can only bin things once!
            if index not in self.binnedIndicies:
                # readability
                px = point[0]
                py = point[1]
                pz = point[2]
                
                # add to the index dict so we can relate the 
                # map back to individual points later
                if (px,py,pz) in self.mappedIndicies:
                    self.mappedIndicies[(px,py,pz)].append(index)
                else:
                    self.mappedIndicies[(px,py,pz)] = [index]
                
                # now increment in the grid
                # for each point we encounter we incrmement
                # it's position + the positions to each side
                # and touching each corner
                multiplier = np.log10(self.DB.contigLengths[index])
                self.incrementAboutPoint(0, px, py, multiplier=multiplier)
                self.incrementAboutPoint(1, self.scaleFactor - pz - 1, py, multiplier=multiplier)
                self.incrementAboutPoint(2, self.scaleFactor - pz - 1, self.scaleFactor - px - 1, multiplier=multiplier)

    def updatePostBin(self, bin):
        """Update data structures after assigning contigs to a new bin"""
        for index in bin.indicies:
            self.binnedIndicies[index] = True
            
            # now update the image map, decrement
            point = np.around(self.DB.transformedCP[index])
            # readability
            px = point[0]
            py = point[1]
            pz = point[2]
            multiplier = np.log10(self.DB.contigLengths[index])
            self.decrementAboutPoint(0, px, py, multiplier=multiplier)
            self.decrementAboutPoint(1, self.scaleFactor - pz - 1, py, multiplier=multiplier)
            self.decrementAboutPoint(2, self.scaleFactor - pz - 1, self.scaleFactor - px - 1, multiplier=multiplier)

    def incrementAboutPoint(self, index, px, py, valP=1, valS=0.6, valC=0.2, multiplier=1):
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
                self.imageMaps[index,px-1,py-1] += valC      # Top left corner
            self.imageMaps[index,px-1,py] += valS            # Top
            if py < self.scaleFactor-1:             
                self.imageMaps[index,px-1,py+1] += valC      # Top right corner

        if py > 0:
            self.imageMaps[index,px,py-1] += valS            # Left side
        self.imageMaps[index,px,py] += valP                  # Point
        if py < self.scaleFactor-1:             
            self.imageMaps[index,px,py+1] += valS            # Right side

        if px < self.scaleFactor-1:
            if py > 0:
                self.imageMaps[index,px+1,py-1] += valC      # Bottom left corner
            self.imageMaps[index,px+1,py] += valS            # Bottom
            if py < self.scaleFactor-1:             
                self.imageMaps[index,px+1,py+1] += valC      # Bottom right corner

    def decrementAboutPoint(self, index, px, py, valP=1, valS=0.6, valC=0.2, multiplier=1):
        """Decrement value at a point in the 2D image maps
        
        multiplier is proportional to the contigs length
        """        
        valP *= multiplier
        valS *= multiplier
        valC *= multiplier
        if px > 0:
            if py > 0:
                self.imageMaps[index,px-1,py-1] -= valC      # Top left corner
                if self.imageMaps[index,px-1,py-1] < np.finfo(float).eps:
                    self.imageMaps[index,px-1,py-1] = 0
                
            self.imageMaps[index,px-1,py] -= valS            # Top
            if self.imageMaps[index,px-1,py] < np.finfo(float).eps:
                self.imageMaps[index,px-1,py] = 0
            if py < self.scaleFactor-1:             
                self.imageMaps[index,px-1,py+1] -= valC      # Top right corner
                if self.imageMaps[index,px-1,py+1] < np.finfo(float).eps:
                    self.imageMaps[index,px-1,py+1] = 0

        if py > 0:
            self.imageMaps[index,px,py-1] -= valS            # Left side
            if self.imageMaps[index,px,py-1] < np.finfo(float).eps:
                self.imageMaps[index,px,py-1] = 0
            
        self.imageMaps[index,px,py] -= valP                  # Point
        if self.imageMaps[index,px,py] < np.finfo(float).eps:
            self.imageMaps[index,px,py] = 0
        if py < self.scaleFactor-1:             
            self.imageMaps[index,px,py+1] -= valS            # Right side
            if self.imageMaps[index,px,py+1] < np.finfo(float).eps:
                self.imageMaps[index,px,py+1] = 0

        if px < self.scaleFactor-1:
            if py > 0:
                self.imageMaps[index,px+1,py-1] -= valC      # Bottom left corner
                if self.imageMaps[index,px+1,py-1] < np.finfo(float).eps:
                    self.imageMaps[index,px+1,py-1] = 0
            self.imageMaps[index,px+1,py] -= valS            # Bottom
            if self.imageMaps[index,px+1,py] < np.finfo(float).eps:
                self.imageMaps[index,px+1,py] = 0
            if py < self.scaleFactor-1:             
                self.imageMaps[index,px+1,py+1] -= valC      # Bottom right corner
                if self.imageMaps[index,px+1,py+1] < np.finfo(float).eps:
                    self.imageMaps[index,px+1,py+1] = 0

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
        if pz < self.scaleFactor-1:
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
        self.blurredMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        #self.maxMaps = np.zeros((3,self.scaleFactor,self.scaleFactor))
        
        for i in range (0,3): # top, front and side
            self.blurredMaps[i,:,:] = ndi.gaussian_filter(self.imageMaps[i,:,:]**0.5, (self.blurRadius,self.blurRadius)) 

        # there's still a lot of background signal to remove
        # we wish to remove 90% of the data, this will leave just the really hot spots
        # Make a histogram of the data (use the top face)
        [vals,points] = np.histogram(np.reshape(self.blurredMaps[0,:,:], (self.scaleFactor, self.scaleFactor,1)), 50)
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
        if(upper >= self.scaleFactor):
            upper = self.scaleFactor - 1
        return (lower, upper)

    def findNewClusterCenter(self):
        """Find a putative cluster"""
        # we work from the top view as this has the base clustering
        max_index = np.argmax(self.blurredMaps[0])
        max_value = self.blurredMaps[0].ravel()[max_index]

        max_x = int(max_index/self.scaleFactor)
        max_y = max_index - self.scaleFactor*max_x
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
        working_block = np.zeros((span_len, span_len, self.scaleFactor))
        
        # go through the entire column
        (x_lower, x_upper) = self.makeCoordRanges(max_x, this_span)
        (y_lower, y_upper) = self.makeCoordRanges(max_y, this_span)
        for z in range(0, self.scaleFactor):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    # check that the point is real and that it has not yet been binned
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if index not in self.binnedIndicies:
                                # this is an unassigned point. 
                                multiplier = np.log10(self.DB.contigLengths[index])
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
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if index not in self.binnedIndicies:
                                center_values = np.append(center_values, self.DB.kmerSigs[index])
                                cv_colours = np.append(cv_colours, self.DB.contigColours[index])
                                c_inc += 1

        # make sure we have something to go on here
        if(np.size(center_values) == 0):
            return (np.array([]), -1)

        # reshape these guys!
        center_values = np.reshape(center_values, (c_inc, np.size(self.DB.kmerSigs[0])))
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
        # get all the real indicies of these points so we can use them in
        # the primary data map
        center_indicies = np.array([])
        for z in range(z_lower, z_upper):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    # check that the point is real and that it has not yet been binned
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if(index not in center_indicies) and (index not in self.binnedIndicies):
                                # make sure the kmer sig is close enough
                                dist = np.linalg.norm(self.DB.kmerSigs[index] - centroid_sig)
                                if(dist < upper_dist):
                                    center_indicies = np.append(center_indicies, index)
        if(np.size(center_indicies) > 0):
            return (center_indicies, max_value)
        
        return (np.array([]), -1)

    def Ablur(self, blur, density, incAtPoint, index, offset, size):
        """AUX: Used when finding the densest point in a small block"""
        point = index + offset;
        if(point >= 0 and point < size):
            blur[point] += incAtPoint[abs(offset)] * density[index]
    
#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def plotBins(self, FNPrefix="BIN_"):
        """Make plots of all the bins"""
        for bid in self.bins:
            self.bins[bid].plotBin(self.DB.transformedCP, self.DB.contigColours, fileName=FNPrefix+str(bid),)
            
    def plotRegion(self, px, py, pz, fileName="", tag="", column=False):
        """Plot the region surrounding a point """
        disp_vals = np.array([])
        disp_cols = np.array([])
        num_points = 0
        # plot all points within span
        (z_lower, z_upper) = self.makeCoordRanges(pz, self.span)
        if(column):
            z_lower = 0
            z_upper = self.scaleFactor - 1

        (x_lower, x_upper) = self.makeCoordRanges(px, self.span)
        (y_lower, y_upper) = self.makeCoordRanges(py, self.span)
        for z in range(z_lower, z_upper):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if index not in self.binnedIndicies:
                                num_points += 1
                                disp_vals = np.append(disp_vals, self.DB.transformedCP[index])
                                disp_cols = np.append(disp_cols, self.DB.contigColours[index])
        
        # make a black mark at the max values
        small_span = self.span/2
        (x_lower, x_upper) = self.makeCoordRanges(px, small_span)
        (y_lower, y_upper) = self.makeCoordRanges(py, small_span)
        (z_lower, z_upper) = self.makeCoordRanges(pz, small_span)
        for z in range(z_lower, z_upper):
            realz = self.scaleFactor - z - 1
            for x in range(x_lower, x_upper):
                for y in range(y_lower, y_upper):
                    if((x,y,realz) in self.mappedIndicies):
                        for index in self.mappedIndicies[(x,y,realz)]:
                            if index not in self.binnedIndicies:
                                num_points += 1
                                disp_vals = np.append(disp_vals, self.DB.transformedCP[index])
                                disp_cols = np.append(disp_cols, colorsys.hsv_to_rgb(0,0,0))
        # reshape
        disp_vals = np.reshape(disp_vals, (num_points, 3))
        disp_cols = np.reshape(disp_cols, (num_points, 3))
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        cm = mpl.colors.LinearSegmentedColormap('my_colormap', disp_cols, 1024)
        result = ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors=disp_cols, c=disp_cols, cmap=cm, marker='.')
        title = str.join(" ", ["Focus at: (",str(px), str(py), str(self.scaleFactor - pz - 1),")\n",tag])
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

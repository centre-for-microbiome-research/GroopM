#!/usr/bin/env python
###############################################################################
#                                                                             #
#    groopmUtils.py                                                           #
#                                                                             #
#    Classes for non-clustering data manipulation and output                  #
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
import math
import colorsys
import random

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

import numpy as np
import scipy.ndimage as ndi
import scipy.spatial.distance as ssdist
from scipy.stats import kstest

import networkx as nx

import time

# GroopM imports
import dataManagers

np.seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BinExplorer:
    """Inspect bins, used for validation"""
    def __init__(self, dbFileName, bids=[]):
        self.PM = dataManagers.ProfileManager(dbFileName)   # based on user specified length
        self.BM = dataManagers.BinManager(dbFileName=dbFileName)   # bins
        if bids is None:
            self.bids = []
        else:
            self.bids = bids

    def plotFlyOver(self, fps=10.0, totalTime=120.0):
        """Plot a flyover of the data with bins being removed"""
        self.BM.loadBins(makeBins=True,silent=True,bids=self.bids)
        all_bids = self.bins.keys()

        # control image form and output
        current_azim = 45.0
        current_elev = 0.0
        current_frame = 0.0
        total_frames = fps * totalTime
        total_azim_shift = 720.0
        total_elev_shift = 360.0
        azim_increment = total_azim_shift / total_frames
        elev_increment = total_elev_shift / total_frames
        
        print "Need",total_frames,"frames:"
        # we need to know when to remove each bin
        bid_remove_rate = total_frames / float(len(all_bids))
        bid_remove_indexer = 1.0
        bid_remove_counter = 0.0
        current_bid_index = 0
        current_bid = all_bids[current_bid_index]
        
        while(current_frame < total_frames):
            print "Frame",int(current_frame)
            file_name = "%04d" % current_frame +".jpg"
            self.PM.renderTransCPData(fileName=file_name,
                                         elev=current_elev,
                                         azim=current_azim,
                                         primaryWidth=6,
                                         dpi=200,
                                         showAxis=True,
                                         format='jpeg'
                                         )
            current_frame += 1
            current_azim += azim_increment
            current_elev += elev_increment

            bid_remove_counter += 1.0
            if(bid_remove_counter >= (bid_remove_rate*bid_remove_indexer)):
                # time to remove a bin!
                self.removeBinAndIndicies(current_bid)
                bid_remove_indexer+=1
                current_bid_index += 1
                if(current_bid_index < len(all_bids)):
                    current_bid = all_bids[current_bid_index]
                else:
                    return

    def plotBinProfiles(self):
        """Plot the distributions of kmer and coverage signatures"""
        self.BM.loadBins(makeBins=True,silent=False,bids=self.bids)
        print "Plotting bin profiles"
        self.BM.plotProfileDistributions()
    
    def plotPoints(self):
        """plot points"""
        self.BM.loadBins(makeBins=True,silent=False,bids=self.bids)
        self.BM.plotBinPoints()
    
    def plotSideBySide(self, coreCut):
        """Plot cores side by side with their contigs"""
        self.PM.loadData(condition="length >= "+str(coreCut))
        self.PM.transformCP()
        self.BM.loadBins(makeBins=True,bids=self.bids)
        print "Creating side by side plots"
        (bin_centroid_points, bin_centroid_colours, bin_ids) = self.BM.findCoreCentres()
        self.plotCoresVsContigs(bin_centroid_points, bin_centroid_colours)

    def plotIds(self):
        """Make a 3d plot of the bins but use IDs instead of points
        
        This function will help users know which bins to merge
        """
        self.BM.loadBins(makeBins=True,silent=False,bids=self.bids)
        self.BM.plotBinIds()

    def plotUnbinned(self, coreCut):
        """Plot all contigs over a certain length which are unbinned"""
        self.PM.plotUnbinned(coreCut)
            
#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def plotCoresVsContigs(self, binCentroidPoints, binCentroidColours):
        """Render the image for validating cores"""
        fig = plt.figure()
        ax1 = fig.add_subplot(121, projection='3d')
        ax1.scatter(self.PM.transformedCP[:,0], self.PM.transformedCP[:,1], self.PM.transformedCP[:,2], edgecolors=self.PM.contigColours, c=self.PM.contigColours, marker='.')
        ax2 = fig.add_subplot(122, projection='3d')
        ax2.scatter(binCentroidPoints[:,0], binCentroidPoints[:,1], binCentroidPoints[:,2], edgecolors=binCentroidColours, c=binCentroidColours)
        try:
            plt.show()
            plt.close(fig)
        except:
            print "Error showing image", sys.exc_info()[0]
            raise
        del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

###############################################################################
###############################################################################
###############################################################################
###############################################################################

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
import os
import sys
import errno

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

import numpy as np

# GroopM imports
import dataManagers
import mstore

np.seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class GMExtractor:
    """Used for extracting reads and contigs based on bin assignments"""
    def __init__(self, dbFileName,
                 bids=[],
                 folder='',
                 ):
        self.dbFileName = dbFileName
        
        if bids is None:
            self.bids = []
        else:
            self.bids = bids
        
        if(folder == ''):
            # write to current working dir
            self.outDir = os.getcwd()
        else:
            self.outDir = folder

        # make the dir if need be
        self.makeSurePathExists(self.outDir)
        
        
    def extractContigs(self, fasta=[], cutoff=0):
        """Extract contigs and write to file"""
        self.BM = dataManagers.BinManager(dbFileName=self.dbFileName)   # bins
        self.BM.loadBins(makeBins=True,silent=False,bids=self.bids)
        self.PM = self.BM.PM
        
############################
############################
        if(False):
            kse = mstore.KmerSigEngine()
            k_weights = kse.getKmerSigWeights()
            print kse.makeKmerColNames()
            sys.exit(0)
            # load all the contigs which have been assigned to bins
            CP = mstore.ContigParser()
            # contigs looks like cid->seq
            contigs = {}
            try:
                for file_name in fasta:
                    with open(file_name, "r") as f:  
                        contigs = CP.getWantedSeqs(f, self.PM.contigNames, storage=contigs)
            except:
                print "Could not parse contig file:",fasta[0],sys.exc_info()[0]
                raise
            
            lengths = np.array([])
            GCs = np.array([])
            dists = np.array([])
            sigs = np.array([])
             
            for bid in self.BM.getBids():
                bin = self.BM.getBin(bid)
                b_lengths = np.array([])
                b_GCs = np.array([])
                b_dists = np.array([])
                b_sigs = np.array([])
                b_names = np.array([])
                for row_index in bin.rowIndicies:
                    dist = np.linalg.norm(self.PM.covProfiles[row_index])
                    b_dists = np.append(b_dists, dist)
                    b_GCs = np.append(b_GCs, kse.getGC(contigs[self.PM.contigNames[row_index]]))
                    b_sigs = np.append(b_sigs, self.PM.kmerSigs[row_index])
                    b_lengths = np.append(b_lengths, self.PM.contigLengths[row_index])
                    b_names = np.append(b_names, self.PM.contigNames[row_index])
                    
                mean_dist = np.mean(b_dists)
                mean_length = np.mean(b_lengths)
                b_sigs = np.reshape(b_sigs, (bin.binSize, len(self.PM.kmerSigs[0])))
                for i in range(len(b_dists)):
                    print "C_"+b_names[i], b_dists[i]/mean_dist, b_GCs[i]
                    print "K_"+b_names[i], b_sigs[i] 
                    dists = np.append(dists, b_dists[i]/mean_dist)
                    GCs = np.append(GCs, b_GCs[i])
                    lengths = np.append(lengths, b_lengths[i])
    
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(dists, GCs, lengths, edgecolors='none', marker='.')
                            
            #fig = plt.figure()
            #plt.subplot(211)
            #plt.plot(dists, lengths, 'b.')
            #plt.subplot(212)
            #plt.plot(dists, GCs, 'r.')
            
            try:
                plt.show()
                plt.close(fig)
            except:
                print "Error showing image", sys.exc_info()[0]
                raise
            del fig
    
            return
#######################
#######################
        
        # load all the contigs which have been assigned to bins
        CP = mstore.ContigParser()
        # contigs looks like cid->seq
        contigs = {}
        try:
            for file_name in fasta:
                with open(file_name, "r") as f:  
                    contigs = CP.getWantedSeqs(f, self.PM.contigNames, storage=contigs)
        except:
            print "Could not parse contig file:",fasta[0],sys.exc_info()[0]
            raise
        
        # now print out the sequences
        print "Writing files"
        for bid in self.BM.getBids():
            file_name = os.path.join(self.outDir, "BIN_%d.fa" % bid)
            try:
                with open(file_name, 'w') as f: 
                    for row_index in self.BM.getBin(bid).rowIndicies:
                        cid = self.PM.contigNames[row_index]
                        if(cid in contigs):
                            f.write(">%s\n%s\n" % (cid, contigs[cid]))
                        else:
                            print "WTF", bid, cid
            except:
                print "Could not open file for writing:",file_name,sys.exc_info()[0]
                raise               
        
    def  extractReads(self, bams=[], shuffled=False):
        """Extract reads from sam files and write to file"""
        print "Soz LOL"
        return
        self.PM = dataManagers.ProfileManager(self.dbFileName)   # based on user specified length
        self.BM = dataManagers.BinManager(dbFileName=self.dbFileName)   # bins
        self.BM.loadBins(makeBins=True,silent=False,bids=self.bids)


    def makeSurePathExists(self, path):
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BinExplorer:
    """Inspect bins, used for validation"""
    def __init__(self, dbFileName, bids=[]):
        self.BM = dataManagers.BinManager(dbFileName=dbFileName)   # bins
        self.PM = self.BM.PM
        self.PM2 = None
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
        
        fig = plt.figure()
        while(current_frame < total_frames):
            print "Frame",int(current_frame)
            file_name = "%04d" % current_frame +".jpg"
            self.PM.renderTransCPData(fileName=file_name,
                                         elev=current_elev,
                                         azim=current_azim,
                                         primaryWidth=6,
                                         dpi=200,
                                         showAxis=True,
                                         format='jpeg',
                                         fig=fig
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
    
    def kmerSig2GC(self, sig_weightings, sig):
        """Calculate the GC content of a contig working from its kmer sig"""
        GC = 0
        outer_index = 0
        for s in sig:
            GC += s * sig_weightings[outer_index]
            outer_index += 1
        return GC
    
    
    def plotSideBySide(self, coreCut):
        """Plot cores side by side with their contigs"""
        self.PM2 = dataManagers.ProfileManager(dbFileName=self.BM.PM.dbFileName)
        self.PM2.loadData(condition="length >= "+str(coreCut))
        (min,max) = self.PM2.transformCP()
        self.BM.loadBins(makeBins=True,bids=self.bids, silent=False, min=min, max=max)
        print "Creating side by side plots"
        (bin_centroid_points, bin_centroid_colours, bin_ids) = self.BM.findCoreCentres()
        self.plotCoresVsContigs(bin_centroid_points, bin_centroid_colours)
        #self.plotCoresVsContigs(bin_centroid_points, bin_centroid_colours, azim=11, elev=43, fileName="forGene")

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

    def plotCoresVsContigs(self, binCentroidPoints, binCentroidColours, azim=0, elev=0, fileName='', dpi=300, format='png'):
        """Render the image for validating cores"""
        if(fileName==""):
            # plot on screen for user
            fig = plt.figure()
            ax1 = fig.add_subplot(121, projection='3d')
            ax1.scatter(self.PM2.transformedCP[:,0], self.PM2.transformedCP[:,1], self.PM2.transformedCP[:,2], edgecolors=self.PM2.contigColours, c=self.PM2.contigColours, marker='.')
            ax2 = fig.add_subplot(122, projection='3d')
            ax2.scatter(binCentroidPoints[:,0], binCentroidPoints[:,1], binCentroidPoints[:,2], edgecolors=binCentroidColours, c=binCentroidColours)
            try:
                plt.show()
                plt.close(fig)
            except:
                print "Error showing image", sys.exc_info()[0]
                raise
            del fig
        else:
            f_name1 = fileName + "_1"
            f_name2 = fileName + "_2"
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.PM2.transformedCP[:,0], self.PM2.transformedCP[:,1], self.PM2.transformedCP[:,2], edgecolors='none', c=self.PM2.contigColours, s=2, marker='.')
            ax.azim = azim
            ax.elev = elev
            ax.set_xlim3d(0,self.PM2.scaleFactor)
            ax.set_ylim3d(0,self.PM2.scaleFactor)
            ax.set_zlim3d(0,self.PM2.scaleFactor)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            try:
                fig.set_size_inches(12,12)            
                plt.savefig(f_name1,dpi=dpi,format=format)
                plt.close(fig)
            except:
                print "Error saving image",f_name1, sys.exc_info()[0]
                raise
            del fig
            
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            outer_index = 0
            for bid in self.BM.getBids():
                ax.text(binCentroidPoints[outer_index,0], 
                        binCentroidPoints[outer_index,1], 
                        binCentroidPoints[outer_index,2], 
                        str(int(bid)), 
                        color=binCentroidColours[outer_index]
                        )
                outer_index += 1
            
            ax.azim = azim
            ax.elev = elev
            ax.set_xlim3d(0,self.PM2.scaleFactor)
            ax.set_ylim3d(0,self.PM2.scaleFactor)
            ax.set_zlim3d(0,self.PM2.scaleFactor)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            try:
                fig.set_size_inches(12,12)            
                plt.savefig(f_name2,dpi=dpi,format=format)
                plt.close(fig)
            except:
                print "Error saving image",f_name1, sys.exc_info()[0]
                raise
            del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

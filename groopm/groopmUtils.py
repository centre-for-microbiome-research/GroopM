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
__copyright__ = "Copyright 2012/2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.2.9"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################
import os
import sys
import errno

import matplotlib.pyplot as plt
from colorsys import hsv_to_rgb as htr

import numpy as np

# GroopM imports
import binManager
import mstore

# other local imports
from bamtyper.utilities import BamParser as BTBP
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
        makeSurePathExists(self.outDir)

    def extractContigs(self, timer, fasta=[], cutoff=0):
        """Extract contigs and write to file"""
        self.BM = binManager.BinManager(dbFileName=self.dbFileName)   # bins
        self.BM.loadBins(timer, makeBins=True,silent=False,bids=self.bids, cutOff=cutoff)
        self.PM = self.BM.PM

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
            if self.BM.PM.isLikelyChimeric[bid]:
                file_name = os.path.join(self.outDir, "BIN_%d.chimeric.fa" % bid)
            else:
                file_name = os.path.join(self.outDir, "BIN_%d.fa" % bid)
            try:
                with open(file_name, 'w') as f:
                    for row_index in self.BM.getBin(bid).rowIndices:
                        cid = self.PM.contigNames[row_index]
                        if(cid in contigs):
                            f.write(">%s\n%s\n" % (cid, contigs[cid]))
                        else:
                            print "WTF", bid, cid
            except:
                print "Could not open file for writing:",file_name,sys.exc_info()[0]
                raise

    def  extractReads(self, timer, bams=[]):
        """Extract reads from sam files and write to file"""
        # load data
        self.BM = binManager.BinManager(dbFileName=self.dbFileName)   # bins
        self.BM.loadBins(timer, makeBins=True,silent=False,bids=self.bids)
        self.PM = self.BM.PM

        print "Extracting reads"

        # work out a set of targets to pass to the parser
        targets = {}
        bids = self.BM.getBids()
        for bid in bids:
            bin = self.BM.getBin(bid)
            for row_index in bin.rowIndices:
                targets[self.PM.contigNames[row_index]] = bid

        # get something to parse the bams with
        bam_parser = BTBP()
        bam_parser.extractReads(bams,
                                '',
                                targets,
                                combineBams=False,
                                headersOnly = True,
                                dontTrustSamFlags=False,
                                folder=self.outDir,
                                verbose=True
                                )

def makeSurePathExists(path):
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
    def __init__(self,
                 dbFileName,
                 bids=[],
                 transform=True,
                 cmstring="HSV",
                 ignoreContigLengths=False,
                 binLabelsFile = "",
                 contigColorsFile = ""):        
        self.ignoreContigLengths = ignoreContigLengths 
        self.transform = transform
        self.cmString = cmstring
        self.BM = binManager.BinManager(dbFileName=dbFileName)   # bins
        self.PM = self.BM.PM
        self.PM2 = None
        if bids is None:
            self.bids = []
        else:
            self.bids = bids

        # labels and colurs
        self.binLabelsFile = binLabelsFile
        self.contigColorsFile = contigColorsFile        
        self.LP = None 

    def plotHighlights(self,
                       timer,
                       elevation,
                       azimuth,
                       file,
                       filetype,
                       dpi,
                       show=False,
                       coreCut=1000):
        """Plot a high def image suitable for publication"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         loadContigLengths=True,
                         loadContigNames=True,
                         transform = self.transform,
                         cutOff=coreCut)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting image"
            if self.bids == []:
                self.bids = self.BM.getBids()

            # bids as labels and randomise colours
            self.LP = LabelParser(self.BM.getBids())
            self.LP.parseBinLabels(self.binLabelsFile)
            if self.contigColorsFile == "":
                self.LP.randomizeCols()
            else:
                self.LP.parseContigColorLabels(self.contigColorsFile, self.PM)

            foreground_disp_vals = []
            foreground_disp_cols = []
            foreground_disp_lens = []
            foreground_num = 0
            background_disp_vals = []
            background_disp_cols = []
            background_disp_lens = []
            background_num = 0

                
            # assign a color to each contig
            for row_index in np.arange(len(self.PM.indices)):
                # plot only binned contigs here
                bid = self.PM.binIds[row_index]
                if bid != 0:
                    if self.LP.loaded[bid]:
                        # this is one we want to show
                        if self.contigColorsFile == "":
                            # set the contigs to the colour of the bins
                            foreground_disp_cols.append(self.LP.bin2Cols[bid])
                        else:
                            # set the contigs to the user defined colour
                            foreground_disp_cols.append(self.LP.contig2Cols[row_index])
                        foreground_disp_vals.append(self.PM.transformedCP[row_index])
                        foreground_disp_lens.append(self.PM.contigLengths[row_index])
                        foreground_num += 1
                    else:
                        # this is one we want to background
                        if self.contigColorsFile == "":
                            # set the contigs to the colour of the bins
                            background_disp_cols.append(self.LP.unbinnedCol)
                            background_disp_vals.append(self.PM.transformedCP[row_index])
                            background_disp_lens.append(self.PM.contigLengths[row_index])
                            background_num += 1
                        elif self.LP.contig2Cols[row_index] != self.LP.unbinnedCol:
                            # set the contigs to the user defined colour, put this in the foreground plot too
                            foreground_disp_cols.append(self.LP.contig2Cols[row_index])
                            foreground_disp_vals.append(self.PM.transformedCP[row_index])
                            foreground_disp_lens.append(self.PM.contigLengths[row_index])
                            foreground_num += 1
                        else:
                            # set the contigs to the colour of the bins
                            background_disp_cols.append(self.LP.unbinnedCol)
                            background_disp_vals.append(self.PM.transformedCP[row_index])
                            background_disp_lens.append(self.PM.contigLengths[row_index])
                            background_num += 1

            foreground_disp_vals = np.reshape(foreground_disp_vals, (foreground_num, 3))
            background_disp_vals = np.reshape(background_disp_vals, (background_num, 3))
            
            fig = plt.figure()
            ax1 = plt.subplot(1,1,1, projection='3d')

            if background_num > 0:
                if self.ignoreContigLengths:
                    sc = ax1.scatter(background_disp_vals[:,0],
                                     background_disp_vals[:,1],
                                     background_disp_vals[:,2],
                                     edgecolors='none',
                                     c=background_disp_cols,
                                     s=10.,
                                     alpha=self.LP.unbinnedAlpha,
                                     marker='.')
                else:
                    sc = ax1.scatter(background_disp_vals[:,0],
                                     background_disp_vals[:,1],
                                     background_disp_vals[:,2],
                                     edgecolors='none',
                                     c=background_disp_cols,
                                     s=np.sqrt(background_disp_lens),
                                     alpha=self.LP.unbinnedAlpha,
                                     marker='.')

            
            if self.ignoreContigLengths:
                sc = ax1.scatter(foreground_disp_vals[:,0],
                                 foreground_disp_vals[:,1],
                                 foreground_disp_vals[:,2],
                                 edgecolors='none',
                                 c=foreground_disp_cols,
                                 s=10.,
                                 marker='.')
            else:
                sc = ax1.scatter(foreground_disp_vals[:,0],
                                 foreground_disp_vals[:,1],
                                 foreground_disp_vals[:,2],
                                 edgecolors='none',
                                 c=foreground_disp_cols,
                                 s=np.sqrt(foreground_disp_lens),
                                 marker='.')

            (bin_centroid_points, bin_centroid_colors, bin_centroid_gc, bin_ids) = self.BM.findCoreCentres(processChimeric=True)
            outer_index = 0
            for bid in self.bids:
                if self.LP.bin2Str[bid] != '':
                    bbox_props = dict(boxstyle="round4", fc="w", ec="none", alpha=0.6)
                    ax1.text(bin_centroid_points[outer_index, 0],
                             bin_centroid_points[outer_index, 1],
                             bin_centroid_points[outer_index, 2],
                             self.LP.bin2Str[bid],
                             color='black',
                             size='large',
                             ha="center", va="center", bbox=bbox_props
                             )
                outer_index += 1
            ax1.set_xticklabels([])
            ax1.set_yticklabels([])
            ax1.set_zticklabels([])
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_zticks([])
            
            if show:
                plt.show()
            elif not file.endswith(filetype):
                file += "." + filetype

                ax1.azim = azimuth
                ax1.elev = elevation
                plt.savefig(file,dpi=dpi,format=filetype)

            plt.close(fig)
            del fig

    def plotFlyOver(self, timer, fps=10.0, totalTime=120.0):
        """Plot a flyover of the data with bins being removed"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         bids=self.bids,
                         loadContigLengths=False,
                         loadContigNames=False,
                         transform = self.transform)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting flyover"
            self.BM.setColorMap(self.cmString)
            all_bids = self.BM.getBids()

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
            bid_remove_rate = float(len(all_bids)) / total_frames
            bids_removed = 0
            current_bid_index = 0
            current_bid = all_bids[current_bid_index]
            restricted_bids = []
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
                                             fig=fig,
                                             restrictedBids = restricted_bids
                                             )
                current_frame += 1
                current_azim += azim_increment
                current_elev += elev_increment
                print bid_remove_rate*current_frame, current_frame, "BR:",bids_removed, int(bid_remove_rate*current_frame)
                while bids_removed < int(bid_remove_rate*current_frame):
                    restricted_bids.append(all_bids[current_bid_index])
                    current_bid_index += 1
                    bids_removed += 1
            del fig

    def plotBinProfiles(self, timer):
        """Plot the distributions of kmer and coverage signatures"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         bids=self.bids,
                         transform=self.transform)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting bin profiles"
            self.BM.setColorMap(self.cmString)
            self.BM.plotProfileDistributions()

    def plotContigs(self, timer, coreCut, all=False):
        """plot contigs"""
        if all:
            print "Plotting all contigs"
            self.PM.plotAll(timer, coreCut, transform=self.transform, ignoreContigLengths=self.ignoreContigLengths)
        else:
            self.BM.loadBins(timer,
                             makeBins=True,
                             silent=False,
                             bids=self.bids,
                             transform=self.transform,
                             cutOff=coreCut)
            if len(self.BM.bins) == 0:
                print "Sorry, no bins to plot"
            else:
                print "Plotting binned contigs"
                self.BM.setColorMap(self.cmString)
                if self.bids == []:
                    self.bids = self.BM.getBids()
                self.BM.plotMultipleBins([self.bids], squash=True, ignoreContigLengths=self.ignoreContigLengths)

    def plotBinAssignents(self, timer, coreCut):
        """visualise bin assignments"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         getUnbinned=True,
                         silent=False,
                         bids=self.bids,
                         cutOff=coreCut,
                         transform=self.transform)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting bin assignments"
            if self.bids == []:
                self.bids = self.BM.getBids()

            # bids as labels and randomise colours
            self.LP = LabelParser(self.BM.getBids(), setBids=True)
            self.LP.randomizeCols()

            # for binned contigs
            b_disp_vals = []
            b_disp_cols = []
            b_disp_lens = []
            b_num = 0
            # for the unwashed masses
            ub_disp_vals = []
            ub_disp_cols = []
            ub_disp_lens = []
            ub_num = 0
                
            # assign a color to each contig
            for row_index in np.arange(len(self.PM.indices)):
                if self.PM.binIds[row_index] == 0:
                    ub_disp_cols.append(self.LP.unbinnedCol)
                    ub_disp_vals.append(self.PM.transformedCP[row_index])
                    ub_disp_lens.append(self.PM.contigLengths[row_index])
                    ub_num += 1
                else:
                    b_disp_cols.append(self.LP.bin2Cols[self.PM.binIds[row_index]])
                    b_disp_vals.append(self.PM.transformedCP[row_index])
                    b_disp_lens.append(self.PM.contigLengths[row_index])
                    b_num += 1

            b_disp_vals = np.reshape(b_disp_vals, (b_num, 3))
            fig = plt.figure()
            ax1 = plt.subplot(1,1,1, projection='3d')
            
            if ub_num != 0:
                ub_disp_vals = np.reshape(ub_disp_vals, (ub_num, 3))
                if self.ignoreContigLengths:
                    sc = ax1.scatter(ub_disp_vals[:,0],
                                     ub_disp_vals[:,1],
                                     ub_disp_vals[:,2],
                                     edgecolors='none',
                                     c=ub_disp_cols,
                                     s=10.,
                                     marker='.',
                                     alpha=self.LP.unbinnedAlpha)
                else:
                    sc = ax1.scatter(ub_disp_vals[:,0],
                                     ub_disp_vals[:,1],
                                     ub_disp_vals[:,2],
                                     edgecolors='none',
                                     c=ub_disp_cols,
                                     s=np.sqrt(ub_disp_lens),
                                     marker='.',
                                     alpha=self.LP.unbinnedAlpha)

            if self.ignoreContigLengths:
                sc = ax1.scatter(b_disp_vals[:,0],
                                 b_disp_vals[:,1],
                                 b_disp_vals[:,2],
                                 edgecolors='none',
                                 c=b_disp_cols,
                                 s=10.,
                                 marker='.')
            else:
                sc = ax1.scatter(b_disp_vals[:,0],
                                 b_disp_vals[:,1],
                                 b_disp_vals[:,2],
                                 edgecolors='none',
                                 c=b_disp_cols,
                                 s=np.sqrt(b_disp_lens),
                                 marker='.')

            (bin_centroid_points, bin_centroid_colors, bin_centroid_gc, bin_ids) = self.BM.findCoreCentres(processChimeric=True)
            outer_index = 0
            for bid in self.bids:
                if self.LP.bin2Str[bid] != '':
                    bbox_props = dict(boxstyle="round4", fc="w", ec=self.LP.bin2Cols[bid], alpha=0.6)
                    ax1.text(bin_centroid_points[outer_index, 0],
                             bin_centroid_points[outer_index, 1],
                             bin_centroid_points[outer_index, 2],
                             self.LP.bin2Str[bid],
                             color='black',
                             size='large',
                             ha="center", va="center", bbox=bbox_props
                             )
                outer_index += 1
            ax1.set_xticklabels([])
            ax1.set_yticklabels([])
            ax1.set_zticklabels([])
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_zticks([])
            
            plt.title(self.PM.dbFileName)

            plt.show()
            plt.close(fig)
            del fig

    def plotPoints(self, timer):
        """plot points"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         bids=self.bids,
                         transform=self.transform,
                         cutOff=coreCut)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting bin points"
            self.BM.setColorMap(self.cmString)
            self.BM.plotBinPoints()

    def plotCompare(self, timer, coreCut):
        """Plot cores side by side with their contigs"""
        self.PM2 = binManager.ProfileManager(dbFileName=self.BM.PM.dbFileName)
        self.PM2.loadData(timer,
                          "length >= "+str(coreCut),
                          bids=self.bids,
                          loadContigNames=False,
                          loadContigLengths=True,
                          )
        if self.transform:
            self.PM2.transformCP(timer)

        self.BM.loadBins(timer,
                         makeBins=True,
                         loadContigNames=False,
                         bids=self.bids,
                         silent=False,
                         transform=self.transform,
                         cutOff=coreCut)

        self.PM2.setColorMap(self.cmString)
        self.BM.setColorMap(self.cmString)
        
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting side by side graphs"
            (bin_centroid_points, bin_centroid_colors, bin_centroid_gc, bin_ids) = self.BM.findCoreCentres()
            self.plotCoresVsContigs(bin_centroid_points, bin_centroid_colors)

    def plotIds(self, timer):
        """Make a 3d plot of the bins but use IDs instead of points

        This function will help users know which bins to merge
        """
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         bids=self.bids,
                         transform=self.transform)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting bin IDs"
            self.BM.setColorMap(self.cmString)
            self.BM.plotBinIds()

    def plotUnbinned(self, timer, coreCut):
        """Plot all contigs over a certain length which are unbinned"""
        print "Plotting unbinned contigs"
        self.PM.plotUnbinned(timer, coreCut, transform=self.transform, ignoreContigLengths=self.ignoreContigLengths)

    def plotSideBySide(self, timer, coreCut):
        """Plot all bins separately in one large image"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         bids=self.bids,
                         transform=self.transform,
                         cutOff=coreCut)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            self.BM.setColorMap(self.cmString)
            self.BM.plotBins(sideBySide=True,
                             plotEllipsoid=True,
                             ignoreContigLengths=self.ignoreContigLengths)

    def plotTogether(self, timer, coreCut, doMers=False):
        """Plot all bins in ellipses on one normal image"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         bids=self.bids,
                         transform=self.transform,
                         cutOff=coreCut)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting all bins together"
            self.BM.setColorMap(self.cmString)
            if self.bids == []:
                p_bids = self.BM.getBids()
            else:
                p_bids = self.bids
            self.BM.plotSelectBins(p_bids,
                                   plotMers=doMers,
                                   plotEllipsoid=True,
                                   ignoreContigLengths=self.ignoreContigLengths)
         

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING

    def plotCoresVsContigs(self, binCentroidPoints, binCentroidColors, azim=0, elev=0, fileName='', dpi=300, format='png'):
        """Render the image for validating cores"""
        if(fileName==""):
            # plot on screen for user
            fig = plt.figure()
            ax1 = fig.add_subplot(121, projection='3d')
            if self.ignoreContigLengths:
                sc = ax1.scatter(self.PM2.transformedCP[:,0], self.PM2.transformedCP[:,1], self.PM2.transformedCP[:,2], edgecolors='none', c=self.PM2.contigGCs, cmap=self.PM2.colorMapGC, vmin=0.0, vmax=1.0, s=10, marker='.')
            else:
                sc = ax1.scatter(self.PM2.transformedCP[:,0], self.PM2.transformedCP[:,1], self.PM2.transformedCP[:,2], edgecolors='k', c=self.PM2.contigGCs, cmap=self.PM2.colorMapGC, vmin=0.0, vmax=1.0, s=np.sqrt(self.PM2.contigLengths), marker='.')                
            sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect

            ax2 = fig.add_subplot(122, projection='3d')
            sc = ax2.scatter(binCentroidPoints[:,0], binCentroidPoints[:,1], binCentroidPoints[:,2], edgecolors=binCentroidColors, c=binCentroidColors)
            sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect

            ax2.set_xlim(ax1.get_xlim())
            ax2.set_ylim(ax1.get_ylim())
            ax2.set_zlim(ax1.get_zlim())


            self.BM.plotStoitNames(ax1)
            self.BM.plotStoitNames(ax2)
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
            ax.scatter(self.PM2.transformedCP[:,0], self.PM2.transformedCP[:,1], self.PM2.transformedCP[:,2], edgecolors='none', c=self.PM2.contigGCs, cmap=self.PM2.colorMapGC, vmin=0.0, vmax=1.0, marker='.')
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
            self.BM.plotStoitNames(ax)

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
                        color=binCentroidColors[outer_index]
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
            self.BM.plotStoitNames(ax)

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

class LabelParser:

    def __init__(self, binIds, setBids=False, unbinnedCol=(0.2,0.2,0.2), unbinnedAlpha=0.2):
        self.unbinnedCol = unbinnedCol          # unassigned color
        self.unbinnedAlpha = unbinnedAlpha      # transparency for unbinned mofos
        self.bin2Cols = {0:self.unbinnedCol}    # map of bin ID to plotting colors
        self.contig2Cols = {}                   # map of rowIndices to contig colours
        self.bin2Str = {}                       # map of bin ID to plot labels
        self.loaded = {}
        if setBids:
            for bid in binIds:
                self.bin2Str[bid] = str(bid)        # labels are easy!
                self.loaded[bid] = True
        else:
            for bid in binIds:
                self.bin2Str[bid] = ''              # labels are easy!
                self.loaded[bid] = False

        # color conversions             
        self._NUMERALS = '0123456789abcdefABCDEF'
        self._HEXDEC = {}
        for v in (x+y for x in self._NUMERALS for y in self._NUMERALS):
            self._HEXDEC[v] = int(v, 16)     

    def rgb(self, triplet):
        """Hex triplet to rgb"""
        triplet = triplet.replace("#", "")
        return (float(self._HEXDEC[triplet[0:2]])/255.,
                float(self._HEXDEC[triplet[2:4]])/255.,
                float(self._HEXDEC[triplet[4:6]])/255.)

    def parseContigColorLabels(self, labelFileName, PM, nullCol=(0.2,0.2,0.2)):
        """parse labels file containing bin text
        
        each line of the label file looks like:
        
        contigID<tab>#colour
        
        If a contig is omitted from this list then it gets coloured nullCol
        """
        # first we need to make a map of contig name to rowIndex
        name_2_row_index = {}
        for row_index in range(len(PM.indices)):
            name_2_row_index[PM.contigNames[row_index]] = row_index

        # first we parse the file
        try:
            with open(labelFileName, "r") as l_fh:
                for line in l_fh:
                    fields = line.rstrip().split("\t")
                    cid = fields[0]
                    try:
                        self.contig2Cols[name_2_row_index[cid]] = self.rgb(fields[1])
                    except KeyError:
                        print "ERROR: contig name %s not recognised" % cid
        except:
            print "ERROR: parsing labels file: %s" % labelFileName
            raise

        # now we parse the reast of the contig names and colour the null colour
        for row_index in range(len(PM.indices)):
            if row_index not in self.contig2Cols:
                self.contig2Cols[row_index] = self.unbinnedCol
            
    def parseBinLabels(self, labelFileName):
        """parse labels file containing bin text
        
        each line of the label file looks like:
        
        binID<tab>label<tab>[BG]
        
        If label is ommitted then the bin gets no label
        If label is '-', then the bin gets it's number as a label
        
        BG is assumed to be False, if it is "BG" the the bin will be
        placed in the background (ie, alpha set...)
        
        lines can be commented out using "#"
        """
        try:
            with open(labelFileName, "r") as l_fh:
                for line in l_fh:
                    if line[0] != "#":
                        fields = line.rstrip().split("\t")
                        bid = int(fields[0])
                        loaded_bin = True
                        # grab the label if any
                        try:
                            if fields[1] == '-':
                                self.bin2Str[bid] = fields[0]
                            else:
                                self.bin2Str[bid] = fields[1]
                            # is backgrounded?
                            try:
                                if fields[2] == 'BG':
                                    loaded_bin = False
                            except IndexError: pass
                        except IndexError: pass
                        self.loaded[bid] = loaded_bin
        except:
            print "ERROR parsing labels file: %s" % labelFileName
            raise
            
    def randomizeCols(self):
        """choose colors randomly"""
        S = 1.0
        V = 1.0
        num_bins = len(self.bin2Str)
        offset = 0.5
        step = 1. /(num_bins-1 + 2. * offset)
        Hs = np.array([step*(i + offset) for i in range(num_bins)])
        cols = [htr(H, S, V) for H in Hs]
        np.random.shuffle(cols)
        i = 0
        for bid in self.bin2Str.keys():
            if self.loaded[bid]:
                # assign the color we picked
                self.bin2Cols[bid] = cols[i]
            else:
                # set to the unbinned color
                self.bin2Cols[bid] = self.unbinnedCol
            i += 1

###############################################################################
###############################################################################
###############################################################################
###############################################################################


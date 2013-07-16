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
__version__ = "0.2.5"
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
        self.BM.loadBins(timer, makeBins=True,silent=False,bids=self.bids)
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
                 ignoreContigLengths=False):
        
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

    def plotHighlights(self, timer, bids, elevation, azimuth, file, filetype, dpi, alpha, invert=False, show=False):
        """Plot a high def image suitable for publication"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         loadContigLengths=True,
                         loadContigNames=False,
                         transform = self.transform)
        if len(self.BM.bins) == 0:
            print "Sorry, no bins to plot"
        else:
            print "Plotting image"
            self.BM.setColorMap(self.cmString)
            fig = plt.figure()
            bins=[]
            if bids is not None:
                for bid in bids:
                    bins.append(self.BM.getBin(bid))

            if show:
                file=""
            elif not file.endswith(filetype):
                file += "." + filetype

            self.PM.renderTransCPData(fileName=file,
                                      elev=elevation,
                                      azim=azimuth,
                                      primaryWidth=6,
                                      dpi=dpi,
                                      showAxis=True,
                                      format=filetype,
                                      fig=fig,
                                      highlight=bins,
                                      alpha=alpha
                                      )
            del fig

            if invert:
                # invert the colors
                from PIL import Image
                import PIL.ImageOps
                image = Image.open(file)
                inverted_image = PIL.ImageOps.invert(image)
                inverted_image.save(file)

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
                             transform=self.transform)
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

            # work out colors per bin
            S = 1.0
            V = 1.0
            num_bins = len(self.bids)
            offset = 0.5
            step = 1. /(num_bins-1 + 2. * offset)
            Hs = np.array([step*(i + offset) for i in range(num_bins)])
            cols = [htr(H, S, V) for H in Hs]
            np.random.shuffle(cols)
            unbinned_col = (0.2,0.2,0.2)
            bin_cols = {0:(0.2,0.2,0.2)}    # unassigned color
            i = 0
            for bid in self.bids:
                bin_cols[bid] = cols[i]
                i += 1

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
                    ub_disp_cols.append(unbinned_col)
                    ub_disp_vals.append(self.PM.transformedCP[row_index])
                    ub_disp_lens.append(self.PM.contigLengths[row_index])
                    ub_num += 1
                else:
                    b_disp_cols.append(bin_cols[self.PM.binIds[row_index]])
                    b_disp_vals.append(self.PM.transformedCP[row_index])
                    b_disp_lens.append(self.PM.contigLengths[row_index])
                    b_num += 1

            b_disp_vals = np.reshape(b_disp_vals, (b_num, 3))
            fig = plt.figure()
            #import matplotlib.gridspec as gridspec
            #gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
            #ax1 = plt.subplot(gs[0], projection='3d')
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
                                     alpha=0.35)
                else:
                    sc = ax1.scatter(ub_disp_vals[:,0],
                                     ub_disp_vals[:,1],
                                     ub_disp_vals[:,2],
                                     edgecolors='none',
                                     c=ub_disp_cols,
                                     s=np.sqrt(ub_disp_lens),
                                     marker='.',
                                     alpha=0.35)

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
                bbox_props = dict(boxstyle="round4", fc="w", ec=bin_cols[bid], alpha=0.6)
                ax1.text(bin_centroid_points[outer_index, 0],
                         bin_centroid_points[outer_index, 1],
                         bin_centroid_points[outer_index, 2],
                         str(int(bid)),
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
            
            if False:
                #plot bin keys
                ax2 = plt.subplot(gs[1])
                
                # work out text placement
                bin_rows = int(np.sqrt(num_bins)) + 1
                bin_columns = int(float(num_bins)/float(bin_rows)) + 2
                bin_rows += 1 
                ax2.set_xlim(0, bin_columns)
                ax2.set_ylim(0, bin_rows)
                bid_coords = zip(np.array([[i]*bin_columns for i in range(bin_rows)]).flatten(), range(bin_columns) * bin_rows)
    
                outer_index = 0
                
                for bid in self.bids:
                    bbox_props = dict(boxstyle="round", fc=bin_cols[bid])
                    ax2.text(bid_coords[outer_index][0] +1,
                             bid_coords[outer_index][1] +1,
                             str(int(bid)),
                             ha="center", va="center",
                             bbox=bbox_props,
                             size='medium' 
                             )
                    outer_index += 1
    
                plt.xticks([], [])
                plt.yticks([], [])

            plt.show()
            plt.close(fig)
            del fig


    def plotPoints(self, timer):
        """plot points"""
        self.BM.loadBins(timer,
                         makeBins=True,
                         silent=False,
                         bids=self.bids,
                         transform=self.transform)
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
                          condition="length >= "+str(coreCut),
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
                         transform=self.transform)

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
                         transform=self.transform)
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
                         transform=self.transform)
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

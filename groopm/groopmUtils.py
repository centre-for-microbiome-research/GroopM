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
__version__ = "0.2.10"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################
import os
import sys
import errno

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, axes3d, proj3d

from colorsys import hsv_to_rgb as htr

import numpy as np
from scipy.spatial.distance import cdist, squareform

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
                       drawRadius=False,
                       show=False,
                       coreCut=1000,
                       testing=False):
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

            if testing:
                # ignore labelling files provided
                self.binLabelsFile = "none"
                raw_input( "****************************************************************\n"
                           " IMAGE MAKING INSTRUCTIONS - PLEASE READ CAREFULLY\n"
                           "****************************************************************\n"
                           " You are using GroopM in highlight mode. Congratulations!\n"
                           " A 3D plot of your contigs coloured by bin assignments will be\n"
                           " displayed shortly. Rotate image until you like what you see \n"
                           " taking note of the 'azimuth' and 'elevation' values displayed in\n"
                           " the lower left hand corner of the plotting window.\n"
                           " When you have decided how your image should be arranged, rerun\n"
                           " highlight mode without the -P parameter, ste the -a and -e\n"
                           " parameters to what you saw here, set bin labels, contig colours...\n\n"
                           " Good Luck!\n\n"
                           " Press return to continue...")
                print "****************************************************************"

            # bids as labels and randomise colours
            self.LP = LabelParser(self.BM.getBids())
            if self.binLabelsFile == "none":
                # force no bin labelling
                draw_circ = False
            elif self.binLabelsFile.lower() == '':
                # set default bin labels
                self.LP.setDefaultBinLabels(self.bids)
                draw_circ = True and drawRadius
            else:
                # an actual file
                self.LP.parseBinLabels(self.binLabelsFile)
                draw_circ = True and drawRadius

            if self.contigColorsFile.lower() == 'none':
                # force no coloring
                self.LP.setAllGrey(self.PM)
            elif self.contigColorsFile == "":
                # assign all contigs belonging to a bin the same randomly selected colour
                self.LP.randomizeCols(setLoaded=True)
            else:
                # use the colors assigned by the user
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

            if self.ignoreContigLengths:
                s_f = 10.
                s_b = 10.
            else:
                s_f = np.sqrt(foreground_disp_lens)
                s_b = np.sqrt(background_disp_lens)

            # get a plotter
            L3P = Labelled3DPlotter(azim=azimuth, elev=elevation, dpi=dpi, drawCirc=draw_circ)
            ax1 = L3P.getInnerAxis()

            if background_num > 0:
                    sc = ax1.scatter(background_disp_vals[:,0],
                                     background_disp_vals[:,1],
                                     background_disp_vals[:,2],
                                     edgecolors='none',
                                     c=background_disp_cols,
                                     s=s_b,
                                     alpha=self.LP.unbinnedAlpha,
                                     marker='.')

            sc = ax1.scatter(foreground_disp_vals[:,0],
                             foreground_disp_vals[:,1],
                             foreground_disp_vals[:,2],
                             edgecolors='none',
                             c=foreground_disp_cols,
                             s=s_f,
                             marker='.')

            # add text labels, make them purdy in a circle and outside of the plot
            num_bins = len(self.bids)
            (bin_centroid_points, bin_centroid_colors, bin_centroid_gc, bin_ids) = self.BM.findCoreCentres(processChimeric=True)
            # keep only those that are 'loaded'
            bin_centroid_points = np.array([bin_centroid_points[i] for i in range(num_bins) if self.LP.loaded[bin_ids[i]]])
            bin_ids = np.array([bin_ids[i] for i in range(num_bins) if self.LP.loaded[bin_ids[i]]])
            num_bins = len(bin_ids)
            # no point going on with labeling crud if there are no bin labels to display
            if num_bins > 0:
                for ii in np.arange(len(bin_ids)):
                    bid = bin_ids[ii]
                    if self.LP.bin2Str[bid] != '':
                        L3P.labelPoint(bin_centroid_points[ii][0], bin_centroid_points[ii][1], bin_centroid_points[ii][2], self.LP.bin2Str[bid])


            if show:
                L3P.renderImage()
            else:
                L3P.renderImage(file, filetype)

    def flat2square(self,idx,side):
        """return a flat index from argmin into square coords"""
        row = int(idx/side)
        col = idx - (row*side)
        return (row,col)

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

            if self.ignoreContigLengths:
                s_ub = 10.
                s_b = 10.
            else:
                s_ub = np.sqrt(ub_disp_lens)
                s_b = np.sqrt(b_disp_lens)

            if ub_num != 0:
                ub_disp_vals = np.reshape(ub_disp_vals, (ub_num, 3))
                sc = ax1.scatter(ub_disp_vals[:,0],
                                 ub_disp_vals[:,1],
                                 ub_disp_vals[:,2],
                                 edgecolors='none',
                                 c=ub_disp_cols,
                                 s=s_ub,
                                 marker='.',
                                 alpha=self.LP.unbinnedAlpha)

            sc = ax1.scatter(b_disp_vals[:,0],
                             b_disp_vals[:,1],
                             b_disp_vals[:,2],
                             edgecolors='none',
                             c=b_disp_cols,
                             s=s_b,
                             marker='.')

            (bin_centroid_points, bin_centroid_colors, bin_centroid_gc, bin_ids) = self.BM.findCoreCentres(processChimeric=True)

            # restrict to the ones we've loaded
            all_bids = self.BM.getBids()
            bin_centroid_points = np.array([bin_centroid_points[i] for i in range(len(all_bids)) if all_bids[i] in self.bids])
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

class Labelled3DPlotter:
    """Class for producing #D plots with labels situated in a circle around the outside"""

    def __init__(self, radiusPercent=0.35, azim=-60, elev=30, dpi=300, side=4., cubeSide=1000., drawCirc=True):
        # the size and center of the cube we're in
        self.cubeSide = cubeSide
        self.cubeCenter = self.cubeSide/2.
        self.radius = radiusPercent
        self.drawCirc = drawCirc

        # 3d size relative to background
        self.threed_size = 0.65
        self.threed_buffer = (1. - self.threed_size)/2.

        self.fig = None
        self.dpi = dpi
        self.side = side

        self.initAxes(azim=azim, elev=elev)

    def initAxes(self, azim=-60, elev=30):
        """make axes for text and line plotting"""
        if self.fig is not None:
            plt.close(self.fig)
            del self.fig
        self.fig = plt.figure(dpi=self.dpi)
        plt.axis('equal')

        # make the background axis which stores lines and text
        self.screenAxis = plt.axes([0.,0.,self.cubeSide,self.cubeSide])
        self.screenAxis.set_xlim(0.,self.cubeSide)
        self.screenAxis.set_ylim(0.,self.cubeSide)
        self.screenAxis.set_xticklabels([])
        self.screenAxis.set_yticklabels([])
        self.screenAxis.set_xticks([])
        self.screenAxis.set_yticks([])

        # create the inner 3D axis and align it's centre to the centre of
        # the background plot
        self.threeDAx = plt.axes([self.threed_buffer,
                                  self.threed_buffer,
                                  self.threed_size,
                                  self.threed_size],
                                 projection='3d')

        self.threeDAx.set_xlim3d(0.,self.cubeSide)
        self.threeDAx.set_ylim3d(0.,self.cubeSide)
        self.threeDAx.set_zlim3d(0.,self.cubeSide)
        self.threeDAx.set_xticklabels([])
        self.threeDAx.set_yticklabels([])
        self.threeDAx.set_zticklabels([])
        self.threeDAx.set_xticks([])
        self.threeDAx.set_yticks([])
        self.threeDAx.set_zticks([])
        self.threeDAx.w_xaxis.line.set_color("white")
        self.threeDAx.w_yaxis.line.set_color("white")
        self.threeDAx.w_zaxis.line.set_color("white")

        # set the angle and azim like a boss
        self.threeDAx.azim = azim
        self.threeDAx.elev = elev

        # make the upper axes invisimable
        self.threeDAx.patch.set_visible(False)

        # get image information to use when converting to screen coords
        self.dpi = self.screenAxis.figure.get_dpi()
        self.twoDImgHeight = self.screenAxis.figure.get_figheight() * self.dpi
        self.twoDImgWidth = self.screenAxis.figure.get_figwidth() * self.dpi


        # get the self.screenAxis coords of the centre of the 3d plot
        [self.flatCubeCenter_x, self.flatCubeCenter_y] = self.flatten3DPoint(self.cubeCenter,
                                                                             self.cubeCenter,
                                                                             self.cubeCenter)
        # plot the radius circle to help things along
        if self.drawCirc:
            rad_circle = plt.Circle((self.flatCubeCenter_x, self.flatCubeCenter_y), radius=self.radius, color="#AAAAAA", fill=False)
            self.screenAxis.add_artist(rad_circle)

    def getInnerAxis(self):
        """return the 3D inner axis"""
        return self.threeDAx

    def makeRandomData(self, numPoints, numSelected=None):
        """make some data points

        useful for testing"""
        x = random.rand(numPoints) * self.cubeSide
        y = random.rand(numPoints) * self.cubeSide
        z = random.rand(numPoints) * self.cubeSide

        if numSelected is not None:
            selected = np.arange(numPoints)
            random.shuffle(selected)
            selected = selected[:numSelected]
            x_selected = x[selected]
            y_selected = y[selected]
            z_selected = z[selected]

            return ((x,y,z),(x_selected, y_selected, z_selected))

        return (x,y,z)

    def flatten3DPoint(self, px, py, pz):
        """project a point in 3D space onto a 2D screen"""
        # convert 3D spave to 2D space
        x2, y2, _ = proj3d.proj_transform(px, py, pz, self.threeDAx.get_proj())
        # convert 2d space to screen space
        [x2,y2] = self.threeDAx.transData.transform((x2, y2))
        # adjust for image dimensions
        x2 = x2/self.twoDImgWidth
        y2 = y2/self.twoDImgHeight
        return [x2, y2]

    def labelPoint(self, px, py, pz, label=None):
        """find the 2d point which is some 2D radial distance from the 2D projected center

        The line from this point to the 2D projected center
        passes through the given point"""

        # flatten the point
        [x2d, y2d] = self.flatten3DPoint(px, py, pz)

        # get the angle around the flattened 2d cube centre
        angle = np.arctan2(y2d-self.flatCubeCenter_y, x2d-self.flatCubeCenter_x)

        # get the extended coords
        xl = self.radius*np.cos(angle) + self.flatCubeCenter_x
        yl = self.radius*np.sin(angle) + self.flatCubeCenter_y

        transFigure = self.fig.transFigure.inverted()
        coord1 = transFigure.transform(self.screenAxis.transData.transform([x2d,y2d]))
        coord2 = transFigure.transform(self.screenAxis.transData.transform([xl,yl]))
        self.fig.lines += self.screenAxis.plot((coord1[0],coord2[0]),(coord1[1],coord2[1]), linestyle=':', c="#999999" )
        if label is None:
            label_text = "(%0.1f,%0.1f,%0.1f)" % (px, py, pz)
        else:
            label_text = label

        if coord2[1] < self.flatCubeCenter_y:
            va = "top"
        elif coord2[1] == self.flatCubeCenter_y:
            va = "center"
        else:
            va = "bottom"

        if coord2[0] > self.flatCubeCenter_x:
            ha = "left"
        elif coord2[0] == self.flatCubeCenter_x:
            ha = "center"
        else:
            ha = "right"


        self.screenAxis.text(coord2[0],
                             coord2[1],
                             label_text,
                             color='black',
                             size='small',
                             ha=ha,
                             va=va)

    def renderImage(self, fileName="",format="svg"):
        """Render the plot and save to file"""
        if fileName == "":
            plt.show()
        else:
            self.fig.set_size_inches(self.side, self.side)
            plt.savefig(fileName+"."+format,dpi=self.dpi,format=format)

        plt.close(self.fig)
        del self.fig



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

    def setAllGrey(self, PM):
        """set all contrigs to the unbinned colour"""
        for row_index in range(len(PM.indices)):
            self.contig2Cols[row_index] = self.unbinnedCol

    def parseBinLabels(self, labelFileName):
        """parse labels file containing bin text

        each line of the label file looks like:

        binID<tab>label

        If label is ommitted then the bin gets no label

        lines can be commented out using "#"
        """
        try:
            with open(labelFileName, "r") as l_fh:
                for line in l_fh:
                    if line[0] != "#":
                        fields = line.rstrip().split("\t")
                        bid = int(fields[0])
                        # grab the label if any
                        try:
                            if fields[1] == '':
                                self.bin2Str[bid] = fields[0]
                            else:
                                self.bin2Str[bid] = fields[1]
                        except IndexError: pass
                        self.loaded[bid] = True
        except:
            print "ERROR parsing labels file: %s" % labelFileName
            raise

    def setDefaultBinLabels(self, bids):
        """Set bin labels to the GM bin ID"""
        for bid in bids:
            self.bin2Str[bid] = str(bid)
            self.loaded[bid] = True

    def randomizeCols(self, setLoaded=False):
        """choose colors randomly"""
        S = 1.0
        V = 1.0
        if setLoaded:
            for bid in self.bin2Str.keys():
                self.loaded[bid] = True
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


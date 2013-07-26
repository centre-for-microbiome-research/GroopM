#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bin.py                                                                   #
#                                                                             #
#    Bins Bins Bins                                                           #
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
__version__ = "0.2.4"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Release"

###############################################################################

import sys

import matplotlib.pyplot as plt
from pylab import show

import numpy as np
from numpy import (around as np_around,
                   array as np_array,
                   mean as np_mean,
                   median as np_median,
                   std as np_std)

from ellipsoid import EllipsoidTool
from groopmExceptions import ModeNotAppropriateException

np.seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Bin:
    """Class for managing collections of contigs

    To (perhaps) simplify things think of a "bin" as an row_index into the
    column names array. The ClusterBlob has a list of bins which it can
    update etc...
    """
    def __init__(self, rowIndices, id, upperCov, covtol=2, mertol=2, gctol=2):
        self.id = id
        self.rowIndices = rowIndices             # all the indices belonging to this bin
        self.binSize = self.rowIndices.shape[0]
        self.upperCov = upperCov
        self.totalBP = 0

        self.covTolerance = covtol
        self.kValTolerance = mertol
        self.gcTolerance = gctol

        # COVERAGE (3D COVERAGE VALUES)
        self.covMedians = np.zeros((3))
        self.covStdevs = np.zeros((3))
        self.covLowerLimits = np.zeros((3)) # lower and upper limits based on tolerance
        self.covUpperLimits = np.zeros((3))

        # AVERAGE COVERAGE
        self.cValMedian = 0.0
        self.cValStdev = 0.0
        self.cValUpperLimit = 0.0
        self.cValLowerLimit = 0.0

        # KMER VALUES for ALL PCs
        self.kMedian = None
        self.kStdev = None

        # NORMALIZED PC1 KMER VALUES
        self.kValMeanNormPC1 = 0.0
        self.kValStdevNormPC1  = 0.0

        # GC VALUES
        self.gcMedian = 0.0
        self.gcStdev = 0.0
        self.gcUpperLimit = 0.0
        self.gcLowerLimit = 0.0

        # contig lengths
        self.lengthMean = 0.0
        self.lengthStd = 0.0

#------------------------------------------------------------------------------
# Tools used for comparing / condensing

    def __cmp__(self, alien):
        """Sort bins based on the normalized first PC of kmer signatures."""
        if self.kValMeanNormPC1 < alien.kValMean:
            return -1
        elif self.kValMeanNormPC1 == alien.kValMean:
            return 0
        else:
            return 1

#------------------------------------------------------------------------------
# Grow and shrink

    def consume(self, transformedCP, averageCoverages, kmerNormPC1, kmerPCs, contigGCs, contigLengths, deadBin, verbose=False):
        """Combine the contigs of another bin with this one"""
        # consume all the other bins rowIndices
        if(verbose):
            print "    BIN:",deadBin.id,"will be consumed by BIN:",self.id
        self.rowIndices = np.concatenate([self.rowIndices, deadBin.rowIndices])
        self.binSize  = self.rowIndices.shape[0]

        # fix the stats on our bin
        self.makeBinDist(transformedCP, averageCoverages, kmerNormPC1, kmerPCs, contigGCs, contigLengths)

    def scoreProfile(self, kmerVal, transformedCP):
        """Determine how similar this profile is to the bin distribution

        This is the norm of the vector containing z distances for both profiles
        """
        #print self.covStdevs, self.binSize
        covZ = np.abs(np.mean(np.abs(transformedCP - self.covMedians)/self.covStdevs))
        merZ = np.abs(kmerVal - self.kValMeanNormPC1)/self.kValStdevNormPC1
        return (covZ,merZ)

    def purge(self, deadIndices, transformedCP, averageCoverages, kmerNormPC1, kmerPCs, contigGCs, contigLengths):
        """Delete some rowIndices and remake stats"""
        old_ri = self.rowIndices
        self.rowIndices = np.array([])
        for i in old_ri:
            if i not in deadIndices:
                self.rowIndices = np.append(self.rowIndices, i)

        # fix the stats on our bin
        self.makeBinDist(transformedCP, averageCoverages, kmerNormPC1, kmerPCs, contigGCs, contigLengths)

#------------------------------------------------------------------------------
# Stats and properties

    def clearBinDist(self):
        """Clear any set distribution statistics"""
        self.totalBP = 0

        self.covMedians = np.zeros((3))
        self.covStdevs = np.zeros((3))
        self.covLowerLimits = np.zeros((3)) # lower and upper limits based on tolerance
        self.covUpperLimits = np.zeros((3))

        self.kMedian = None
        self.kStdevs = None

        self.kValMeanNormPC1 = 0.0
        self.kValStdevNormPC1 = 0.0

        self.gcMedian = 0.0
        self.gcStdev = 0.0
        self.gcUpperLimit = 0.0
        self.gcLowerLimit = 0.0

    def makeBinDist(self, transformedCP, averageCoverages, kmerNormPC1, kmerPCs, contigGCs, contigLengths):
        """Determine the distribution of the points in this bin

        The distribution is largely normal, except at the boundaries.
        """
        #print "MBD", self.id, self.binSize
        self.binSize = self.rowIndices.shape[0]
        if(0 == np.size(self.rowIndices)):
            return

        # get the centroids
        (self.covMedians, self.covStdevs) = self.getCentroidStats(transformedCP)
        (self.lengthMean, self.lengthStd) = self.getCentroidStats(contigLengths)

        self.kValMeanNormPC1 = np_median(kmerPCs[self.rowIndices])
        self.kValStdevNormPC1 = np_std(kmerPCs[self.rowIndices])

        self.kMedian = np_median(kmerPCs[self.rowIndices], axis=0)
        self.kStdevs = np_std(kmerPCs[self.rowIndices], axis=0)

        cvals = self.getAverageCoverageDist(averageCoverages)
        self.cValMedian = np_around(np_median(cvals), decimals=3)
        self.cValStdev = np_around(np_std(cvals), decimals=3)

        self.gcMedian = np_median(contigGCs[self.rowIndices])
        self.gcStdev = np_std(contigGCs[self.rowIndices])

        # work out the total size
        self.totalBP = sum([contigLengths[i] for i in self.rowIndices])

        # set the acceptance ranges
        self.makeLimits()

    def makeLimits(self):
        """Set inclusion limits based on mean, variance and tolerance settings"""
        covTol=self.covTolerance
        gcTol=self.gcTolerance

        for i in range(0,3):
            self.covLowerLimits[i] = int(self.covMedians[i] - covTol * self.covStdevs[i])
            if(self.covLowerLimits[i] < 0):
                self.covLowerLimits[i] = 0.0
            self.covUpperLimits[i] = int(self.covMedians[i] + covTol * self.covStdevs[i])
            if(self.covUpperLimits[i] > self.upperCov):
                self.covUpperLimits[i] = self.upperCov

        self.gcLowerLimit = self.gcMedian - gcTol * self.gcStdev
        if(self.gcLowerLimit < 0):
            self.gcLowerLimit = 0
        self.gcUpperLimit = self.gcMedian + gcTol * self.gcStdev

        self.cValLowerLimit = self.cValMedian - covTol * self.cValStdev
        if(self.cValLowerLimit < 0):
            self.cValLowerLimit = 0
        self.cValUpperLimit = self.cValMedian + covTol * self.cValStdev

    def getCentroidStats(self, profile):
        """Calculate the centroids of the profile"""
        working_list = profile[self.rowIndices]
        
        # return the mean and stdev
        # we divide by std so we need to make sure it's never 0
        tmp_stds = np_std(working_list, axis=0)
        mean_std = np_mean(tmp_stds)
        try:
            std = np_array([x if x != 0 else mean_std for x in tmp_stds])
        except:
            std = mean_std
        return (np_median(working_list,axis=0), std)
    
    def getkmerValDist(self, kmerNormPC1):
        """Return an array of kmer vals for this bin"""
        return np.array([kmerNormPC1[i] for i in self.rowIndices])

    def getGC_Dist(self, GCs):
        """Return an array of GCs for this bin"""
        return np.array([GCs[i] for i in self.rowIndices])

    def getAverageCoverageDist(self, averageCoverages):
        """Return the average coverage for all contigs in this bin"""
        return np.array([averageCoverages[i] for i in self.rowIndices])

    def getAverageTransformedCoverageDist(self, coverages):
        """Return the average transformed coverage for all contigs in this bin"""
        return np.array([np.mean(coverages[i]) for i in self.rowIndices])

    def getInnerVariance(self, profile, mode="kmer"):
        """Work out the variance for the coverage/kmer/gc profile"""
        dists = []
        if(mode == "kmer"):
            dists = [np.abs(self.kValMeanNormPC1 - profile[i]) for i in self.rowIndices]
        elif(mode =="cov"):
            dists = [self.getCDist(profile[i]) for i in self.rowIndices]
        elif(mode =="gc"):
            dists = [np.abs(self.gcMedian - profile[i]) for i in self.rowIndices]
        else:
            raise ModeNotAppropriateException("Mode",mode,"unknown")
        
        dist_range = np.max(np.array(dists)) - np.min(np.array(dists))
        return (np.mean(np.array(dists)), np.std(np.array(dists)), dist_range)

    def getCDist(self, Csig, centroid=None):
        """Get the distance of this contig from the coverage centroid"""
        # z-norm and then distance!
        if centroid is None:
            centroid = self.covMedians
        return np.linalg.norm(Csig-centroid)

    def getBoundingEllipsoid(self, transformedCP, ET=None, retA=False):
        """Return the minimum bounding ellipsoid

        returns (center, radii, rotation) or (A, center, radii, rotation)
        """
        bin_points = np.array([transformedCP[i] for i in self.rowIndices])
        if len(bin_points) > 1:
            if ET is None:
                ET = EllipsoidTool()
            try:
                return ET.getMinVolEllipse(bin_points, retA=retA)
            except:
                print bin_points
                raise
        else: # minimum bounding ellipse of a point is 0
            if retA:
                return (np.zeros((3,3)), transformedCP[self.rowIndices[0]], np.zeros((3)), np.eye(3))
            else:
                return (transformedCP[self.rowIndices[0]], np.zeros((3)), np.eye(3))

    def getBoundingCEllipsoidVol(self, transformedCP, ET=None, retA=False):
        """Return the volume of the minimum bounding coverage ellipsoid"""
        if ET is None:
            ET = EllipsoidTool()
        (A, center, radii, _rotation) = self.getBoundingEllipsoid(transformedCP, ET=ET, retA=True)
        if retA:
            return ((A, center), ET.getEllipsoidVolume(radii))
        else:
            return ET.getEllipsoidVolume(radii)

    def getBoundingKEllipseArea(self, KPCAs, ET=None, retA=False):
        """Return the area of the minimum bounding kmer PCA ellipse"""
        if len(KPCAs) > 1:
            if ET is None:
                ET = EllipsoidTool()
            (A, center, radii, _rotation) = ET.getMinVolEllipse(KPCAs, retA=True)
            if retA:
                return ((A, center), ET.getEllipsoidVolume(radii))
            else:
                return ET.getEllipsoidVolume(radii)
        else: # minimum bounding ellipse of a point is 0
            if retA:
                return ((np.zeros((2,2)), KPCAs[0]), 0)
            else:
                return 0

#------------------------------------------------------------------------------
# Grow the bin

    def makeRanges(self, pos, span, limit):
        """Make search ranges which won't go out of bounds"""
        lower = pos-span
        upper = pos+span+1
        if(lower < 0):
            lower = 0
        if(upper > limit):
            upper = limit
        return np.arange(lower, upper)

    def recruit(self,
                PM,
                GT,
                im2RowIndices,
                inclusivity=1):
        """Recruit more contigs into the bin, used during coring only"""
        num_recruited = 0

        # make the distribution
        self.makeBinDist(PM.transformedCP, PM.averageCoverages, PM.kmerNormPC1, PM.kmerPCs, PM.contigGCs, PM.contigLengths)
        c_lens = PM.contigLengths[self.rowIndices]

        for x in self.makeRanges(self.covMedians[0], inclusivity*self.covStdevs[0], PM.scaleFactor):
            for y in self.makeRanges(self.covMedians[1], inclusivity*self.covStdevs[1], PM.scaleFactor):
                for z in self.makeRanges(self.covMedians[2], inclusivity*self.covStdevs[2], PM.scaleFactor):
                    # make sure it's a legit point
                    try:
                        for row_index in im2RowIndices[(x,y,z)]:
                            if (row_index not in PM.binnedRowIndices) and (row_index not in self.rowIndices) and (row_index not in PM.restrictedRowIndices):
                                # check the length
                                length_wrong = GT.isMaxOutlier(PM.contigLengths[row_index],
                                                               c_lens)
                                if not length_wrong:
                                    # fits length cutoff
                                    (covZ,merZ) = self.scoreProfile(PM.kmerNormPC1[row_index], PM.transformedCP[row_index])
                                    if covZ <= inclusivity and merZ <= inclusivity:
                                        # we can recruit
                                        self.rowIndices = np.append(self.rowIndices,row_index)
                                        num_recruited += 1
                    except KeyError: pass

    def shuffleMembers(self, adds, removes):
        """add some guys, take some guys away"""
        for row_index in self.rowIndices:
            if(row_index not in removes):
                adds.append(row_index)
        self.rowIndices = np.array(adds)
        self.binSize = self.rowIndices.shape[0]


#------------------------------------------------------------------------------
# IO and IMAGE RENDERING
#
    def plotProfileDistributions(self, transformedCP, kmerSigs, fileName=""):
        """plot the profile distibutions for this bin"""
        cov_working_array = np.zeros((self.binSize,3))
        mer_working_array = np.zeros((self.binSize,np.size(kmerSigs[0])))
        outer_index = 0
        for row_index in self.rowIndices:
            for i in range(0,3):
                cov_working_array[outer_index][i] = transformedCP[row_index][i]
            mer_working_array[outer_index] = kmerSigs[row_index]
            outer_index += 1

        # calculate the mean and stdev
        covMeans = np.mean(cov_working_array,axis=0)
        covStdevs = np.std(cov_working_array,axis=0)
        merMeans = np.mean(mer_working_array, axis=0)
        merStdevs = np.std(mer_working_array, axis=0)

        # z-normalise each column in each working array
        for index in range(0,np.size(self.rowIndices)):
            mer_working_array[index] = (mer_working_array[index]-merMeans)/merStdevs
            cov_working_array[index] = (cov_working_array[index]-covMeans)/covStdevs

        # work out the distribution of distances from z-normed sigs to the centroid
        k_dists = np.array([])
        c_dists = np.array([])
        merZeros = np.zeros((np.size(kmerSigs[0])))
        covZeros = np.zeros((3))

        for i in range(0,self.binSize):
            k_dists = np.append(k_dists, np.linalg.norm(mer_working_array[i]-merZeros))
            c_dists = np.append(c_dists, np.linalg.norm(cov_working_array[i]-covZeros))

        k_dists = np.sort(k_dists)
        c_dists = np.sort(c_dists)

        kDistMean = np.mean(k_dists)
        kDistStdev = np.std(k_dists)
        cDistMean = np.mean(c_dists)
        cDistStdev = np.std(c_dists)

        for i in range(0,self.binSize):
            k_dists[i] = (k_dists[i] - kDistMean)/kDistStdev
            c_dists[i] = (c_dists[i] - cDistMean)/cDistStdev

        B = np.arange(0, self.binSize, 1)

        fig = plt.figure()
        plt.subplot(211)
        plt.plot(B, k_dists, 'r-')
        plt.xlabel("kmer distribution")
        plt.subplot(212)
        plt.plot(B, c_dists, 'b-')
        plt.xlabel("coverage distribution")
        if(fileName != ""):
            try:
                fig.set_size_inches(10,4)
                plt.savefig(fileName,dpi=300)
            except:
                print "Error saving image:", fileName, sys.exc_info()[0]
                raise
        else:
            try:
                plt.show()
            except:
                print "Error showing image:", sys.exc_info()[0]
                raise
        del fig


    def plotBin(self, transformedCP, contigGCs, kmerNormPC1, contigLengths, colorMapGC, isLikelyChimeric, fileName="", ignoreContigLengths=False, ET=None):
        """Plot a single bin"""
        fig = plt.figure()
        title = self.plotOnFig(fig, 1, 1, 1,
                               transformedCP,
                               contigGCs,
                               contigLengths,
                               colorMapGC,
                               isLikelyChimeric,
                               fileName=fileName,
                               ignoreContigLengths=ignoreContigLengths,
                               ET=ET)

        plt.title(title)
        if(fileName != ""):
            try:
                fig.set_size_inches(6,6)
                plt.savefig(fileName+".png",dpi=300)
            except:
                print "Error saving image:", fileName, sys.exc_info()[0]
                raise
        elif(show):
            try:
                plt.show()
            except:
                print "Error showing image:", sys.exc_info()[0]
                raise
        plt.close(fig)
        del fig

    def plotOnFig(self,
                  fig,
                  plot_rows,
                  plot_cols,
                  plot_num,
                  transformedCP,
                  contigGCs,
                  contigLengths,
                  colorMapGC,
                  isLikelyChimeric,
                  fileName="",
                  ignoreContigLengths=False,
                  ET=None,
                  plotColorbar=True,
                  extents=None):
        ax = fig.add_subplot(plot_rows, plot_cols, plot_num, projection='3d')
        return self.plotOnAx(ax,
                             transformedCP,
                             contigGCs,
                             contigLengths,
                             colorMapGC,
                             isLikelyChimeric,
                             fileName=fileName,
                             ignoreContigLengths=ignoreContigLengths,
                             ET=ET,
                             plotColorbar=plotColorbar,
                             extents=extents)

    def plotOnAx(self,
                 ax,
                 transformedCP,
                 contigGCs,
                 contigLengths,
                 colorMapGC,
                 isLikelyChimeric,
                 fileName="",
                 ignoreContigLengths=False,
                 plotCentroid=True,
                 ET=None,
                 printID=False,
                 plotColorbar=True,
                 extents=None):
        """Plot a bin in a given subplot

        If you pass through an EllipsoidTool then it will plot the minimum bounding ellipsoid as well!
        """

        disp_vals = np.array([])
        disp_lens = np.array([])
        num_points = 0
        for row_index in self.rowIndices:
            num_points += 1
            disp_vals = np.append(disp_vals, transformedCP[row_index])
            disp_lens = np.append(disp_lens, np.sqrt(contigLengths[row_index]))

        # make a black mark at the max values
        cc_string = ""
        if plotCentroid and printID == False:
            self.makeLimits()
            px = self.covMedians[0]
            py = self.covMedians[1]
            #pz = self.covMedians[2]
            #num_points += 1
            #disp_vals = np.append(disp_vals, [px,py,pz])
            #disp_lens = np.append(disp_lens, 100)
            cc_string = "Coverage centroid: %d %d [%d -> %d]\n" % (px,py,self.covLowerLimits[2],self.covUpperLimits[2])

        # fix these
        self.makeLimits()

        # reshape
        disp_vals = np.reshape(disp_vals, (num_points, 3))

        if ignoreContigLengths:
            sc = ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors='none', c=contigGCs[self.rowIndices], cmap=colorMapGC, vmin=0.0, vmax=1.0, s=10, marker='.')
        else:
            sc = ax.scatter(disp_vals[:,0], disp_vals[:,1], disp_vals[:,2], edgecolors='k', c=contigGCs[self.rowIndices], cmap=colorMapGC, vmin=0.0, vmax=1.0, s=disp_lens, marker='.')            
        sc.set_edgecolors = sc.set_facecolors = lambda *args:None # disable depth transparency effect

        ax.set_xlabel('x coverage')
        ax.set_ylabel('y coverage')
        ax.set_zlabel('z coverage')
        if plotColorbar:
            cbar = plt.colorbar(sc, shrink=0.5)
            cbar.ax.tick_params()
            cbar.ax.set_title("% GC", size=10)
            cbar.set_ticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
            cbar.ax.set_ylim([0.15, 0.85])
            cbar.outline.set_ydata([0.15] * 2 + [0.85] * 4 + [0.15] * 3)

        if ET != None:
            (center, radii, rotation) = self.getBoundingEllipsoid(transformedCP, ET=ET)
            centroid_gc = np.mean(contigGCs[self.rowIndices])
            centroid_color = colorMapGC(centroid_gc)
            if printID:
                ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=False, cageColor=centroid_color, label=self.id)
            else:
                ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=False, cageColor=centroid_color)

        if extents:
            ax.set_xlim([extents[0], extents[1]])
            ax.set_ylim([extents[2], extents[3]])
            ax.set_zlim([extents[4], extents[5]])

        from locale import format, setlocale, LC_ALL # purdy commas
        setlocale(LC_ALL, "")
        title = str.join(" ", ["Bin: %d : %d contigs : %s BP\n" %(self.id,self.binSize,format('%d', self.totalBP, True)),
                               cc_string,
                               "GC: median: %.4f stdev: %.4f\n" % (self.gcMedian, self.gcStdev)]
                         )

        if isLikelyChimeric[self.id]:
            title += "Likely Chimeric"

        return title

    def plotMersOnAx(self, ax, kPCA1, kPCA2, contigGCs, contigLengths, colorMapGC, fileName="", ET=None, printID=False, plotColorbar=True):
        """Plot a bins kmer sig PCAs in a given subplot

        If you pass through an EllipsoidTool then it will plot the minimum bounding ellipse as well!
        """
        disp_vals = np.array(zip([kPCA1[i] for i in self.rowIndices],
                                 [kPCA2[i] for i in self.rowIndices]))
        disp_lens = np.array([np.sqrt(contigLengths[i]) for i in self.rowIndices])

        # reshape
        disp_vals = np.reshape(disp_vals, (len(self.rowIndices), 2))

        sc = ax.scatter(disp_vals[:,0], disp_vals[:,1], edgecolors='k', c=contigGCs[self.rowIndices], cmap=colorMapGC, vmin=0.0, vmax=1.0, s=disp_lens, marker='.')

        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        if plotColorbar:
            cbar = plt.colorbar(sc, shrink=0.5)
            cbar.ax.tick_params()
            cbar.ax.set_title("% GC", size=10)
            cbar.set_ticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
            cbar.ax.set_ylim([0.15, 0.85])
            cbar.outline.set_ydata([0.15] * 2 + [0.85] * 4 + [0.15] * 3)

        if ET != None:
            (center, radii, rotation) = ET.getMinVolEllipse(disp_vals)
            centroid_gc = np.mean(contigGCs[self.rowIndices])
            centroid_color = colorMapGC(centroid_gc)
            if printID:
                ET.plotEllipse(center, radii, rotation, ax=ax, plotAxes=False, cageColor=centroid_color, label=self.id)
            else:
                ET.plotEllipse(center, radii, rotation, ax=ax, plotAxes=False, cageColor=centroid_color)

    def printBin(self, contigNames, covProfiles, contigGCs, contigLengths, isLikelyChimeric, outFormat="summary", separator="\t", stream=sys.stdout):
        """print this bin info in csvformat"""
        kvm_str = "%.4f" % self.kValMeanNormPC1
        kvs_str = "%.4f" % self.kValStdevNormPC1
        cvm_str = "%.4f" % self.cValMedian
        cvs_str = "%.4f" % self.cValStdev
        gcm_str = "%.4f" % self.gcMedian
        gcs_str = "%.4f" % self.gcStdev

        if(outFormat == 'summary'):
            stream.write(separator.join([str(self.id), str(isLikelyChimeric[self.id]), str(self.totalBP), str(self.binSize), cvm_str, cvs_str, gcm_str, gcs_str])+"\n")
        elif(outFormat == 'full'):
            stream.write("#bid_"+str(self.id)+
                  "_likelyChimeric_"+str(isLikelyChimeric[self.id])+
                  "_totalBP_"+str(self.totalBP)+
                  "_numCons_"+str(self.binSize)+
                  "_gcMean_"+gcm_str+
                  "_gcStdev_"+gcs_str+
                  "_kMean_"+kvm_str+
                  "_kStdev_"+kvs_str+
                  "\n")
            stream.write(separator.join(["#\"bid\"""\"cid\"","\"length\""])+"\n")
            for row_index in self.rowIndices:
                stream.write(separator.join([str(self.id), contigNames[row_index], str(contigLengths[row_index])])+"\n")
        elif(outFormat == 'contigs'):
            for row_index in self.rowIndices:
                stream.write(separator.join([str(self.id), contigNames[row_index], str(contigLengths[row_index]), '%.4f' % contigGCs[row_index]])+"\n")
        elif(outFormat == 'bins'):
            data = [str(self.id), str(isLikelyChimeric[self.id]), str(self.totalBP), str(self.binSize), gcm_str, gcs_str]
            cov_mean = np.mean(covProfiles[self.rowIndices], axis=0)
            cov_std = np.std(covProfiles[self.rowIndices], axis=0)
            for i in xrange(0, len(cov_mean)):
                data.append('%.4f' % cov_mean[i])
                data.append('%.4f' % cov_std[i])
            stream.write(separator.join(data)+"\n")
        else:
            stream.write("--------------------------------------\n")
            stream.write("Bin:", self.id,"\n")
            stream.write("Bin size:", self.binSize,"\n")
            stream.write("Total BP:", self.totalBP,"\n")
            stream.write("GC Mean:", gcm_str,"\n")
            stream.write("GC Stdev:", gcs_str,"\n")
            stream.write("KMean:", kvm_str,"\n")
            stream.write("KStdev:", kvs_str,"\n")
            stream.write("--------------------------------------\n")

###############################################################################
###############################################################################
###############################################################################
###############################################################################

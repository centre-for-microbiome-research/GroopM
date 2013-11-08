#!/usr/bin/env python
###############################################################################
#                                                                             #
#    groopm.py                                                                #
#                                                                             #
#    Wraps coarse workflows                                                   #
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
__version__ = "0.2.10.11"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Release"

###############################################################################

import matplotlib as mpl

# GroopM imports
import mstore
import cluster
import refine
import binManager
import groopmUtils
import groopmTimekeeper as gtime
from groopmExceptions import ExtractModeNotAppropriateException
from mstore import GMDataManager

###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Track rogue print statements
#from groopmExceptions import Tracer
#import sys
#sys.stdout = Tracer(sys.stdout)
#sys.stderr = Tracer(sys.stderr)
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class GroopMOptionsParser():
    def __init__(self, version):
        # set default value for matplotlib
        mpl.rcParams['lines.linewidth'] = 1

        mpl.rcParams['xtick.labelsize'] = 8
        mpl.rcParams['ytick.labelsize'] = 8
        mpl.rcParams['legend.fontsize'] = 10
        mpl.rcParams['axes.labelsize'] = 10
        mpl.rcParams['axes.titlesize'] = 12
        mpl.rcParams['axes.linewidth'] = 0.25

        mpl.rcParams['axes3d.grid'] = True

        mpl.rcParams['savefig.dpi'] = 300

        mpl.rcParams['figure.figsize'] = [6.5, 6.5]
        mpl.rcParams['figure.facecolor'] = '1.0'

        self.GMVersion = version

    def parseOptions(self, options ):
        timer = gtime.TimeKeeper()
        if(options.subparser_name == 'parse'):
            # parse raw input
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in data parsing mode..." % self.GMVersion
            print "*******************************************************************************"
            GMdata = mstore.GMDataManager()
            success = GMdata.createDB(options.bamfiles,
                                      options.reference,
                                      options.dbname,
                                      timer,
                                      force=options.force)
            if not success:
                print options.dbname,"not updated"

        elif(options.subparser_name == 'core'):
            # make bin cores
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in core creation mode..." % self.GMVersion
            print "*******************************************************************************"
            CE = cluster.ClusterEngine(options.dbname,
                                       timer,
                                       force=options.force,
                                       finalPlot=options.plot,
                                       plot=options.multiplot,
                                       minSize=options.size,
                                       minVol=options.bp)
            if options.graphfile is None:
                gf = ""
            else:
                gf=options.graphfile
            CE.makeCores(coreCut=options.cutoff,
                         gf=gf)

        elif(options.subparser_name == 'refine'):
            # refine bin cores
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in core refining mode..." % self.GMVersion
            print "*******************************************************************************"
            bids = []
            #if options.bids is not None:
            #    bids = options.bids
            auto = options.auto
            transform=True^options.no_transform

            RE = refine.RefineEngine(timer,
                                     dbFileName=options.dbname,
                                     transform=transform,
                                     bids=bids,
                                     loadContigNames=True)

            if options.plot:
                pfx="REFINED"
            else:
                pfx=""
            print "Refine bins"

            RE.refineBins(timer,
                          auto=auto,
                          saveBins=True,
                          plotFinal=pfx)

        elif(options.subparser_name == 'recruit'):
            # make bin cores
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in bin expansion mode..." % self.GMVersion
            print "*******************************************************************************"
            RE = refine.RefineEngine(timer,
                                     dbFileName=options.dbname,
                                     getUnbinned=True,
                                     loadContigNames=False,
                                     cutOff=options.cutoff)

            RE.recruitWrapper(timer,
                              inclusivity=options.inclusivity,
                              step=options.step,
                              saveBins=True)

        elif(options.subparser_name == 'extract'):
            # Extract data
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in '%s' extraction mode..." % (self.GMVersion, options.mode)
            print "*******************************************************************************"
            bids = []
            if options.bids is not None:
                bids = options.bids
            BX = groopmUtils.GMExtractor(options.dbname,
                                          bids=bids,
                                          folder=options.outfolder
                                          )
            if(options.mode=='contigs'):
                BX.extractContigs(timer, fasta=options.data, cutoff=options.cutoff)
            elif(options.mode=='reads'):
                BX.extractReads(timer, bams=options.data)
            else:
                raise ExtractModeNotAppropriateException("mode: "+ options.mode + " is unknown")

        elif(options.subparser_name == 'merge'):
            # make bin cores
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in bin merging mode..." % self.GMVersion
            print "*******************************************************************************"
            BM = binManager.BinManager(dbFileName=options.dbname)
            BM.loadBins(timer, makeBins=True, silent=False)
            BM.merge(options.bids, options.force, saveBins=True)

        elif(options.subparser_name == 'split'):
            # make bin cores
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in bin splitting mode..." % self.GMVersion
            print "*******************************************************************************"
            BM = binManager.BinManager(dbFileName=options.dbname)
            BM.loadBins(timer, makeBins=True, silent=False)
            BM.split(options.bid, options.parts, mode=options.mode, saveBins=True, auto=options.force)

        elif(options.subparser_name == 'delete'):
            # make bin cores
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in bin deleting mode..." % self.GMVersion
            print "*******************************************************************************"
            BM = binManager.BinManager(dbFileName=options.dbname)
            BM.loadBins(timer, makeBins=True, silent=True)#, bids=options.bids)
            BM.deleteBins(options.bids, force=options.force, saveBins=True, freeBinnedRowIndices=True)

        elif(options.subparser_name == 'plot'):
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in bin plotting mode..." % self.GMVersion
            print "*******************************************************************************"
            BM = binManager.BinManager(dbFileName=options.dbname)

            if options.bids is None:
                bids = []
            else:
                bids = options.bids
            BM.loadBins(timer, makeBins=True, silent=False, bids=bids, loadContigNames=False)

            BM.setColorMap(options.cm)

            BM.plotBins(FNPrefix=options.tag,
                        plotEllipsoid=True,
                        ignoreContigLengths=options.points,
                        folder=options.folder)

        elif(options.subparser_name == 'explore'):
            # make bin cores
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in bin '%s' explorer mode..." % (self.GMVersion, options.mode)
            print "*******************************************************************************"
            transform=True^options.no_transform
            bids = []
            if options.bids is not None:
                bids = options.bids
            BE = groopmUtils.BinExplorer(options.dbname,
                                         bids=bids,
                                         transform=transform,
                                         cmstring=options.cm,
                                         ignoreContigLengths=options.points)
            if(options.mode == 'binpoints'):
                BE.plotPoints(timer)
            elif(options.mode == 'binids'):
                BE.plotIds(timer)
            elif(options.mode == 'allcontigs'):
                BE.plotContigs(timer, coreCut=options.cutoff, all=True)
            elif(options.mode == 'unbinnedcontigs'):
                BE.plotUnbinned(timer, coreCut=options.cutoff)
            elif(options.mode == 'binnedcontigs'):
                BE.plotContigs(timer, coreCut=options.cutoff)
            elif(options.mode == 'flyover'):
                BE.plotFlyOver(timer)
            elif(options.mode == 'binassignments'):
                BE.plotBinAssignents(timer, coreCut=options.cutoff)
            elif(options.mode == 'compare'):
                BE.plotCompare(timer, coreCut=options.cutoff)
            elif (options.mode == 'together'):
                BE.plotTogether(timer, coreCut=options.cutoff, doMers=options.kmers)
            elif (options.mode == 'sidebyside'):
                BE.plotSideBySide(timer, coreCut=options.cutoff)
            else:
                print "**Error: unknown mode:",options.mode

        elif(options.subparser_name == 'highlight'):
            # make bin cores
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in highlighter mode..." % self.GMVersion
            print "*******************************************************************************"
            BE = groopmUtils.BinExplorer(options.dbname,
                                         binLabelsFile = options.binlabels,
                                         contigColorsFile = options.contigcolors,
                                         ignoreContigLengths=options.points)
            BE.plotHighlights(timer,
                              options.elevation,
                              options.azimuth,
                              options.file,
                              options.filetype,
                              options.dpi,
                              drawRadius=options.radius,
                              show=options.show,
                              coreCut=options.cutoff,
                              testing=options.place
                              )

        elif(options.subparser_name == 'print'):
            BM = binManager.BinManager(dbFileName=options.dbname)
            bids = []
            if options.bids is not None:
                bids = options.bids
            BM.loadBins(timer, getUnbinned=options.unbinned, makeBins=True, silent=True, bids=bids)
            BM.printBins(options.format, fileName=options.outfile)

        elif(options.subparser_name == 'dump'):
            print "*******************************************************************************"
            print " [[GroopM %s]] Running in data dumping mode..." % self.GMVersion
            print "*******************************************************************************"

            # prep fields. Do this first cause users are mot likely to
            # mess this part up!
            allowable_fields = ['names', 'mers', 'gc', 'coverage', 'tcoverage', 'ncoverage', 'lengths', 'bins', 'all']
            fields = options.fields.split(',')
            for field in fields:
                if field not in allowable_fields:
                    print "ERROR: field '%s' not recognised. Allowable fields are:" % field
                    print '\t',",".join(allowable_fields)
                    return
            if options.separator == '\\t':
                separator = '\t'
            else:
                separator = options.separator

            DM = GMDataManager()
            DM.dumpData(options.dbname,
                        fields,
                        options.outfile,
                        separator,
                        not options.no_headers)

        return 0


###############################################################################
###############################################################################
###############################################################################
###############################################################################

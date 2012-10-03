#!/usr/bin/env python
###############################################################################
#                                                                             #
#    groopm.py                                                                #
#                                                                             #
#    Implements groopm shell and wraps coarse workflows                       #
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

import argparse
import sys
import os 
#import code
#import readline 
#import rlcompleter 
#import atexit

# GroopM imports
import mstore
import cluster
import bin
import dataManagers
import groopmUtils
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Track rogue print statements
#from groopmExceptions import Tracer
#sys.stdout = Tracer(sys.stdout)
#sys.stderr = Tracer(sys.stderr)
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class GroopMOptionsParser():
    def __init__(self):
        return
    
    def parseOptions(self, options ):

        if(options.subparser_name == 'parse'):
            # parse raw input
            print "****************************************************************"
            print " [[GroopM]] Running in data parsing mode..."
            print "****************************************************************"
            GMdata = mstore.GMDataManager()
            success = GMdata.createDB(options.bamfiles,
                                      options.reference,
                                      options.dbname,
                                      dumpAll=options.dump,
                                      force=options.force
                                      )
            if not success:
                print options.dbname,"not updated" 
                            
        elif(options.subparser_name == 'core'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in core creation mode..."
            print "****************************************************************"
            CE = cluster.ClusterEngine(options.dbname,
                                       force=options.force,
                                       plot=options.plot
                                       )
            CE.makeCores(coreCut=options.cutoff, minSize=options.size, minVol=options.bp)

        elif(options.subparser_name == 'refine'):
            # refine bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in core condensing mode..."
            print "****************************************************************"
            plotter = False
            chimera_check = False
            outlier_check = False
            if(options.mode == 'plot'):
                plotter = True
            #elif(options.mode == 'chimera'):
            #    chimera_check = True
            elif(options.mode == 'outlier'):
                outlier_check = True
            else:
                print "**Error: unknown mode:",options.mode
                return
            BM = dataManagers.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=False)
            BM.condenseWrapper(save=True,
                               chimeraCheck=chimera_check,
                               outlierCheck=outlier_check,
                               plotter=plotter
                               )

        elif(options.subparser_name == 'makesoms'):
            # make SOMs
            print "****************************************************************"
            print " [[GroopM]] Running in SOM training mode..."
            print "****************************************************************"
            BM = dataManagers.BinManager(dbFileName=options.dbname)
            SM = dataManagers.SOMManager(BM, somSide=options.side, somIterations=options.iterations)
            do_merge = not options.no_merge
            SM.DoSOMPipeline(merge=do_merge, force=options.force, tag=options.tag)

        elif(options.subparser_name == 'expand'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin expansion mode..."
            print "****************************************************************"
            CE = cluster.ClusterEngine(options.dbname)
            CE.expandBins(force=options.force)
        
        elif(options.subparser_name == 'extract'):
            # Extract data
            print "****************************************************************"
            print " [[GroopM]] Running in extraction mode(",options.mode,")..."
            print "****************************************************************"
            bids = []
            if options.bids is not None:
                bids = options.bids
            BX = groopmUtils.GMExtractor(options.dbname,
                                          bids=bids,
                                          folder=options.outfolder
                                          )
            if(options.mode=='contigs'):
                BX.extractContigs(fasta=options.data, cutoff=options.cutoff)
            elif(options.mode=='reads'):
                BX.extractReads(bams=options.data, shuffle=options.shuffle)
            else:
                raise ge.ExtractModeNotAppropriateException("mode: "+mode+" is unknown")

        elif(options.subparser_name == 'merge'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin merging mode..."
            print "****************************************************************"
            BM = dataManagers.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=False)
            BM.merge(options.bids, options.force, saveBins=True)

        elif(options.subparser_name == 'split'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin splitting mode..."
            print "****************************************************************"
            BM = dataManagers.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=False)
            BM.split(options.bid, options.parts, mode=options.mode, saveBins=True, auto=options.force)

        elif(options.subparser_name == 'delete'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin deleting mode..."
            print "****************************************************************"
            BM = dataManagers.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=True, bids=options.bids)
            BM.deleteBins(options.bids, force=options.force, saveBins=True, freeBinnedRowIndicies=True)

        elif(options.subparser_name == 'explore'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin explorer mode (",options.mode,")..."
            print "****************************************************************"
            bids = []
            if options.bids is not None:
                bids = options.bids
            BE = groopmUtils.BinExplorer(options.dbname, bids=bids)
            if(options.mode == 'points'):
                BE.plotPoints()
            elif(options.mode == 'ids'):
                BE.plotIds()
            elif(options.mode == 'flyover'):
                BE.plotFlyOver()
            elif(options.mode == 'profile'):
                BE.plotBinProfiles()
            elif(options.mode == 'compare'):
                BE.plotSideBySide(coreCut=options.cutoff)
            elif(options.mode == 'unbinned'):
                BE.plotUnbinned(coreCut=options.cutoff)
            else:
                print "**Error: unknown mode:",options.mode
            
        elif(options.subparser_name == 'print'):
            BM = dataManagers.BinManager(dbFileName=options.dbname)
            bids = []
            if options.bids is not None:
                bids = options.bids
            BM.loadBins(getUnbinned=options.unbinned, makeBins=True, silent=True, bids=bids)
            BM.printBins(options.format, fileName=options.outfile)

        elif(options.subparser_name == 'plot'):
            print "****************************************************************"
            print " [[GroopM]] Plot bins..."
            print "****************************************************************"
            BM = dataManagers.BinManager(dbFileName=options.dbname)
            bids = []
            if options.bids is not None:
                bids = options.bids
            BM.loadBins(makeBins=True, silent=False, bids=bids)
            BM.plotBins(FNPrefix=options.tag, sideBySide=options.sidebyside)
            
        return 0
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################

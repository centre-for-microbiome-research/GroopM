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
import code
import readline 
import rlcompleter 
import os 
import atexit

# GroopM imports
import mstore
import cluster
import binUtils

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class GroopMOptionsParser():
    def __init__(self):
        return
    
    def parseOptions(self, options ):
        if(options.subparser_name == 'prompt'):
            # drop to shell prompt
            vars = globals().copy()
            vars.update(locals())
            interactive_shell = GroopMInteractiveConsole(vars)
            interactive_shell.startConsole()
        
        elif(options.subparser_name == 'batch'):
            # run batch commands            
            print "****************************************************************"
            print " [[GroopM]] Running in batch mode..."
            print "****************************************************************"
            #batch_shell = cmd.GroopMBatchShell()
            #batch_shell.run()
        
        elif(options.subparser_name == 'parse'):
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

        elif(options.subparser_name == 'expand'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin expansion mode..."
            print "****************************************************************"
            CE = cluster.ClusterEngine(options.dbname,
                                       force=options.force,
                                       plot=options.plot
                                       )
            CE.expandBins()
        
        elif(options.subparser_name == 'refine'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in core condensing mode..."
            print "****************************************************************"
            BM = binUtils.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=False)
            BM.condenseWrapper(save=True,auto=options.auto)

        elif(options.subparser_name == 'separate'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in core separating mode..."
            print "****************************************************************"
            BM = binUtils.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=False)
            BM.chimeraWrapper(save=True,auto=options.auto)

        elif(options.subparser_name == 'merge'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin merging mode..."
            print "****************************************************************"
            BM = binUtils.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=False)
            BM.merge(options.bids, options.auto, saveBins=True)

        elif(options.subparser_name == 'split'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin splitting mode..."
            print "****************************************************************"
            BM = binUtils.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=False)
            BM.split(options.bid, options.parts, mode=options.mode, saveBins=True, test=options.auto, auto=False)

        elif(options.subparser_name == 'delete'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin deleting mode..."
            print "****************************************************************"
            BM = binUtils.BinManager(dbFileName=options.dbname)
            BM.loadBins(makeBins=True, silent=True, bids=options.bids)
            BM.deleteBins(options.bids, force=options.force, saveBins=True, freeBinnedRowIndicies=True)

        elif(options.subparser_name == 'explore'):
            # make bin cores
            print "****************************************************************"
            print " [[GroopM]] Running in bin explorer mode..."
            print "****************************************************************"
            bids = []
            if options.bids is not None:
                bids = options.bids
            BE = binUtils.BinExplorer(options.dbname, bids=bids)
            if(options.flyover):
                BE.plotFlyOver()
            elif(options.profiles):
                BE.plotBinProfiles()
            else:
                BE.plotSideBySide(coreCut=options.cutoff)
            
        elif(options.subparser_name == 'print'):
            BM = binUtils.BinManager(dbFileName=options.dbname)
            bids = []
            if options.bids is not None:
                bids = options.bids
            BM.loadBins(getUnbinned=options.unbinned, makeBins=True, silent=True, bids=bids)
            BM.printBins(options.format, fileName=options.outfile)

        elif(options.subparser_name == 'plot'):
            print "****************************************************************"
            print " [[GroopM]] Plot bins..."
            print "****************************************************************"
            BM = binUtils.BinManager(dbFileName=options.dbname)
            bids = []
            if options.bids is not None:
                bids = options.bids
            BM.loadBins(makeBins=True, silent=False, bids=bids)
            BM.plotBins(FNPrefix=options.tag, sideBySide=options.sidebyside)
            
        else:
            print "****************************************************************"
            print " [[GroopM]] - Use -h for help"
            print "****************************************************************"
        return 0
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class GroopMInteractiveConsole(code.InteractiveConsole):
    "Wrapper around Python that can filter input/output to the shell"
    def __init__(self, locals=None, filename="<console>"):
        
        self.shell_banner="""----------------------------------------------------------------
 Welcome to the GroopM interactive console (%s)
 Copyright (C) Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
 
 Type "gmhelp" for help
 
 This is a python interactive shell so
 you can use python commands here too
 To exit, press [CTRL]-D or type "exit()"

 Saving shell history to %s
----------------------------------------------------------------"""
        
        code.InteractiveConsole.__init__(self, locals)
        return

    def push(self,line):
        """Intercept user input and react accordingly"""
        # split the line on space
        broken_cmd = line.split(" ")
        if (broken_cmd[0] == "gmhelp"):
            if(len(broken_cmd) == 1):
                self.printHelp()
            else:
                self.printHelp(subHelp=broken_cmd[1])
            return
        
        return code.InteractiveConsole.push(self,line)

    def startConsole(self):
        """Start the python console"""
        # tab completion 
        readline.parse_and_bind('tab: complete') 
        
        # history file 
        histfile = os.path.join(os.environ['HOME'], '.GroopMhistory') 
        try: 
            readline.read_history_file(histfile) 
        except IOError: 
            pass 
        atexit.register(readline.write_history_file, histfile) 

        self.interact(banner=self.shell_banner % (__version__, os.path.join(os.environ['HOME'], '.GroopMhistory')))
        
    def printHelp(self, subHelp=""):
        """Print user help - quite helpful"""
        subcomands = [["loadsams", "\nHelp on parsing samfiles", ""],
                      ["loadcons", "\nHelp on parsing contig files", ""],
                      ["loadsec", "\nHelp on parsing secondary data", ""],
                      ["transform", "\nHelp on data transformations", ""],
                      ["cluster", "\nHelp on initial clustering", ""],
                      ["refine", "\nHelp on cluster refinement", ""],
                      ["som", "\nHelp on GroopM's SOMs", ""],
                      ["output", "\nHelp on printing / saving output", ""],
                      ] 
        
        if("" == subHelp):
            print """For help on a subtopic type "gmhelp subtopic"
    
    Subtopics:
"""
            for sub in subcomands:
                print "\t",sub[0]
        else:
            for sub in subcomands:
                if(subHelp == sub[0]):
                    print sub[1]

###############################################################################
###############################################################################
###############################################################################
###############################################################################
                    
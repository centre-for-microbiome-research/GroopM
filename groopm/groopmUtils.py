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
import numpy as np

# GroopM imports
import mstore
import cluster

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class PrintEngine:
    """Print bin information"""
    def __init__(self, dbFileName, outFormat, fileName=""):
        self.dbFileName = dbFileName
        self.dataManager = mstore.GMDataManager()
        self.fileName = fileName
        self.format = outFormat

        # primary data storage
        self.indicies = np.array([])
        self.contigNames = np.array([])
        self.contigLengths = np.array([])
        self.bins = np.array([])
        self.numBins = 0
        
        # munged data 
        self.bin_sizes = {}
        self.bin_members = {}
        
    def loadData(self, getUnbinned=False, bins=[]):
        """load information from the DB"""
        if(getUnbinned):
            # get the lot!
            self.indicies = self.dataManager.getConditionalIndicies(self.dbFileName)
        else:
            # only get those with a non-zero bin ID
            if(len(bins) == 0):
                self.indicies = self.dataManager.getConditionalIndicies(self.dbFileName,
                                                                        condition='bin != 0')
            else:
                # get only those we're told to get
                condition = "((bin == "+str(bins[0])+")"
                for index in range (1,len(bins)):
                    condition += " | (bin == "+str(bins[index])+")"
                
                condition += ")"
                self.indicies = self.dataManager.getConditionalIndicies(self.dbFileName,
                                                                        condition=condition)
        self.contigNames = self.dataManager.getContigNames(self.dbFileName,
                                                           indicies=self.indicies)
        self.contigLengths = self.dataManager.getContigLengths(self.dbFileName,
                                                           indicies=self.indicies)
        self.bins = self.dataManager.getBins(self.dbFileName, indicies=self.indicies)
        self.numBins = self.dataManager.getNumBins(self.dbFileName)
        
        self.initialiseContainers()

    def initialiseContainers(self):
        """Munge the raw data into something more usable"""
        # intialise these containers
        for index in range(0,self.numBins+1):
            self.bin_sizes[index] = 0;
            self.bin_members[index] = []
        
        # fill them up
        for index in range(0, np.size(self.indicies)):
            self.bin_members[self.bins[index]].append(index)
            self.bin_sizes[self.bins[index]] += self.contigLengths[index]

    def printBins(self):
        """Wrapper for print handles piping to file or stdout"""
        if("" != self.fileName):
            try:
                # redirect stdout to a file
                sys.stdout = open(self.fileName, 'w')
                self.printInner()
            except:
                print "Error diverting stout to file:", fileName, sys.exc_info()[0]
                raise
        else:
            self.printInner()           
        
    def printInner(self):
        """Print bin information to STDOUT"""
        if(self.format == 'summary'):
            print "#\"bid\"\t\"totalBP\"\t\"numCons\""
            for bid in self.bin_members:
                if(np.size(self.bin_members[bid]) > 0):
                    print str(bid)+"\t"+str(self.bin_sizes[bid])+"\t"+str(np.size(self.bin_members[bid]))
        elif(self.format == 'full'):
            for bid in self.bin_members:
                if(np.size(self.bin_members[bid]) > 0):
                    print "#bid_"+str(bid)+"_totalBP_"+str(self.bin_sizes[bid])+"_numCons_"+str(np.size(self.bin_members[bid]))
                    print "#\"bid\"\t\"cid\"\t\"length\""            
                    for member in self.bin_members[bid]:
                        print bid, self.contigNames[member], self.contigLengths[member]
        elif(self.format == 'minimal'):
            print "#\"bid\"\t\"cid\"\t\"length\""            
            for bid in self.bin_members:
                if(np.size(self.bin_members[bid]) > 0):
                    for member in self.bin_members[bid]:
                        print bid, self.contigNames[member], self.contigLengths[member]
            pass
        else:
            print "Error: Unrecognised format:", self.format
        
###############################################################################
###############################################################################
###############################################################################
###############################################################################

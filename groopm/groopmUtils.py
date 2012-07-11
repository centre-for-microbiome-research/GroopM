#!/usr/bin/env python
###############################################################################
#                                                                             #
#    groopm_utils.py                                                          #
#                                                                             #
#    Classes for heavy lifting of raw data + output methods                   #
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
import argparse
import sys
from pprint import pprint

import tables
import sys
import pysam
import csv
import os

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################
# CLASSES
###############################################################################

class GMProj:
    """main object to poke when doing gm work"""
    def __init__(self):
        self.conParser = None
        self.bamParser = None
        self.dbName = "unset"
        self.contigsFile = "unset"
        self.bamFiles = None
        
    def loadDB(self, fileName):
        pass
    
    def createDB(self, bamFiles, contigs, profileVals, dbName):
        """Main wrapper for parsing all the files in"""
        # load all the passed vars
        self.dbName = dbName
        self.contigsFile = contigs
        self.bamFiles = str.split(bamFiles)
        self.conParser = ContigParser(self.contigsFile)
        self.bamParser = BamParser(self.bamFiles)
        
        # parse the files
        # TODO Thread this!
        do_continue = self.bamParser.parse(self.dbName)
        if(do_continue):
            self.conParser.parse(self.dbName)
        

class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self, contigsFile):
        self.contigsFile = contigsFile

    def parse(self, dbName):
        try:
            with open(dbName) as f: pass
        except IOError as e:
            print "Error: database:",dbName,"does not exist - perhaps you should call the bamParser first"

class BamParser:
    """Parse multiple bam files and write the output to hdf5 """

    def __init__(self, bamFiles):
        
        self.bamFiles = split

    def parse(self, dbName):
        """ parse multiple bam files and store the results in the main DB"""
        try:
            with open('filename') as f:
                option = raw_input("Warning: database:",dbName,"exists. Overwrite? (y,n)")
                if(option.upper() != "Y"):
                    print "Operation cancelled"
                    return False
                else:
                    print "Overwriting database",dbName
        except IOError as e:
            print "Creating new database", dbName
            
        # build a table template based on the number of bamfiles we have
        db_desc = { 'cid' : StringCol(64) }
        for bf in self.bamFiles:
            db_desc[bf] = tables.FloatCol()
        
        # create the db
        with tables.openFile(filename, mode = "w", title = "GroopM") as h5file:
            # Create a new group under "/" (root)
            group = h5file.createGroup("/", 'profile', 'Assembly profiles')
            # Create one table on it
            table = h5file.createTable(group, 'coverage', db_desc, "Bam based coverage")
            # open the SAM/BAM
            for bam_file in args:
                
                if (self.isBam):
                    try:
                        bamFile = pysam.Samfile(bam_file, 'rb')
                    except:
                        print "Unable to open BAM file -- did you supply a SAM file instead?"
                        return False
                else:
                    try:
                        bamFile = pysam.Samfile(bam_file, 'r')
                    except:
                        print "Unable to open SAM file -- did you supply a BAM file instead?"
                        return False                

                for reference, length in zip(bamFile.references, bamFile.lengths):
                    num_reads = bamFile.count(reference, 0, length)
                    print num_reads, length
                    tid = bamFile.gettid(reference)
            table.flush()


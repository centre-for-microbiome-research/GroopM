#!/usr/bin/env python
###############################################################################
#                                                                             #
#    groopmUtils.py                                                           #
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
import numpy as np

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
        self.numContigs = -1
        self.numBamFiles = -1
        
    def loadDB(self, fileName):
        pass
    
    def createDB(self, bamFiles, contigs, profileVals, dbName):
        """Main wrapper for parsing all the files in"""
        # load all the passed vars
        self.dbName = dbName
        self.contigsFile = contigs
        self.bamFiles = str.split(bamFiles)
        self.numBamFiles = len(self.bamFiles)
        self.profileVals = profileVals
        
        try:
            self.conParser = ContigParser(self.contigsFile)
        except:
            print "Could not create ContigParser:", sys.exc_info()[0]
        try:
            self.bamParser = BamParser(self.bamFiles)
        except:
            print "Could not create BamParser:", sys.exc_info()[0]
        
        # parse bam files
        do_continue = False
        try:
            do_continue = self.bamParser.parse(self.dbName)
        except:
            print "Error parsing bam files:", sys.exc_info()[0]
        
        # parse contigs 
        if(do_continue):
            try:
                self.numContigs = self.conParser.parse(self.dbName)
            except:
                print "Error parsing contigs:", sys.exc_info()[0]
        
        if(self.numContigs == -1):
            print "Something went wrong while parsing contigs"
            return False
        
        self.bamParser.dumpCovTable(self.dbName)

class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self, contigsFile):
        self.contigsFile = contigsFile

    def parse(self, dbName):
        try:
            with open(dbName) as f: pass
        except IOError as e:
            print "Error: database:",dbName,"does not exist - perhaps you should call the bamParser first"
        
        return 3

class Counter:
    """Callback for counting aligned reads"""
    counts = 0
    def __call__(self, alignment):
        self.counts += 1
           
class BamParser:
    """Parse multiple bam files and write the output to hdf5 """

    def __init__(self, bamFiles):
        self.bamFiles = bamFiles[0].split(',')
        self.numBamFiles = len(self.bamFiles)
        self.bamColNames = []
        
    def getBamDescriptor(self, fullPath):
        """reduce a full path to just the file name minus extension"""
        return os.path.splitext(os.path.basename(fullPath))[0]
    
    def parse(self, dbName):
        """ parse multiple bam files and store the results in the main DB"""
        # make sure we'r eonly overwriting existing DBs with the users consent
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
            
        # parse the BAMs
        tmp_storage = {}
        bam_count = 0
        try:
            for bf in self.bamFiles:
                bam_file = None
                try:
                    bam_file = pysam.Samfile(bf, 'rb')
                except:
                    print "Unable to open BAM file",bf,"-- did you supply a SAM file instead?"
                    return False

                try:
                    for reference, length in zip(bam_file.references, bam_file.lengths):
                        c = Counter()
                        try:
                            bam_file.fetch(reference, 0, length, callback = c )
                            num_reads = c.counts
                        except ValueError as e:
                            print "Could not calculate num reads for:",reference,"in",bf,"\t",e
                            return False
                        except:
                            print "Could not calculate num reads for:",reference,"in",bf, sys.exc_info()[0]
                            return False
                        
                        this_bf = self.getBamDescriptor(bf)
                        
                        # save this guy into the tmp_storage
                        if(reference not in tmp_storage):
                            tmp_storage[reference] = np.zeros((self.numBamFiles+1))
                            tmp_storage[reference][self.numBamFiles] = length
                        # we need to divide the count by the length if we are going
                        # to use a normalised coverage
                        tmp_storage[reference][bam_count] = float(num_reads)/float(length)
                except:
                    print "Error parsing bam file:",bf, sys.exc_info()[0]
                    return False
                bam_count += 1
        except:
            print "Error parsing bam files:", sys.exc_info()[0]
            return False

        # create the db and write results
        # build a table template based on the number of bamfiles we have
        db_desc = { 'cid' : tables.StringCol(64) }
        for bf in self.bamFiles:
            # assume the file is called somethign like "fred.bam"
            # we want to rip off the ".bam" part
            bam_desc = self.getBamDescriptor(bf)
            db_desc[bam_desc] = tables.FloatCol()
            self.bamColNames.append(bam_desc)
        db_desc['length'] = tables.FloatCol()
        
        try:        
            with tables.openFile(dbName, mode = "w", title = "GroopM") as h5file:
                # Create a new group under "/" (root)
                group = h5file.createGroup("/", 'profile', 'Assembly profiles')
                # Create one table on it
                table = h5file.createTable(group, 'coverage', db_desc, "Bam based coverage")
                # go through all the contigs sorted by name
                for cid in sorted(tmp_storage.keys()):
                    # make a new row
                    cov_row = table.row
                    # punch in the data
                    cov_row['cid'] = cid
                    for i in range(0,len(self.bamColNames)):
                        cov_row[self.bamColNames[i]] = tmp_storage[cid][i]
                    cov_row['length'] = tmp_storage[cid][self.numBamFiles]
                    cov_row.append()
                table.flush()
        except:
            print "Error saving results to DB"
            return False

        return True

    def dumpCovTable(self, dbName):
        with tables.openFile(dbName, mode = "r") as h5file:
            table = h5file.root.profile.coverage
            for row in table:
                print row['cid'],",",
                for i in range(0,len(self.bamColNames)):
                    print row[self.bamColNames[i]],",",
                print row['length']


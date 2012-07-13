#!/usr/bin/env python
###############################################################################
#                                                                             #
#    mstore.py                                                                #
#                                                                             #
#    GroopM data management                                                   #
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

import tables
import sys

# GroopM imports
from groopmUtils import ContigParser,BamParser,AuxParser,KmerSigEngine,getBamDescriptor

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class GMDataManager:
    """Top level class for manipulating GroopM data
    
    Use this class for parsing in raw data into a hdf DB and 
    for reading from and updating same DB
    """
    def __init__(self):
        self.conParser = None
        self.bamParser = None
        self.auxParser = None
        self.contigsFile = "unset"
        self.bamFiles = []
        self.dbName = "unset"
        self.bamColNames = []
        self.contigNames = []
        
    def loadDB(self, fileName): pass
    
    def createDB(self, bamFiles, contigs, aux, dbName, kmerSize=4, dumpAll=False):
        """Main wrapper for parsing all input files"""
        # load all the passed vars
        self.dbName = dbName
        self.contigsFile = contigs
        self.bamFiles = bamFiles.split(",")
        self.auxFile = aux
        
        self.kse = KmerSigEngine(kmerSize)
        self.conParser = ContigParser()
        self.bamParser = BamParser()
        self.auxParser = AuxParser()

        # make sure we're only overwriting existing DBs with the users consent
        try:
            with open(dbName) as f:
                option = raw_input("****************************************************************\n"\
                                   "    !!!WARNING!!! Database: "+dbName+" exists.\n" \
                                   "    If you continue you *WILL* delete any previous analyses!\n" \
                                   "****************************************************************\n"\
                                   " Overwrite? (y,n)")
                if(option.upper() != "Y"):
                    print "Operation cancelled"
                    return False
                else:
                    print "Overwriting database",dbName
        except IOError as e:
            print "Creating new database", dbName
            raise
        
        # create the db
        try:        
            with tables.openFile(dbName, mode = "w", title = "GroopM") as h5file:
                # Create a new group under "/" (root)
                group = h5file.createGroup("/", 'profile', 'Assembly profiles')

                #------------------------
                # parse bam files
                #------------------------
                # build a table template based on the number of bamfiles we have
                db_desc = { 'cid' : tables.StringCol(64) }
                for bf in self.bamFiles:
                    # assume the file is called somethign like "fred.bam"
                    # we want to rip off the ".bam" part
                    bam_desc = getBamDescriptor(bf)
                    db_desc[bam_desc] = tables.FloatCol()
                    self.bamColNames.append(bam_desc)
                db_desc['length'] = tables.Int32Col()
                try:
                    COV_table = h5file.createTable(group, 'coverage', db_desc, "Bam based coverage")
                    self.bamParser.parse(COV_table, self.bamFiles, self.bamColNames)
                except:
                    print "Error creating coverage table:", sys.exc_info()[0]
                    raise
                
                #------------------------
                # parse contigs 
                #------------------------
                db_desc = { 'cid' : tables.StringCol(64) }
                for mer in self.kse.kmerCols:
                     db_desc[mer] = tables.FloatCol()
                try:
                    TN_table = h5file.createTable(group, 'tns', db_desc, "TetraNucl signature")
                except:
                    print "Error creating TNS table:", sys.exc_info()[0]
                    raise
                try:
                    f = open(self.contigsFile, "r")
                    self.contigNames = self.conParser.parse(f, self.kse, TN_table)
                    f.close()
                except:
                    print "Could not parse contig file:",self.contigsFile,sys.exc_info()[0]
                    raise

                #------------------------
                # parse aux profile 
                #------------------------
                db_desc = { 'cid' : tables.StringCol(64), 'aux' :  tables.StringCol(64) }
                try:
                    AUX_table = h5file.createTable(group, 'aux', db_desc, "Secondary profile")
                except:
                    print "Error creating AUX table:", sys.exc_info()[0]
                    raise
                try:
                    f = open(self.auxFile, "r")
                    self.auxParser.parse(f, AUX_table)
                except:
                    print "Could not parse the auxilary profile file:",self.contigsFile,sys.exc_info()[0]
                    raise
                
                #------------------------
                # Add a table for the bins
                #------------------------
                db_desc = { 'cid' : tables.StringCol(64), 'bin' :  tables.Int32Col() }
                try:
                    BIN_table = h5file.createTable(group, 'bin', db_desc, "Bin IDs")
                except:
                    print "Error creating BIN table:", sys.exc_info()[0]
                    raise
                self.initBins(BIN_table)
        except:
            print "Error creating database:", dbName, sys.exc_info()[0]
            raise

        if(dumpAll):
            with tables.openFile(dbName, mode = "r") as h5file:
                self.bamParser.dumpCovTable(h5file.root.profile.coverage, self.bamColNames)
                self.conParser.dumpTnTable(h5file.root.profile.tns, self.kse.kmerCols)
                self.auxParser.dumpAUXTable(h5file.root.profile.aux)
                
    def initBins(self, table):
        """Initialise the bins table
        
        set to 0 for no bin assignment
        """
        for cid in self.contigNames:
            BIN_row = table.row
            BIN_row['cid'] = cid
            BIN_row['bin'] = 0
        table.flush()
    
    def dumpBins(self, table):
        print "-----------------------------------"
        print "Bins table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",row['bin']

###############################################################################
###############################################################################
###############################################################################
###############################################################################

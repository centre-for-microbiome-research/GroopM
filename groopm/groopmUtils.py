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
import string
import re

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################
# MAIN CLASS
###############################################################################
class GMProj:
    """main object to poke when doing gm work"""
    def __init__(self):
        self.conParser = None
        self.bamParser = None
        self.auxParser = None
        self.contigsFile = "unset"
        self.bamFiles = []
        self.dbName = "unset"
        self.bamColNames = []
        self.contigNames = []
        
    def loadDB(self, fileName):
        pass
    
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
# AUX PARSING CLASSES
###############################################################################
class AuxParser:
    """Class for parsing auxilary profile values"""
    def __init__(self): pass
    
    def parse(self, auxFile, table):
        """Do the heavy lifting of parsing"""
        print "Parsing auxilary profile"
        search=re.compile(r"^['\"]").search # line starts with a quote char
        for line in auxFile:
            if(not search(line)):
                fields = line.rstrip().split(",")
                AUX_row = table.row
                # punch in the data
                AUX_row['cid'] = fields[0]
                AUX_row['aux'] = fields[1]
                AUX_row.append()
        table.flush()

    def dumpAUXTable(self, table):
        """dump the guts of the AUX table"""
        print "-----------------------------------"
        print "Aux profile table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",row['aux']
    
###############################################################################
# CONTIG PARSING CLASSES
###############################################################################
class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readfq(self, fp): # this is a generator function
        """https://github.com/lh3"""
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].split()[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break
                
    def parse(self, contigFile, kse, table):
        """Do the heavy lifting of parsing"""
        print "Parsing contigs"        
        con_names = []
        for cid,seq,qual in self.readfq(contigFile):
            con_names.append(cid)
            sig = kse.getKSig(seq.upper())
            # make a new row
            TN_row = table.row
            # punch in the data
            TN_row['cid'] = cid
            for mer in sig.keys():
                TN_row[mer] = sig[mer]
            TN_row.append()
        table.flush()
        return con_names

    def dumpTnTable(self, table, kmerCols):
        """dump the guts of the TN table"""

        print "-----------------------------------"
        print "TNuclSig table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],
            for mer in kmerCols:
                print ",",row[mer],
            print ""

class KmerSigEngine:
    """Simple class for determining kmer signatures"""

    def __init__(self, kLen):
        self.kLen = kLen
        self.compl = string.maketrans('ACGT', 'TGCA')
        self.kmerCols = self.makeKmerColNames()
        self.numMers = len(self.kmerCols)
        
    def makeKmerColNames(self):
        """work out the range of kmers required based on kmer length"""
        # build up the big list
        base_words = ("A","C","G","T")
        out_list = ["A","C","G","T"]
        for i in range(1,self.kLen):
            working_list = []
            for mer in out_list:
                for char in base_words:
                    working_list.append(mer+char)
            out_list = working_list
        
        # pear it down based on lexicographical ordering
        ret_list = []
        for mer in out_list:
            lmer = self.shiftLowLexi(mer)
            if lmer not in ret_list:
                ret_list.append(lmer)
        return ret_list
    
    def shiftLowLexi(self, seq):
        """return the lexicographically lowest form of this sequence"""
        rseq = self.revComp(seq)
        if(seq < rseq):
            return seq
        return rseq
        
    def revComp(self, seq):
        """return the reverse complement of a sequence"""
        return seq.translate(self.compl)[::-1]
    
    def getKSig(self, seq):
        """Work out kmer signature for a nucleotide sequence"""
        sig = dict(zip(self.kmerCols, [0.0] * self.numMers))
        ll = len(seq)
        for i in range(0,ll-self.kLen+1):
            sig[self.shiftLowLexi(seq[i:i+self.kLen])] += 1
        # normalise by length and return
        return dict(zip(self.kmerCols, [ X / ll for X in sig.values()]))


###############################################################################
# BAM PARSING CLASSES
###############################################################################
class BamParser:
    """Parse multiple bam files and write the output to hdf5 """

    def __init__(self): pass
    
    def parse(self, table, bamFiles, bamColNames):
        """ parse multiple bam files and store the results in the main DB
        
        table: a table in an open h5 file like "CID,COV_1,...,COV_n,length"
        bamColNames: names of the COV_x columns
        """
        
        # parse the BAMs
        tmp_storage = {}
        bam_count = 0
        for bf in bamFiles:
            bam_file = None
            try:
                bam_file = pysam.Samfile(bf, 'rb')
                print "Parsing",bamColNames[bam_count],"(",(bam_count+1),"of",len(bamColNames),")"
                self.parseBam(bam_file, bam_count, len(bamColNames), tmp_storage)                
                bam_count += 1
            except:
                print "Unable to open BAM file",bf,"-- did you supply a SAM file instead?"
                raise

        # go through all the contigs sorted by name and write to the DB
        try:
            for cid in sorted(tmp_storage.keys()):
                # make a new row
                cov_row = table.row
                # punch in the data
                cov_row['cid'] = cid
                for i in range(0,len(bamColNames)):
                    cov_row[bamColNames[i]] = tmp_storage[cid][i]
                cov_row['length'] = tmp_storage[cid][len(bamColNames)]
                cov_row.append()
            table.flush()
        except:
            print "Error saving results to DB"
            raise

    def parseBam(self, bamFile, bamCount, numBams, storage):
        """Parse the bam file handle and store the number of reads mapped"""
        for reference, length in zip(bamFile.references, bamFile.lengths):
            c = Counter()
            try:
                bamFile.fetch(reference, 0, length, callback = c )
                num_reads = c.counts
            except ValueError as e:
                print "Could not calculate num reads for:",reference,"in",bf,"\t",e
                raise
            except:
                print "Could not calculate num reads for:",reference,"in",bf, sys.exc_info()[0]
                raise
            
            # save this guy into the tmp_storage
            if(reference not in storage):
                storage[reference] = np.zeros((numBams+1))
                storage[reference][numBams] = length

            # we need to divide the count by the length if we are going
            # to use a normalised coverage
            storage[reference][bamCount] = float(num_reads)/float(length)
        
    def dumpCovTable(self, table, bamColNames):
        """dump the guts of the coverage table"""

        print "-----------------------------------"
        print "Coverage table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",
            for colName in bamColNames:
                print row[colName],",",
            print row['length']

class Counter:
    """Call back for counting aligned reads
    
    Used in conjunction with pysam.fetch 
    """
    counts = 0
    def __call__(self, alignment):
        self.counts += 1

def getBamDescriptor(fullPath):
    """reduce a full path to just the file name minus extension"""
    return os.path.splitext(os.path.basename(fullPath))[0]

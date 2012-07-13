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
import pysam
import os
import numpy as np
import string
import re

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class GMDataManager:
    """Top level class for manipulating GroopM data
    
    Use this class for parsing in raw data into a hdf DB and 
    for reading from and updating same DB
    """
    def __init__(self): pass
        
    def loadData(self, dbFileName, covProfiles=None, auxProfiles=None, bins=None, contigNames=None, kmerSigs=None, condition=''):
        """Load data from a pre parsed DB
        
        Data loaded depends on which containers are passed in
        """
        if('' == condition):
            condition = 'length > 0' # no condition breaks everything!
            
        try:
            with tables.openFile(dbFileName, mode = "r") as h5file:
                # get the metadata
                META_row = h5file.root.meta.meta.read()
                bamColNames = META_row['bamColNames'].split(",")
                merColNames = META_row['merColNames'].split(",")
                numCons = META_row['numCons']
                numBins = META_row['numBins']
                
                for row1 in h5file.root.profile.coverage.where(condition):
                    row2 = h5file.root.profile.kms[row1.nrow]
                    print row1['cid'],row1['length'],row2['AAAA']
                    
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise
                           
    def createDB(self, bamFiles, contigs, aux, dbName, kmerSize=4, dumpAll=False):
        """Main wrapper for parsing all input files"""
        # load all the passed vars
        dbName = dbName
        contigsFile = contigs
        bamFiles = bamFiles.split(",")
        bamColNames = []
        auxFile = aux
        
        kse = KmerSigEngine(kmerSize)
        conParser = ContigParser()
        bamParser = BamParser()
        auxParser = AuxParser()

        # make sure we're only overwriting existing DBs with the users consent
        try:
            with open(dbName) as f:
                option = raw_input(" !!!WARNING!!! Database: "+dbName+" exists.\n" \
                                   " If you continue you *WILL* delete any previous analyses!\n" \
                                   " Overwrite? (y,n) : ")
                print "****************************************************************"
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
                for bf in bamFiles:
                    # assume the file is called somethign like "fred.bam"
                    # we want to rip off the ".bam" part
                    bam_desc = getBamDescriptor(bf)
                    db_desc[bam_desc] = tables.FloatCol()
                    bamColNames.append(bam_desc)
                db_desc['length'] = tables.Int32Col()
                try:
                    COV_table = h5file.createTable(group, 'coverage', db_desc, "Bam based coverage")
                    rows_created = bamParser.parse(COV_table, bamFiles, bamColNames)
                except:
                    print "Error creating coverage table:", sys.exc_info()[0]
                    raise
                
                #------------------------
                # parse contigs 
                #------------------------
                db_desc = { 'cid' : tables.StringCol(64) }
                for mer in kse.kmerCols:
                     db_desc[mer] = tables.FloatCol()
                try:
                    KMER_table = h5file.createTable(group, 'kms', db_desc, "kmer signature")
                except:
                    print "Error creating KMERSIG table:", sys.exc_info()[0]
                    raise
                try:
                    f = open(contigsFile, "r")
                    contigNames = conParser.parse(f, kse, KMER_table)
                    f.close()
                except:
                    print "Could not parse contig file:",contigsFile,sys.exc_info()[0]
                    raise

                # make sure everything is the same size
                if(rows_created != len(contigNames)):
                    raise Exception("Contig file contains: "+str(len(contigNames))+" but BAM files reference "+str(rows_created)+" contigs!")
                
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
                    f = open(auxFile, "r")
                    rows_created = auxParser.parse(f, AUX_table)
                except:
                    print "Could not parse the auxilary profile file:",contigsFile,sys.exc_info()[0]
                    raise

                # make sure everything is the same size
                if(rows_created != len(contigNames)):
                    raise Exception("Contig file contains: "+str(len(contigNames))+" but AUX file references "+str(rows_created)+" contigs!")
                
                #------------------------
                # Add a table for the bins
                #------------------------
                db_desc = { 'cid' : tables.StringCol(64), 'bin' :  tables.Int32Col() }
                try:
                    BIN_table = h5file.createTable(group, 'bin', db_desc, "Bin IDs")
                    self.initBins(BIN_table, contigNames)
                except:
                    print "Error creating BIN table:", sys.exc_info()[0]
                    raise
                
                #------------------------
                # Add metadata
                #------------------------
                group = h5file.createGroup("/", 'meta', 'Project meta data')
                db_desc = {'bamColNames' : tables.StringCol(500),
                           'merColNames' : tables.StringCol(5000),
                           'numCons' : tables.Int32Col(),
                           'numBins' : tables.Int32Col(),
                           }
                try:
                    META_table = h5file.createTable(group, 'meta', db_desc, "Descriptive data")
                    self.initMeta(META_table, str.join(',',bamColNames), str.join(',',kse.kmerCols), len(contigNames))
                except:
                    print "Error creating META table:", sys.exc_info()[0]
                    raise
        except:
            print "Error creating database:", dbName, sys.exc_info()[0]
            raise
        
        print "****************************************************************"
        print "Data loaded successfully!"
        print " ->",len(contigNames),"contigs"
        print " ->",len(bamColNames),"BAM files"
        print "Written to: '"+dbName+"'"
        print "****************************************************************"

        if(dumpAll):
            with tables.openFile(dbName, mode = "r") as h5file:
                bamParser.dumpCovTable(h5file.root.profile.coverage, bamColNames)
                conParser.dumpTnTable(h5file.root.profile.kms, kse.kmerCols)
                auxParser.dumpAUXTable(h5file.root.profile.aux)
                self.dumpBins(h5file.root.profile.bin)
                self.dumpMeta(h5file.root.meta.meta)
                
    def initBins(self, table, contigNames):
        """Initialise the bins table
        
        set to 0 for no bin assignment
        """
        for cid in sorted(contigNames):
            BIN_row = table.row
            BIN_row['cid'] = cid
            BIN_row['bin'] = 0
            BIN_row.append()
        table.flush()
    
    def dumpBins(self, table):
        print "-----------------------------------"
        print "Bins table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",row['bin']

    def initMeta(self, table, bamColNames, merColNames, numCons):
        """Initialise the meta-data table"""
        META_row = table.row
        META_row['bamColNames'] = bamColNames
        META_row['merColNames'] = merColNames
        META_row['numCons'] = numCons
        META_row['numBins'] = 0
        META_row.append()
        table.flush()

    def dumpMeta(self, table):
        print "-----------------------------------"
        print "MetaData table"
        print "-----------------------------------"
        for row in table:
            print row['bamColNames']
            print row['merColNames']
            print row['numCons']
            print row['numBins']

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class AuxParser:
    """Class for parsing auxilary profile values"""
    def __init__(self): pass
    
    def parse(self, auxFile, table):
        """Do the heavy lifting of parsing"""
        print "Parsing auxilary profile"
        search=re.compile(r"^['\"]").search # line starts with a quote char
        tmp_storage = {} # so we can have all tables sorted accordingly
        for line in auxFile:
            if(not search(line)):
                fields = line.rstrip().split(",")
                tmp_storage[fields[0]] = fields[1]

        rows_created = 0
        for cid in sorted(tmp_storage.keys()):     
            AUX_row = table.row
            # punch in the data
            AUX_row['cid'] = cid
            AUX_row['aux'] = tmp_storage[cid]
            AUX_row.append()
            rows_created += 1
        table.flush()
        return rows_created

    def dumpAUXTable(self, table):
        """Dump the guts of the AUX table"""
        print "-----------------------------------"
        print "Aux profile table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",row['aux']
    
###############################################################################
###############################################################################
###############################################################################
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
        tmp_storage = {} # save everything here first so we can sort accordingly
        for cid,seq,qual in self.readfq(contigFile):
            tmp_storage[cid] = kse.getKSig(seq.upper()) 
        
        for cid in sorted(tmp_storage.keys()):     
            con_names.append(cid)
            # make a new row
            KMER_row = table.row
            # punch in the data
            KMER_row['cid'] = cid
            for mer in tmp_storage[cid].keys():
                KMER_row[mer] = tmp_storage[cid][mer]
            KMER_row.append()
        table.flush()
        return con_names

    def dumpTnTable(self, table, kmerCols):
        """Dump the guts of the TN table"""
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
        """Work out the range of kmers required based on kmer length"""
        # build up the big list
        base_words = ("A","C","G","T")
        out_list = ["A","C","G","T"]
        for i in range(1,self.kLen):
            working_list = []
            for mer in out_list:
                for char in base_words:
                    working_list.append(mer+char)
            out_list = working_list
        
        # pare it down based on lexicographical ordering
        ret_list = []
        for mer in out_list:
            lmer = self.shiftLowLexi(mer)
            if lmer not in ret_list:
                ret_list.append(lmer)
        return ret_list
    
    def shiftLowLexi(self, seq):
        """Return the lexicographically lowest form of this sequence"""
        rseq = self.revComp(seq)
        if(seq < rseq):
            return seq
        return rseq
        
    def revComp(self, seq):
        """Return the reverse complement of a sequence"""
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
###############################################################################
###############################################################################
###############################################################################
class BamParser:
    """Parse multiple bam files and write the output to hdf5 """

    def __init__(self): pass
    
    def parse(self, table, bamFiles, bamColNames):
        """Parse multiple bam files and store the results in the main DB
        
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
        rows_created = 0
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
                rows_created += 1
            table.flush()
        except:
            print "Error saving results to DB"
            raise
        return rows_created

    def parseBam(self, bamFile, bamCount, numBams, storage):
        """Parse a bam file (handle) and store the number of reads mapped"""
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
        """Dump the guts of the coverage table"""
        print "-----------------------------------"
        print "Coverage table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",
            for colName in bamColNames:
                print row[colName],",",
            print row['length']

class Counter:
    """AUX: Call back for counting aligned reads
    
    Used in conjunction with pysam.fetch 
    """
    counts = 0
    def __call__(self, alignment):
        self.counts += 1

def getBamDescriptor(fullPath):
    """AUX: Reduce a full path to just the file name minus extension"""
    return os.path.splitext(os.path.basename(fullPath))[0]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

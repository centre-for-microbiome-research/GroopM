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

    NOTE: All tables are kept in the same order indexed by the contig ID
    Tables managed by this class are listed below

    ------------------------
     PROFILES
    group = '/profile'
    ------------------------
    **Kmer Signature**
    table = 'kms'
    'mer1' : tables.FloatCol(pos=1)
    'mer2' : tables.FloatCol(pos=2)
    'mer3' : tables.FloatCol(pos=3)
    ...
    
    **Bam files (coverage profile)**
    table = 'coverage'
    'bam1' : tables.FloatCol(pos=1)
    'bam2' : tables.FloatCol(pos=2)
    'bam3' : tables.FloatCol(pos=3)
    ...
    
    **Aux profile**
    table = 'aux'
    'aux' : tables.FloatCol(pos=1)
    
    ------------------------
     METADATA
    group = '/meta'
    ------------------------
    ** Metadata **
    table = 'meta'
    'bamColNames' : tables.StringCol(500, pos=0)
    'numBams'     : tables.Int32Col(pos=1)
    'merColNames' : tables.StringCol(5000,pos=2)
    'merSize'     : tables.Int32Col(pos=3)
    'numMers'     : tables.Int32Col(pos=4)
    'numCons'     : tables.Int32Col(pos=5)
    'numBins'     : tables.Int32Col(pos=6)
    'clustered'   : tables.BoolCol(pos=7)           # set to true after clustering is complete
    'complete'    : tables.BoolCol(pos=8)           # set to true after clustering finishing is complete
    
    ** Bins **
    table = 'bin'
    'cid'    : tables.StringCol(500, pos=0)
    'bin'    : tables.Int32Col(pos=1)
    'length' : tables.Int32Col(pos=2)
    """
    def __init__(self): pass

#------------------------------------------------------------------------------
# DB CREATION / INITIALISATION 

    def createDB(self, bamFiles, contigs, aux, dbFileName, kmerSize=4, dumpAll=False, force=False):
        """Main wrapper for parsing all input files"""
        # load all the passed vars
        dbFileName = dbFileName
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
            with open(dbFileName) as f:
                if(not force):
                    option = raw_input(" ****WARNING**** Database: '"+dbFileName+"' exists.\n" \
                                       " If you continue you *WILL* delete any previous analyses!\n" \
                                       " Overwrite? (y,n) : ")
                    print "****************************************************************"
                    if(option.upper() != "Y"):
                        print "Operation cancelled"
                        return False
                    else:
                        print "Overwriting database",dbFileName
        except IOError as e:
            print "Creating new database", dbFileName
        
        # create the db
        try:        
            with tables.openFile(dbFileName, mode = "w", title = "GroopM") as h5file:
                # Create groups under "/" (root) for storing profile information and metadata
                profile_group = h5file.createGroup("/", 'profile', 'Assembly profiles')
                meta_group = h5file.createGroup("/", 'meta', 'Project meta data')
                #------------------------
                # parse contigs and make kmer sigs
                #
                # Contig IDs are key. Any keys existing in other files but not in this file will be
                # ignored. Any missing keys in other files will be given the default profile value 
                # (typically 0). Ironically, we don't store the CIDs here, these are saved one time
                # only in the bin table 
                #------------------------
                db_desc = {}
                ppos = 0
                for mer in kse.kmerCols:
                     db_desc[mer] = tables.FloatCol(pos=ppos)
                     ppos += 1
                try:
                    KMER_table = h5file.createTable(profile_group, 'kms', db_desc, "kmer signature")
                except:
                    print "Error creating KMERSIG table:", sys.exc_info()[0]
                    raise
                try:
                    f = open(contigsFile, "r")
                    # keep all the contig names so we can check other tables
                    # contigNames is a dict of type ID -> Length
                    contigNames = conParser.parse(f, kse, KMER_table)
                    f.close()
                except:
                    print "Could not parse contig file:",contigsFile,sys.exc_info()[0]
                    raise
                #------------------------
                # Add a table for the bins
                #------------------------
                db_desc = { 'cid' : tables.StringCol(64, pos=0), 'bin' : tables.Int32Col(pos=1), 'length' : tables.Int32Col(pos=2) }
                try:
                    BIN_table = h5file.createTable(meta_group, 'bin', db_desc, "Bin IDs")
                    self.initBins(BIN_table, contigNames)
                except:
                    print "Error creating BIN table:", sys.exc_info()[0]
                    raise
                #------------------------
                # Add metadata
                #------------------------
                # Create a new group under "/" (root) for storing profile information
                db_desc = {'bamColNames' : tables.StringCol(500, pos=0),
                           'numBams' : tables.Int32Col(pos=1),
                           'merColNames' : tables.StringCol(5000,pos=2),
                           'merSize' : tables.Int32Col(pos=2),
                           'numMers' : tables.Int32Col(pos=4),
                           'numCons' : tables.Int32Col(pos=5),
                           'numBins' : tables.Int32Col(pos=6),
                           'clustered' : tables.BoolCol(pos=7),                  # set to true after clustering is complete
                           'complete' : tables.BoolCol(pos=8)                    # set to true after clustering finishing is complete
                           }
                try:
                    META_table = h5file.createTable(meta_group, 'meta', db_desc, "Descriptive data")
                    self.initMeta(META_table, str.join(',',bamColNames), len(bamColNames), str.join(',',kse.kmerCols), kmerSize, len(kse.kmerCols), len(contigNames))
                except:
                    print "Error creating META table:", sys.exc_info()[0]
                    raise
                #------------------------
                # parse bam files
                #------------------------
                # build a table template based on the number of bamfiles we have
                db_desc = {}
                ppos = 0
                for bf in bamFiles:
                    # assume the file is called something like "fred.bam"
                    # we want to rip off the ".bam" part
                    bam_desc = getBamDescriptor(bf)
                    db_desc[bam_desc] = tables.FloatCol(pos=ppos)
                    bamColNames.append(bam_desc)
                    ppos += 1
                try:
                    COV_table = h5file.createTable(profile_group, 'coverage', db_desc, "Bam based coverage")
                    bamParser.parse(bamFiles, bamColNames, COV_table, contigNames)
                except:
                    print "Error creating coverage table:", sys.exc_info()[0]
                    raise
                #------------------------
                # parse aux profile 
                #------------------------
                db_desc = { 'aux' : tables.FloatCol(pos=0) }
                try:
                    AUX_table = h5file.createTable(profile_group, 'aux', db_desc, "Secondary profile")
                except:
                    print "Error creating AUX table:", sys.exc_info()[0]
                    raise
                try:
                    f = open(auxFile, "r")
                    auxParser.parse(f, AUX_table, contigNames)
                except:
                    print "Could not parse the auxilary profile file:",contigsFile,sys.exc_info()[0]
                    raise
        except:
            print "Error creating database:", dbFileName, sys.exc_info()[0]
            raise
        
        print "****************************************************************"
        print "Data loaded successfully!"
        print " ->",len(contigNames),"contigs"
        print " ->",len(bamColNames),"BAM files"
        print "Written to: '"+dbFileName+"'"
        print "****************************************************************"

        if(dumpAll):
            self.dumpAll(dbFileName)
            
        # all good!
        return True
                
    def initBins(self, table, contigNames):
        """Initialise the bins table
        
        set to 0 for no bin assignment
        """
        for cid in sorted(contigNames):
            BIN_row = table.row
            BIN_row['cid'] = cid
            BIN_row['length'] = contigNames[cid]
            BIN_row['bin'] = 0
            BIN_row.append()
        table.flush()
    
    def initMeta(self, table, bamColNames, numBams, merColNames, merSize, numMers, numCons):
        """Initialise the meta-data table"""
        META_row = table.row
        META_row['bamColNames'] = bamColNames
        META_row['numBams'] = numBams
        META_row['merColNames'] = merColNames
        META_row['merSize'] = merSize
        META_row['numMers'] = numMers
        META_row['numCons'] = numCons
        META_row['numBins'] = 0
        META_row['clustered'] = False
        META_row['complete'] = False
        META_row.append()
        table.flush()
        
#------------------------------------------------------------------------------
# GET / SET DATA TABLES 

    def getConditionalIndicies(self, dbFileName, condition=''):
        """return the indicies into the db which meet the condition"""
        if('' == condition):
            condition = "cid != ''" # no condition breaks everything!
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return np.array([x.nrow for x in h5file.root.meta.bin.where(condition)])
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def getCoverageProfiles(self, dbFileName, condition='', indicies=np.array([])):
        """Load coverage profiles"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([list(h5file.root.profile.coverage[x]) for x in indicies])
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(h5file.root.profile.coverage[x.nrow]) for x in h5file.root.meta.bin.where(condition)])
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise
        
    def getAuxProfiles(self, dbFileName, condition='', indicies=np.array([])):
        """Load aux profiles"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([list(h5file.root.profile.aux[x])[0] for x in indicies]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(h5file.root.profile.aux[x.nrow])[0] for x in h5file.root.meta.bin.where(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def getBins(self, dbFileName, condition='', indicies=np.array([])):
        """Load bins"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([h5file.root.meta.bin[x][1] for x in indicies]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[1] for x in h5file.root.meta.bin.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def setBins(self, dbFileName, updates):
        """Set bin stats
        
        updates is a dictionary which looks like:
        { tableRow : binValue }
        """
        row_nums = updates.keys()
        new_bins = updates.values()
        try:
            with tables.openFile(dbFileName, mode='a') as h5file:
                table = h5file.root.meta.bin
                for row_num in updates.keys():
                    new_row = np.zeros((1,),dtype=('S64,i4,i4'))
                    new_row[:] = [(table[row_num][0],updates[row_num],table[row_num][2])]
                    table.modifyRows(start=row_num, rows=new_row)
                table.flush()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def getContigNames(self, dbFileName, condition='', indicies=np.array([])):
        """Load contig names"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([h5file.root.meta.bin[x][0] for x in indicies]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[0] for x in h5file.root.meta.bin.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def getContigLengths(self, dbFileName, condition='', indicies=np.array([])):
        """Load contig lengths"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([h5file.root.meta.bin[x][2] for x in indicies]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[2] for x in h5file.root.meta.bin.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def getKmerSigs(self, dbFileName, condition='', indicies=np.array([])):
        """Load kmer sigs"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([list(h5file.root.profile.kms[x]) for x in indicies])
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(h5file.root.profile.kms[x.nrow]) for x in h5file.root.meta.bin.where(condition)])
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

#------------------------------------------------------------------------------
# GET / SET WORKFLOW FLAGS 

    def isClustered(self, dbFileName):
        """Has this data set been clustered?"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return h5file.root.meta.meta.read()['clustered']
        except:
            print "Error opening database:", dbFileName, sys.exc_info()[0]
            raise
            
    def setClustered(self, dbFileName, state=True):
        """Set the state of clustering"""
        try:
            with tables.openFile(dbFileName, mode='a') as h5file:
                META_table = h5file.root.meta.meta
                for META_row in META_table: # there is only one meta row
                    META_row['clustered'] = state
                    META_row.update()
                META_table.flush()
        except:
            print "Error opening database:", dbFileName, sys.exc_info()[0]
            raise
            
    def isComplete(self, dbFileName):
        """Has this data set been clustered?"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return h5file.root.meta.meta.read()['complete']
        except:
            print "Error opening database:", dbFileName, sys.exc_info()[0]
            raise

    def setComplete(self, dbFileName, state=True):
        """Set the state of completion"""
        try:
            with tables.openFile(dbFileName, mode='a') as h5file:
                META_table = h5file.root.meta.meta
                for META_row in META_table: # there is only one meta row
                    META_row['complete'] = state
                    META_row.update()
                META_table.flush()
        except:
            print "Error opening database:", dbFileName, sys.exc_info()[0]
            raise

#------------------------------------------------------------------------------
# FILE / IO 

    def dumpBins(self, table):
        print "-----------------------------------"
        print "Bins table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",row['length'],",",row['bin']

    def dumpMeta(self, table):
        print "-----------------------------------"
        print "MetaData table"
        print "-----------------------------------"
        for row in table:
            print row['bamColNames']
            print row['merColNames']
            print row['merSize']
            print row['numMers']
            print row['numCons']
            print row['numBins']
            print row['clustered']
            print row['complete']

    def dumpAll(self, dbFileName):
        """Dump all contents of all DBs to screen"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                
                # get the metadata
                META_row = h5file.root.meta.meta.read()
                print META_row['bamColNames']
    
                bamColNames = META_row['bamColNames'][0].split(",")
                merSize = META_row['merSize']
     
                kse = KmerSigEngine(merSize)
                conParser = ContigParser()
                bamParser = BamParser()
                auxParser = AuxParser()
    
                bamParser.dumpCovTable(h5file.root.profile.coverage, bamColNames)
                conParser.dumpTnTable(h5file.root.profile.kms, kse.kmerCols)
                auxParser.dumpAUXTable(h5file.root.profile.aux)
                self.dumpBins(h5file.root.meta.bin)
                self.dumpMeta(h5file.root.meta.meta)
        except:
            print "Error opening database:", dbFileName, sys.exc_info()[0]
            raise

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class AuxParser:
    """Class for parsing auxilary profile values"""
    def __init__(self): pass
    
    def parse(self, auxFile, table, contigNames):
        """Do the heavy lifting of parsing AUX data"""
        print "Parsing auxilary profile"
        search=re.compile(r"^['\"]").search # line starts with a quote char
        tmp_storage = {} # so we can have all tables sorted accordingly
        for cid in contigNames.keys():
            tmp_storage[cid] = 0.0
        
        for line in auxFile:
            if(not search(line)):
                fields = line.rstrip().split(",")
                cid = fields[0]
                if cid in contigNames:
                    tmp_storage[fields[0]] = fields[1]

        rows_created = 0
        for cid in sorted(tmp_storage.keys()):     
            AUX_row = table.row
            # punch in the data
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
            print row['aux']
    
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
        tmp_storage = {} # save everything here first so we can sort accordingly
        for cid,seq,qual in self.readfq(contigFile):
            tmp_storage[cid] = (kse.getKSig(seq.upper()), len(seq)) 
        
        con_names = {}
        for cid in sorted(tmp_storage.keys()):     
            con_names[cid] = tmp_storage[cid][1]
            # make a new row
            KMER_row = table.row
            # punch in the data
            for mer in tmp_storage[cid][0].keys():
                KMER_row[mer] = tmp_storage[cid][0][mer]
            KMER_row.append()
        table.flush()
        return con_names

    def dumpTnTable(self, table, kmerCols):
        """Dump the guts of the TN table"""
        print "-----------------------------------"
        print "TNuclSig table"
        print "-----------------------------------"
        for row in table:
            for mer in kmerCols:
                print ",",row[mer],
            print ""

###############################################################################
###############################################################################
###############################################################################
###############################################################################
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
            this_mer = self.shiftLowLexi(seq[i:i+self.kLen])
            if this_mer in sig:
                sig[this_mer] += 1
        # normalise by length and return
        return dict(zip(self.kmerCols, [ X / ll for X in sig.values()]))


###############################################################################
###############################################################################
###############################################################################
###############################################################################
class BamParser:
    """Parse multiple bam files and write the output to hdf5 """

    def __init__(self): pass
    
    def parse(self, bamFiles, bamColNames, table, contigNames):
        """Parse multiple bam files and store the results in the main DB
        
        table: a table in an open h5 file like "CID,COV_1,...,COV_n,length"
        bamColNames: names of the COV_x columns
        """
        # parse the BAMs
        # we need to have some type of entry for each contig
        # so start by putting 0's here
        tmp_storage = {}
        num_bams = len(bamColNames)
        for cid in contigNames.keys():
            tmp_storage[cid] = np.zeros((num_bams))

        bam_count = 0
        for bf in bamFiles:
            bam_file = None
            try:
                bam_file = pysam.Samfile(bf, 'rb')
                print "Parsing",bamColNames[bam_count],"(",(bam_count+1),"of",num_bams,")"
                self.parseBam(bam_file, bam_count, tmp_storage, contigNames)                
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
                for i in range(0,len(bamColNames)):
                    cov_row[bamColNames[i]] = tmp_storage[cid][i]
                cov_row.append()
                rows_created += 1
            table.flush()
        except:
            print "Error saving results to DB"
            raise
        return rows_created

    def parseBam(self, bamFile, bamCount, storage, contigNames):
        """Parse a bam file (handle) and store the number of reads mapped"""
        for reference, length in zip(bamFile.references, bamFile.lengths):
            if(reference in contigNames): # we only care about contigs we have seen IN
                c = Counter()             # the fasta file during contig parsing
                try:
                    bamFile.fetch(reference, 0, length, callback = c )
                    num_reads = c.counts
                except ValueError as e:
                    print "Could not calculate num reads for:",reference,"in",bf,"\t",e
                    raise
                except:
                    print "Could not calculate num reads for:",reference,"in",bf, sys.exc_info()[0]
                    raise
                
                # we have already made storage for this guy above so we can gaurantee
                # there is space to save it!
    
                # we need to divide the count by the length if we are going
                # to use a normalised coverage
                storage[reference][bamCount] = float(num_reads)/float(length)
        
    def dumpCovTable(self, table, bamColNames):
        """Dump the guts of the coverage table"""
        print "-----------------------------------"
        print "Coverage table"
        print "-----------------------------------"
        for row in table:
            for colName in bamColNames:
                print ",",row[colName],
            print ""

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

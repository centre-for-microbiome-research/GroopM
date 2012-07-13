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
import pysam
import os
import numpy as np
import string
import re

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

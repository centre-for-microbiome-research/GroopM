#!/usr/bin/env python
###############################################################################
#                                                                             #
#    mstore.py                                                                #
#                                                                             #
#    GroopM - Low level data management and file parsing                      #
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
__version__ = "0.2.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Alpha"

###############################################################################

from sys import exc_info, exit
from os.path import splitext as op_splitext, basename as op_basename
from string import maketrans as s_maketrans

import tables
import numpy as np
import pysam

# GroopM imports
import groopmExceptions as ge
import groopmTimekeeper as gtime

np.seterr(all='raise')     

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
    
    **Coverage profile**
    table = 'coverage'
    'stoit1' : tables.FloatCol(pos=1)
    'stoit2' : tables.FloatCol(pos=2)
    'stoit3' : tables.FloatCol(pos=3)
    ...

    ------------------------
     LINKS
    group = '/links'
    ------------------------
    ** Links **
    table = 'links'
    'contig1'    : tables.Int32Col(pos=0)            # reference to index in meta/contigs
    'contig2'    : tables.Int32Col(pos=1)            # reference to index in meta/contigs
    'numReads'   : tables.Int32Col(pos=2)            # number of reads supporting this link 
    'linkType'   : tables.Int32Col(pos=3)            # the type of the link (SS, SE, ES, EE)
    'gap'        : tables.Int32Col(pos=4)            # the estimated gap between the contigs
    
    ------------------------
     METADATA
    group = '/meta'
    ------------------------
    ** Metadata **
    table = 'meta'
    'stoitColNames' : tables.StringCol(512, pos=0)
    'numStoits'     : tables.Int32Col(pos=1)
    'merColNames'   : tables.StringCol(4096,pos=2)
    'merSize'       : tables.Int32Col(pos=3)
    'numMers'       : tables.Int32Col(pos=4)
    'numCons'       : tables.Int32Col(pos=5)
    'numBins'       : tables.Int32Col(pos=6)
    'clustered'     : tables.BoolCol(pos=7)           # set to true after clustering is complete
    'complete'      : tables.BoolCol(pos=8)           # set to true after clustering finishing is complete

    ** Contigs **
    table = 'contigs'
    'cid'    : tables.StringCol(512, pos=0)
    'bid'    : tables.Int32Col(pos=1)
    'length' : tables.Int32Col(pos=2)
    
    ** Bins **
    table = 'bin'
    'bid'        : tables.Int32Col(pos=0)
    'numMembers' : tables.Int32Col(pos=1)

    """
    def __init__(self): pass

#------------------------------------------------------------------------------
# DB CREATION / INITIALISATION  - PROFILES

    def createDB(self, bamFiles, contigs, dbFileName, kmerSize=4, dumpAll=False, force=False):
        """Main wrapper for parsing all input files"""
        # load all the passed vars
        dbFileName = dbFileName
        contigsFile = contigs
        stoitColNames = []
        
        kse = KmerSigEngine(kmerSize)
        conParser = ContigParser()
        bamParser = BamParser()

        cid_2_indices ={}
        
        # make sure we're only overwriting existing DBs with the users consent
        try:
            with open(dbFileName) as f:
                if(not force):
                    user_option = self.promptOnOverwrite(dbFileName)
                    if(user_option != "Y"):
                        print "Operation cancelled"
                        return False
                    else:
                        print "Overwriting database",dbFileName
        except IOError as e:
            print "Creating new database", dbFileName
        
        # create the db
        timer = gtime.TimeKeeper()
        try:        
            with tables.openFile(dbFileName, mode = "w", title = "GroopM") as h5file:
                # Create groups under "/" (root) for storing profile information and metadata
                profile_group = h5file.createGroup("/", 'profile', 'Assembly profiles')
                meta_group = h5file.createGroup("/", 'meta', 'Associated metadata')
                links_group = h5file.createGroup("/", 'links', 'Paired read link information')
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
                ksig_data = None
                for mer in kse.kmerCols:
                     db_desc[mer] = tables.FloatCol(pos=ppos)
                     ppos += 1
                try:
                    f = open(contigsFile, "r")
                    # keep all the contig names so we can check other tables
                    # contigNames is a dict of type ID -> Length
                    (ksig_data, contigNames) = conParser.parse(f, kse)
                    f.close()
                except:
                    print "Could not parse contig file:",contigsFile,exc_info()[0]
                    raise

                try:
                    KMER_table = h5file.createTable(profile_group, 'kms', db_desc, "Kmer signature", expectedrows=len(contigNames))
                except:
                    print "Error creating KMERSIG table:", exc_info()[0]
                    raise
                try:
                    conParser.storeSigs(ksig_data, KMER_table)
                except:
                    print "Could not load kmer sigs:",contigsFile,exc_info()[0]
                    raise

                #------------------------
                # Add a table for the contigs
                #------------------------
                db_desc = {'cid' : tables.StringCol(512, pos=0),
                           'bid' : tables.Int32Col(dflt=0,pos=1),
                           'length' : tables.Int32Col(pos=2) }
                try:
                    CONTIG_table = h5file.createTable(meta_group, 'contigs', db_desc, "Contig information", expectedrows=len(contigNames))
                    cid_2_indices = self.initContigs(CONTIG_table, contigNames)
                except:
                    print "Error creating CONTIG table:", exc_info()[0]
                    raise

                #------------------------
                # Add a table for the bins
                #------------------------
                db_desc = {'bid' : tables.Int32Col(pos=0),
                           'numMembers' : tables.Int32Col(dflt=0,pos=1) }
                try:
                    BIN_table = h5file.createTable(meta_group, 'bins', db_desc, "Bin information")
                    BIN_table.flush()
                except:
                    print "Error creating BIN table:", exc_info()[0]
                    raise
                print "    %s" % timer.getTimeStamp()

                #------------------------
                # parse bam files
                #------------------------
                # build a table template based on the number of bamfiles we have
                db_desc = {}
                ppos = 0
                links = {}
                for bf in bamFiles:
                    # assume the file is called something like "fred.bam"
                    # we want to rip off the ".bam" part
                    bam_desc = getBamDescriptor(bf)
                    db_desc[bam_desc] = tables.FloatCol(pos=ppos)
                    stoitColNames.append(bam_desc)
                    ppos += 1
                try:
                    COV_table = h5file.createTable(profile_group, 'coverage', db_desc, "Bam based coverage", expectedrows=len(contigNames))
                except:
                    print "Error creating coverage table:", exc_info()[0]
                    raise

                # now fill in coverage details
                rowwise_links = bamParser.parse(bamFiles, stoitColNames, COV_table, contigNames, cid_2_indices)
                
                #------------------------
                # contig links
                #------------------------
                # set table size according to the number of links returned from
                # the previous call
                db_desc = {'contig1' : tables.Int32Col(pos=0),
                           'contig2' : tables.Int32Col(pos=1),
                           'numReads' : tables.Int32Col(pos=2),
                           'linkType' : tables.Int32Col(pos=3),
                           'gap' : tables.Int32Col(pos=4) }
                try:
                    LINKS_table = h5file.createTable(links_group, 'links', db_desc, "ContigLinks", expectedrows=len(rowwise_links))
                    bamParser.initLinks(rowwise_links, LINKS_table)
                except:
                    print "Error creating links table:", exc_info()[0]
                    raise
                print "    %s" % timer.getTimeStamp()
                
                #------------------------
                # Add metadata
                #------------------------
                # Create a new group under "/" (root) for storing profile information
                db_desc = {'stoitColNames' : tables.StringCol(512, pos=0),
                           'numStoits' : tables.Int32Col(pos=1),
                           'merColNames' : tables.StringCol(4096,pos=2),
                           'merSize' : tables.Int32Col(pos=2),
                           'numMers' : tables.Int32Col(pos=4),
                           'numCons' : tables.Int32Col(pos=5),
                           'numBins' : tables.Int32Col(dflt=0, pos=6),
                           'clustered' : tables.BoolCol(dflt=False, pos=7),                  # set to true after clustering is complete
                           'complete' : tables.BoolCol(dflt=False, pos=8)                    # set to true after clustering finishing is complete
                           }
                try:
                    META_table = h5file.createTable(meta_group, 'meta', db_desc, "Descriptive data", expectedrows=1)
                    self.initMeta(META_table, str.join(',',stoitColNames), len(stoitColNames), str.join(',',kse.kmerCols), kmerSize, len(kse.kmerCols), len(contigNames))
                except:
                    print "Error creating META table:", exc_info()[0]
                    raise
        except:
            print "Error creating database:", dbFileName, exc_info()[0]
            raise
        
        print "****************************************************************"
        print "Data loaded successfully!"
        print " ->",len(contigNames),"contigs"
        print " ->",len(stoitColNames),"BAM files"
        print "Written to: '"+dbFileName+"'"
        print "****************************************************************"
        print "    %s" % timer.getTimeStamp()

        if(dumpAll):
            self.dumpAll(dbFileName)
            
        # all good!
        return True
                
    def initContigs(self, table, contigNames):
        """Initialise the contigs table
        
        set to 0 for no bin assignment
        """
        cid_2_indices = {}
        outer_index = 0
        for cid in sorted(contigNames):
            CONTIG_row = table.row
            CONTIG_row['cid'] = cid
            CONTIG_row['length'] = contigNames[cid]
            CONTIG_row.append()
            
            cid_2_indices[cid] = outer_index
            outer_index += 1
        table.flush()
        return cid_2_indices 
    
    def initMeta(self, table, stoitColNames, numStoits, merColNames, merSize, numMers, numCons):
        """Initialise the meta-data table"""
        META_row = table.row
        META_row['stoitColNames'] = stoitColNames
        META_row['numStoits'] = numStoits
        META_row['merColNames'] = merColNames
        META_row['merSize'] = merSize
        META_row['numMers'] = numMers
        META_row['numCons'] = numCons
        META_row.append()
        table.flush()

    def promptOnOverwrite(self, dbFileName, minimal=False):
        """Check that the user is ok with overwriting the db"""
        input_not_ok = True
        valid_responses = ['Y','N']
        vrs = ",".join([str.lower(str(x)) for x in valid_responses])
        while(input_not_ok):
            if(minimal):
                option = raw_input(" Overwrite? ("+vrs+") : ")
            else: 
                
                option = raw_input(" ****WARNING**** Database: '"+dbFileName+"' exists.\n" \
                                   " If you continue you *WILL* delete any previous analyses!\n" \
                                   " Overwrite? ("+vrs+") : ")
            if(option.upper() in valid_responses):
                print "****************************************************************"
                return option.upper()
            else:
                print "Error, unrecognised choice '"+option.upper()+"'"
                minimal = True

#------------------------------------------------------------------------------
# GET LINKS 

    def restoreLinks(self, dbFileName, indices=[]):
        """Restore the links hash for a given set of indicies"""
        full_record = []
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                full_record = [list(x) for x in h5file.root.links.links.readWhere("contig1 >= 0")]
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise
        
        if indices == []:
            # get all!
            indices = self.getConditionalIndicies(dbFileName)
        
        links_hash = {}
        if full_record != []:
            for record in full_record:
                # make sure we have storage
                if record[0] in indices and record[1] in indices:
                    try:
                        links_hash[record[0]].append(record[1:])
                    except KeyError:
                        links_hash[record[0]] = [record[1:]]
        return links_hash
    
#------------------------------------------------------------------------------
# GET / SET DATA TABLES - PROFILES 

    def getConditionalIndicies(self, dbFileName, condition=''):
        """return the indices into the db which meet the condition"""
        if('' == condition):
            condition = "cid != ''" # no condition breaks everything!
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return np.array([x.nrow for x in h5file.root.meta.contigs.where(condition)])
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getCoverageProfiles(self, dbFileName, condition='', indices=np.array([])):
        """Load coverage profiles"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([list(h5file.root.profile.coverage[x]) for x in indices])
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(h5file.root.profile.coverage[x.nrow]) for x in h5file.root.meta.contigs.where(condition)])
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def nukeBins(self, dbFileName):
        """Reset all bin information, completely"""
        print "    Clearing all old bin information from",dbFileName
        contig_names = {}
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                # get a dictionary of contigIDs to lengths
                contig_names = dict(zip(
                                        [list(x)[0] for x in h5file.root.meta.contigs.readWhere("cid != ''")],
                                        [list(x)[2] for x in h5file.root.meta.contigs.readWhere("cid != ''")]
                                        )
                                    )
                
            num_stoits = self.getNumStoits(dbFileName)
            mer_col_names = self.getMerColNames(dbFileName)
            mer_size = self.getMerSize(dbFileName)
            num_mers = self.getNumMers(dbFileName)
            num_cons = self.getNumCons(dbFileName)
            stoit_col_names = self.getStoitColNames(dbFileName)
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise
        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/meta") as meta_group:
                # clear bin stats
                # try remove any older failed attempts
                try:
                    meta_group.removeNode('/', 'tmp_bins')
                except:
                    pass
                # make a new tmp table
                db_desc = {'bid' : tables.Int32Col(pos=0), 'numMembers' : tables.Int32Col(dflt=0,pos=1) }
                BIN_table = meta_group.createTable('/', 'tmp_bins', db_desc, "Bin information")
                # rename as the bins table
                meta_group.renameNode('/', 'bins', 'tmp_bins', overwrite=True)       

                # clear contig bin ids
                # try remove any older failed attempts
                try:
                    meta_group.removeNode('/', 'tmp_contigs')
                except:
                    pass
                # make a new tmp table
                db_desc = {'cid' : tables.StringCol(512, pos=0),
                           'bid' : tables.Int32Col(dflt=0,pos=1),
                           'length' : tables.Int32Col(pos=2) }
                CONTIG_table = meta_group.createTable('/', 'tmp_contigs', db_desc, "Contig information", expectedrows=len(contig_names))
                cid_2_indices = self.initContigs(CONTIG_table, contig_names)
                # do the rename
                meta_group.renameNode('/', 'contigs', 'tmp_contigs', overwrite=True)
                
                # we need to reset the num_bins region of meta
                try:
                    meta_group.removeNode('/', 'tmp_meta')
                except:
                    pass
                # make a new tmp table
                # Create a new group under "/" (root) for storing profile information
                db_desc = {'stoitColNames' : tables.StringCol(512, pos=0),
                           'numStoits' : tables.Int32Col(pos=1),
                           'merColNames' : tables.StringCol(4096,pos=2),
                           'merSize' : tables.Int32Col(pos=2),
                           'numMers' : tables.Int32Col(pos=4),
                           'numCons' : tables.Int32Col(pos=5),
                           'numBins' : tables.Int32Col(dflt=0, pos=6),
                           'clustered' : tables.BoolCol(dflt=False, pos=7),                  # set to true after clustering is complete
                           'complete' : tables.BoolCol(dflt=False, pos=8)                    # set to true after clustering finishing is complete
                           }

                META_table = meta_group.createTable('/', 'tmp_meta', db_desc, "Descriptive data", expectedrows=1)
                self.initMeta(META_table, stoit_col_names, num_stoits, mer_col_names, mer_size, num_mers, num_cons)
                meta_group.renameNode('/', 'meta', 'tmp_meta', overwrite=True)
                                
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise


    def setBinStats(self, dbFileName, updates):
        """Set bins table 
        
        updates is a dictionary which looks like:
        { bid : numMembers }
        """
        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/meta") as meta_group:
                # pytables is a little "dumb" thus it's easier to 
                # make a new table and then copy everything over
                
                # try remove any older failed attempts
                try:
                    meta_group.removeNode('/', 'tmp_bins')
                except:
                    pass
                # make a new tmp table
                db_desc = {'bid' : tables.Int32Col(pos=0), 'numMembers' : tables.Int32Col(dflt=0,pos=1) }
                BIN_table = meta_group.createTable('/', 'tmp_bins', db_desc, "Bin information")
                # add in the new stuff
                for bid in updates.keys(): 
                    BIN_row = BIN_table.row
                    BIN_row['bid'] = bid
                    BIN_row['numMembers'] = updates[bid]
                    BIN_row.append()
                BIN_table.flush()
                
                # do the rename
                meta_group.renameNode('/', 'bins', 'tmp_bins', overwrite=True)
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getBinStats(self, dbFileName):
        """Load data from bins table
        
        Returns a dict of type:
        { bid : numMembers }
        """
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                ret_dict = {}
                all_rows = h5file.root.meta.bins.read()
                for row in all_rows:
                    ret_dict[row[0]] = row[1]
                return ret_dict
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise
        return {}
        
                
    def getBins(self, dbFileName, condition='', indices=np.array([])):
        """Load per-contig bins"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([h5file.root.meta.contigs[x][1] for x in indices]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[1] for x in h5file.root.meta.contigs.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def setBinAssignments(self, dbFileName, assignments):
        """Set per-contig bins
        
        updates is a dictionary which looks like:
        { tableRow : binValue }
        """
        row_nums = assignments.keys()
        try:
            with tables.openFile(dbFileName, mode='a') as h5file:
                table = h5file.root.meta.contigs
                for row_num in row_nums:
                    new_row = np.zeros((1,),dtype=('S512,i4,i4'))
                    new_row[:] = [(table[row_num][0],assignments[row_num],table[row_num][2])]
                    table.modifyRows(start=row_num, rows=new_row)
                table.flush()
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getContigNames(self, dbFileName, condition='', indices=np.array([])):
        """Load contig names"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([h5file.root.meta.contigs[x][0] for x in indices]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[0] for x in h5file.root.meta.contigs.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getContigLengths(self, dbFileName, condition='', indices=np.array([])):
        """Load contig lengths"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([h5file.root.meta.contigs[x][2] for x in indices]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[2] for x in h5file.root.meta.contigs.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getKmerSigs(self, dbFileName, condition='', indices=np.array([])):
        """Load kmer sigs"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([list(h5file.root.profile.kms[x]) for x in indices])
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(h5file.root.profile.kms[x.nrow]) for x in h5file.root.meta.contigs.where(condition)])
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getMetaField(self, dbFileName, fieldName):
        """return the value of fieldName in the metadata tables"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                # theres only one value
                return h5file.root.meta.meta.read()[fieldName][0]
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getNumStoits(self, dbFileName):
        """return the value of numStoits in the metadata tables"""
        return self.getMetaField(dbFileName, 'numStoits')
            
    def getMerColNames(self, dbFileName):
        """return the value of merColNames in the metadata tables"""
        return self.getMetaField(dbFileName, 'merColNames')
            
    def getMerSize(self, dbFileName):
        """return the value of merSize in the metadata tables"""
        return self.getMetaField(dbFileName, 'merSize')

    def getNumMers(self, dbFileName):
        """return the value of numMers in the metadata tables"""
        return self.getMetaField(dbFileName, 'numMers')

    def getNumCons(self, dbFileName):
        """return the value of numCons in the metadata tables"""
        return self.getMetaField(dbFileName, 'numCons')

    def getNumBins(self, dbFileName):
        """return the value of numBins in the metadata tables"""
        return self.getMetaField(dbFileName, 'numBins')
        
    def setNumBins(self, dbFileName, numBins):
        """set the number of bins"""
        try:
            with tables.openFile(dbFileName, mode='a') as h5file:
                META_table = h5file.root.meta.meta
                for META_row in META_table: # there is only one meta row
                    META_row['numBins'] = numBins
                    META_row.update()
                META_table.flush()
        except:
            print "Error opening database:", dbFileName, exc_info()[0]
            raise
        
    def getStoitColNames(self, dbFileName):
        """return the value of stoitColNames in the metadata tables"""
        return self.getMetaField(dbFileName, 'stoitColNames')

#------------------------------------------------------------------------------
# GET / SET WORKFLOW FLAGS 

    def isClustered(self, dbFileName):
        """Has this data set been clustered?"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return h5file.root.meta.meta.read()['clustered']
        except:
            print "Error opening database:", dbFileName, exc_info()[0]
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
            print "Error opening database:", dbFileName, exc_info()[0]
            raise
            
    def isComplete(self, dbFileName):
        """Has this data set been *completely* clustered?"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return h5file.root.meta.meta.read()['complete']
        except:
            print "Error opening database:", dbFileName, exc_info()[0]
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
            print "Error opening database:", dbFileName, exc_info()[0]
            raise

#------------------------------------------------------------------------------
# FILE / IO 

    def dumpContigs(self, table):
        """Raw dump of contig information"""
        print "-----------------------------------"
        print "Contigs table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",row['length'],",",row['bid']

    def dumpbins(self, table):
        """Raw dump of bin information"""
        print "-----------------------------------"
        print "Bins table"
        print "-----------------------------------"

    def dumpMeta(self, table):
        """Raw dump of metadata"""
        print "-----------------------------------"
        print "MetaData table"
        print "-----------------------------------"
        for row in table:
            print row['stoitColNames']
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
                print META_row['stoitColNames']
    
                stoitColNames = META_row['stoitColNames'][0].split(",")
                merSize = META_row['merSize']
     
                kse = KmerSigEngine(merSize)
                conParser = ContigParser()
                bamParser = BamParser()
    
                bamParser.dumpCovTable(h5file.root.profile.coverage, stoitColNames)
                conParser.dumpTnTable(h5file.root.profile.kms, kse.kmerCols)
                self.dumpContigs(h5file.root.meta.contigs)
                self.dumpMeta(h5file.root.meta.meta)
        except:
            print "Error opening database:", dbFileName, exc_info()[0]
            raise

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
                
    def parse(self, contigFile, kse):
        """Do the heavy lifting of parsing"""
        print "Parsing contigs"        
        tmp_storage = {} # save everything here first so we can sort accordingly
        for cid,seq,qual in self.readfq(contigFile):
            tmp_storage[cid] = (kse.getKSig(seq.upper()), len(seq))

        # get a list of contig names
        # This array is used like everywhere dude...
        con_names = {}
        for cid in sorted(tmp_storage.keys()):     
            con_names[cid] = tmp_storage[cid][1]
        return (tmp_storage, con_names)
        
    def storeSigs(self, data, table):
        for cid in sorted(data.keys()):     
            # make a new row
            KMER_row = table.row
            # punch in the data
            for mer in data[cid][0].keys():
                KMER_row[mer] = data[cid][0][mer]
            KMER_row.append()
        table.flush()

    def getWantedSeqs(self, contigFile, wanted, storage={}):
        """Do the heavy lifting of parsing"""
        print "Parsing contigs"        
        for cid,seq,qual in self.readfq(contigFile):
            if(cid in wanted):
                storage[cid] = seq
        return storage 

    def dumpTnTable(self, table, kmerCols):
        """Dump the guts of the TN table"""
        print "-----------------------------------"
        print "KmerSig table"
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
    def __init__(self, kLen=4):
        self.kLen = kLen
        self.compl = s_maketrans('ACGT', 'TGCA')
        (self.kmerCols, self.llDict) = self.makeKmerColNames(makeLL=True)
        self.numMers = len(self.kmerCols)
        
    def makeKmerColNames(self, makeLL=False):
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
        ll_dict = {}
        for mer in out_list:
            lmer = self.shiftLowLexi(mer)
            ll_dict[mer] = lmer 
            if lmer not in ret_list:
                ret_list.append(lmer)
        if makeLL:
            return (ret_list, ll_dict)
        else:
            return ret_list

    def getKmerSigWeights(self):
        """Return a hash of index into kmer sig -> GC %"""
        kmer_names = self.makeKmerColNames()
        weights = []
        compl = s_maketrans('ACGTacgt', '01100110')
        for kmer in kmer_names:
            weights.append(sum([float(x) for x in list(kmer.translate(compl))])/float(self.kLen))
        return weights

    def getGC(self, seq):
        """Get the GC of a sequence"""
        Ns = seq.count('N') + seq.count('n')
        compl = s_maketrans('ACGTacgtnN', '0110011000')
        return sum([float(x) for x in list(seq.translate(compl))])/float(len(seq) - Ns)

    def shiftLowLexi(self, seq):
        """Return the lexicographically lowest form of this sequence"""
        rseq = self.revComp(seq)
        if(seq < rseq):
            return seq
        return rseq

    def shiftLowLexiMer(self, seq):
        """Return the lexicographically lowest form of this kmer"""
        try:
            return self.llDict[seq]
        except KeyError:
            return seq
        
    def revComp(self, seq):
        """Return the reverse complement of a sequence"""
        # build a dictionary to know what letter to switch to
        return seq.translate(self.compl)[::-1]
    
    def getKSig(self, seq):
        """Work out kmer signature for a nucleotide sequence"""
        sig = dict(zip(self.kmerCols, [0.0] * self.numMers))
        ll = len(seq)
        for i in range(0,ll-self.kLen+1):
            try:
                this_mer = self.llDict[seq[i:i+self.kLen]]
                try:
                    sig[this_mer] += 1.0
                except KeyError:
                    sig[this_mer] = 1.0
            except KeyError:
                pass

        # normalise by length and return
        return dict(zip(self.kmerCols, [ X / ll for X in sig.values()]))


###############################################################################
###############################################################################
###############################################################################
###############################################################################
class BamParser:
    """Parse multiple bam files and write the output to hdf5 """

    def __init__(self): pass
    
    def parse(self, bamFiles, stoitColNames, covTable, contigNames, cid2Indicies):
        """Parse multiple bam files and store the results in the main DB
        
        table: a table in an open h5 file like "CID,COV_1,...,COV_n,length"
        stoitColNames: names of the COV_x columns
        """
        print "Importing BAM files"
        from bamtyper.utilities import BamParser as BTBP
        BP = BTBP()
        (links, ref_lengths, coverages) = BP.getLinks(bamFiles, full=False, verbose=True, doCoverage=True, minJoin=5)

        # go through all the contigs sorted by name and write to the DB
        try:
            for cid in sorted(contigNames.keys()):
                # make a new row
                cov_row = covTable.row
                # punch in the data
                for i in range(len(stoitColNames)):
                    try:
                        cov = coverages[i][cid]
                    except KeyError:
                        # may be no coverage for this contig
                        cov = 0.0
                    cov_row[stoitColNames[i]] = cov 
                cov_row.append()
            covTable.flush()
        except:
            print "Error saving results to DB"
            raise
        
        # transform the links into something a little easier to parse later
        rowwise_links = []
        for cid in links:
            for link in links[cid]:
                try:
                    rowwise_links.append([cid2Indicies[cid],          # contig 1 
                                          cid2Indicies[link[0]],      # contig 2
                                          int(link[1]),               # numReads
                                          int(link[2]),               # linkType
                                          int(link[3])                # gap
                                          ])
                except KeyError:
                    pass
        return rowwise_links
    
    def dumpCovTable(self, table, stoitColNames):
        """Dump the guts of the coverage table"""
        print "-----------------------------------"
        print "Coverage table"
        print "-----------------------------------"
        for row in table:
            for colName in stoitColNames:
                print ",",row[colName],
            print ""

    def initLinks(self, rowwiseLinks, linksTable):
        # go through all the contigs sorted by name and write to the DB
        try:
            for link in rowwiseLinks:
                # make a new row
                link_row = linksTable.row
                # punch in the data
                link_row['contig1'] = link[0]
                link_row['contig2'] = link[1]
                link_row['numReads'] = link[2]
                link_row['linkType'] = link[3]
                link_row['gap'] = link[4]
                link_row.append()
            linksTable.flush()
        except:
            print "Error saving results to DB"
            raise

def getBamDescriptor(fullPath):
    """AUX: Reduce a full path to just the file name minus extension"""
    return op_splitext(op_basename(fullPath))[0]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

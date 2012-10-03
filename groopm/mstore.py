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
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import sys
import os
import string

import tables
import numpy as np
import pysam

# GroopM imports
import groopmExceptions as ge

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
    'core'   : tables.BoolCol(pos=3)                  # is this contig part of a bin's core?
    
    ** Bins **
    table = 'bin'
    'bid'        : tables.Int32Col(pos=0)
    'numMembers' : tables.Int32Col(pos=1)

    ------------------------
     SOM
    group = '/som'
    ------------------------
    ** Metadata **
    table = 'meta'
    'side'            : tables.Int32Col(pos=1)      # number of rows in the som
    'covDimension'    : tables.Int32Col(pos=2)      # number of stoits
    'merDimension'    : tables.Int32Col(pos=3)      # number of mers being used
    'covWeightIds'    : tables.StringCol(64,pos=4)
    'merWeightIds'    : tables.StringCol(64,pos=5)
    'covRegionIds'    : tables.StringCol(64,pos=6)
    'merRegionIds'    : tables.StringCol(64,pos=7)
    
    ** Weights **                                   # raw weights formed from training
    table = 'covWeights[1,2,3]'                     # 3 SOMS for coverage
    'col1' : tables.FloatCol(pos=1)
    'col2' : tables.FloatCol(pos=2)
    'col3' : tables.FloatCol(pos=3)
    ...
    
    table = 'merWeights[1,2,3]'                     # 3 SOMS for mers too
    'col1' : tables.FloatCol(pos=1)
    'col2' : tables.FloatCol(pos=2)
    'col3' : tables.FloatCol(pos=3)
    ...
    
    ** Regions **                                   # maps points in the map to bin ids
    table = 'covRegions[1,2,3]'
    'col1' : tables.Int32Col(pos=1)
    'col2' : tables.Int32Col(pos=2)
    'col3' : tables.Int32Col(pos=3)
    ...
    
    table = 'merRegions[1,2,3]'
    'col1' : tables.Int32Col(pos=1)
    'col2' : tables.Int32Col(pos=2)
    'col3' : tables.Int32Col(pos=3)
    ...

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
        try:        
            with tables.openFile(dbFileName, mode = "w", title = "GroopM") as h5file:
                # Create groups under "/" (root) for storing profile information and metadata
                profile_group = h5file.createGroup("/", 'profile', 'Assembly profiles')
                meta_group = h5file.createGroup("/", 'meta', 'Associated metadata')
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
                    KMER_table = h5file.createTable(profile_group, 'kms', db_desc, "Kmer signature")
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
                # Add a table for the contigs
                #------------------------
                db_desc = {'cid' : tables.StringCol(512, pos=0),
                           'bid' : tables.Int32Col(dflt=0,pos=1),
                           'length' : tables.Int32Col(pos=2),
                           'core' : tables.BoolCol(dflt=False, pos=3) }
                try:
                    CONTIG_table = h5file.createTable(meta_group, 'contigs', db_desc, "Contig information")
                    self.initContigs(CONTIG_table, contigNames)
                except:
                    print "Error creating CONTIG table:", sys.exc_info()[0]
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
                    print "Error creating BIN table:", sys.exc_info()[0]
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
                    stoitColNames.append(bam_desc)
                    ppos += 1
                try:
                    COV_table = h5file.createTable(profile_group, 'coverage', db_desc, "Bam based coverage")
                    bamParser.parse(bamFiles, stoitColNames, COV_table, contigNames)
                except:
                    print "Error creating coverage table:", sys.exc_info()[0]
                    raise
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
                    META_table = h5file.createTable(meta_group, 'meta', db_desc, "Descriptive data")
                    self.initMeta(META_table, str.join(',',stoitColNames), len(stoitColNames), str.join(',',kse.kmerCols), kmerSize, len(kse.kmerCols), len(contigNames))
                except:
                    print "Error creating META table:", sys.exc_info()[0]
                    raise
        except:
            print "Error creating database:", dbFileName, sys.exc_info()[0]
            raise
        
        print "****************************************************************"
        print "Data loaded successfully!"
        print " ->",len(contigNames),"contigs"
        print " ->",len(stoitColNames),"BAM files"
        print "Written to: '"+dbFileName+"'"
        print "****************************************************************"

        if(dumpAll):
            self.dumpAll(dbFileName)
            
        # all good!
        return True
                
    def initContigs(self, table, contigNames):
        """Initialise the contigs table
        
        set to 0 for no bin assignment
        """
        for cid in sorted(contigNames):
            CONTIG_row = table.row
            CONTIG_row['cid'] = cid
            CONTIG_row['length'] = contigNames[cid]
            CONTIG_row.append()
        table.flush()
    
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
# DB CREATION / INITIALISATION  - SOMS

    def updateSOMTables(self,
                        dbFileName,
                        side,
                        covDim,
                        merDim,
                        covWeights={},          # use these to do updates
                        merWeights={},
                        covRegions={},
                        merRegions={}
                        ):
        try:
            # first, make sure we've got a som group`
            ids_in_use = self.getSOMDataInfo(dbFileName)
            
            # now we need to fix the ids_in_use structure to suit our updates
            for i in covWeights.keys():
                if i not in ids_in_use["weights"]["cov"]:
                     ids_in_use["weights"]["cov"].append(i)
            for i in merWeights.keys():
                if i not in ids_in_use["weights"]["mer"]:
                     ids_in_use["weights"]["mer"].append(i)
            for i in covRegions.keys():
                if i not in ids_in_use["regions"]["cov"]:
                     ids_in_use["regions"]["cov"].append(i)
            for i in merRegions.keys():
                if i not in ids_in_use["regions"]["mer"]:
                     ids_in_use["regions"]["mer"].append(i)
            
            with tables.openFile(dbFileName, mode='a', rootUEP="/som") as som_group:
                #------------------------
                # Metadata
                #------------------------
                db_desc = {'side'            : tables.Int32Col(pos=1),
                           'covDimension'    : tables.Int32Col(pos=2),
                           'merDimension'    : tables.Int32Col(pos=3),
                           'covWeightIds'    : tables.StringCol(64,pos=4),
                           'merWeightIds'    : tables.StringCol(64,pos=5),
                           'covRegionIds'    : tables.StringCol(64,pos=6),
                           'merRegionIds'    : tables.StringCol(64,pos=7)
                           }

                # clear previous metadata
                # try remove any older failed attempts
                try:
                    som_group.removeNode('/', 'tmp_meta')
                except:
                    pass

                # make a new tmp table
                try:
                    META_table = som_group.createTable('/', 'tmp_meta', db_desc, "SOM metadata")
                    self.initSOMMeta(META_table, side, covDim, merDim, ids_in_use)
                except:
                    print "Error creating META table:", sys.exc_info()[0]
                    raise
                # rename as the meta table
                som_group.renameNode('/', 'meta', 'tmp_meta', overwrite=True)
                
                #------------------------
                # Weights matricies
                #------------------------
                for i in covWeights.keys():
                    self.updateSOMData(som_group, "covWeights"+str(i), covDim, "COVERAGE SOM weights", covWeights[i], type="float")
                for i in merWeights.keys():
                    self.updateSOMData(som_group, "merWeights"+str(i), merDim, "KMER SOM weights", merWeights[i], type="float")
                    
                #------------------------
                # Regions
                #------------------------
                for i in covRegions.keys():
                    self.updateSOMData(som_group, "covRegions"+str(i), 1, "COVERAGE SOM regions", covRegions[i], type="int")
                for i in merRegions.keys():
                    self.updateSOMData(som_group, "merRegions"+str(i), 1, "KMER SOM regions", merRegions[i], type="int")
                    
        except:
            print "Error creating SOM database:", dbFileName, sys.exc_info()[0]
            raise
        
    def initSOMMeta(self, table, side, covDim, merDim, idsInUse):
        """Initialise the meta-data table"""
        META_row = table.row
        META_row['side'] = side
        META_row['covDimension'] = covDim
        META_row['merDimension'] = merDim
        
        if(len(idsInUse["weights"]["cov"]) != 0):
            META_row['covWeightIds'] = ",".join([str(x) for x in idsInUse["weights"]["cov"]])
        else:
            META_row['covWeightIds'] = ""
        
        if(len(idsInUse["weights"]["mer"]) != 0):
            META_row['merWeightIds'] = ",".join([str(x) for x in idsInUse["weights"]["mer"]])
        else:
            META_row['merWeightIds'] = ""
        
        if(len(idsInUse["regions"]["cov"]) != 0):
            META_row['covRegionIds'] = ",".join([str(x) for x in idsInUse["regions"]["cov"]])
        else:
            META_row['covRegionIds'] = ""
        
        if(len(idsInUse["regions"]["mer"]) != 0):
            META_row['merRegionIds'] = ",".join([str(x) for x in idsInUse["regions"]["mer"]])
        else:
            META_row['merRegionIds'] = ""
            
        META_row.append()
        table.flush()

    def updateSOMData(self, group, tableName, dimension, description, data, type):
        """Make a table and update data"""
        db_desc = {}                # make the db desc dynamically
        for j in range(dimension):
            if(type == "float"):
                db_desc["dim"+str(j)] = tables.FloatCol(pos=j)
            else:
                db_desc["dim"+str(j)] = tables.Int32Col(pos=j)
        try:
            group.removeNode('/', 'tmp_'+tableName)
        except:
            pass
        try:
            NEW_table = group.createTable('/', 'tmp_'+tableName, db_desc, description)
            self.injectSOMData(NEW_table, data, dimension)
        except:
            print "Error creating table:", sys.exc_info()[0]
            raise
        group.renameNode('/', tableName, 'tmp_'+tableName, overwrite=True)

    def injectSOMData(self, table, data, dimension):
        """Load SOM data into the database"""
        for row in data:
            for val in row:
                DATA_row = table.row
                for i in range(dimension):
                    DATA_row["dim"+str(i)] = val[i]
                DATA_row.append()
        table.flush()

#------------------------------------------------------------------------------
# GET / SET DATA TABLES - PROFILES 

    def getSOMDataInfo(self, dbFileName):
        """Return the numbers and types of soms stored in the database"""
        ids_in_use = { "weights" : { "mer":[], "cov":[] }, "regions": { "mer" : [], "cov": [] } }
        # open the database
        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                # open the som table                  
                try:
                    som_group = h5file.createGroup("/", 'som', 'SOM data')
                except:
                    # we already have a group, so we should assume that there is 
                    # a meta table associated with this group. Read that to get the
                    # Ids currently in use
                    try:
                        tmp_table = h5file.root.som._f_getChild('meta')
                        meta_row = h5file.root.som.meta.read(start=0, stop=1, step=1)[0]
                        if(len(meta_row[3]) > 0):
                            ids_in_use["weights"]["cov"] = [int(x) for x in meta_row[3].split(",")]
                        if(len(meta_row[4]) > 0):
                            ids_in_use["weights"]["mer"] = [int(x) for x in meta_row[4].split(",")]
                        if(len(meta_row[5]) > 0):
                            ids_in_use["regions"]["cov"] = [int(x) for x in meta_row[5].split(",")]
                        if(len(meta_row[6]) > 0):
                            ids_in_use["regions"]["mer"] = [int(x) for x in meta_row[6].split(",")]
                    except: pass
        except:
            print "Error creating SOM database:", dbFileName, sys.exc_info()[0]
            raise
        return ids_in_use

    def getSOMMetaFields(self, dbFileName):
        """return the value of fieldName in the SOM metadata tables"""
        ret_hash = {'side': -1,
                    'covDimension' : -1,
                    'merDimension' : -1
                    }
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                # theres only one value
                meta_row = h5file.root.som.meta.read(start=0, stop=1, step=1)[0]
                ret_hash['side'] = meta_row[0]
                ret_hash['covDimension'] = meta_row[1]
                ret_hash['merDimension'] = meta_row[2]
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise
        return ret_hash
        
    def getSOMData(self, dbFileName, index, type="weights", flavour="mer"):
        """Load SOM data from the database"""
        # make sure that this table exists
        ids_in_use = self.getSOMDataInfo(dbFileName)
        if index not in ids_in_use[type][flavour]:
            raise ge.SOMDataNotFoundException("No such fish:",flavour,type,index)

        # work out the dimensions of the data
        meta_fields = self.getSOMMetaFields(dbFileName)
        dimension = meta_fields['merDimension']
        if(flavour == "cov"):
             dimension = meta_fields['covDimension']
        
        # now get the data
        table_name = flavour+type.title()+str(index) 
        with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
            table = h5file.root.som._f_getChild(table_name)
            if(type=="weights"):
                data = np.reshape(np.array([list(x) for x in table.read()]), (meta_fields['side'],meta_fields['side'],dimension))
            elif(type == "regions"):
                data = np.reshape(np.array([list(x) for x in table.read()]), (meta_fields['side'],meta_fields['side'],1))
            else:
                raise ge.SOMTypeException("Unknown SOM type: "+type)
            return data
        return np.array([])
    
#------------------------------------------------------------------------------
# GET / SET DATA TABLES - PROFILES 

    def getConditionalIndicies(self, dbFileName, condition=''):
        """return the indicies into the db which meet the condition"""
        if('' == condition):
            condition = "cid != ''" # no condition breaks everything!
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return np.array([x.nrow for x in h5file.root.meta.contigs.where(condition)])
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
                    return np.array([list(h5file.root.profile.coverage[x.nrow]) for x in h5file.root.meta.contigs.where(condition)])
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
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
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
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
                           'length' : tables.Int32Col(pos=2),
                           'core' : tables.BoolCol(dflt=False, pos=3) }
                CONTIG_table = meta_group.createTable('/', 'tmp_contigs', db_desc, "Contig information", expectedrows=len(contig_names))
                self.initContigs(CONTIG_table, contig_names)
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
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
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
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
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
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise
        return {}
        
                
    def getBins(self, dbFileName, condition='', indicies=np.array([])):
        """Load per-contig bins"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([h5file.root.meta.contigs[x][1] for x in indicies]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[1] for x in h5file.root.meta.contigs.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def setBins(self, dbFileName, updates):
        """Set per-contig bins
        
        updates is a dictionary which looks like:
        { tableRow : binValue }
        """
        row_nums = updates.keys()
        try:
            with tables.openFile(dbFileName, mode='a') as h5file:
                table = h5file.root.meta.contigs
                for row_num in updates.keys():
                    new_row = np.zeros((1,),dtype=('S512,i4,i4,b1'))
                    new_row[:] = [(table[row_num][0],updates[row_num],table[row_num][2],table[row_num][3])]
                    table.modifyRows(start=row_num, rows=new_row)
                table.flush()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def getCores(self, dbFileName, condition='', indicies=np.array([])):
        """Load bin core info"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([h5file.root.meta.contigs[x][3] for x in indicies]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[3] for x in h5file.root.meta.contigs.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def setCores(self, dbFileName, updates):
        """Set bin cores
        
        updates is a dictionary which looks like:
        { tableRow : coreValue }
        """
        row_nums = updates.keys()
        try:
            with tables.openFile(dbFileName, mode='a') as h5file:
                table = h5file.root.meta.contigs
                for row_num in updates.keys():
                    new_row = np.zeros((1,),dtype=('S512,i4,i4,b1'))
                    new_row[:] = [(table[row_num][0],table[row_num][1],table[row_num][2],updates[row_num])]
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
                    return np.array([h5file.root.meta.contigs[x][0] for x in indicies]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[0] for x in h5file.root.meta.contigs.readWhere(condition)]).ravel()
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def getContigLengths(self, dbFileName, condition='', indicies=np.array([])):
        """Load contig lengths"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indicies) != 0):
                    return np.array([h5file.root.meta.contigs[x][2] for x in indicies]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[2] for x in h5file.root.meta.contigs.readWhere(condition)]).ravel()
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
                    return np.array([list(h5file.root.profile.kms[x.nrow]) for x in h5file.root.meta.contigs.where(condition)])
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
            raise

    def getMetaField(self, dbFileName, fieldName):
        """return the value of fieldName in the metadata tables"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                # theres only one value
                return h5file.root.meta.meta.read()[fieldName][0]
        except:
            print "Error opening DB:",dbFileName, sys.exc_info()[0]
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
            print "Error opening database:", dbFileName, sys.exc_info()[0]
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
        """Has this data set been *completely* clustered?"""
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

    def dumpContigs(self, table):
        """Raw dump of contig information"""
        print "-----------------------------------"
        print "Contigs table"
        print "-----------------------------------"
        for row in table:
            print row['cid'],",",row['length'],",",row['bid'],",",row['core']

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
            print "Error opening database:", dbFileName, sys.exc_info()[0]
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
    def __init__(self, kLen=4):
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

    def getKmerSigWeights(self):
        """Return a hash of index into kmer sig -> GC %"""
        kmer_names = self.makeKmerColNames()
        weights = []
        compl = string.maketrans('ACGTacgt', '01100110')
        for kmer in kmer_names:
            weights.append(sum([float(x) for x in list(kmer.translate(compl))])/float(self.kLen))
        return weights

    def getGC(self, seq):
        """Get the GC of a sequence"""
        Ns = seq.count('N') + seq.count('n')
        compl = string.maketrans('ACGTacgtnN', '0110011000')
        return sum([float(x) for x in list(seq.translate(compl))])/float(len(seq) - Ns)

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
    
    def parse(self, bamFiles, stoitColNames, table, contigNames):
        """Parse multiple bam files and store the results in the main DB
        
        table: a table in an open h5 file like "CID,COV_1,...,COV_n,length"
        stoitColNames: names of the COV_x columns
        """
        # parse the BAMs
        # we need to have some type of entry for each contig
        # so start by putting 0's here
        tmp_storage = {}
        num_bams = len(stoitColNames)
        for cid in contigNames.keys():
            tmp_storage[cid] = np.zeros((num_bams))

        bam_count = 0
        for bf in bamFiles:
            bam_file = None
            try:
                bam_file = pysam.Samfile(bf, 'rb')
                print "Parsing",stoitColNames[bam_count],"(",(bam_count+1),"of",num_bams,")"
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
                for i in range(0,len(stoitColNames)):
                    cov_row[stoitColNames[i]] = tmp_storage[cid][i]
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
        
    def dumpCovTable(self, table, stoitColNames):
        """Dump the guts of the coverage table"""
        print "-----------------------------------"
        print "Coverage table"
        print "-----------------------------------"
        for row in table:
            for colName in stoitColNames:
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

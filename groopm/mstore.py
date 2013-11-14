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
__copyright__ = "Copyright 2012/2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.3.6.2"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Release"
__current_GMDB_version__ = 5

###############################################################################

from sys import exc_info
from os.path import splitext as op_splitext, basename as op_basename
from string import maketrans as s_maketrans

import tables
import numpy as np
from scipy.spatial.distance import cdist, squareform
import pysam

# GroopM imports
from PCA import PCA, Center

np.seterr(all='raise')

# shut up pytables!
import warnings
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

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
    'mer1' : tables.FloatCol(pos=0)
    'mer2' : tables.FloatCol(pos=1)
    'mer3' : tables.FloatCol(pos=2)
    ...

    **Kmer Vals**
    table = 'kpca'
    'pc1' : tables.FloatCol(pos=0)
    'pc2' : tables.FloatCol(pos=1)
    'pc3' : tables.FloatCol(pos=2)
    ...

    **Coverage profile**
    table = 'coverage'
    'stoit1' : tables.FloatCol(pos=0)
    'stoit2' : tables.FloatCol(pos=1)
    'stoit3' : tables.FloatCol(pos=2)
    ...

    **Transformed coverage profile**
    table = 'transCoverage'
    'x' : tables.FloatCol(pos=0)
    'y' : tables.FloatCol(pos=1)
    'z' : tables.FloatCol(pos=2)

    **Normalised coverage profile**
    table = 'normCoverage'
    'normCov' : tables.FloatCol(pos=0)

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
    'formatVersion' : tables.Int32Col(pos=9)          # groopm file version

    **PC variance**
    table = 'kpca_variance'
    'pc1_var' : tables.FloatCol(pos=0)
    'pc2_var' : tables.FloatCol(pos=1)
    'pc3_var' : tables.FloatCol(pos=2)
    ...

    ** Contigs **
    table = 'contigs'
    'cid'    : tables.StringCol(512, pos=0)
    'bid'    : tables.Int32Col(pos=1)
    'length' : tables.Int32Col(pos=2)
    'gc'     : tables.FloatCol(pos=3)

    ** Bins **
    table = 'bins'
    'bid'        : tables.Int32Col(pos=0)
    'numMembers' : tables.Int32Col(pos=1)
    'isLikelyChimeric' : tables.BoolCol(pos=2)

    **Transformed coverage corners**
    table = 'transCoverageCorners'
    'x' : tables.FloatCol(pos=0)
    'y' : tables.FloatCol(pos=1)
    'z' : tables.FloatCol(pos=2)

    """
    def __init__(self): pass

#------------------------------------------------------------------------------
# DB CREATION / INITIALISATION  - PROFILES

    def createDB(self, bamFiles, contigs, dbFileName, timer, kmerSize=4, force=False):
        """Main wrapper for parsing all input files"""
        # load all the passed vars
        dbFileName = dbFileName
        contigsFile = contigs
        stoitColNames = []

        kse = KmerSigEngine(kmerSize)
        conParser = ContigParser()
        bamParser = BamParser()

        cid_2_indices = {}

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
                links_group = h5file.createGroup("/", 'links', 'Paired read link information')
                #------------------------
                # parse contigs
                #
                # Contig IDs are key. Any keys existing in other files but not in this file will be
                # ignored. Any missing keys in other files will be given the default profile value
                # (typically 0). Ironically, we don't store the CIDs here, these are saved one time
                # only in the bin table
                #
                # Before writing to the database we need to make sure that none of them have
                # 0 coverage @ all stoits.
                #------------------------
                try:
                    with open(contigsFile, "r") as f:
                        try:
                            (con_names, con_gcs, con_lengths, con_ksigs) = conParser.parse(f, kse)
                            num_cons = len(con_names)
                            cid_2_indices = dict(zip(con_names, range(num_cons)))
                        except:
                            print "Error parsing contigs"
                            raise
                except:
                    print "Could not parse contig file:",contigsFile,exc_info()[0]
                    raise

                #------------------------
                # parse bam files
                #------------------------
                (rowwise_links, cov_profiles) = bamParser.parse(bamFiles,
                                                                con_names,
                                                                cid_2_indices)

                tot_cov = [sum(profile) for profile in cov_profiles]
                good_indices = np.nonzero(tot_cov)[0]
                bad_indices = [i for i in range(num_cons) if i not in good_indices]

                if len(bad_indices) > 0:
                    # report the bad contigs to the user
                    # and strip them before writing to the DB
                    print "****************************************************************"
                    print " IMPORTANT! - there are %d contigs with 0 coverage" % len(bad_indices)
                    print " across all stoits. They will be ignored:"
                    print "****************************************************************"
                    for i in xrange(0, min(5, len(bad_indices))):
                        print con_names[bad_indices[i]]
                    if len(bad_indices) > 5:
                      print '(+ %d additional contigs)' % (len(bad_indices)-5)
                    print "****************************************************************"

                    con_names = con_names[good_indices]
                    con_lengths = con_lengths[good_indices]
                    con_gcs = con_gcs[good_indices]
                    cov_profiles = cov_profiles[good_indices]

                # these will need to be tupalized regardless...
                con_ksigs = [tuple(i) for i in con_ksigs[good_indices]]
                num_cons = len(con_names)

                #------------------------
                # calculate PCAs and write kmer sigs
                #------------------------
                # store the raw calculated kmer sigs in one table
                db_desc = []
                for mer in kse.kmerCols:
                     db_desc.append((mer, float))
                try:
                    h5file.createTable(profile_group,
                                       'kms',
                                       np.array(con_ksigs, dtype=db_desc),
                                       title='Kmer signatures',
                                       expectedrows=num_cons
                                       )
                except:
                    print "Error creating KMERSIG table:", exc_info()[0]
                    raise

                # compute the PCA of the ksigs and store these too
                pc_ksigs, sumvariance = conParser.PCAKSigs(con_ksigs)

                db_desc = []
                for i in xrange(0, len(pc_ksigs[0])):
                  db_desc.append(('pc' + str(i+1), float))

                try:
                    h5file.createTable(profile_group,
                                       'kpca',
                                       np.array(pc_ksigs, dtype=db_desc),
                                       title='Kmer signature PCAs',
                                       expectedrows=num_cons
                                       )
                except:
                    print "Error creating KMERVALS table:", exc_info()[0]
                    raise

                #------------------------
                # write cov profiles
                #------------------------
                # build a table template based on the number of bamfiles we have
                for i, bf in enumerate(bamFiles):
                    # assume the file is called something like "fred.bam"
                    # we want to rip off the ".bam" part
                    bam_desc = getBamDescriptor(bf, i + 1)
                    stoitColNames.append(bam_desc)

                stoitColNames = np.array(stoitColNames)

                # determine normalised cov profiles and shuffle the BAMs
                norm_coverages = np.array([np.linalg.norm(cov_profiles[i]) for i in range(num_cons)])
                CT = CoverageTransformer(num_cons,
                                         len(stoitColNames),
                                         norm_coverages,
                                         np.array(pc_ksigs)[:,0],
                                         cov_profiles,
                                         stoitColNames)

                CT.transformCP()
                # these will need to be tupalized regardless...
                cov_profiles = [tuple(i) for i in cov_profiles]

                CT.transformedCP = [tuple(i) for i in CT.transformedCP]
                CT.corners = [tuple(i) for i in CT.corners]
                # now CT stores the transformed coverages and other important information
                # the ordering of stoitColNames and cov_profiles should be fixed
                # so we will write this to the database without further modification

                # raw coverages
                db_desc = []
                for scn in CT.stoitColNames:
                    db_desc.append((scn, float))

                try:
                    h5file.createTable(profile_group,
                                       'coverage',
                                       np.array(cov_profiles, dtype=db_desc),
                                       title="Bam based coverage",
                                       expectedrows=num_cons)
                except:
                    print "Error creating coverage table:", exc_info()[0]
                    raise

                # transformed coverages
                db_desc = [('x', float),
                           ('y', float),
                           ('z', float)]
                try:
                    h5file.createTable(profile_group,
                                       'transCoverage',
                                       np.array(CT.transformedCP , dtype=db_desc),
                                       title="Transformed coverage",
                                       expectedrows=num_cons)
                except:
                    print "Error creating transformed coverage table:", exc_info()[0]
                    raise

                # transformed coverage corners
                db_desc = [('x', float),
                           ('y', float),
                           ('z', float)]
                try:
                    h5file.createTable(meta_group,
                                       'transCoverageCorners',
                                       np.array(CT.corners , dtype=db_desc),
                                       title="Transformed coverage corners",
                                       expectedrows=len(stoitColNames))
                except:
                    print "Error creating transformed coverage corner table:", exc_info()[0]
                    raise

                # normalised coverages
                db_desc = [('normCov', float)]
                try:
                    h5file.createTable(profile_group,
                                       'normCoverage',
                                       np.array(CT.normCoverages , dtype=db_desc),
                                       title="Normalised coverage",
                                       expectedrows=num_cons)
                except:
                    print "Error creating normalised coverage table:", exc_info()[0]
                    raise

                #------------------------
                # Add a table for the contigs
                #------------------------
                self.setBinAssignments((h5file, meta_group),
                                       image=zip(con_names,
                                                 [0]*num_cons,
                                                 con_lengths, con_gcs)
                                       )

                #------------------------
                # Add a table for the bins
                #------------------------
                self.setBinStats(dbFileName, [], firstWrite=True)

                print "    %s" % timer.getTimeStamp()

                #------------------------
                # contig links
                #------------------------
                # set table size according to the number of links returned from
                # the previous call
                db_desc = [('contig1', int),
                           ('contig2', int),
                           ('numReads', int),
                           ('linkType', int),
                           ('gap', int)]
                try:
                    h5file.createTable(links_group,
                                       'links',
                                       np.array(rowwise_links, dtype=db_desc),
                                       title="ContigLinks",
                                       expectedrows=len(rowwise_links))
                except:
                    print "Error creating links table:", exc_info()[0]
                    raise
                print "    %s" % timer.getTimeStamp()

                #------------------------
                # Add metadata
                #------------------------
                meta_data = (str.join(',',stoitColNames),
                             len(stoitColNames),
                             str.join(',',kse.kmerCols),
                             kmerSize,
                             len(kse.kmerCols),
                             num_cons,
                             0,
                             False,
                             False,
                             __current_GMDB_version__)
                self.setMeta(h5file, meta_data)

                # kmer signature variance table
                pc_var = [sumvariance[0]]
                for i in xrange(1, len(sumvariance)):
                  pc_var.append(sumvariance[i]-sumvariance[i-1])
                pc_var = tuple(pc_var)

                db_desc = []
                for i in xrange(0, len(pc_var)):
                    db_desc.append(('pc' + str(i+1) + '_var', float))

                try:
                    h5file.createTable(meta_group,
                                            'kpca_variance',
                                            np.array([pc_var], dtype=db_desc),
                                            title='Variance of kmer signature PCAs',
                                            expectedrows=1
                                            )
                except:
                    print "Error creating tmp_kpca_variance table:", exc_info()[0]
                    raise

        except:
            print "Error creating database:", dbFileName, exc_info()[0]
            raise

        print "****************************************************************"
        print "Data loaded successfully!"
        print " ->",num_cons,"contigs"
        print " ->",len(stoitColNames),"BAM files"
        print "Written to: '"+dbFileName+"'"
        print "****************************************************************"
        print "    %s" % timer.getTimeStamp()

        # all good!
        return True

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
# DB UPGRADE

    def checkAndUpgradeDB(self, dbFileName, silent=False):
        """Check the DB and upgrade if necessary"""
        # get the DB format version
        this_DB_version = self.getGMDBFormat(dbFileName)
        if __current_GMDB_version__ == this_DB_version:
            if not silent:
                print "    GroopM DB version (%s) up to date" % this_DB_version
            return

        # now, if we get here then we need to do some work
        upgrade_tasks = {}
        upgrade_tasks[(0,1)] = self.upgradeDB_0_to_1
        upgrade_tasks[(1,2)] = self.upgradeDB_1_to_2
        upgrade_tasks[(2,3)] = self.upgradeDB_2_to_3
        upgrade_tasks[(3,4)] = self.upgradeDB_3_to_4
        upgrade_tasks[(4,5)] = self.upgradeDB_4_to_5

        # we need to apply upgrades in order!
        # keep applying the upgrades as long as we need to
        while this_DB_version < __current_GMDB_version__:
            task = (this_DB_version, this_DB_version+1)
            upgrade_tasks[task](dbFileName)
            this_DB_version += 1

    def upgradeDB_0_to_1(self, dbFileName):
        """Upgrade a GM db from version 0 to version 1"""
        print "*******************************************************************************\n"
        print "              *** Upgrading GM DB from version 0 to version 1 ***"
        print ""
        print "                            please be patient..."
        print ""
        # the change in this version is that we'll be saving the first
        # two kmerSig PCA's in a separate table
        print "    Calculating and storing the kmerSig PCAs"

        # compute the PCA of the ksigs
        ksigs = self.getKmerSigs(dbFileName)
        CP = ContigParser()
        pc_ksigs, sumvariance = CP.PCAKSigs(ksigs)
        num_cons = len(pc_ksigs)

        db_desc = [('pc1', float),
                   ('pc2', float)]
        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/profile") as profile_group:
                try:
                    profile_group.createTable('/',
                                              'kpca',
                                              np.array(pc_ksigs, dtype=db_desc),
                                              title='Kmer signature PCAs',
                                              expectedrows=num_cons
                                              )
                except:
                    print "Error creating KMERVALS table:", exc_info()[0]
                    raise
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

        # update the formatVersion field and we're done
        self.setGMDBFormat(dbFileName, 1)
        print "*******************************************************************************"

    def upgradeDB_1_to_2(self, dbFileName):
        """Upgrade a GM db from version 1 to version 2"""
        print "*******************************************************************************\n"
        print "              *** Upgrading GM DB from version 1 to version 2 ***"
        print ""
        print "                            please be patient..."
        print ""
        # the change in this version is that we'll be saving a variable number of kmerSig PCA's
        # and GC information for each contig
        print "    Calculating and storing the kmer signature PCAs"

        # grab any data needed from database before opening if for modification
        bin_ids = self.getBins(dbFileName)
        orig_con_names = self.getContigNames(dbFileName)
        ksigs = self.getKmerSigs(dbFileName)

        # compute the PCA of the ksigs
        conParser = ContigParser()
        pc_ksigs, sumvariance = conParser.PCAKSigs(ksigs)
        num_cons = len(pc_ksigs)

        db_desc = []
        for i in xrange(0, len(pc_ksigs[0])):
          db_desc.append(('pc' + str(i+1), float))

        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                pg = h5file.getNode('/', name='profile')
                try:
                    try:
                        h5file.removeNode(pg, 'tmp_kpca')
                    except:
                        pass

                    h5file.createTable(pg,
                                       'tmp_kpca',
                                       np.array(pc_ksigs, dtype=db_desc),
                                       title='Kmer signature PCAs',
                                       expectedrows=num_cons
                                      )

                    h5file.renameNode(pg, 'kpca', 'tmp_kpca', overwrite=True)

                except:
                    print "Error creating kpca table:", exc_info()[0]
                    raise

                # Add GC
                contigFile = raw_input('\nPlease specify fasta file containing the bam reference sequences: ')
                with open(contigFile, "r") as f:
                    try:
                        contigInfo = {}
                        for cid,seq,qual in conParser.readfq(f):
                            contigInfo[cid] = (len(seq), conParser.calculateGC(seq))

                        # sort the contig names here once!
                        con_names = np.array(sorted(contigInfo.keys()))

                        # keep everything in order...
                        con_gcs = np.array([contigInfo[cid][1] for cid in con_names])
                        con_lengths = np.array([contigInfo[cid][0] for cid in con_names])
                    except:
                        print "Error parsing contigs"
                        raise

                # remove any contigs not in the current DB (these were removed due to having zero coverage)
                good_indices = [i for i in range(len(orig_con_names)) if orig_con_names[i] in con_names]

                con_names = con_names[good_indices]
                con_lengths = con_lengths[good_indices]
                con_gcs = con_gcs[good_indices]
                bin_ids = bin_ids[good_indices]

                mg = h5file.getNode('/', name='meta')
                self.setBinAssignments((h5file, mg),
                               image=zip(con_names,
                                         bin_ids,
                                         con_lengths, con_gcs)
                               )
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

        # update the formatVersion field and we're done
        self.setGMDBFormat(dbFileName, 2)
        print "*******************************************************************************"

    def upgradeDB_2_to_3(self, dbFileName):
        """Upgrade a GM db from version 2 to version 3"""
        print "*******************************************************************************\n"
        print "              *** Upgrading GM DB from version 2 to version 3 ***"
        print ""
        print "                            please be patient..."
        print ""
        # the change in this version is that we'll be saving the variance for each kmerSig PCA
        print "    Calculating and storing variance of kmer signature PCAs"

        # compute the PCA of the ksigs
        conParser = ContigParser()
        ksigs = self.getKmerSigs(dbFileName)
        pc_ksigs, sumvariance = conParser.PCAKSigs(ksigs)

        # calcualte variance of each PC
        pc_var = [sumvariance[0]]
        for i in xrange(1, len(sumvariance)):
          pc_var.append(sumvariance[i]-sumvariance[i-1])
        pc_var = tuple(pc_var)

        db_desc = []
        for i in xrange(0, len(pc_var)):
          db_desc.append(('pc' + str(i+1) + '_var', float))

        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                meta = h5file.getNode('/', name='meta')
                try:
                    try:
                        h5file.removeNode(meta, 'tmp_kpca_variance')
                    except:
                        pass

                    h5file.createTable(meta,
                                              'tmp_kpca_variance',
                                              np.array([pc_var], dtype=db_desc),
                                              title='Variance of kmer signature PCAs',
                                              expectedrows=1
                                              )

                    h5file.renameNode(meta, 'kpca_variance', 'tmp_kpca_variance', overwrite=True)

                except:
                    print "Error creating kpca_variance table:", exc_info()[0]
                    raise
        except:
            print "Error opening DB:", dbFileName, exc_info()[0]
            raise

        # update the formatVersion field and we're done
        self.setGMDBFormat(dbFileName, 3)
        print "*******************************************************************************"

    def upgradeDB_3_to_4(self, dbFileName):
        """Upgrade a GM db from version 3 to version 4"""
        print "*******************************************************************************\n"
        print "              *** Upgrading GM DB from version 3 to version 4 ***"
        print ""
        print "                            please be patient..."
        print ""
        # the change in this version is that we'll be saving the variance for each kmerSig PCA
        print "    Adding chimeric flag for each bin."
        print "    !!! Groopm core must be run again for this flag to be properly set. !!!"

        # read existing data in 'bins' table
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                ret_dict = {}
                all_rows = h5file.root.meta.bins.read()
                for row in all_rows:
                    ret_dict[row[0]] = row[1]
        except:
            print "Error opening DB:", dbFileName, exc_info()[0]
            raise

        # write new table with chimeric flag set to False by default
        db_desc = [('bid', int),
                   ('numMembers', int),
                   ('isLikelyChimeric', bool)]

        data = []
        for bid in ret_dict:
          data.append((bid, ret_dict[bid], False))

        bd = np.array(data, dtype=db_desc)

        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                mg = h5file.getNode('/', name='meta')

                try:
                    h5file.removeNode(mg, 'tmp_bins')
                except:
                    pass

                try:
                    h5file.createTable(mg,
                                       'tmp_bins',
                                       bd,
                                       title="Bin information",
                                       expectedrows=1)
                except:
                    print "Error creating META table:", exc_info()[0]
                    raise

                h5file.renameNode(mg, 'bins', 'tmp_bins', overwrite=True)
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

        # update the formatVersion field and we're done
        self.setGMDBFormat(dbFileName, 4)
        print "*******************************************************************************"

    def upgradeDB_4_to_5(self, dbFileName):
        """Upgrade a GM db from version 4 to version 5"""
        print "*******************************************************************************\n"
        print "              *** Upgrading GM DB from version 4 to version 5 ***"
        print ""
        print "                            please be patient..."
        print ""
        # the change in this version is that we'll be saving the transformed coverage coords
        print "    Saving transformed coverage profiles"
        print "    You will not need to re-run parse or core due to this change"

        # we need to get the raw coverage profiles and the kmerPCA1 data
        indices = self.getConditionalIndices(dbFileName, silent=False, checkUpgrade=False)
        kPCA_1 = self.getKmerPCAs(dbFileName, indices=indices)[:,0]
        raw_coverages = self.getCoverageProfiles(dbFileName, indices=indices)
        norm_coverages = np.array([np.linalg.norm(raw_coverages[i]) for i in range(len(indices))])

        CT = CoverageTransformer(len(indices),
                                 self.getNumStoits(dbFileName),
                                 norm_coverages,
                                 kPCA_1,
                                 raw_coverages,
                                 np.array(self.getStoitColNames(dbFileName).split(",")))

        CT.transformCP()
        CT.transformedCP = [tuple(i) for i in CT.transformedCP]
        CT.covProfiles = [tuple(i) for i in CT.covProfiles]
        CT.corners = [tuple(i) for i in CT.corners]

        # now CT stores the transformed coverages and other important information
        # we will write this to the database
        with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
            meta_group = h5file.getNode('/', name='meta')
            profile_group = h5file.getNode('/', name='profile')

            # raw coverages - we may have reordered rows, so we should fix this now!
            db_desc = []
            for scn in CT.stoitColNames:
                db_desc.append((scn, float))

            try:
                h5file.removeNode(mg, 'tmp_coverages')
            except:
                pass

            try:
                h5file.createTable(profile_group,
                                   'tmp_coverages',
                                   np.array(CT.covProfiles, dtype=db_desc),
                                   title="Bam based coverage",
                                   expectedrows=CT.numContigs)
            except:
                print "Error creating coverage table:", exc_info()[0]
                raise

            h5file.renameNode(profile_group, 'coverage', 'tmp_coverages', overwrite=True)

            # transformed coverages
            db_desc = [('x', float),
                       ('y', float),
                       ('z', float)]
            try:
                h5file.createTable(profile_group,
                                   'transCoverage',
                                   np.array(CT.transformedCP , dtype=db_desc),
                                   title="Transformed coverage",
                                   expectedrows=CT.numContigs)
            except:
                print "Error creating transformed coverage table:", exc_info()[0]
                raise

            # transformed coverage corners
            db_desc = [('x', float),
                       ('y', float),
                       ('z', float)]
            try:
                h5file.createTable(meta_group,
                                   'transCoverageCorners',
                                   np.array(CT.corners , dtype=db_desc),
                                   title="Transformed coverage corners",
                                   expectedrows=CT.numStoits)
            except:
                print "Error creating transformed coverage corner table:", exc_info()[0]
                raise


            # normalised coverages
            db_desc = [('normCov', float)]
            try:
                h5file.createTable(profile_group,
                                   'normCoverage',
                                   np.array(CT.normCoverages , dtype=db_desc),
                                   title="Normalised coverage",
                                   expectedrows=CT.numContigs)
            except:
                print "Error creating normalised coverage table:", exc_info()[0]
                raise

        # stoit col names may have been shuffled
        meta_data = (",".join([str(i) for i in CT.stoitColNames]),
                    CT.numStoits,
                    self.getMerColNames(dbFileName),
                    self.getMerSize(dbFileName),
                    self.getNumMers(dbFileName),
                    self.getNumCons(dbFileName),
                    self.getNumBins(dbFileName),
                    self.isClustered(dbFileName),
                    self.isComplete(dbFileName),
                    self.getGMDBFormat(dbFileName))

        with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
            self.setMeta(h5file, meta_data, overwrite=True)

        # update the formatVersion field and we're done
        self.setGMDBFormat(dbFileName, 5)
        print "*******************************************************************************"


#------------------------------------------------------------------------------
# GET LINKS

    def restoreLinks(self, dbFileName, indices=[], silent=False):
        """Restore the links hash for a given set of indices"""
        full_record = []
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                full_record = [list(x) for x in h5file.root.links.links.readWhere("contig1 >= 0")]
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

        if indices == []:
            # get all!
            indices = self.getConditionalIndices(dbFileName, silent=silent)

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

    def getConditionalIndices(self, dbFileName, condition='', silent=False, checkUpgrade=True):
        """return the indices into the db which meet the condition"""
        # check the DB out and see if we need to change anything about it
        if checkUpgrade:
            self.checkAndUpgradeDB(dbFileName, silent=silent)

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

    def getTransformedCoverageProfiles(self, dbFileName, condition='', indices=np.array([])):
        """Load transformed coverage profiles"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([list(h5file.root.profile.transCoverage[x]) for x in indices])
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(h5file.root.profile.transCoverage[x.nrow]) for x in h5file.root.meta.contigs.where(condition)])
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getNormalisedCoverageProfiles(self, dbFileName, condition='', indices=np.array([])):
        """Load normalised coverage profiles"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([list(h5file.root.profile.normCoverage[x]) for x in indices])
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(h5file.root.profile.normCoverage[x.nrow]) for x in h5file.root.meta.contigs.where(condition)])
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def nukeBins(self, dbFileName):
        """Reset all bin information, completely"""
        print "    Clearing all old bin information from",dbFileName
        self.setBinStats(dbFileName, [])
        self.setNumBins(dbFileName, 0)
        self.setBinAssignments(dbFileName, updates={}, nuke=True)

    def setBinStats(self, dbFileName, updates, firstWrite=False):
        """Set bins table

        updates is a list of tuples which looks like:
        [ (bid, numMembers, isLikelyChimeric) ]
        """

        db_desc = [('bid', int),
                   ('numMembers', int),
                   ('isLikelyChimeric', bool)]
        bd = np.array(updates, dtype=db_desc)

        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                mg = h5file.getNode('/', name='meta')
                if not firstWrite:
                    t_name = 'tmp_bins'
                    # nuke any previous failed attempts
                    try:
                        h5file.removeNode(mg, 'tmp_bins')
                    except:
                        pass
                else:
                    t_name = 'bins'

                try:
                    h5file.createTable(mg,
                                       t_name,
                                       bd,
                                       title="Bin information",
                                       expectedrows=1)
                except:
                    print "Error creating META table:", exc_info()[0]
                    raise

                if not firstWrite:
                    # rename the tmp table to overwrite
                    h5file.renameNode(mg, 'bins', 'tmp_bins', overwrite=True)
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getBinStats(self, dbFileName):
        """Load data from bins table

        Returns a dict of type:
        { bid : [numMembers, isLikelyChimeric] }
        """
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                ret_dict = {}
                all_rows = h5file.root.meta.bins.read()
                for row in all_rows:
                    ret_dict[row[0]] = [row[1], row[2]]

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

    def setBinAssignments(self, storage, updates=None, image=None, nuke=False):
        """Set per-contig bins

        updates is a dictionary which looks like:
        { tableRow : binValue }
        if updates is set then storage is the
        path to the hdf file

        image is a list of tuples which look like:
        [(cid, bid, len, gc)]
        if image is set then storage is a tuple of type:
        (h5file, group)
        """
        db_desc = [('cid', '|S512'),
                   ('bid', int),
                   ('length', int),
                   ('gc', float)]
        closeh5 = False
        if updates is not None:
            # we need to build the image
            dbFileName = storage
            contig_names = self.getContigNames(dbFileName)
            contig_lengths = self.getContigLengths(dbFileName)
            contig_gcs = self.getContigGCs(dbFileName)
            num_cons = len(contig_lengths)
            if nuke:
                # clear all bin assignments
                bins = [0]*num_cons
            else:
                bins = self.getBins(dbFileName)

            # now apply the updates
            for tr in updates.keys():
                bins[tr] = updates[tr]

            # and build the image
            image = np.array(zip(contig_names, bins, contig_lengths, contig_gcs),
                             dtype=db_desc)

            try:
                h5file = tables.openFile(dbFileName, mode='a')
            except:
                print "Error opening DB:",dbFileName, exc_info()[0]
                raise
            meta_group = h5file.getNode('/', name='meta')
            closeh5 = True

        elif image is not None:
            h5file = storage[0]
            meta_group = storage[1]
            num_cons = len(image)
            image = np.array(image,
                             dtype=db_desc)
        else:
            print "get with the program dude"
            return

        # now we write the data
        try:
            # get rid of any failed attempts
            h5file.removeNode(meta_group, 'tmp_contigs')
        except:
            pass

        try:
            h5file.createTable(meta_group,
                               'tmp_contigs',
                               image,
                               title="Contig information",
                               expectedrows=num_cons)
        except:
            print "Error creating CONTIG table:", exc_info()[0]
            raise

        # rename the tmp table to overwrite
        h5file.renameNode(meta_group, 'contigs', 'tmp_contigs', overwrite=True)
        if closeh5:
            h5file.close()

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

    def getContigGCs(self, dbFileName, condition='', indices=np.array([])):
        """Load contig gcs"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([h5file.root.meta.contigs[x][3] for x in indices]).ravel()
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(x)[3] for x in h5file.root.meta.contigs.readWhere(condition)]).ravel()
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

    def getKmerPCAs(self, dbFileName, condition='', indices=np.array([])):
        """Load kmer sig PCAs"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                if(np.size(indices) != 0):
                    return np.array([list(h5file.root.profile.kpca[x]) for x in indices])
                else:
                    if('' == condition):
                        condition = "cid != ''" # no condition breaks everything!
                    return np.array([list(h5file.root.profile.kpca[x.nrow]) for x in h5file.root.meta.contigs.where(condition)])
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

#------------------------------------------------------------------------------
# GET / SET METADATA

    def getKmerVarPC(self, dbFileName, condition='', indices=np.array([])):
        """Load variance of kmer sig PCAs"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
              return np.array(h5file.root.meta.kpca_variance)[0]
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getTransformedCoverageCorners(self, dbFileName):
        """Load transformed coverage corners"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return np.array([list(x) for x in h5file.root.meta.transCoverageCorners.read()])
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def setMeta(self, h5file, metaData, overwrite=False):
        """Write metadata into the table

        metaData should be a tuple of values
        """
        db_desc = [('stoitColNames', '|S512'),
                   ('numStoits', int),
                   ('merColNames', '|S4096'),
                   ('merSize', int),
                   ('numMers', int),
                   ('numCons', int),
                   ('numBins', int),
                   ('clustered', bool),     # set to true after clustering is complete
                   ('complete', bool),      # set to true after clustering finishing is complete
                   ('formatVersion', int)]
        md = np.array([metaData], dtype=db_desc)

        # get hold of the group
        mg = h5file.getNode('/', name='meta')

        if overwrite:
            t_name = 'tmp_meta'
            # nuke any previous failed attempts
            try:
                h5file.removeNode(mg, 'tmp_meta')
            except:
                pass
        else:
            t_name = 'meta'

        try:
            h5file.createTable(mg,
                               t_name,
                               md,
                               "Descriptive data",
                               expectedrows=1)
        except:
            print "Error creating META table:", exc_info()[0]
            raise

        if overwrite:
            # rename the tmp table to overwrite
            h5file.renameNode(mg, 'meta', 'tmp_meta', overwrite=True)

    def getMetaField(self, dbFileName, fieldName):
        """return the value of fieldName in the metadata tables"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                # theres only one value
                return h5file.root.meta.meta.read()[fieldName][0]
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def setGMDBFormat(self, dbFileName, version):
        """Update the GMDB format version"""
        stoit_col_names = self.getStoitColNames(dbFileName)
        meta_data = (stoit_col_names,
                    len(stoit_col_names.split(',')),
                    self.getMerColNames(dbFileName),
                    self.getMerSize(dbFileName),
                    self.getNumMers(dbFileName),
                    self.getNumCons(dbFileName),
                    self.getNumBins(dbFileName),
                    self.isClustered(dbFileName),
                    self.isComplete(dbFileName),
                    version)
        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                self.setMeta(h5file, meta_data, overwrite=True)
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getGMDBFormat(self, dbFileName):
        """return the format version of this GM file"""
        # this guy needs to be a bit different to the other meta methods
        # becuase earlier versions of GM didn't include a format parameter
        with tables.openFile(dbFileName, mode='r') as h5file:
            # theres only one value
            try:
                this_DB_version = h5file.root.meta.meta.read()['formatVersion'][0]
            except ValueError:
                # this happens when an oldskool formatless DB is loaded
                this_DB_version = 0
        return this_DB_version

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

    def setNumBins(self, dbFileName, numBins):
        """set the number of bins"""
        stoit_col_names = self.getStoitColNames(dbFileName)
        meta_data = (stoit_col_names,
                    len(stoit_col_names.split(',')),
                    self.getMerColNames(dbFileName),
                    self.getMerSize(dbFileName),
                    self.getNumMers(dbFileName),
                    self.getNumCons(dbFileName),
                    numBins,
                    self.isClustered(dbFileName),
                    self.isComplete(dbFileName),
                    self.getGMDBFormat(dbFileName))
        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                self.setMeta(h5file, meta_data, overwrite=True)
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def getNumBins(self, dbFileName):
        """return the value of numBins in the metadata tables"""
        return self.getMetaField(dbFileName, 'numBins')

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

    def setClustered(self, dbFileName, state):
        """Set the state of clustering"""
        stoit_col_names = self.getStoitColNames(dbFileName)
        meta_data = (stoit_col_names,
                    len(stoit_col_names.split(',')),
                    self.getMerColNames(dbFileName),
                    self.getMerSize(dbFileName),
                    self.getNumMers(dbFileName),
                    self.getNumCons(dbFileName),
                    self.getNumBins(dbFileName),
                    state,
                    self.isComplete(dbFileName),
                    self.getGMDBFormat(dbFileName))
        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                self.setMeta(h5file, meta_data, overwrite=True)
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

    def isComplete(self, dbFileName):
        """Has this data set been *completely* clustered?"""
        try:
            with tables.openFile(dbFileName, mode='r') as h5file:
                return h5file.root.meta.meta.read()['complete']
        except:
            print "Error opening database:", dbFileName, exc_info()[0]
            raise

    def setComplete(self, dbFileName, state):
        """Set the state of completion"""
        stoit_col_names = self.getStoitColNames(dbFileName)
        meta_data = (stoit_col_names,
                    len(stoit_col_names.split(',')),
                    self.getMerColNames(dbFileName),
                    self.getMerSize(dbFileName),
                    self.getNumMers(dbFileName),
                    self.getNumCons(dbFileName),
                    self.getNumBins(dbFileName),
                    self.isClustered(dbFileName),
                    state,
                    self.getGMDBFormat(dbFileName))
        try:
            with tables.openFile(dbFileName, mode='a', rootUEP="/") as h5file:
                self.setMeta(h5file, meta_data, overwrite=True)
        except:
            print "Error opening DB:",dbFileName, exc_info()[0]
            raise

#------------------------------------------------------------------------------
# FILE / IO

    def dumpData(self, dbFileName, fields, outFile, separator, useHeaders):
        """Dump data to file"""
        header_strings = []
        data_arrays = []

        if fields == ['all']:
            fields = ['names', 'lengths', 'gc', 'bins', 'coverage', 'tcoverage', 'ncoverage', 'mers']

        num_fields = len(fields)
        data_converters = []

        for field in fields:
            if field == 'names':
                header_strings.append('cid')
                data_arrays.append(self.getContigNames(dbFileName))
                data_converters.append(lambda x : x)

            elif field == 'lengths':
                header_strings.append('length')
                data_arrays.append(self.getContigLengths(dbFileName))
                data_converters.append(lambda x : str(x))

            elif field == 'gc':
                header_strings.append('GCs')
                data_arrays.append(self.getContigGCs(dbFileName))
                data_converters.append(lambda x : str(x))

            elif field == 'bins':
                header_strings.append('bid')
                data_arrays.append(self.getBins(dbFileName))
                data_converters.append(lambda x : str(x))

            elif field == 'coverage':
                stoits = self.getStoitColNames(dbFileName).split(',')
                for stoit in stoits:
                    header_strings.append(stoit)
                data_arrays.append(self.getCoverageProfiles(dbFileName))
                data_converters.append(lambda x : separator.join(["%0.4f" % i for i in x]))

            elif field == 'tcoverage':
                header_strings.append('transformedCoverageX')
                header_strings.append('transformedCoverageY')
                header_strings.append('transformedCoverageZ')
                data_arrays.append(self.getTransformedCoverageProfiles(dbFileName))
                data_converters.append(lambda x : separator.join(["%0.4f" % i for i in x]))

            elif field == 'ncoverage':
                header_strings.append('normalisedCoverage')
                data_arrays.append(self.getNormalisedCoverageProfiles(dbFileName))
                data_converters.append(lambda x : separator.join(["%0.4f" % i for i in x]))

            elif field == 'mers':
                mers = self.getMerColNames(dbFileName).split(',')
                for mer in mers:
                    header_strings.append(mer)
                data_arrays.append(self.getKmerSigs(dbFileName))
                data_converters.append(lambda x : separator.join(["%0.4f" % i for i in x]))

        try:
            with open(outFile, 'w') as fh:
                if useHeaders:
                    header = separator.join(header_strings) + "\n"
                    fh.write(header)

                num_rows = len(data_arrays[0])
                for i in range(num_rows):
                    fh.write(data_converters[0](data_arrays[0][i]))
                    for j in range(1, num_fields):
                        fh.write(separator+data_converters[j](data_arrays[j][i]))
                    fh.write('\n')
        except:
            print "Error opening output file %s for writing" % outFile
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
        contigInfo = {} # save everything here first so we can sort accordingly
        for cid,seq,qual in self.readfq(contigFile):
            contigInfo[cid] = (kse.getKSig(seq.upper()), len(seq), self.calculateGC(seq))

        # sort the contig names here once!
        con_names = np.array(sorted(contigInfo.keys()))

        # keep everything in order...
        con_gcs = np.array([contigInfo[cid][2] for cid in con_names])
        con_lengths = np.array([contigInfo[cid][1] for cid in con_names])
        con_ksigs = np.array([contigInfo[cid][0] for cid in con_names])

        return (con_names, con_gcs, con_lengths, con_ksigs)

        # store the PCA'd kmersigs
        k_PCA_data = np.reshape(k_PCA_data, (rows,cols))
        self.storeSigPCAs(k_PCA_data, kPCATable)

    def calculateGC(self, seq):
      """Calculate fraction of nucleotides that are G or C."""
      testSeq = seq.upper()
      gc = testSeq.count('G') + testSeq.count('C')
      at = testSeq.count('A') + testSeq.count('T')

      return float(gc) / (gc + at)

    def PCAKSigs(self, kSigs, variance = 0.8):
        """PCA kmer sig data. All PCs require to capture the specified variance are returned.

        returns an array of tuples [(pc11, pc21, ..., pcN1), (pc12, pc22, ..., pnN2), ...]
        """

        # make a copy
        data = np.copy(kSigs)

        Center(data,verbose=0)
        p = PCA(data, fraction=variance)
        components = p.pc()

        return [tuple(i) for i in components], p.sumvariance[0:len(components[0])]

    def getWantedSeqs(self, contigFile, wanted, storage={}):
        """Do the heavy lifting of parsing"""
        print "Parsing contigs"
        for cid,seq,qual in self.readfq(contigFile):
            if(cid in wanted):
                storage[cid] = seq
        return storage

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
        """Work out the range of kmers required based on kmer length

        returns a list of sorted kmers and optionally a llo dict
        """
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
            return (sorted(ret_list), ll_dict)
        else:
            return sorted(ret_list)

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
        """Work out kmer signature for a nucleotide sequence

        returns a tuple of floats which is the kmer sig
        """
        # tmp storage
        sig = dict(zip(self.kmerCols, [0.0] * self.numMers))
        # the number fo kmers in this sequence
        num_mers = len(seq)-self.kLen+1
        for i in range(0,num_mers):
            try:
                sig[self.llDict[seq[i:i+self.kLen]]] += 1.0
            except KeyError:
                # typically due to an N in the sequence. Reduce the number of mers we've seen
                num_mers -= 1

        # normalise by length and return
        try:
            return tuple([sig[x] / num_mers for x in self.kmerCols])
        except ZeroDivisionError:
            print "***WARNING*** Sequence '%s' is not playing well with the kmer signature engine " % seq
            return tuple([0.0] * self.numMers)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class BamParser:
    """Parse multiple bam files and write the output to hdf5 """

    def __init__(self): pass

    def parse(self, bamFiles, contigNames, cid2Indices):
        """Parse multiple bam files and store the results in the main DB"""
        print "Importing BAM files"
        from bamtyper.utilities import BamParser as BTBP
        BP = BTBP()
        (links, ref_lengths, coverages) = BP.getLinks(bamFiles, full=False, verbose=True, doCoverage=True, minJoin=5)

        # go through all the contigs sorted by name.
        # we want to make an array of tuples of coverage sigs
        cov_sigs = []
        bam_range = range(len(bamFiles))
        for cid in contigNames:
            tmp_cov = []
            for i in bam_range:
                try:
                    tmp_cov.append(coverages[i][cid])
                except KeyError:
                    # may be no coverage for this contig
                    tmp_cov.append(0.0)
            cov_sigs.append(tuple(tmp_cov))

        # transform the links into something a little easier to parse later
        rowwise_links = []
        for cid in links:
            for link in links[cid]:
                try:
                    rowwise_links.append((cid2Indices[cid],          # contig 1
                                          cid2Indices[link[0]],      # contig 2
                                          int(link[1]),               # numReads
                                          int(link[2]),               # linkType
                                          int(link[3])                # gap
                                          ))
                except KeyError:
                    pass
        return (rowwise_links, np.array(cov_sigs))

def getBamDescriptor(fullPath, index_num):
    """AUX: Reduce a full path to just the file name minus extension"""
    return str(index_num) + '_' + op_splitext(op_basename(fullPath))[0]

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CoverageTransformer:
    """Transform coverage profiles GroopM-style"""

    def __init__(self,
                 numContigs,
                 numStoits,
                 normCoverages,
                 kmerNormPC1,
                 coverageProfiles,
                 stoitColNames,
                 scaleFactor=1000):
        self.numContigs = numContigs
        self.numStoits = numStoits
        self.normCoverages = normCoverages
        self.kmerNormPC1 = kmerNormPC1
        self.covProfiles = coverageProfiles
        self.stoitColNames = stoitColNames
        self.indices = range(self.numContigs)
        self.scaleFactor = scaleFactor

        # things we care about!
        self.TCentre = None
        self.transformedCP = np.zeros((self.numContigs,3))
        self.corners = np.zeros((self.numStoits,3))

    def transformCP(self, silent=False, nolog=False):
        """Do the main transformation on the coverage profile data"""
        shrinkFn = np.log10
        if(nolog):
            shrinkFn = lambda x:x

        if(not silent):
            print "    Reticulating splines"
            print "    Dimensionality reduction"

        unit_vectors = [(np.cos(i*2*np.pi/self.numStoits),np.sin(i*2*np.pi/self.numStoits)) for i in range(self.numStoits)]

        # make sure the bams are ordered consistently
        if self.numStoits > 3:
            self.shuffleBAMs()

        for i in range(len(self.indices)):
            shifted_vector = np.array([0.0,0.0])
            try:
                flat_vector = (self.covProfiles[i] / sum(self.covProfiles[i]))
            except FloatingPointError:
                flat_vector = self.covProfiles[i]

            for j in range(self.numStoits):
                shifted_vector[0] += unit_vectors[j][0] * flat_vector[j]
                shifted_vector[1] += unit_vectors[j][1] * flat_vector[j]

            # log scale it towards the centre
            scaling_vector = shifted_vector * self.scaleFactor
            sv_size = np.linalg.norm(scaling_vector)
            if sv_size > 1:
                shifted_vector /= shrinkFn(sv_size)

            self.transformedCP[i,0] = shifted_vector[0]
            self.transformedCP[i,1] = shifted_vector[1]
            # should always work cause we nuked
            # all 0 coverage vecs in parse
            self.transformedCP[i,2] = shrinkFn(self.normCoverages[i])

        # finally scale the matrix to make it equal in all dimensions
        min = np.amin(self.transformedCP, axis=0)
        self.transformedCP -= min
        max = np.amax(self.transformedCP, axis=0)
        max = max / (self.scaleFactor-1)
        self.transformedCP /= max

        # get the corner points
        XYcorners = np.reshape([i for i in np.array(unit_vectors)],
                               (self.numStoits, 2))

        for i in range(self.numStoits):
            self.corners[i,0] = XYcorners[i,0]
            self.corners[i,1] = XYcorners[i,1]

        # shift the corners to match the space
        self.corners -= min
        self.corners /= max

        # scale the corners to fit the plot
        cmin = np.amin(self.corners, axis=0)
        self.corners -= cmin
        cmax = np.amax(self.corners, axis=0)
        cmax = cmax / (self.scaleFactor-1)
        self.corners[:,0] /= cmax[0]
        self.corners[:,1] /= cmax[1]
        for i in range(self.numStoits):
            self.corners[i,2] = self.scaleFactor + 100 # only affect the z axis

        self.TCentre = np.mean(self.corners, axis=0)

    def small2indices(self, index, side):
        """Return the indices of the comparative items
        when given an index into a condensed distance matrix
        """
        step = 0
        while index >= (side-step):
            index = index - side + step
            step += 1
        return (step, step + index + 1)

    def shuffleBAMs(self, ordering=None):
        """Make the data transformation deterministic by reordering the bams"""
        if ordering is None:
            # we will need to deduce the ordering of the contigs
            # first we should make a subset of the total data
            # we'd like to take it down to about 1500 or so RI's
            # but we'd like to do this in a repeatable way
            ideal_contig_num = 1500
            sub_cons = np.arange(self.numContigs)
            while len(sub_cons) > ideal_contig_num:
                # select every second contig when sorted by norm cov
                cov_sorted = np.argsort(self.normCoverages[sub_cons])
                sub_cons = np.array([sub_cons[cov_sorted[i*2]] for i in np.arange(int(len(sub_cons)/2))])

                if len(sub_cons) > ideal_contig_num:
                    # select every second contig when sorted by mer PC1
                    mer_sorted = np.argsort(self.kmerNormPC1[sub_cons])
                    sub_cons = np.array([sub_cons[mer_sorted[i*2]] for i in np.arange(int(len(sub_cons)/2))])

            # now that we have a subset, calculate the distance between each of the untransformed vectors
            num_sc = len(sub_cons)

            # log shift the coverages towards the origin
            sub_covs = np.transpose([self.covProfiles[i]*(np.log10(self.normCoverages[i])/self.normCoverages[i]) for i in sub_cons])
            sq_dists = cdist(sub_covs,sub_covs,'cityblock')
            dists = squareform(sq_dists)

            # initialise a list of left, right neighbours
            lr_dict = {}
            for i in range(self.numStoits):
                lr_dict[i] = []

            too_big = 10000
            while True:
                closest = np.argmin(dists)
                if dists[closest] == too_big:
                    break
                (i,j) = self.small2indices(closest, self.numStoits-1)
                lr_dict[j].append(i)
                lr_dict[i].append(j)

                # mark these guys as neighbours
                if len(lr_dict[i]) == 2:
                    # no more than 2 neighbours
                    sq_dists[i,:] = too_big
                    sq_dists[:,i] = too_big
                    sq_dists[i,i] = 0.0
                if len(lr_dict[j]) == 2:
                    # no more than 2 neighbours
                    sq_dists[j,:] = too_big
                    sq_dists[:,j] = too_big
                    sq_dists[j,j] = 0.0

                # fix the dist matrix
                sq_dists[j,i] = too_big
                sq_dists[i,j] = too_big
                dists = squareform(sq_dists)

            # now make the ordering
            ordering = [0, lr_dict[0][0]]
            done = 2
            while done < self.numStoits:
                last = ordering[done-1]
                if lr_dict[last][0] == ordering[done-2]:
                    ordering.append(lr_dict[last][1])
                    last = lr_dict[last][1]
                else:
                    ordering.append(lr_dict[last][0])
                    last = lr_dict[last][0]
                done+=1

        # reshuffle the contig order!
        # yay for bubble sort!
        working = np.arange(self.numStoits)
        for i in range(1, self.numStoits):
            # where is this guy in the list
            loc = list(working).index(ordering[i])
            if loc != i:
                # swap the columns
                self.covProfiles[:,[i,loc]] = self.covProfiles[:,[loc,i]]
                self.stoitColNames[[i,loc]] = self.stoitColNames[[loc,i]]
                working[[i,loc]] = working[[loc,i]]

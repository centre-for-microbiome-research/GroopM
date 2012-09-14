#!/usr/bin/env python
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

from unittest import TestCase, main
import numpy as np
from os import getcwd, path
from groopm import mstore

###############################################################################

class MyTests(TestCase):
    def setUp(self):
        self.dataManager = mstore.GMDataManager()

    def test_writing(self):
        """Test methods for writing to a DB"""
        this_dir = getcwd()
        bam_files = [path.join(this_dir,'data','test1.bam'),
                     path.join(this_dir,'data','test2.bam'),
                     path.join(this_dir,'data','test3.bam')
                     ]
        reference = path.join(this_dir,'data','test.fa')
        db_file_name = path.join(this_dir,'data','test.h5')
        self.assertEqual(self.dataManager.createDB(bam_files,reference,db_file_name,dumpAll=False,force=True), True, "Test failed loading data")
    
    def test_reading(self):
        """Test methods for reading from a stored DB"""
        # test whether the stuff got saved correctly, use the slower way of retrieval
        db_file_name = path.join(getcwd(),'data','test.h5')
        
        # this is what the data should look like
        true_cov = np.array([[ 0.10017046,  0.0722952,   0.00030081],[ 0.03642941,  0.03653902,  0.0366121 ],[ 0.01494815,  0.01494815,  0.00845416]])
        true_bins = np.array([0,0,0])
        true_contig_names = np.array(['Contig_1_Cov_98.26','Contig_2_Cov_38.45','Contig_3_Cov_1195.29'])
        true_contig_lengths = np.array([9973,27368,66831])
        true_kmer_sigs = np.array([[0.00180487,0.00190514,0.00952572,0.00260704,0.04381831,0.00050135,
                                0.00381029,0.00862328,0.00260704,0.04913266,0.0089241,0.01353655,
                                0.01133059,0.00230623,0.01423844,0.01052843,0.03298907,0.00942545,
                                0.01383736,0.00812193,0.00310839,0.0157425,0.00421137,0.03399178,
                                0.00772085,0.00230623,0.00350948,0.00401083,0.00130352,0.00080217,
                                0.00842274,0.0007019,0.0123333,0.01463953,0.00360975,0.00230623,
                                0.0007019,0.00812193,0.00110298,0.00030081,0.00330893,0.00220596,
                                0.00521408,0.00190514,0.01123032,0.004813,0.00511381,0.00180487,
                                0.01514088,0.00531435,0.00421137,0.00150406,0.00451218,0.0,
                                0.00040108,0.00310839,0.0007019,0.00320866,0.00691868,0.01484007,
                                0.00742003,0.00220596,0.00501354,0.00541462,0.00300812,0.0024065,
                                0.00150406,0.00090244,0.00050135,0.00110298,0.01393763,0.004813,
                                0.00050135,0.00270731,0.0106287,0.01022761,0.0065176,0.00431164,
                                0.01433871,0.00631706,0.01143086,0.01032789,0.03359069,0.01463953,
                                0.00040108,0.00601624,0.00150406,0.00310839,0.01804873,0.00661787,
                                0.01614359,0.01113005,0.01363682,0.03679936,0.0164444,0.0123333,
                                0.00461245,0.00100271,0.00210569,0.00551489,0.00260704,0.00421137,
                                0.0075203,0.0082222,0.00220596,0.00320866,0.00531435,0.00661787,
                                0.00661787,0.01203249,0.00381029,0.00150406,0.00591597,0.00501354,
                                0.00671814,0.0058157,0.00902437,0.01443899,0.0089241,0.00010027,
                                0.00862328,0.00220596,0.00160433,0.00421137,0.00451218,0.0041111,
                                0.00290785,0.00350948,0.0116314,0.01353655,0.00862328,0.00200541,
                                0.00030081,0.0007019,0.00080217,0.00742003],
                                [0.00186349,0.01088863,0.0016808,0.00990208,0.00566355,0.01209442,
                                0.00274043,0.00175387,0.00292312,0.00847705,0.00551739,0.00588278,
                                0.00047501,0.00928091,0.00950015,0.00204619,0.00248465,0.00361736,
                                0.00336159,0.00570009,0.00887898,0.00175387,0.01388483,0.00186349,
                                0.00102309,0.00303274,0.00507892,0.00884244,0.0001827,0.01074247,
                                0.00562701,0.00179041,0.00650395,0.00021923,0.00570009,0.00836744,
                                0.00672318,0.00255773,0.00979246,0.01055978,0.01571178,0.00855013,
                                0.00829436,0.01070593,0.00774627,0.00511546,0.0059924,0.01289828,
                                0.00456738,0.00712511,0.01001169,0.01476177,0.01516369,0.01322713,
                                0.02093686,0.00855013,0.01227711,0.01319059,0.00851359,0.00968284,
                                0.00723473,0.00869629,0.01443291,0.00204619,0.01487138,0.00723473,
                                0.00675972,0.00723473,0.01286174,0.00690588,0.00208272,0.00544431,
                                0.03043701,0.00423853,0.00489623,0.0053347,0.00325197,0.00511546,
                                0.00632125,0.00445776,0.00500585,0.0061751,0.00146156,0.00668664,
                                0.01541947,0.00803859,0.00621163,0.00818474,0.00570009,0.0088059,
                                0.00164426,0.00091348,0.00624817,0.00354429,0.00522508,0.00752704,
                                0.01512716,0.01150979,0.00350775,0.00120579,0.00756358,0.00719819,
                                0.00650395,0.00906168,0.00584624,0.00825782,0.00690588,0.00427507,
                                0.00369044,0.00986554,0.0088059,0.00884244,0.00675972,0.00920783,
                                0.00080386,0.00190003,0.01008477,0.00931745,0.01026747,0.00672318,
                                0.00858667,0.01191172,0.00705203,0.0074905,0.00401929,0.00409237,
                                0.00869629,0.01344636,0.00332505,0.00303274,0.01344636,0.0066501,
                                0.00493277,0.02543116,0.00887898,0.00701549],
                                [0.00390537,0.00695785,0.0060451,0.01074352,0.0054765,0.00653888,
                                0.0049079,0.00567102,0.00582065,0.00655384,0.00510242,0.0094118,
                                0.00143646,0.00948662,0.01212012,0.00231928,0.00674837,0.00848409,
                                0.00395026,0.00451886,0.00670348,0.00550643,0.01219494,0.00417471,
                                0.00207987,0.00466849,0.00767608,0.01068366,0.00010474,0.01565142,
                                0.00938187,0.00258862,0.01605542,0.00118209,0.00811001,0.00972603,
                                0.00498272,0.00746659,0.00894794,0.00381559,0.00939684,0.00432434,
                                0.01238946,0.00656881,0.00649399,0.00620969,0.00694289,0.0087235,
                                0.01077344,0.00917239,0.0143197,0.00447397,0.00897787,0.00813994,
                                0.00932202,0.00767608,0.01101285,0.01351169,0.01396059,0.01098293,
                                0.0126588,0.00676333,0.00914246,0.00209484,0.00995047,0.0091275,
                                0.00549146,0.0070177,0.0061648,0.01013003,0.00347144,0.00763119,
                                0.01122234,0.00520716,0.00480316,0.00719726,0.0046236,0.00837934,
                                0.01674373,0.00873846,0.00439916,0.00923224,0.00350137,0.01174605,
                                0.01252413,0.00812497,0.00586554,0.00504257,0.00369589,0.02061917,
                                0.00202002,0.00284299,0.0057608,0.00667355,0.00613488,0.00965121,
                                0.01042929,0.00882824,0.00712244,0.00374078,0.00772097,0.00561117,
                                0.00640421,0.00854394,0.00453382,0.00691296,0.00620969,0.00549146,
                                0.00824468,0.0081549,0.00617977,0.0061648,0.01342191,0.00617977,
                                0.00543161,0.00326196,0.00942676,0.01078841,0.00966617,0.00248388,
                                0.00544657,0.00603014,0.00231928,0.00543161,0.00398019,0.00556628,
                                0.00671844,0.00923224,0.00592539,0.00706259,0.01340695,0.00523709,
                                0.00230432,0.00803519,0.00650896,0.01013003]
                                ])
        
        # get the data
        cov_profiles = self.dataManager.getCoverageProfiles(db_file_name)
        bins = self.dataManager.getBins(db_file_name)
        contig_names = self.dataManager.getContigNames(db_file_name)
        contig_lengths = self.dataManager.getContigLengths(db_file_name)
        kmer_sigs = self.dataManager.getKmerSigs(db_file_name)

        # and check it's good
        print np.around(true_cov.ravel(), decimals=3)
        print cov_profiles.ravel()
        self.assertTrue(np.array_equal(np.around(true_cov.ravel(), decimals=3), np.around(cov_profiles.ravel(), decimals=3)))
        self.assertTrue(np.equal(true_bins.all(), bins.all()))
        self.assertTrue(np.equal(true_contig_names, contig_names))
        self.assertTrue(np.equal(true_contig_lengths.all(), contig_lengths.all()))
        self.assertTrue(np.array_equal(np.around(true_kmer_sigs.ravel(), decimals=3), np.around(kmer_sigs.ravel(), decimals=3)))
        
        # test the pre-computed indicies method
        indicies = self.dataManager.getConditionalIndicies(db_file_name)
        
        self.assertTrue(np.array_equal(np.around(cov_profiles.ravel(), decimals=3), np.around(self.dataManager.getCoverageProfiles(db_file_name, indicies=indicies).ravel(), decimals=3)), "Discrepancy comparing cov profile retrieval method")
        self.assertTrue(np.equal(bins.all(), self.dataManager.getBins(db_file_name, indicies=indicies).all()), "Discrepancy comparing bins retrieval method")
        self.assertTrue(np.equal(contig_names, self.dataManager.getContigNames(db_file_name, indicies=indicies)), "Discrepancy comparing contig names retrieval method")
        self.assertTrue(np.equal(contig_lengths.all(), self.dataManager.getContigLengths(db_file_name, indicies=indicies).all()), "Discrepancy comparing contig lengths retrieval method")
        self.assertTrue(np.array_equal(np.around(kmer_sigs.ravel(), decimals=3), np.around(self.dataManager.getKmerSigs(db_file_name, indicies=indicies).ravel(), decimals=3)), "Discrepancy comparing kmer sigs profile retrieval method")

    def test_updateBin(self):
        """Test the ability to update the bins"""
        db_file_name = path.join(getcwd(),'data','test.h5')
        self.dataManager.setBins(db_file_name, {0:7,1:0,2:3})
        self.assertTrue(np.equal(np.array([7,0,3]).all(), self.dataManager.getBins(db_file_name).all()), "Error updating bins")
    
if __name__ == '__main__':
    main()


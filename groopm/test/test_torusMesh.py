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
from groopm import torusMesh

class MyTests(TestCase):
    
    def test_euclidean_dist(self):
        """Test euclidean distance"""
        tm = torusMesh.TorusMesh(49, columns=99)
        self.assertEqual(tm.dist((1,0),(1,1)),1.0,"dist error 1")
        self.assertEqual(tm.dist((1,0),(2,1)),np.sqrt(2),"dist error 2")
        self.assertEqual(tm.dist((0,0),(25,25)),np.sqrt(24**2+25**2),"dist error 3")
        self.assertEqual(tm.dist((0,0),(24,24)),np.sqrt(2*24**2),"dist error 4")
        self.assertEqual(tm.dist((0,0),(25,50)),np.sqrt(24**2+49**2),"dist error 5")
        
    def test_manhattan_dist(self):
        """Test manhattan distance"""
        tm = torusMesh.TorusMesh(49, columns=99)
        self.assertEqual(tm.dist((1,0),(1,1),type="man"),1.0,"dist error 1")
        self.assertEqual(tm.dist((1,0),(2,1),type="man"),2.0,"dist error 2")
        self.assertEqual(tm.dist((0,0),(25,25),type="man"),24.0+25.0,"dist error 3")
        self.assertEqual(tm.dist((0,0),(24,24),type="man"),48.0,"dist error 4")
        self.assertEqual(tm.dist((0,0),(25,50),type="man"),24.0+49.0,"dist error 5")
        
    def test_findNeighborhood(self):
        tm = torusMesh.TorusMesh(49, columns=99)
        fred = tm.findNeighborhood((5,5), 15)
        self.assertEqual(len(fred),749)
        
if __name__ == '__main__':
    main()

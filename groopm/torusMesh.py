#!/usr/bin/env python
###############################################################################
#                                                                             #
#    torusMesh.py                                                             #
#                                                                             #
#    Implements structure and methods for working with a torus shaped mesh    #
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
import numpy as np
from random import random
from PIL import Image, ImageDraw

np.seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class UnknownDistanceType(BaseException):
    pass

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TorusMesh:
    """A closed mesh, in the shape of a torus"""
    
    def __init__(self, rows, columns=0, dimension=1, randomize=False):
        """ init
        
        Set columns if you'd like something other than a square
        By default the torus is a scalar field. 
        Increase dimension if you'd like a vector field 
        """
        self.rows = rows
        if(columns == 0): # make it square
            self.columns = rows
        else:
            self.columns = columns
        self.dimension = dimension
        
        # we use these values in dist many many times
        self.halfRow = float(self.rows)/2
        self.halfColumn = float(self.columns)/2
        
        # the first two dimensions repesent points on the surface
        # the remainder represent values
        if(randomize):
            self.nodes = np.array([[ [random() for i in range(int(self.dimension))] for c in range(self.columns)] for r in range(self.rows)])
        else:
            self.nodes = np.array([[ [0.0 for i in range(int(self.dimension))] for c in range(self.columns)] for r in range(self.rows)])

        # make the colour lookup map
        if(self.dimension > 2):
            (self.colorLookup, self.cVec, self.maxAngle) = self.makeColorScheme()

#------------------------------------------------------------------------------
# WORKING WITH THE DATA

    def dist(self,A,B,type="euc"):
        """Work out the distance between points A and B
        
        Each point is a tuple of the form (R,C)
        """
        # we need to work out the deltas, not so straight forward on
        # the surface of a torus
        dr = np.abs(A[0] - B[0])
        dc = np.abs(A[1] - B[1])
        if(dr > self.halfRow):
            dr = self.rows - dr
        if(dc > self.halfColumn):
            dc = self.columns - dc
        if(type == "euc"):
            return np.sqrt(dr**2+dc**2)
        elif(type == "man"):
            return (dr + dc)
        else:
            raise UnknownDistanceType(type)

    def bestMatch(self, targetVector):
        """Returns location of the best match to an existing vector
        
        uses Euclidean distance
        """
        loc = np.argmin((((self.nodes - targetVector)**2).sum(axis=2))**0.5)
        col = np.mod(loc,self.columns)
        row = (loc - col)/self.columns
        col -= 1 
        return (row, col)

#------------------------------------------------------------------------------
# COLORING
    
    def makeColorScheme(self, dimension=0.0):
        """Work out how we are going to change a vector into a colour
        
        Returns a list of starts and lengths
        """
        if(dimension == 0.0):
            dimension = self.dimension
         
        step = np.ceil(dimension/3)
        start3 = dimension - step
        start2 = int(start3/2)

        # work out the central vector and corresponding max angle        
        c_vec = np.ones((step))
        c_vec = c_vec / np.linalg.norm(c_vec)
        top_vec = np.zeros_like(c_vec)
        top_vec[0] = 1
        max_angle = self.getAngBetweenNormed(top_vec, c_vec)

        return ([0,start2,start3,step], c_vec, max_angle)

    def getColor(self, vector):
        """return a colour for a given weight vector"""
        col = [0,0,0]
        if(self.dimension == 3):
            return [int(x * 255) for x in vector]
        for l in range(3):
            # grab the subset of the vector
            sub_vec = vector[self.colorLookup[l]:(self.colorLookup[l]+self.colorLookup[3]):1]
            # average and the turn into an rgb value
            sub_vec_size = np.linalg.norm(sub_vec) + 0.001
            col[l] = int(self.getAngBetweenNormed(sub_vec/sub_vec_size, self.cVec)/self.maxAngle*255)
        return col

    def getAngBetweenNormed(self, P1, P2):
        """Return the angle between two points (in radians)"""
        # find the existing angle between them theta
        c = np.dot(P1,P2) 
        # rounding errors hurt everyone...
        if(c > 1.0):
            c = 1.0
        elif(c < -1.0):
            c = -1.0
        return np.arccos(c) # in radians

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING 

    def __str__(self):
        """string method"""
        ret_array = []
        for r in range(self.rows):
            for c in range(self.columns):
                ret_array.append("[ ")
                for v in range(self.dimension):
                    ret_array.append(str(self.nodes[r,c,v])+" ")
                ret_array.append("],")
            ret_array.append("\n")
        return "".join(ret_array)

    def renderSurface(self, fileName):
        """make an image of the weights in the som"""
        try:
            img = Image.new("RGB", (self.columns, self.rows))
            for r in range(self.rows):
                # build a color value for a vector value
                for c in range(self.columns):
                    col = self.getColor(self.nodes[r,c])
                    img.putpixel((c,r), (col[0], col[1], col[2]))
            img = img.resize((self.columns*10, self.rows*10),Image.NEAREST)
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise
        
###############################################################################
###############################################################################
###############################################################################
###############################################################################

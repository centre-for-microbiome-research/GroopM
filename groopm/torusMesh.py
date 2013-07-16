#!/usr/bin/env python
###############################################################################
#                                                                             #
# torusMesh.py                                                                #
#                                                                             #
# Implements structure and methods for working with a torus shaped mesh       #
#                                                                             #
# Copyright (C) Michael Imelfort                                              #
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
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012/2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.1.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################
import sys
import numpy as np
from PIL import Image, ImageDraw
from scipy.spatial.distance import cdist
from colorsys import hsv_to_rgb as htr

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
        self.shape = (self.rows,self.columns,self.dimension)
        self.flatShape = (self.rows*self.columns,self.dimension)
        self.size = self.dimension*self.rows*self.columns
        
        # we use these values in dist many many times
        self.halfRow = float(self.rows)/2
        self.halfColumn = float(self.columns)/2
        
        # the first two dimensions repesent points on the surface
        # the remainder represent values
        if(randomize):
            self.nodes = np.random.random(self.size).reshape(self.shape)
        else:
            self.nodes = np.zeros(self.size).reshape(self.shape)

        # make an array of flattened nodes
        self.flatNodes = self.nodes.reshape(self.flatShape)
        
        # work out the central vector and corresponding max angle
        c_vec = np.ones((self.dimension))
        self.largestMag = np.linalg.norm(c_vec)
        self.cVec = c_vec / self.largestMag 
        top_vec = np.zeros_like(self.cVec)
        top_vec[0] = 1
        self.maxAngle = self.getAngBetweenNormed(top_vec, self.cVec)

    def fixFlatNodes(self, weights=None):
        """Make sure everything is in sync"""
        if weights is not None:
            self.nodes = weights
        self.flatNodes = self.nodes.reshape(self.flatShape)
        return self.flatNodes     
            
#------------------------------------------------------------------------------
# WORKING WITH THE DATA

    def bestMatch(self, targetVector):
        """Returns location of the best match to an existing vector
        uses Euclidean distance
        """
        loc = np.argmin(cdist(self.flatNodes, [targetVector]))
        row = int(loc/self.columns)
        col = loc-(row*self.rows)
        return [row, col]

    def buildVarianceSurface(self):
        """Work out the difference between each point and it's eight neighbours"""
        diff_array = np.zeros(self.shape)
        shift_array = np.zeros(self.shape)
        shift_array2 = np.zeros(self.shape)
        shift_diff = np.zeros(self.shape)

        # shift horizontal
        # ---
        # BA-
        # ---
        shift_array[:,:-1,:] = self.nodes[:,1:self.columns,:]   # shift left
        shift_array[:,-1,:] = self.nodes[:,0,:]                 # first node col to the end
        tmp_diff = np.abs(shift_array-self.nodes)           # get the difference
        diff_array += tmp_diff                              # add it on
        shift_diff[:,0,:] = tmp_diff[:,-1,:]                    # shift diff right
        shift_diff[:,1:self.columns,:] = tmp_diff[:,:-1,:]      # last to the first
        diff_array += shift_diff                            # add it on

        # shift horizontal vertical
        # B--
        # -A-
        # ---
        shift_array2[:-1,:,:] = shift_array[1:self.columns,:,:]
        shift_array2[-1,:,:] = shift_array[0,:,:]                
        tmp_diff = np.abs(shift_array2-self.nodes)          
        diff_array += tmp_diff                             
        shift_diff[0,:,:] = tmp_diff[-1,:,:]                   
        shift_diff[1:self.columns,:,:] = tmp_diff[:-1,:,:]     
        tmp_diff[:,0,:] = shift_diff[:,-1,:]
        tmp_diff[:,1:self.columns,:] = shift_diff[:,:-1,:]
        diff_array += tmp_diff

        # shift vertical
        # -B-
        # -A-
        # ---
        shift_array[:-1,:,:] = self.nodes[1:self.columns,:,:]
        shift_array[-1,:,:] = self.nodes[0,:,:]                
        tmp_diff = np.abs(shift_array-self.nodes)          
        diff_array += tmp_diff                             
        shift_diff[0,:,:] = tmp_diff[-1,:,:]                   
        shift_diff[1:self.columns,:,:] = tmp_diff[:-1,:,:]     
        diff_array += shift_diff                           

        # shift vertical horizontal
        # --B
        # -A-
        # ---
        shift_array2[:,0,:] = shift_array[:,-1,:]
        shift_array2[:,1:self.columns,:] = shift_array[:,:-1,:]
        tmp_diff = np.abs(shift_array2-self.nodes)
        diff_array += tmp_diff
        shift_diff[:,:-1,:] = tmp_diff[:,1:self.columns,:]
        shift_diff[:,-1,:] = tmp_diff[:,0,:]
        tmp_diff[0,:,:] = shift_diff[-1,:,:]                   
        tmp_diff[1:self.columns,:,:] = shift_diff[:-1,:,:]     
        diff_array += tmp_diff

        return diff_array
    
#------------------------------------------------------------------------------
# COLORING
    
    def getColor(self, vector):
        """return a colour for a given weight vector"""
        sn = np.linalg.norm(vector)
        if sn > 0:
            vv = vector / sn
            ang_perc = self.getAngBetweenNormed(vv, self.cVec)/self.maxAngle
            mag_perc = sn / self.largestMag
        else:
            ang_perc = 0.0
            mag_perc = 0.0
        V = 1       # VAL remain fixed at 1. Reduce to make pastels if that's your preference...
        col = [int(i*255) for i in htr(ang_perc, mag_perc, V)]
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

    def renderSurface(self, fileName, nodes=None):
        """make an image of the weights in the som"""
        if nodes is None:
            nodes = self.nodes
        ns = np.shape(nodes)
        rows = ns[0]
        columns = ns[1]
        try:
            img = Image.new("RGB", (columns, rows))
            for r in range(rows):
                # build a color value for a vector value
                for c in range(columns):
                    col = self.getColor(nodes[r,c])
                    img.putpixel((c,r), (col[0], col[1], col[2]))
            img = img.resize((columns*10, rows*10),Image.NEAREST)
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise
        
###############################################################################
###############################################################################
###############################################################################
###############################################################################

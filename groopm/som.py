#!/usr/bin/env python
from __future__ import division
###############################################################################
#                                                                             #
# som.py                                                                      #
#                                                                             #
# GroopM self organising map engine                                           #
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
# This script uses / is inspired by chunks of code and ideas originally
# developed by Kyle Dickerson and others
# Here: http://www.singingtree.com/~kyle/blog/archives/00000251.htm
# and Here: http://paraschopra.com/sourcecode/SOM/index.php
###############################################################################
## Kyle Dickerson
## kyle.dickerson@gmail.com
## Jan 15, 2008
##
## Self-organizing map using scipy
## This code is licensed and released under the GNU GPL
##
## This code uses a square grid rather than hexagonal grid, as scipy allows for
## fast square grid computation. I designed sompy for speed, so attempting to
## read the code may not be very intuitive. If you're trying to learn how SOMs
## work, I would suggest starting with Paras Chopras SOMPython code:
## http://www.paraschopra.com/sourcecode/SOM/index.php
## It has a more intuitive structure for those unfamiliar with scipy,
## however it is much slower. If you do use this code for something,
## please let me know, I'd like to know if has been useful to anyone.
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
import time
from random import *
from math import *
import sys
import numpy as np
import os
from scipy.spatial.distance import cdist
from PIL import Image, ImageDraw
import string
from numpy.random import random, randint
from random import randrange
np.seterr(all='raise')

# GroopM imports
from torusMesh import TorusMesh as TM
import rainbow

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class SOM:
    """A single instance of a self organising map"""
    def __init__(self, side, dimension, lc=None, uc=None):
        self.side = side # side length of the grid
        self.dimension = dimension # size of our grid vectors
        self.bestMatchCoords = [] # x/y coords of classified vectors

        self.radius = float(side)/2 # the radius of neighbour nodes which will be influenced by each new training vector
        
        # initialise the nodes to random values between 0 -> 1
        self.weights = TM(self.side, dimension=self.dimension, randomize=True)
        
        # usable area
        self.boundaryMask = np.zeros((self.side,self.side)) 
        
        # bin assignments
        self.binAssignments = np.zeros((self.side,self.side)) # the 0 bin is 'not assigned'
        
        if lc is not None:
            diff = uc - lc
            self.weights.nodes *= uc
            self.weights.nodes += lc
        self.regions = None
        
        # we'd like to know who is next to whom
        self.regionNeighbours = {}
        
    def getWeights(self):
        """Get the weights nodes"""
        return self.weights.nodes

    def getRegions(self):
        """Get the regions nodes"""
        return self.regions.nodes
    
    def loadWeights(self, nodes):
        """Use externally supplied data"""
        self.weights.nodes = nodes
        
    def loadRegions(self, nodes):
        """Use externally supplied regions"""
        self.regions = TM(self.side, dimension=1)
        self.regions.nodes = nodes
    
#------------------------------------------------------------------------------
# CLASSIFICATION

    # Classify!
    # run this after training to get an X/Y for each vector
    # you'd like to classify
    def classify(self, point):
        """Classify an individual point"""
        (row,col) = self.weights.bestMatch(point)
        return self.regions.nodes[row,col][0]

    def regionalise(self, bids, trainVector):
        """Create regions on the torus based on matches with the training vector
        Train vector is a list of numpy arrays
        """
        # build a regions structure
        self.regions = TM(self.side, dimension=1)
        for row in range(self.side):
            for col in range(self.side):
                self.regions.nodes[row,col] = self.classifyPoint(self.weights.nodes[row,col],
                                                                 trainVector,
                                                                 bids)
             
    def classifyPoint(self, point, trainVector, bids):
        """Returns the bid if of the best match to the trainVector
        trainVector and bids must be in sync
        """
        return bids[np.argmin((((trainVector - point)**2).sum(axis=1))**0.5)]

    def findRegionNeighbours(self):
        """Find out which regions neighbour which other regions"""
        neighbours = {}
        self.regionNeighbours = {}
        for row in range(self.side-1):
            for col in range(self.side-1):
                s_bid = self.regions.nodes[row,col][0]
                # test against right, down and diagonally down
                q_bids = [self.regions.nodes[row,col+1][0],
                          self.regions.nodes[row+1,col][0],
                          self.regions.nodes[row+1,col+1][0]
                         ]
                for q_bid in q_bids:
                    if(s_bid != q_bid):
                        # make storage for this region
                        if(s_bid not in self.regionNeighbours):
                            self.regionNeighbours[s_bid] = []
                        if(q_bid not in self.regionNeighbours):
                            self.regionNeighbours[q_bid] = []
                        # add in the neighbours
                        if(q_bid not in self.regionNeighbours[s_bid]):
                            self.regionNeighbours[s_bid].append(q_bid)
                        if(s_bid not in self.regionNeighbours[q_bid]):
                            self.regionNeighbours[q_bid].append(s_bid)

                        # we only need to return a tuple
                        nt = self.makeNTuple(s_bid,q_bid)
                        neighbours[nt] = True
        return neighbours.keys()
        
    def makeNTuple(self, bid1, bid2):
        """A way for making standard tuples from bids"""
        if(bid1 < bid2): return (bid1, bid2)
        return (bid2, bid1)

    def getNeighbours(self, bids):
        """return the neighbours of these bids"""
        ret_list = []
        for bid in bids:
            if(bid in self.regionNeighbours):
                ret_list.extend([i for i in self.regionNeighbours[bid] if i not in ret_list])
        return ret_list
        
#------------------------------------------------------------------------------
# TRAINING
        
    def train(self, trainVector, iterations=1000, vectorSubSet=1000, weightImgFileName="", epsilom=0.001):
        """Train the SOM
        Train vector is a list of numpy arrays
        """
        
        print " Start training. Max:", iterations, "iterations"

        # over time we'll shrink the radius of nodes which
        # are influenced by the current training node
        time_constant = iterations/log(self.radius)
        
        # we would ideally like to select guys from the training set at random
        if(len(trainVector) <= vectorSubSet):
            index_array = np.arange(len(trainVector))
            cut_off = len(trainVector) # if less than 1000 training vectors, set this to suit
        else:
            rand_index_array = np.arange(len(trainVector))
            cut_off = vectorSubSet
        
        for i in range(1, iterations+1):
            sys.stdout.write("\r Iteration: % 4d of % 4d" % (i, iterations))
            sys.stdout.flush()                

#--------
# Make stamp
            # gaussian decay on radius and amount of influence
            radius_decaying=self.radius*exp(-1.0*i/time_constant)
            if(radius_decaying < 2):
                return
            
            #phi = radius_decaying/3
            #rad_div_val = (phi)**2
            #front_multiplier = 1/np.sqrt(2 * np.pi * phi)
            front_multiplier = 0.4
            grad = -1 * front_multiplier / radius_decaying
            # we will make a "stamp" to speed things up
            max_radius = int(radius_decaying)
            q_stamp = np.zeros((max_radius+1,max_radius+1))

            for row in range(0, max_radius+1):
                for col in range(0, max_radius + 1):
                    # now we check to see that the euclidean distance is less than
                    # the specified distance.
                    true_dist = np.sqrt( row**2 + col**2 )
                    #if true_dist > 0.0: 
                    check_dist = np.round(true_dist+0.00001)
                    if(check_dist <= radius_decaying):
                        # influence is propotional to distance
                        influence = true_dist*grad + front_multiplier
                        q_stamp[row, col] = influence
                        q_stamp[col, row] = influence

            # divide by 2, so we don't mess up the stamp
            q_stamp[:,0] /= 2
            q_stamp[0,:] /= 2
            stamp = np.zeros((2*max_radius+1,2*max_radius+1))
            # bottom right
            stamp[max_radius:,max_radius:] += q_stamp
            # top right
            stamp[:max_radius+1,max_radius:] += np.rot90(q_stamp,1)
            # top left
            stamp[:max_radius+1,:max_radius+1] += np.rot90(q_stamp,2)
            # bottom left
            stamp[max_radius:,:max_radius+1] += np.rot90(q_stamp,3)
            # center
            stamp[max_radius, max_radius] = front_multiplier
            
            
            # now find where the useless info is and cull it from the stamp
            max_vals = np.max(stamp, axis=0)
            k = 0
            while k < len(max_vals):
                if max_vals[k] > 0.0:
                    break
                k += 1
            if k < len(max_vals):
                stamp = stamp[k:2*max_radius+1-k,k:2*max_radius+1-k]
            
            # keep track of how big the stamp is now
            stamp_side = len(stamp)
            stamp_radius = int((stamp_side-1)/2) 
            
            # if there are more than vectorSubSet training vecs
            # take a random selection
            if(len(trainVector) > vectorSubSet):
                np.random.shuffle(rand_index_array)
                index_array = rand_index_array[:cut_off]
#--------
# Make worksheet
            worksheet = np.zeros(self.dimension*self.side*self.side*9).reshape((self.side*3,
                                                                                self.side*3,
                                                                                self.dimension))
            worksheet[0:self.side,0:self.side] = self.weights.nodes
            worksheet[0:self.side,self.side:self.side*2] = self.weights.nodes
            worksheet[0:self.side,self.side*2:self.side*3] = self.weights.nodes
            worksheet[self.side:self.side*2,0:self.side] = self.weights.nodes
            worksheet[self.side:self.side*2,self.side:self.side*2] = self.weights.nodes
            worksheet[self.side:self.side*2,self.side*2:self.side*3] = self.weights.nodes
            worksheet[self.side*2:self.side*3,0:self.side] = self.weights.nodes
            worksheet[self.side*2:self.side*3,self.side:self.side*2] = self.weights.nodes
            worksheet[self.side*2:self.side*3,self.side*2:self.side*3] = self.weights.nodes

            # make a set of "delta nodes"
            # these contain the changes to the set of grid nodes
            # and we will add their values to the grid nodes
            # once we have input all the training nodes
            deltasheet = np.zeros_like(worksheet)

            for j in index_array:
                # find the best match between then training vector and the
                # current grid, inlined for greater speed
                loc = np.argmin(cdist(self.weights.flatNodes, [trainVector[j]]))
                row = int(loc/self.side)
                col = loc-(row*self.side)
                
                # row col represent the center of the stamp
                weights_patch = worksheet[self.side+row-stamp_radius:self.side+row+stamp_radius+1,
                                          self.side+col-stamp_radius:self.side+col+stamp_radius+1]
                weights_patch = -1*(weights_patch - trainVector[j])
                weights_patch[:,:,0] *= stamp
                weights_patch[:,:,1] *= stamp
                weights_patch[:,:,2] *= stamp
                weights_patch[:,:,3] *= stamp

                deltasheet[self.side+row-stamp_radius:self.side+row+stamp_radius+1,
                           self.side+col-stamp_radius:self.side+col+stamp_radius+1] += weights_patch
                           
            # now fold the deltas and update the weights
            deltasheet[:,self.side:2*self.side] += deltasheet[:,0:self.side] 
            deltasheet[:,self.side:2*self.side] += deltasheet[:,2*self.side:3*self.side] 
            deltasheet[self.side:2*self.side,self.side:2*self.side] += deltasheet[0:self.side,self.side:2*self.side]                 
            deltasheet[self.side:2*self.side,self.side:2*self.side] += deltasheet[2*self.side:3*self.side,self.side:2*self.side]

            # add the deltas to the grid nodes and clip to keep between 0 and 1
            self.weights.nodes = np.clip(self.weights.nodes + deltasheet[self.side:2*self.side,self.side:2*self.side], 0, 1)
            
            # make a tmp image, perhaps
            if(weightImgFileName != ""):
                filename = "%s_%04d.jpg" % (weightImgFileName, i)
                print " writing: %s" % filename
                self.weights.renderSurface(filename)

    def makeBoundaryMask(self, plotMaskFile=""):
        """Make a mask for cutting out boudaries"""
        # First create the mask
        VS = self.weights.buildVarianceSurface()
        VS_flat = np.array([[int(j) for j in i] for i in np.array(VS[:,:,0] + VS[:,:,1] + VS[:,:,2] + VS[:,:,3])*250]).reshape((self.side, self.side))
        # 10 is a magic number, but it seems to work OK
        self.boundaryMask = np.where(VS_flat > 10, 1., 0.)

        if plotMaskFile != "":
            self.renderBoundaryMask(plotMaskFile)

    def maskBoundaries(self, addNoise=False):
        """mask boundaries and add some random noise to
        some non-masked areas if asked to"""
        max_noise = 0.1
        noise_targets = 3
        for r in range(self.side):
            for c in range(self.side):
                if self.boundaryMask[r,c] == 1:
                    # on the boundary, mask as -1's
                    self.weights.nodes[r,c] = [-1.]*self.dimension
                elif addNoise:
                    if randint(10) <= noise_targets:
                        # add some noise
                        noise_amount = random() * max_noise + 1.0
                        self.weights.nodes[r,c] *= noise_amount 
        self.weights.flatNodes = self.weights.nodes.reshape((self.weights.rows*self.weights.columns,self.weights.dimension))

    def defineBinRegions(self, bids, binProfiles):
        """Work out which bins go where"""
        rcols = {}
        rand_col_lower = 15
        rand_col_upper = 200
        
        # use a flood fill algorithm to color in adjacent spots
        # and assign bins to unmasked points
        for i in range(len(bids)):
            # find out where this bin matches bestest
            [row, col] = self.weights.bestMatch(binProfiles[i])
            bid = bids[i]
            self.expandAssign(row, col, bid)
            rcols[bid] = (randrange(rand_col_lower, rand_col_upper),
                          randrange(rand_col_lower, rand_col_upper),
                          randrange(rand_col_lower, rand_col_upper)
                         )
            
        self.renderBoundaryMask("parry.png", rcols)

        # now clean up the mask
        for r in range(self.side):
            for c in range(self.side):
                if self.boundaryMask[r,c] == 0 and self.binAssignments[r,c] == 0:
                    # unmasked AND unassigned
                    self.boundaryMask[r,c] = 1
        
        self.renderBoundaryMask("parry_clean.png", rcols)

    def expandAssign(self, startR, startC, bid):
        """A floodfill algorithm, used to assign more and more
        points to a bin"""
        toFill = set()
        toFill.add((startR, startC))
        seen = {(startR, startC) : True}
        while len(toFill) != 0:
            (r,c) = toFill.pop()
            if self.boundaryMask[r,c] == 1:
                # we are at the boundary of a region
                continue
            if self.binAssignments[r,c] != 0:
                if self.binAssignments[r,c] != bid:
                    # we have already assigned this point to a bin
                    # oops...
                    print (r,c), "already assigned to bin %d, trying to reassign to bin %d" % (self.binAssignments[r,c], bid)
                    continue
            self.binAssignments[r,c] = bid
            
            # don't forget we're on a torus
            if r == 0: rm1 = self.side - 1
            else: rm1 = r - 1
            if c == 0: cm1 = self.side - 1
            else: cm1 = c - 1
            if r == self.side - 1: rp1 = 0
            else: rp1 = r + 1
            if c == self.side - 1:cp1 = 0
            else: cp1 = c + 1
            if (rm1,c) not in seen: toFill.add((rm1,c)); seen[(rm1,c)] = True
            if (rp1,c) not in seen: toFill.add((rp1,c)); seen[(rp1,c)] = True
            if (r,cm1) not in seen: toFill.add((r,cm1)); seen[(r,cm1)] = True
            if (r,cp1) not in seen: toFill.add((r,cp1)); seen[(r,cp1)] = True
                    
    def secondsToStr(self, t):
        rediv = lambda ll,b : list(divmod(ll[0],b)) + ll[1:]
        return "%d:%02d:%02d.%03d" % tuple(reduce(rediv,[[t*1000,],1000,60,60]))

#------------------------------------------------------------------------------
# IO and IMAGE RENDERING

    def renderWeights(self, tag):
        """Render the surface weights"""
        filename = tag+".png"
        self.weights.renderSurface(filename)

    def renderRegions(self, tag, palette):
        """Render the regions
        palette is a hash of bid -> color
        """
        filename = tag+".png"
        if(self.regions is None):
            raise ge.RegionsDontExistException
        try:
            img = Image.new("RGB", (self.weights.rows,self.regions.columns))
            for row in range(self.side):
                for col in range(self.side):
                    img.putpixel((col,row), palette[self.regions.nodes[row,col][0]])
            img = img.resize((self.weights.columns*10, self.weights.rows*10),Image.NEAREST)
            img.save(filename)
        except:
            print sys.exc_info()[0]
            raise
        
    def renderBoundaryMask(self, fileName, colMap=None):
        """Plot the boundary mask"""
        try:
            img = Image.new("RGB", (self.side, self.side))
            for r in range(self.side):
                for c in range(self.side):
                    if self.boundaryMask[r,c] == 0:
                        if colMap is not None:
                            try:
                                col = colMap[self.binAssignments[r,c]]
                            except KeyError:
                                col = (255,255,255)  
                        else:
                            col = (255,255,255)
                        img.putpixel((c,r), col)
                    else:
                        img.putpixel((c,r), (0,0,0))
            img = img.resize((self.side*10, self.side*10),Image.NEAREST)
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise
        
    def transColour(self, val):
        """Transform color value"""
        return 10 * log(val)

    def renderBestMatches(self, fileName, weighted=False):
        """make an image of where best matches lie
        set weighted to use a heatmap view of where they map
        """
        img_points = np.zeros((self.weights.rows,self.weights.columns))
        try:
            img = Image.new("RGB", (self.weights.columns, self.weights.rows))
            if(weighted): # color points by bestmatch density
                max = 0
                for point in self.bestMatchCoords:
                    img_points[point[0],point[1]] += 1
                    if(max < img_points[point[0],point[1]]):
                        max = img_points[point[0],point[1]]
                max += 1
                resolution = 200
                if(max < resolution):
                    resolution = max - 1
                max = self.transColour(max)
                rainbow = Rainbow.rainbow(0, max, resolution, "gbr")
                for point in self.bestMatchCoords:
                    img.putpixel((point[1],point[0]), rainbow.getColour(self.transColour(img_points[point[0],point[1]])))
            else: # make all best match points white
                for point in self.bestMatchCoords:
                    img.putpixel((point[1],point[0]), (255,255,255))
            img = img.resize((self.weights.columns*10, self.weights.rows*10),Image.NEAREST)
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise
    
###############################################################################
###############################################################################
###############################################################################
############################################################################### 

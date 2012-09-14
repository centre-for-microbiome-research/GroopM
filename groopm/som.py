#!/usr/bin/env python
from __future__ import division
###############################################################################
#                                                                             #
#    som.py                                                                   #
#                                                                             #
#    GroopM self organising map engine                                        #
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
##  http://www.paraschopra.com/sourcecode/SOM/index.php
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
from PIL import Image, ImageDraw
import string

# GroopM imports
from torusMesh import TorusMesh as tm
import rainbow

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class SOM:
    """A single instance of a self organising map"""
    def __init__(self, side, dimension):
        self.side = side                            # side length of the grid
        self.dimension = dimension                  # size of our grid vectors
        self.bestMatchCoords = []                   # x/y coords of classified vectors

        self.radius = float(side)/2.0               # the radius of neighbour nodes which will be influenced by each new training vector
        
        # initialise the nodes to random values between 0 -> 1
        self.weights = tm(self.side, dimension=self.dimension, randomize=True)
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
        self.regions = tm(self.side, dimension=1)
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
        self.regions = tm(self.side, dimension=1)
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
        
    def train(self, trainVector, iterations=1000, vectorSubSet = 1000, weightImgFileName="", epsilom=0.001):
        """Train the SOM
        
        Train vector is a list of numpy arrays
        """
        
        print "    Start training. Max:", iterations, "iterations"
        t0 = time.time()
        
        # over time we'll shrink the radius of nodes which
        # are influenced by the current training node
        time_constant = iterations/log(self.radius)
        
        # make a set of "delta nodes"
        # these contain the changes to the set of grid nodes
        # and we will add their values to the grid nodes
        # once we have input all the training nodes
        delta_nodes = np.zeros_like(self.weights.nodes)
        
        current_disp_time = "0" # sed for display
        current_total_time = "0"
        
        for i in range(1, iterations+1):
            # reset the deltas
            delta_nodes.fill(0.0)

            # gaussian decay on radius and amount of influence
            radius_decaying=self.radius*exp(-1.0*i/time_constant)
            if(radius_decaying < 2):
                return
            
            radius_ratio = radius_decaying/ self.side
             
            # calculate once!
            phi = radius_decaying/3
            rad_div_val = (phi)**2#2 * radius_decaying * i
            front_multiplier = 1/np.sqrt(2 * np.pi * phi)

            # we would ideally like to select guys from the training set at random
            
            index_array = np.arange(len(trainVector))
            np.random.shuffle(index_array)
            cut_off = int(len(index_array)/vectorSubSet)
            if(len(trainVector) < vectorSubSet):
                cut_off = len(trainVector) # if less than 1000 training vectors, set this to suit 
            elif(cut_off < vectorSubSet): # else, try to make sure there are at least 1000
                cut_off = vectorSubSet
            counter = 0
            t1 = time.time()
            for j in index_array:

                counter += 1
                if(counter > cut_off):
                    break 

                sys.stdout.write("\r    Iteration: % 4d of % 4d, Comparison: % 4d of % 4d, Radius: %04f / %04f " % (i, iterations, counter, cut_off, radius_decaying, radius_ratio) + "- (PREV: "+current_disp_time+", TOTAL: "+current_total_time+")")
                sys.stdout.flush()

                # find the best match between then training vector and the
                # current grid
                best = self.weights.bestMatch(trainVector[j])

                # apply the learning to the neighbours on the grid
                # the origin of the best match is alwys within the distance
                if(front_multiplier > epsilom):
                    delta_nodes[best[0],best[1]] += front_multiplier*(trainVector[j]-self.weights.nodes[best[0],best[1]])
                else:
                    return
                
                # now do the crosshairs
                max_radius = 0
                for dist in range(1, int(radius_decaying) + 1):
                    influence = front_multiplier * exp(-1*(dist**2)/rad_div_val)
                    if(influence > epsilom):
                        max_radius = dist 
                        # Up                            
                        Urow = np.mod((best[0]+dist),self.weights.rows) 
                        Ucol = best[1]
                        
                        # Down
                        Drow = np.mod((best[0]-dist),self.weights.rows) 
                        Dcol = best[1]
                        
                        # Left
                        Lrow = best[0]
                        Lcol = np.mod((best[1]-dist),self.weights.columns)
                        
                        # Right
                        Rrow = best[0]
                        Rcol = np.mod((best[1]+dist),self.weights.columns)
                                            
                        delta_nodes[Urow,Ucol] += influence*(trainVector[j]-self.weights.nodes[Urow,Ucol])
                        delta_nodes[Drow,Dcol] += influence*(trainVector[j]-self.weights.nodes[Drow,Dcol])
                        delta_nodes[Lrow,Lcol] += influence*(trainVector[j]-self.weights.nodes[Lrow,Lcol])
                        delta_nodes[Rrow,Rcol] += influence*(trainVector[j]-self.weights.nodes[Rrow,Rcol])
                    else:
                        break
                                
                # now do the main body
                for row in range(1, int(max_radius) + 1):
                    for col in range(1, int(max_radius) + 1):
                        # now we check to see that the euclidean distance is less than
                        # the specified distance.
                        true_dist = np.sqrt( row**2 + col**2 )
                        check_dist = np.round(true_dist+0.00001)
                        if(check_dist <= radius_decaying):

                            # influence is propotional to distance
                            influence = front_multiplier * exp(-1*(true_dist**2)/rad_div_val)
                            if(influence > epsilom):
                                #True                            
                                Trow = np.mod((best[0]+row),self.weights.rows) 
                                Tcol = np.mod((best[1]+col),self.weights.columns)
                                
                                # Mirrored down
                                Drow = np.mod((best[0]-row),self.weights.rows) 
                                Dcol = Tcol
                                
                                # Mirrored Left
                                Lrow = Trow
                                Lcol = np.mod((best[1]-col),self.weights.columns)
                                
                                # Mirrored diagonal
                                Grow = Drow
                                Gcol = Lcol
    
                                delta_nodes[Trow,Tcol] += influence*(trainVector[j]-self.weights.nodes[Trow,Tcol])
                                delta_nodes[Drow,Dcol] += influence*(trainVector[j]-self.weights.nodes[Drow,Dcol])
                                delta_nodes[Lrow,Lcol] += influence*(trainVector[j]-self.weights.nodes[Lrow,Lcol])
                                delta_nodes[Grow,Gcol] += influence*(trainVector[j]-self.weights.nodes[Grow,Gcol])
                            
            # add the deltas to the grid nodes and clip to keep between 0 and 1
            self.weights.nodes = np.clip(self.weights.nodes + delta_nodes, 0, 1)

            t2 = time.time()
            current_disp_time = self.secondsToStr(t2-t1)
            current_total_time = self.secondsToStr(t2-t0)
            
            # make a tmp image, perhaps
            if(weightImgFileName != ""):
                filename = weightImgFileName+"_%04d" % i+".png"
                print "   writing:",filename
                self.weights.renderSurface(filename)

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
            if(weighted):   # color points by bestmatch density 
                max = 0
                for point in self.bestMatchCoords:
                    img_points[point[0],point[1]] += 1
                    if(max < img_points[point[0],point[1]]):
                        max  = img_points[point[0],point[1]]
                max += 1
                resolution = 200
                if(max < resolution):
                    resolution = max - 1
                max  = self.transColour(max)
                rainbow = Rainbow.rainbow(0, max, resolution, "gbr")
                for point in self.bestMatchCoords:
                    img.putpixel((point[1],point[0]), rainbow.getColour(self.transColour(img_points[point[0],point[1]])))
            else:       # make all best match points white
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

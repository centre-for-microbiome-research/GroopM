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
from random import *
from math import *
import sys
import numpy as np
import os
from PIL import Image, ImageDraw
import string

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class SOM:
    def __init__(self, side, V_size, learning_rate=0.05):
        self.side = side                            # side length of the grid
        self.V_size = V_size                      # size of our grid vectors
        self.radius = side/2                        # the radius of neighbour nodes which will be influenced by each new training vector
        self.learning_rate = learning_rate          # the amount of influence
        self.plots = []                             # x/y coords of classified vectors
        self.ref_positions = {}                     # map of contig refs to positions in the SOM         
        
        # initialise the nodes to random values between 0 -> 1
        self.nodes = np.array([[ [random() for i in range(V_size)] for x in range(self.side)] for y in range(self.side)])

        # for when we want to autocolor the background
        self.colorLookup = self.getColorLookup()

#------------------------------------------------------------------------------
# CLASSIFICATION 

    # Classify!
    # run this after training to get an X/Y for each vector 
    # you'd like to classify
    def classify(self, names_vector=[], data_vector=[[]]):
        print "Start classification..."
        index_array = np.arange(len(data_vector))
        for j in index_array:
            # find the best match between then data vector and the
            # current grid
            best = self.best_match(data_vector[j])
            self.plots.append(best)
        print "Done"

#------------------------------------------------------------------------------
# TRAINING 
        
    def train(self, iterations=1000, train_vector=[[]], weightImgFileName=""):
        """Train the SOM
        
        Train vector is a list of numpy arrays
        """
        
        print "Start training -->", iterations, "iterations"
        
        # over time we'll shrink the radius of nodes which
        # are influenced by the current training node
        time_constant = iterations/log(self.radius)
        
        # make a set of "delta nodes"
        # these contain the changes to the set of grid nodes
        # and we will add their values to the grid nodes
        # once we have input all the training nodes
        delta_nodes = np.array([[[0.0 for i in range(self.V_size)] for x in range(self.side)] for y in range(self.side)])
        
        for i in range(1, iterations+1):
            delta_nodes.fill(0.0)
            # gaussian decay on radius and amount of influence
            radius_decaying=self.radius*exp(-1.0*i/time_constant)
            learning_rate_decaying=self.learning_rate*exp(-1.0*i/time_constant)

            # calculate once!
            rad_div_val = 2 * radius_decaying * i
            
            # we would ideally like to select guys from the training set at random
            index_array = np.arange(len(train_vector))
            np.random.shuffle(index_array)
            cut_off = int(len(index_array)/1000)
            if(cut_off < 1000):
                cut_off = 1000
            counter = 0
            for j in index_array:
                counter += 1
                if(counter > cut_off):
                    break 
                # find the best match between then training vector and the
                # current grid
                best = self.best_match(train_vector[j])
                
                # apply the learning to the neighbors on the grid
                for loc in self.find_neighborhood(best, radius_decaying):
                    influence = exp( (-1.0 * (loc[2]**2)) / rad_div_val)
                    inf_lrd = influence*learning_rate_decaying
                    delta_nodes[loc[0],loc[1]] += inf_lrd*(train_vector[j]-self.nodes[loc[0],loc[1]])
            
            # add the deltas to the grid nodes
            self.nodes += delta_nodes
            
            # renormalise the vector
            self.reNorm()
            
            # make a tmp image, perhaps
            if(weightImgFileName != ""):
                filename = "%04d%s" % (i, "." + weightImgFileName)
                self.renderWeightsImage(filename)

#------------------------------------------------------------------------------
# SOM STUFF 

    def reNorm(self):
        """set all weight vectors to values between 0 and 1
        
        some of the values will want to creep above 1 or below 0
        this function renormalises the grid
        """
        for r in range(self.side):
            for c in range(self.side):
                for v in range(self.V_size):
                    if(self.nodes[r,c,v] < 0):
                        self.nodes[r,c,v] = 0
                    elif(self.nodes[r,c,v] > 1):
                        self.nodes[r,c,v] = 1
   
    def find_neighborhood(self, pt, dist):
        """Returns a list of points which live within 'dist' of 'pt'
        
        mapped onto a torus! pt = [row,col]
        """  
        neighbors = []
        # first we bound ourselves to the oriented square surrounding 
        # the point of interest 
        for y in range(int(pt[0] - dist), int(pt[0] + dist)):
            for x in range(int(pt[1] - dist), int(pt[1] + dist)):
                # now we check to see that the euclidean distance is less than
                # the specified distance.
                this_dist = sqrt( pow(pt[1] - x , 2) + pow(pt[0] - y , 2) )
                if(this_dist <= dist):
                    x_out = x
                    y_out = y
                    if(x >= self.side):
                        x_out -= self.side
                    elif(x < 0):
                        x_out += self.side
                    if(y >= self.side):
                        y_out -= self.side
                    elif(y < 0):
                        y_out += self.side
                    neighbors.append((y_out,x_out,this_dist))
        return neighbors
    
    def best_match(self, target_FV):
        """Returns location of best match, uses Euclidean distance
        
        target_FV is a numpy array
        """
        loc = np.argmin((((self.nodes - target_FV)**2).sum(axis=2))**0.5)
        r = 0
        while loc > self.side:
            loc -= self.side
            r += 1
        c = loc - 1
        return (r, c)

#------------------------------------------------------------------------------
# IMAGE RENDERING 

    def renderWeightsImage(self, fileName):
        """make an image of the weights in the som"""
        try:
            img = Image.new("RGB", (self.side, self.side))
            for r in range(self.side):
                # build a color value for a vector value
                for c in range(self.side):
                    col = self.getColor(self.nodes[r,c])
                    img.putpixel((c,r), (col[0], col[1], col[2]))
            img = img.resize((self.side*10, self.side*10),Image.NEAREST)
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise

    def renderClassImage(self, fileName):
        """make an image of raw classifications""" 
        try:
            img = Image.new("RGB", (self.side, self.side))
            for point in self.plots:
                img.putpixel((point[1],point[0]), (255,0,0))
            img = img.resize((self.side*10, self.side*10),Image.NEAREST)
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise

    def renderWeightedClassImage(self, fileName):
        """make an image of weighted classifications"""
        weighted_nodes = np.array([[0 for x in range(self.side)] for y in range(self.side)])
        max = 0
        for point in self.plots:
            weighted_nodes[point[1],point[0]] += 1
            if(max < weighted_nodes[point[1],point[0]]):
                max  = weighted_nodes[point[1],point[0]]

        max += 1
        res = 200
        if(max < res):
            res = max - 1
        max  = self.transColour(max)

        rainbow = Rainbow.rainbow(0, max, res, "gbr")
        
        try:
            img = Image.new("RGB", (self.side, self.side))
            for point in self.plots:
                img.putpixel((point[1],point[0]), rainbow.getColour(self.transColour(weighted_nodes[point[1],point[0]])))
            img = img.resize((self.side*10, self.side*10),Image.NEAREST)
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise
        
    def renderContigImage(self, fileName):
        """make an colour image of raw classifications"""
        sep_size = 50
        try:
            img = Image.new("RGB", (self.side, self.side))
            img2 = Image.new("RGB", (self.side, self.side))
            draw = ImageDraw.Draw(img)
            for key,val in self.ref_positions.items():
                lrp = len(val)
                if(lrp > 1):
                    for index in range(lrp-1):
                        if((abs(val[index][0] - val[index+1][0]) < sep_size) and (abs(val[index][1] - val[index+1][1]) < sep_size)): 
                            draw.line((val[index][1], val[index][0], val[index+1][1], val[index+1][0]), fill=128)
                        else:
                            img2.putpixel((val[index][1],val[index][0]), (255,0,0))
                            img2.putpixel((val[index+1][1],val[index+1][0]), (255,0,0))
            del draw 
            img = img.resize((self.side*10, self.side*10),Image.NEAREST)
            img2 = img2.resize((self.side*10, self.side*10),Image.NEAREST)
            img.save(fileName)
            img2.save(fileName+".png")
        except:
            print sys.exc_info()[0]
            raise
    
    def transColour(self, val):
        """Transform color value"""
        return 10 * log(val)
    
    def getColor(self, target_FV):
        """return a colour for a given weight vector"""
        col = [0,0,0]
        for l in range(3):
            # grab the subset of the vector we care about
            sub_vec = target_FV[self.colorLookup[l]:(self.colorLookup[l]+self.colorLookup[3]):1]

            # multiply each value in the sub vector by it's index in the vector
            vec_sum = 0
            divisor = 1
            for s in range(self.colorLookup[3]):
                m = 1 + (s + 1)/10
                vec_sum += (sub_vec[s]*m) 
                divisor *= m
            
            # renormalise
            vec_sum /= divisor
            
            # average and the turn into an rgb value
            col[l] = int(vec_sum/self.colorLookup[3]*255)
            
        return col
        
    def getColorLookup(self):
        """Work out how we are going to change a vector into a colour
        
        Returns a list of starts and lengths
        """
        if(self.V_size == 3):
            return [0,1,2,1]
        elif(self.V_size == 4):
            return [0,1,2,2]
        elif(self.V_size == 5):
            return [0,1,2,3]
        elif(self.V_size == 6):
            return [0,2,4,2]
        elif(self.V_size == 7):
            return [0,2,4,3]
        elif(self.V_size == 8):
            return [0,2,4,4]
        elif(self.V_size == 9):
            return [0,2,4,5]
        elif(self.V_size == 10):
            return [0,3,6,4]
        elif(self.V_size == 11):
            return [0,3,6,5]
        elif(self.V_size == 12):
            return [0,3,6,6]
        elif(self.V_size == 13):
            return [0,3,6,7]
        elif(self.V_size == 14):
            return [0,4,8,6]
        elif(self.V_size == 15):
            return [0,4,8,7]
        else:
            print "***ERROR: Max vector size of 15! Time to be a haxxor!"
            sys.exit(1)
    
#------------------------------------------------------------------------------
# FILE IO 
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Rainbow:
    def __init__(self, lb, ub, res, type="rb"):
        """Simple dimple heatmap
        
        Specify the upper and lower bounds for your data. 
        resolution refers to the number of bins which are available in this space
        Supports four heatmap types is red-blue, blue-red, red-green-blue and blue-green-red
        """ 

        # constants
        self.RB_lower_offset = 0.5
        self.RB_divisor = (2.0/3.0)
        self.RB_ERROR_COLOUR = [0,0,0]

        # set the limits
        self.lowerBound = lb
        self.upperBound = ub
        self.resolution = res
        self.tickSize = (self.upperBound - self.lowerBound)/(self.resolution - 1)
        
        # set the type, red-blue by default
        self.type = type
        self.redOffset = 0.0
        self.greenOffset = self.RB_divisor * math.pi * 2.0
        self.blueOffset = self.RB_divisor * math.pi
        
        self.ignoreRed = False
        self.ignoreGreen = True
        self.ignoreBlue = False
        
        self.lowerScale = 0.0
        self.upperScale = self.RB_divisor * math.pi
        
        if(self.type == "rbg"): # red-blue-green
            self.redOffset = 0.0
            self.greenOffset = self.RB_divisor * math.pi * 2.0
            self.blueOffset = self.RB_divisor * math.pi
            
            self.ignoreRed = False
            self.ignoreGreen = False
            self.ignoreBlue = False
            
            self.lowerScale = 0.0
            self.upperScale = (self.RB_divisor * math.pi * 2.0)

        elif(self.type == "gbr"): # green-blue-red
            self.redOffset = self.RB_divisor * math.pi * 2.0
            self.greenOffset = 0.0
            self.blueOffset = self.RB_divisor * math.pi

            self.ignoreRed = False
            self.ignoreGreen = False
            self.ignoreBlue = False
            
            self.lowerScale = 0.0
            self.upperScale = (self.RB_divisor * math.pi * 2.0)

        elif(self.type == "br"): # blue-red
            self.redOffset = self.RB_divisor * math.pi
            self.greenOffset = self.RB_divisor * math.pi * 2.0
            self.blueOffset = 0.0
            
            self.ignoreRed = False
            self.ignoreGreen = True
            self.ignoreBlue = False

            self.lowerScale = 0.0
            self.upperScale = (self.RB_divisor * math.pi)

        self.scaleMultiplier = (self.upperScale - self.lowerScale)/(self.upperBound - self.lowerBound)
    
    def getValue(self, val):
        """Get a raw value, not a colour"""
        return (math.cos(val) + self.RB_lower_offset) * self.RB_divisor

    def getColour(self, val):
        """Return a colour for the given value.
        
        If nothing makes sense. return black
        """
        if(val > self.upperBound or val < self.lowerBound):
            return self.RB_ERROR_COLOUR
        
        # normalise the value to suit the ticks
        normalised_value = round(val/self.tickSize) * self.tickSize
    
        # map the normalised value onto the horizontal scale
        scaled_value = (normalised_value - self.lowerBound) * self.scaleMultiplier + self.lowerScale
            
        red = 0
        green = 0
        blue = 0
        
        if(not self.ignoreRed):
            red = int(round(self.getValue(scaled_value - self.redOffset) * 255))
        if(not self.ignoreGreen):
            green = int(round(self.getValue(scaled_value - self.greenOffset) * 255))
        if(not self.ignoreBlue):
            blue = int(round(self.getValue(scaled_value - self.blueOffset) * 255))
    
        return (red, green, blue)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

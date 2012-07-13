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

import argparse
import sys
from random import *
from math import *
import sys
import scipy
import numpy
import csv
import os
import Rainbow
from PIL import Image, ImageDraw
import string

###############################################################################
###############################################################################
###############################################################################
###############################################################################
def doWork( options ):
    # get a working dir
    if not os.path.exists(options.working_dir):
        os.makedirs(options.working_dir)

    # load the data file
    data = []
    names = []
    with open(options.data_file, 'rb') as f:
        reader = csv.reader(f, delimiter=options.delimiter)
        headerline = reader.next()
        for row in reader:
            names.append(row.pop(0))
            for t in range(len(row)):
                row[t] = float(row[t])
            data.append(scipy.array(row))

    # get us a som
    som = SOM(options.working_dir, options.grid_size, len(data[0]), options.learning_rate, options.weight_img)

    # we only bother doing training if we're told to
    if(not options.classify_only):
        # start training
        som.train(options.reps, data, options.animate)
        som.renderWeightsImage(os.path.join(options.working_dir, options.weight_img))
        som.printWeights(os.path.join(options.working_dir, options.weight_csv))
    else:
        # check to see if the weights file is there:
        if(not os.path.isfile(os.path.join(options.working_dir, options.weight_csv))):
           print "Error: no weights file found at classification step"
           sys.exit(1)
        
        # override the weights in the som, no need if we've just finished training
        print "Loading weights..."
        som.overrideGrid(os.path.join(options.working_dir, options.weight_csv))
        
    # now classify this mofo!
    som.classify(names, data)
    #som.renderWeightedClassImage(os.path.join(options.working_dir, options.class_img))
    
    # only do contig validation if there are refs
    if("" != options.refs):
        som.groupContigs(options.refs, names)
        som.renderContigImage(os.path.join(options.working_dir, options.contig_img))
    
    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################
class SOM:
    def __init__(self, working_dir, side=10, FV_size=10, learning_rate=0.05, img="messom.png"):
        self.side = side                            # side length of the grid
        self.FV_size = FV_size                      # size of our grid vectors
        self.radius = side/2                        # the radius of neighbour nodes which will be influenced by each new training vector
        self.learning_rate = learning_rate          # the amount of influence
        self.weight_img_name = img                  # name of the image we produce to visualise the weights
        self.working_dir = working_dir              # where we will save any files etc...
        self.plots = []                             # x/y coords of classified vectors
        self.ref_positions = {}                     # map of contig refs to positions in the SOM         
        
        # initialise the nodes to random values between 0 -> 1
        self.nodes = scipy.array([[ [random() for i in range(FV_size)] for x in range(self.side)] for y in range(self.side)])
        
        print "Initiated grid with:", (self.side * self.side), "points. (", self.side , "X", self.side, ")" 

        self.colorLookup = self.getColorLookup()

#------------------------------------------------------------------------------
# CLASSIFICATION 

    # Classify!
    # run this after training to get an X/Y for each vector 
    # you'd like to classify
    def classify(self, names_vector=[], data_vector=[[]]):
        print "Start classification..."
        index_array = numpy.arange(len(data_vector))
        for j in index_array:
            # find the best match between then data vector and the
            # current grid
            best = self.best_match(data_vector[j])
            self.plots.append(best)
        print "Done"

#------------------------------------------------------------------------------
# TRAINING 
        
    def train(self, iterations=1000, train_vector=[[]], animate=False):
        """Train the SOM
        
        Train vector is a list of scipy arrays
        """
        
        print "Start training -->", iterations, "iterations"
        
        # over time we'll shrink the radius of nodes which
        # are influenced by the current training node
        time_constant = iterations/log(self.radius)
        
        # make a set of "delta nodes"
        # these contain the changes to the set of grid nodes
        # and we will add their values to the grid nodes
        # once we have input all the training nodes
        delta_nodes = scipy.array([[[0.0 for i in range(self.FV_size)] for x in range(self.side)] for y in range(self.side)])
        
        for i in range(1, iterations+1):
            print "Iteration:", i
            delta_nodes.fill(0.0)
            # gaussian decay on radius and amount of influence
            radius_decaying=self.radius*exp(-1.0*i/time_constant)
            learning_rate_decaying=self.learning_rate*exp(-1.0*i/time_constant)

            # calculate once!
            rad_div_val = 2 * radius_decaying * i
            
            # we would ideally like to select guys from the training set at random
            index_array = numpy.arange(len(train_vector))
            numpy.random.shuffle(index_array)
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
            if(animate):
                filename = "%04d%s" % (i, "." + self.weight_img_name)
                full_file_name = os.path.join(self.working_dir, filename)
                self.renderWeightsImage(full_file_name)
            
        print "Done"

#------------------------------------------------------------------------------
# EXTRA 

    # Group contigs together so that we can see if they actually do group together
    def groupContigs(self, refFiles, names):
        # first check if it is a dir or a list of files
        refs = []
        if os.path.isdir(refFiles):
            short_refs = [file for file in os.listdir(refFiles) if file.lower().endswith(".ref")]
            for ref in short_refs:
                refs.append(os.path.join(refFiles, ref))
        else:
            refs = string.split(refFiles, ",")

        # now that we have consistent interface to refs,
        # it's time to make a dict of all guys from each contig
        for ref in refs:
            with open(ref, 'rb') as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    self.ref_positions[row[0]] = []
                    
        # make a dict of all frag names and link to their positions
        pos_lookup = {}
        for dat_counter in range(len(self.plots)):
            pos_lookup[names[dat_counter]] = self.plots[dat_counter]
             
        # now go through all the names
        # we have and clump them together
        for name in names:
            name_fields = string.split(name, "_")
            keepable = name_fields[0:len(name_fields)-4:1]
            short_name = string.join(keepable,"_");
            self.ref_positions[short_name].append(pos_lookup[name])
        
#------------------------------------------------------------------------------
# SOM STUFF 
        
    # override the values in the weights grid
    # fileName is a vanilla tsv file of weights
    def overrideGrid(self, fileName):
        with open(fileName, 'rb') as f:
            r = 0
            c = 0
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                self.nodes[r,c] = scipy.array(row)
                c += 1
                if(c == self.side):
                    c = 0
                    r += 1
        # check to see thatwe are sane!
        if(r != self.side or c != 0):
            print "Error: Do the dimesions of the supplied weight data match those of the som?"
        return 0

    # some of the values will want to creep above 1 or below 0
    # this function renormalises the grid
    def reNorm(self):
        for r in range(self.side):
            for c in range(self.side):
                for v in range(self.FV_size):
                    if(self.nodes[r,c,v] < 0):
                        self.nodes[r,c,v] = 0
                    elif(self.nodes[r,c,v] > 1):
                        self.nodes[r,c,v] = 1
   

    # Returns a list of points which live within 'dist' of 'pt'
    # mapped onto a torus! pt = [row,col]
    def find_neighborhood(self, pt, dist):
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
    
    # Returns location of best match, uses Euclidean distance
    # target_FV is a scipy array
    def best_match(self, target_FV):
        loc = scipy.argmin((((self.nodes - target_FV)**2).sum(axis=2))**0.5)
        r = 0
        while loc > self.side:
            loc -= self.side
            r += 1
        c = loc - 1
        return (r, c)

#------------------------------------------------------------------------------
# IMAGE RENDERING 

    # make an image of the weights in the som
    def renderWeightsImage(self, fileName):
        print "Rendering weights image"
        try:
            img = Image.new("RGB", (self.side, self.side))
            for r in range(self.side):
                # build a color value for a vector value
                for c in range(self.side):
                    col = self.getColor(self.nodes[r,c])
                    img.putpixel((c,r), (col[0], col[1], col[2]))
            img = img.resize((self.side*10, self.side*10),Image.NEAREST)
            print "Saving image to:", fileName
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise

    # make an image of raw classifications 
    def renderClassImage(self, fileName):
        print "Rendering classification image"
        try:
            img = Image.new("RGB", (self.side, self.side))
            for point in self.plots:
                img.putpixel((point[1],point[0]), (255,0,0))
            img = img.resize((self.side*10, self.side*10),Image.NEAREST)
            print "Saving image to:", fileName
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise

    # make an image of weighted classifications
    def renderWeightedClassImage(self, fileName):
        print "Rendering weighted classification image"
        
        weighted_nodes = scipy.array([[0 for x in range(self.side)] for y in range(self.side)])
        max = 0
        for point in self.plots:
            weighted_nodes[point[1],point[0]] += 1
            if(max < weighted_nodes[point[1],point[0]]):
                max  = weighted_nodes[point[1],point[0]]

        max += 1
        res = 100
        if(max < res):
            res = max - 1
        max  = self.transColour(max)

        rainbow = Rainbow.rainbow(0, max, res, "gbr")
        
        try:
            img = Image.new("RGB", (self.side, self.side))
            for point in self.plots:
                img.putpixel((point[1],point[0]), rainbow.getColour(self.transColour(weighted_nodes[point[1],point[0]])))
            img = img.resize((self.side*10, self.side*10),Image.NEAREST)
            print "Saving image to:", fileName
            img.save(fileName)
        except:
            print sys.exc_info()[0]
            raise
        
    # make an image of raw classifications 
    def renderContigImage(self, fileName):
        print "Rendering contig validation image"
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
            print "Saving image to:", fileName
            img.save(fileName)
            img2.save(fileName+".png")
        except:
            print sys.exc_info()[0]
            raise
    
    def transColour(self, val):
        return 10 * log(val)
    
    # return a color for a given vector
    def getColor(self, target_FV):
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
        
    # Work out how we are going to change our vector into a color
    # returns a list of starts and lengths
    def getColorLookup(self):
        if(self.FV_size == 3):
            return [0,1,2,1]
        elif(self.FV_size == 4):
            return [0,1,2,2]
        elif(self.FV_size == 5):
            return [0,1,2,3]
        elif(self.FV_size == 6):
            return [0,2,4,2]
        elif(self.FV_size == 7):
            return [0,2,4,3]
        elif(self.FV_size == 8):
            return [0,2,4,4]
        elif(self.FV_size == 9):
            return [0,2,4,5]
        elif(self.FV_size == 10):
            return [0,3,6,4]
        elif(self.FV_size == 11):
            return [0,3,6,5]
        elif(self.FV_size == 12):
            return [0,3,6,6]
        elif(self.FV_size == 13):
            return [0,3,6,7]
        elif(self.FV_size == 14):
            return [0,4,8,6]
        elif(self.FV_size == 15):
            return [0,4,8,7]
        else:
            print "***ERROR: Max vector size of 15! Time to be a haxxor!"
            sys.exit(1)
    
#------------------------------------------------------------------------------
# FILE IO 
    
    # print final weights to disk
    # save the data
    def printWeights(self, fileName):
        print "Saving weights"
        weights_writer = csv.writer(open(fileName, 'wb'), delimiter='\t')
        for r in range(self.side):
            for c in range(self.side):
                weights_writer.writerow(self.nodes[r,c])

#------------------------------------------------------------------------------
# 

###############################################################################
# TEMPLATE SUBS
###############################################################################
#
# Entry point, parse command line args and call out to doWork
#
if __name__ == '__main__':
    # intialise the options parser
    parser = argparse.ArgumentParser(description='Self organising map engine for Messy.')

    # add options here:
    parser.add_argument("working_dir", default=".", help="A place to write all the files etc...")
    parser.add_argument("data_file", help="Full path to the data file")
    parser.add_argument("--classify_only", action="store_true", default=False, help="Do not repeat the training, assumes that there is a valid weights file which corresponds to the data file")
    
    datagroup = parser.add_argument_group('Data options')
    datagroup.add_argument("--no_header", action="store_true", default=False, help="Use this only if the data file has NO header [default: has header]")
    datagroup.add_argument("--delimiter", "-l", default="\t", help="Delimeter in the data file [default: \"\\t\"]")    
    
    alggroup = parser.add_argument_group('Algorithm options')    
    alggroup.add_argument("--grid_size", "-g", type=int, default=32, help="Size of the data grid [default: 32]")
    alggroup.add_argument("--reps", "-r", type=int, default=50, help="Number of reps to use on the training data [default: 50]")
    alggroup.add_argument("--learning_rate", type=float, default=0.05, help="The diffusion speed [default: 0.05]")
    
    classgroup = parser.add_argument_group('Classification and validation')
    classgroup.add_argument("--class_img", "-L", default="messom_class.png", help="Full path to the classification image file [default: \"messom_class.png\"]")  
    classgroup.add_argument("--contig_img", "-G", default="messom_contigs.png", help="Full path to the contig validation image file [default: \"messom_contigs.png\"]")  
    classgroup.add_argument("--refs", "-R", default="", help="Comma separated list of refs files created by cleanAssem.pl OR directory containing refs files")  
 
    iogroup = parser.add_argument_group('IO options')
    iogroup.add_argument("--weight_img", "-w", default="messom_weights.png", help="Full path to the weights image file [default: \"messom_weights.png\"]")
    iogroup.add_argument("--weight_csv", "-W", default="weights.txt", help="specify a name for the weights file [default: \"weights.txt\"]")  
    iogroup.add_argument("--animate", "-a", action="store_true", default=False, help="Create a separate image for each iteration [default: False]")


    # get and check options
    args = parser.parse_args()

    # 
    # do what we came here to do
    #
    doWork(args)


###############################################################################
###############################################################################
###############################################################################
###############################################################################

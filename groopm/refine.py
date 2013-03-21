#!/usr/bin/env python
###############################################################################
#                                                                             #
#    refine.py                                                                #
#                                                                             #
#    A collection of classes / methods used when bin refinment and expansion  #
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
__version__ = "0.2.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Alpha"

###############################################################################

from sys import exc_info, exit, stdout as sys_stdout
from Queue import Queue

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import plot,subplot,axis,stem,show,figure

from numpy import arange as numpy_arange, copy as np_copy, arange as np_arange, ravel as np_ravel, ones as np_ones, eye as np_eye, shape as np_shape, around as np_around, argmax as np_argmax, arccos as np_cos, dot as np_dot, sum as np_sum, abs as np_abs, amax as np_amax, amin as np_amin, append as np_append, arccos as np_arccos, argmin as np_argmin, argsort as np_argsort, array as np_array, ceil as np_ceil, concatenate as np_concatenate, delete as np_delete, log10 as np_log10, max as np_max, mean as np_mean, median as np_median, min as np_min, pi as np_pi, reshape as np_reshape, seterr as np_seterr, size as np_size, sort as np_sort, sqrt as np_sqrt, std as np_std, where as np_where, zeros as np_zeros, cos as np_cos, sin as np_sin
from numpy.linalg import norm as np_norm 
from numpy.random import shuffle as shuffle

from scipy.spatial import KDTree as kdt
from scipy.cluster.vq import kmeans,vq,whiten,kmeans2
from scipy.spatial.distance import cdist, pdist, squareform

# GroopM imports
from profileManager import ProfileManager
from binManager import BinManager
from bin import Bin
from ellipsoid import EllipsoidTool
from PCA import PCA, Center
import groopmTimekeeper as gtime
import groopmExceptions as ge

np_seterr(all='raise')     

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class RefineEngine:
    """Workhorse wrapper for bin refinement"""
    def __init__(self): pass
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

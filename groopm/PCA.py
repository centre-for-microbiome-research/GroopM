#!/usr/bin/env python
###############################################################################
#                                                                             #
#    PCA.py                                                                   #
#                                                                             #
#    Tools for performing PCA.                                                #
#                                                                             #
#    Original Author:                                                         #
#    denis-bz-py at t-online dot de                                           #
#                                                                             #
#    See                                                                      #
#           http://stackoverflow.com/questions/1730600/                       #
#             ...principal-component-analysis-in-python                       #
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
import numpy as np

""" a small class for Principal Component Analysis

Notes:
    PCA does not center or scale A; you usually want to first
        A -= A.mean(A, axis=0)
        A /= A.std(A, axis=0)
    with the little class Center or the like, below.

See also:
    http://en.wikipedia.org/wiki/Principal_component_analysis
    http://en.wikipedia.org/wiki/Singular_value_decomposition
    Press et al., Numerical Recipes (2 or 3 ed), SVD
    PCA micro-tutorial
    iris-pca .py .png

"""

class PCA:
    def __init__( self, A, fraction=0.80):
        assert 0 <= fraction <= 1
            # A = U . diag(d) . Vt, O( m n^2 ), lapack_lite --
        self.U, self.d, self.Vt = np.linalg.svd( A, full_matrices=False )
        # different versions of numpy can return U and Vt such that
        # U * Vt is constant but the signs may be switched. Gross...
        # numpy, you owe me one day buster!
        # force a check here...
        if self.Vt[0,0] < 0:
            self.Vt *= -1.
            self.U *= -1. 
        assert np.all( self.d[:-1] >= self.d[1:] )  # sorted
        self.eigen = self.d**2
        self.sumvariance = np.cumsum(self.eigen)
        try:
            self.sumvariance /= self.sumvariance[-1]
        except:
            print len(A), len(self.sumvariance), len(self.eigen)
            raise

        self.npc = np.searchsorted( self.sumvariance, fraction ) + 1
        while(self.npc == 1):   # prevents less than 2 pcs being found
            fraction *= 1.1
            self.npc = np.searchsorted( self.sumvariance, fraction ) + 1
        self.dinv = np.array([ 1/d if d > self.d[0] * 1e-6  else 0
                                for d in self.d ])

    def pc( self ):
        """ e.g. 1000 x 2 U[:, :npc] * d[:npc], to plot etc. """
        n = self.npc
        return self.U[:, :n] * self.d[:n]

    # These 1-line methods may not be worth the bother;
    # then use U d Vt directly --

    def vars_pc( self, x ):
        n = self.npc
        return self.d[:n] * np.dot( self.Vt[:n], x.T ).T  # 20 vars -> 2 principal

    def pc_vars( self, p ):
        n = self.npc
        return np.dot( self.Vt[:n].T, (self.dinv[:n] * p).T ) .T  # 2 PC -> 20 vars

    def pc_obs( self, p ):
        n = self.npc
        return np.dot( self.U[:, :n], p.T )  # 2 principal -> 1000 obs

    def obs_pc( self, obs ):
        n = self.npc
        return np.dot( self.U[:, :n].T, obs ) .T  # 1000 obs -> 2 principal

    def obs( self, x ):
        return self.pc_obs( self.vars_pc(x) )  # 20 vars -> 2 principal -> 1000 obs

    def vars( self, obs ):
        return self.pc_vars( self.obs_pc(obs) )  # 1000 obs -> 2 principal -> 20 vars

class Center:
    """ A -= A.mean() /= A.std(), inplace -- use A.copy() if need be
        uncenter(x) == original A . x
    """
        # mttiw
    def __init__( self, A, axis=0, scale=True, verbose=1 ):
        self.mean = A.mean(axis=axis)
        if verbose:
            print "Center -= A.mean:", self.mean
        A -= self.mean
        if scale:
            std = A.std(axis=axis)
            self.std = np.where( std, std, 1. )
            if verbose:
                print "Center /= A.std:", self.std
            A /= self.std
        else:
            self.std = np.ones( A.shape[-1] )
        self.A = A

    def uncenter( self, x ):
        return np.dot( self.A, x * self.std ) + np.dot( x, self.mean )

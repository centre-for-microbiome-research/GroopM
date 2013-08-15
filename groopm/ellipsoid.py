#!/usr/bin/python
###############################################################################
#                                                                             #
#    ellipsoid.py                                                             #
#                                                                             #
#    Playing with ellipses!                                                   #
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
__copyright__ = "Copyright 2012/2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.1.2"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Beta"

###############################################################################

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg

np.seterr(all='raise')

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class EllipsoidTool:
    """Some stuff for playing with ellipsoids"""
    def __init__(self): pass

    def getMinVolEllipse(self, P, tolerance=0.01, retA=False):
        """ Find the minimum volume ellipsoid which holds all the points

        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!

        Here, P is a numpy array of 3D points like this:
        P = [[x,y,z],
             [x,y,z],
             [x,y,z]]

        Returns:
        (center, radii, rotation)

        """
        (N, d) = np.shape(P)

        # Q will be out working array
        Q = np.copy(P.T)
        Q = np.vstack([Q, np.ones(N)])
        QT = Q.T

        # initializations
        err = 1 + tolerance
        u = np.array([1.0 / N for i in range(N)]) # first iteration

        # Khachiyan Algorithm
        singular = False
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            try:
                M = np.diag(np.dot(QT , np.dot(linalg.inv(V), Q)))    # M the diagonal vector of an NxN matrix
            except linalg.linalg.LinAlgError:
                # most likely a singular matrix
                # permute the values a little and then we'll try again
                from random import random, randint
                PP = np.copy(P)
                for i in range(N):
                    if randint(0,3) == 0:
                        j = randint(0,2)
                        if randint(0,1) != 0:
                            PP[i,j] += random()
                        else:
                            PP[i,j] -= random()
                (A, center, radii, rotation) = self.getMinVolEllipse(PP, retA=True)
                singular = True
                break

            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u

        if not singular:
            # center of the ellipse
            center = np.dot(P.T, u)

            # the A matrix for the ellipse
            try:
                A = linalg.inv(
                               np.dot(P.T, np.dot(np.diag(u), P)) -
                               np.array([[a * b for b in center] for a in center])
                               ) / d
            except linalg.linalg.LinAlgError:
                # the matrix is singular so we need to return a degenerate ellipse
                #print '[Notice] Degenerate ellipse constructed indicating a bin with extremely small coverage divergence.'
                center = np.mean(P, axis=0)
                radii = np.max(P,axis=0) - np.min(P, axis=0)

                if len(P[0]) == 3:
                    rotation = [[0,0,0],[0,0,0],[0,0,0]]
                else:
                    rotation = [[0,0],[0,0]]

                if retA:
                    return (None, center, radii, rotation)
                else:
                    return (center, radii, rotation)
            
            # Get the values we'd like to return
            try:
                U, s, rotation = linalg.svd(A)
                radii = 1.0/np.sqrt(s)
            except np.linalg.linalg.LinAlgError:
                # hack -> better than crashing...
                rotation = np.eye(3)
                radii = np.ones(3)
        else:
            # hack -> better than crashing...
            rotation = np.eye(3)
            radii = np.ones(3)            
        if retA:
            return (A, center, radii, rotation)
        else:
            return (center, radii, rotation)

    def getEllipsoidVolume(self, radii):
        """Calculate the volume of the blob"""
        if len(radii) == 2:
            return np.pi*radii[0]*radii[1]
        else:
            return (4.0/3.0)*np.pi*radii[0]*radii[1]*radii[2]

    def doesIntersect3D(self, A, cA, B, cB):
        """Rough test to see if ellipsoids A and B intersect

        Not perfect, should work for "well overlapping" ones
        We assume that the volume of B is less than (or =) volume of A
        """
        #To make things simple, we just check if the points on a wire frame of
        #B lie within A

        # Quick check if the centre of B is within ellipse A. This deals with
        # degenerate cases where B is only a single point or an otherwise
        # degenerate ellipse.
        p_c = cB - cA
        try:
            if np.dot(p_c.T, np.dot(A, p_c)) <= 1:
                return True
        except TypeError:
                return False

        if A == None or B == None: # degenerate ellipse that can't be processed
            return False

        U, s, rotation = linalg.svd(B)
        try:
            radii_B = 1.0/np.sqrt(s)
        except FloatingPointError:
            # the given matrix B was made on a group of only one point
            # we need only check if the one point (the center)
            # in in A
            p_c = cB - cA
            return np.dot(p_c.T, np.dot(A, p_c)) <= 1

        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)

        # cartesian coordinates that correspond to the spherical angles:
        x = radii_B[0] * np.outer(np.cos(u), np.sin(v))
        y = radii_B[1] * np.outer(np.sin(u), np.sin(v))
        z = radii_B[2] * np.outer(np.ones_like(u), np.cos(v))

        # rotate accordingly
        for i in range(len(x)):
            for j in range(len(x)):
                # make a point on the wireFrame
                wire_point = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + cB
                # test if it's inside
                # work out (p-c)'A(p-c) and see if it's <= 1
                p_c = wire_point - cA
                if np.dot(p_c.T, np.dot(A, p_c)) <= 1:
                    return True

        return False

    def doesIntersect2D(self, A, cA, B, cB):
        """Rough test to see if ellipsoids A and B intersect

        Not perfect, should work for "well overlapping" ones
        We assume that the volume of B is less than (or =) volume of A
        """
        #To make things simple, we just check if the points on a wire frame of
        #B lie within A

        # Quick check if the centre of B is within ellipse A. This deals with
        # degenerate cases where B is only a single point or an otherwise
        # degenerate ellipse.
        p_c = cB - cA
        if np.dot(p_c.T, np.dot(A, p_c)) <= 1:
            return True

        if A == None or B == None:  # degenerate ellipse that can't be processed
            return False

        U, s, rotation = linalg.svd(B)
        try:
            radii_B = 1.0/np.sqrt(s)
        except FloatingPointError:
            # the given matrix B was made on a group of only one point
            # we need only check if the one point (the center)
            # in in A
            p_c = cB - cA
            return np.dot(p_c.T, np.dot(A, p_c)) <= 1

        u = np.linspace(0.0, 2.0 * np.pi, 100)

        # cartesian coordinates that correspond to the spherical angles:
        x = radii_B[0] * np.cos(u)
        y = radii_B[1] * np.sin(u)
        # rotate accordingly
        for i in range(len(x)):
            # make a point on the wireFrame
            edge_point = np.dot([x[i],y[i]], rotation) + cB
            # test if it's inside
            # work out (p-c)'A(p-c) and see if it's <= 1
            p_c = edge_point - cA
            if np.dot(p_c.T, np.dot(A, p_c)) <= 1:
                return True
        return False

    def plotEllipsoid(self, center, radii, rotation, ax=None, plotAxes=False, cageColor='b', cageAlpha=0.2, label=None):
        """Plot an ellipsoid"""
        make_ax = ax == None
        if make_ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)

        # cartesian coordinates that correspond to the spherical angles:
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
        # rotate accordingly
        for i in range(len(x)):
            for j in range(len(x)):
                [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center

        if plotAxes:
            # make some purdy axes
            axes = np.array([[radii[0],0.0,0.0],
                             [0.0,radii[1],0.0],
                             [0.0,0.0,radii[2]]])
            # rotate accordingly
            for i in range(len(axes)):
                axes[i] = np.dot(axes[i], rotation)


            # plot axes
            for p in axes:
                X3 = np.linspace(-p[0], p[0], 100) + center[0]
                Y3 = np.linspace(-p[1], p[1], 100) + center[1]
                Z3 = np.linspace(-p[2], p[2], 100) + center[2]
                ax.plot(X3, Y3, Z3, color=cageColor)

        # plot ellipsoid
        ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=cageColor, alpha=cageAlpha)

        if label is not None:
            ax.text(center[0],
                    center[1],
                    center[2],
                    label,
                    color=[0,0,0],
                    weight='bold'
                    )

        if make_ax:
            plt.show()
            plt.close(fig)
            del fig

    def plotEllipse(self, center, radii, rotation, ax=None, plotAxes=False, cageColor='b', cageAlpha=0.2, label=None, linewidth=-1):
        """plot an ellipse"""
        make_ax = ax == None
        if make_ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        u = np.linspace(0.0, 2.0 * np.pi, 100)

        # cartesian coordinates that correspond to the spherical angles:
        x = radii[0] * np.cos(u)
        y = radii[1] * np.sin(u)

        # rotate accordingly
        for i in range(len(x)):
            [x[i],y[i]] = np.dot([x[i],y[i]], rotation) + center

        if plotAxes:
            # make some purdy axes
            axes = np.array([[radii[0],0.0],[0.0,radii[1]]])
            # rotate accordingly
            for i in range(len(axes)):
                axes[i] = np.dot(axes[i], rotation)

            # plot axes
            for p in axes:
                X3 = np.linspace(-p[0], p[0], 100) + center[0]
                Y3 = np.linspace(-p[1], p[1], 100) + center[1]
                ax.plot(X3, Y3, color=cageColor)

        # plot ellipsoid
        if linewidth == -1:
            ax.plot(x, y, color=cageColor, alpha=cageAlpha)
        else:
            ax.plot(x, y, color=cageColor, alpha=cageAlpha, linewidth=linewidth, zorder = 10)

        if label is not None:
            ax.text(center[0],
                    center[1],
                    label,
                    color=[0,0,0],
                    weight='bold'
                    )

        if make_ax:
            plt.show()
            plt.close(fig)
            del fig

###############################################################################
###############################################################################
###############################################################################
###############################################################################

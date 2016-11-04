#!/usr/bin/env python

import os
import sys
import warnings
import numpy as np

from elements import ELEMENTS

warnings.filterwarnings("error")
verbosity = False

# Dictionary for energy conversion


energy_conversion = {'au'       : 1,
                     'eV'       : 27.21138505,
                     'wn'       : 219474.63,
                     'kj/mol'   : 2625.5,
                     'kcal/mol' : 627.503}

# au to eV = 27.21128505 / 1
# eV to au = 1 / 27.21128505
# u1 to u2 = dict[u2]/dict[u1]


# Dictionaries comparison:
# Options are stored in dictionaries. dictA is the default options dictionary
# and dictB is the one passed to the function.
# We want to compare dictA and dictB. We want to add to dictB all the missing 
# keys in dictA with their value.

def checkfile(filename):

    if not os.path.isfile(filename):
        print(banner(text='ERROR', ch='#', length=80))
        print("File %s not found!" % filename)
        sys.exit()


def skiplines(openfile, nlines=0):
    '''Skips nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.'''

    for i in range(nlines):
        next(openfile)

    return next(openfile)


def refframe(A, B, C):
    '''Returns a reference frame where the x axis goes from A to B, the y axis
    passes through C and the z axis is built accordingly.'''

    x = (B - A) / np.linalg.norm(B - A)

    # Define the point P on x whose perpendicular to x passes through C
    P = A + np.dot((C - A), x) * x
    y = (C - P) / np.linalg.norm(C - P)

    z = np.cross(x, y)

    ref = np.array([x, y, z])

    return ref


def symm_mat(M):
    '''Symmetrize an upper- or lower diagonal matrix.'''
    return M + M.T - np.diag(M.diagonal())


def v1v2_angle(v1, v2):
    '''Returns the angle between two vectors.'''
    # Remember that the angle between a plane and a vector equals
    # 90 - alpha, where alpha is the angle between the vector and
    # the normal to the plane. To obtain such an angle, you could
    # do angle = 90 - v1v2_angle(v1, np.cross(x, y)), where x, y
    # are the two vectors that define the plane.

    dotprod = np.dot(v1, v2)
    try:
        theta = np.degrees(np.arccos(dotprod / (np.linalg.norm(v1) * np.linalg.norm(v2))))
    except:
        theta = 0.0

    return theta


def kabsch(struct1, struct2, matrix=False):
    '''Returns the RMSD calculated with Kabsch's algorithm.'''

    # Modify structures to get rid of the atomic symbol or number and convert
    # to np.array
    struct1 = np.array([ [atom[0], atom[1], atom[2]] for atom in struct1 ])    
    struct2 = np.array([ [atom[0], atom[1], atom[2]] for atom in struct2 ])    

    # check for consistency in number of atoms
    try:
        assert len(struct1) == len(struct2)
        L = len(struct1)
        assert L > 0

    except AssertionError:
        print(" The selected coordinates sets have a different number of atoms!")
        sys.exit()

    # Center the two fragments to their center of coordinates
    com1 = np.sum(struct1, axis=0) / float(L)
    com2 = np.sum(struct2, axis=0) / float(L)
    struct1 -= com1
    struct2 -= com2

    # Initial residual, see Kabsch.
    E0 = np.sum(np.sum(struct1 * struct1, axis=0), axis=0) + \
         np.sum(np.sum(struct2 * struct2, axis=0), axis=0)

    # This beautiful step provides the answer. V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    V, S, Wt = np.linalg.svd(np.dot(np.transpose(struct2), struct1))

    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation. V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))

    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))

    # The rotation matrix U is simply V*Wt
    U = np.dot(V, Wt)
 
    # rotate and translate the molecule
    struct2 = np.dot((struct2), U)
    struct2 = struct2 + com1

    if matrix:
        return U, com2, com1

    else:
        return struct2, U


if __name__ == '__main__':
    pass

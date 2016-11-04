#!/usr/bin/env python

import os
import re
import sys
import numpy as np
from chrom import *
from ..util import util as u
from ..util.elements import ELEMENTS

#
# XYZ file
#
class XYZ(Chrom):
    '''A class to parse xyz files.'''

    def __init__(self, outfile):

        # Inherit from Chrom
        Chrom.__init__(self, outfile)
        # Fill in attributes
        self.type = self.__class__.__name__
        self.get_data()


    def get_data(self):
        '''Fills in class attributes with the data from td logfile.'''

        #
        # Molecular Properties
        #
        parsed_data = self.extract()

        # Geometrical properties
        self.Z_atoms = parsed_data[0]
        self.atoms = parsed_data[1]
        self.coords = parsed_data[2]
        self.natoms = len(self.atoms)
        self.com = self.calc_com(self.atoms, self.coords)

        # Electronic properties
        self.energies = None
        self.ntran = None
        self.f_osc = None
        self.r_vel = None
        self.r_len = None
        self.mu_vel = None
        self.mu_len = None
        self.mag = None

        pass


    def extract(self):
        '''Parses XYZ file for geometric properties.'''

        with open(self.file) as f:

            structure = []
            line = u.skiplines(f, 2)

            data = line.split()

            while len(data) == 4:

                try:
                    atom = int(data[0])
                except ValueError:
                    atom = data[0]

                atom_x = float(data[1])
                atom_y = float(data[2])
                atom_z = float(data[3])
                structure.append([atom, atom_x, atom_y, atom_z])

                try:
                    data = next(f).split()
                except:
                    break


            try:
                Z_atoms = [ int(x[0]) for x in structure ]
                atoms = [ ELEMENTS[x[0]].symbol for x in structure ]
            except ValueError:
                Z_atoms = [ ELEMENTS[x[0]].number for x in structure ]
                atoms = [ x[0] for x in structure ]

            coords = np.array([ x[1:] for x in structure ])

            return Z_atoms, atoms, coords


if __name__ == '__main__':
    pass

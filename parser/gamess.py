#!/usr/bin/env python

import os
import re
import sys
import numpy as np
from chrom import *
from ..util import util as u
from ..util.elements import ELEMENTS

au2ang = 0.529177249

#
# GAMESS US output file
#
class Gamess(Chrom):
    '''A class to parse GAMESS US td logfiles.'''

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
        self.energies = parsed_data[3]
        self.ntran = len(self.energies)
        self.f_osc = parsed_data[4]
        self.mu_len = parsed_data[5]
        self.r_vel = None
        self.r_len = None
        self.mu_vel = None
        self.mag = None
        self.trchgs = np.zeros((self.natoms, self.ntran))

        pass


    def extract(self):
        '''Parses Q-Chem td logfile for geometric and electronic properties.'''

        with open(self.file) as f:

            structure = []
            energies = []
            oscillators = []
            mu_len = []

            for line in f:

                #
                # Geometry
                #
                if "COORDINATES (BOHR)" in line:
                    line = u.skiplines(f, 1)
                    data = line.split()
                    
                    while len(data) == 5:

                        atom = data[0]
                        atom_x = float(data[2]) * au2ang
                        atom_y = float(data[3]) * au2ang
                        atom_z = float(data[4]) * au2ang
                        structure.append([atom, atom_x, atom_y, atom_z])

                        data = next(f).split()

                #
                # Electronic Properties
                #
                if "SUMMARY OF TDDFT RESULTS" in line:

                    line = u.skiplines(f, 4)
                    data = line.split()

                    while len(data) == 8:

                        #
                        # Excitation Energy
                        #
                        energy = float(data[3])
                        energies.append(energy)

                        #
                        # Transition Electric Dipole Moment
                        #
                        mu_x = float(data[4])
                        mu_y = float(data[5])
                        mu_z = float(data[6])
                        mu_len.append([mu_x, mu_y, mu_z])

                        #
                        # Oscillator Strength
                        #
                        oscillator = float(data[-1])
                        oscillators.append(oscillator)

                        # Update variables for while cycle
                        line = u.skiplines(f)
                        data = line.split()

            Z_atoms = [ ELEMENTS[x[0]].number for x in structure ]
            atoms = [ x[0] for x in structure ]
            coords = np.array([ x[1:] for x in structure ])
            energies = np.array(energies)
            oscillators = np.array(oscillators)
            mu_len = np.array(mu_len)

            return Z_atoms, atoms, coords, energies, oscillators, mu_len


if __name__ == '__main__':
    pass

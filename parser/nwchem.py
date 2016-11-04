#!/usr/bin/env python

import os
import re
import sys
import numpy as np
from chrom import *
from ..util import util as u
from ..util.elements import ELEMENTS

#
# NWChem output file
#
class NWChem(Chrom):
    '''A class to parse NWChem td logfiles.'''

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
        '''Parses G09 td logfile for geometric and electronic properties.'''

        with open(self.file) as f:

            structure = []
            energies = []
            oscillators = []
            mu_len = []

            for line in f:

                #
                # Reoriented geometry, overwrite previous
                #
                if "XYZ format geometry" in line:
                    structure = []
                    line = u.skiplines(f, 3)
                    data = line.split()

                    while len(data) == 4:

                        atom = data[0]
                        atom_x = float(data[1])
                        atom_y = float(data[2])
                        atom_z = float(data[3])
                        structure.append([atom, atom_x, atom_y, atom_z])

                        data = next(f).split()

                #
                # Excitation Energies
                #
                pattern = re.compile('\s+Root\s+\d+\s+singlet')
                if pattern.search(line):

                    energy = float(line.split()[-2])
                    energies.append(energy)
                    line = u.skiplines(f, 1)

                    while '-------------' not in line:

                        #
                        # Transition Electric Dipole Moment
                        #
                        if 'Transition Moments    X' in line:

                            mu_x = float(line.split()[3])
                            mu_y = float(line.split()[5])
                            mu_z = float(line.split()[7])
                            mu_len.append([mu_x, mu_y, mu_z])

                        #
                        # Oscillator Strength
                        #
                        if 'Dipole Oscillator Strength' in line:
                            oscillators.append(float(line.split()[-1]))

                        # Update line value for while cycle
                        line = u.skiplines(f)

            Z_atoms = [ ELEMENTS[x[0]].number for x in structure ]
            atoms = [ x[0] for x in structure ]
            coords = np.array([ x[1:] for x in structure ])
            energies = np.array(energies)
            oscillators = np.array(oscillators)
            mu_len = np.array(mu_len)

            return Z_atoms, atoms, coords, energies, oscillators, mu_len


if __name__ == '__main__':

    a = NWChem("../data/NWChem/nwchem_td.out")
    a.save_es_analysis(2, 25, 31)
    a.save_visdip(center=1)
    pass

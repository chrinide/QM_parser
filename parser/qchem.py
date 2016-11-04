#!/usr/bin/env python

import os
import re
import sys
import numpy as np
from chrom import *
from ..util import util as u
from ..util.elements import ELEMENTS

#
# Q-Chem output file
#
class QChem(Chrom):
    '''A class to parse Q-Chem td logfiles.'''

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
                if "Standard Nuclear Orientation" in line:
                    line = u.skiplines(f, 2)
                    data = line.split()

                    while len(data) == 5:

                        atom = data[1]
                        atom_x = float(data[2])
                        atom_y = float(data[3])
                        atom_z = float(data[4])
                        structure.append([atom, atom_x, atom_y, atom_z])

                        data = next(f).split()

                #
                # Electronic Properties
                #
                if "Excitation Energies" in line:

                    # Set empty lists, so that if, for example, both TDA
                    # and TDDFT are performed, only the last one will be
                    # stored
                    energies = []
                    oscillators = []
                    mu_len = []

                    line = u.skiplines(f, 2)

                    while '---------------------------------------------------' not in line:

                        #
                        # Excitation Energies
                        #
                        if 'excitation energy' in line:
                            energy = float(line.split()[-1])
                            energies.append(energy)
                      
                        #
                        # Discard Triplet States
                        #
                        if 'Multiplicity' in line:

                            if line.split()[-1] == 'Triplet':
                                energies = energies[:-1]
                                line = u.skiplines(f, 5)
                                continue

                        #
                        # Transition Electric Dipole Moment
                        #
                        if 'Trans. Mom.' in line:
                            mu_x = float(line.split()[2])
                            mu_y = float(line.split()[4])
                            mu_z = float(line.split()[6])
                            mu_len.append([mu_x, mu_y, mu_z])

                        #
                        # Oscillator Strengths
                        #
                        if 'Strength' in line:
                            oscillator = float(line.split()[-1])
                            oscillators.append(oscillator)

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

    a = QChem("../data/QChem/qchem_td.out")
    a.save_es_analysis(2,4,3)
    a.save_visdip(center=[2,4])


    b = QChem("../data/QChem/qchem_td2.out")
    b.save_es_analysis(4, 1, [5,6])
    b.save_visdip(center=[1,4])

    c = QChem("opt_qmmm.out")
    c.save_xyz()
    pass

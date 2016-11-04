#!/usr/bin/env python

import os
import re
import sys
import numpy as np
from chrom import *
from ..util import util as u
from ..util.elements import ELEMENTS

wn2eV = u.energy_conversion['eV'] / u.energy_conversion['wn']

#
# ORCA output file
#
class Orca(Chrom):
    '''A class to parse Orca td logfiles.'''

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
        self.r_vel = parsed_data[5]
        self.r_len = parsed_data[6]
        self.mu_vel = parsed_data[7]
        self.mu_len = parsed_data[8]
        self.mag = parsed_data[9]
        self.trchgs = np.zeros((self.natoms, self.ntran))

        pass


    def extract(self):
        '''Parses G09 td logfile for geometric and electronic properties.'''

        with open(self.file) as f:

            structure = []
            energies = []
            oscillators = []
            rot_vel = []
            rot_len = []
            mu_vel = []
            mu_len = []
            mag = []

            for line in f:

                #
                # Geometry
                #
                if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                    structure = []
                    line = u.skiplines(f, 1)
                    data = line.split()

                    while len(data) == 4:

                        atom = data[0]
                        atom_x = float(data[1])
                        atom_y = float(data[2])
                        atom_z = float(data[3])
                        structure.append([atom, atom_x, atom_y, atom_z])

                        data = next(f).split()

                #
                # Electronic Properties
                #
                if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:

                    line = u.skiplines(f, 4)
                    data = line.split()

                    while len(data) == 8:

                        #
                        # Excitation Energies
                        #
                        energy = float(data[1]) * wn2eV
                        energies.append(energy)

                        #
                        # Oscillator Strengths
                        #
                        oscillator = float(data[3])
                        oscillators.append(oscillator)

                        #
                        # Length transition electric dipole moment
                        #
                        mu_x = float(data[-3])
                        mu_y = float(data[-2])
                        mu_z = float(data[-1])
                        mu_len.append([mu_x, mu_y, mu_z])

                        data = next(f).split()

                #
                # Electronic Properties
                #
                if "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" in line:

                    line = u.skiplines(f, 4)
                    data = line.split()

                    while len(data) == 8:

                        #
                        # Excitation Energies and Oscillator Strengths are
                        # the same as in the already parsed table
                        # Velocity transition electric dipole moment
                        #
                        mu_x = float(data[-3])
                        mu_y = float(data[-2])
                        mu_z = float(data[-1])
                        mu_vel.append([mu_x, mu_y, mu_z])

                        data = next(f).split()

                #
                # Electronic Properties
                #
                if "CD SPECTRUM" in line:

                    line = u.skiplines(f, 4)
                    data = line.split()

                    while len(data) == 7:

                        #
                        # Rotatory Strength, presumably in Length formalism
                        #
                        r = float(data[3])
                        rot_len.append(r)

                        #
                        # Transition Magnetic Dipole Moment
                        #
                        m_x = float(data[-3])
                        m_y = float(data[-2])
                        m_z = float(data[-1])
                        mag.append([m_x, m_y, m_z])

                        data = next(f).split()

            Z_atoms = [ ELEMENTS[x[0]].number for x in structure ]
            atoms = [ x[0] for x in structure ]
            coords = np.array([ x[1:] for x in structure ])
            energies = np.array(energies)
            oscillators = np.array(oscillators)

            # To do: calculate rot_vel from the transition dipole moments
            # parsed. Look for the formula in the literature and if possible
            # compare the result with G09.
            rot_vel = np.array(rot_vel)
            rot_vel = None # TEMPORARY
            rot_len = np.array(rot_len)
            mu_vel = np.array(mu_vel)
            mu_len = np.array(mu_len)
            mag = np.array(mag)

            return Z_atoms, atoms, coords, energies, oscillators, rot_vel, rot_len, mu_vel, mu_len, mag


if __name__ == '__main__':

    a = Orca("../data/ORCA/orca_td.out")
    a.save_es_analysis(1, 4, [2, 3])
    a.save_visdip(center=[1,4])
    pass

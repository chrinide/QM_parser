#!/usr/bin/env python

import os
import re
import sys
import numpy as np
from chrom import *
from ..util import util as u
from ..util.elements import ELEMENTS

#
# G03 output file
#
class G03(Chrom):
    '''A class to parse G03 td logfiles.'''

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
            struct_done = False
            energies = []
            oscillators = []
            rot_vel = []
            rot_len = []
            mu_vel = []
            mu_len = []
            mag = []

            for line in f:

                #
                # Non reoriented geometry
                #
                if "Input orientation" in line and not struct_done:
                    line = u.skiplines(f, 4)
                    data = line.split()

                    while len(data) == 6:

                        Z_atom = int(data[1])
                        atom_x = float(data[3])
                        atom_y = float(data[4])
                        atom_z = float(data[5])
                        structure.append([Z_atom, atom_x, atom_y, atom_z])

                        data = next(f).split()

                #
                # Reoriented geometry, overwrite previous
                #
                if "Standard orientation" in line:
                    structure = []
                    line = u.skiplines(f, 4)
                    data = line.split()

                    while len(data) == 6:

                        Z_atom = int(data[1])
                        atom_x = float(data[3])
                        atom_y = float(data[4])
                        atom_z = float(data[5])
                        structure.append([Z_atom, atom_x, atom_y, atom_z])

                        data = next(f).split()
                        struct_done = True

                #
                # Length transition electric dipole
                #
                if "excited state Transition electric dipole" in line:

                    line = u.skiplines(f, 1)
                    data = line.split()
                    while len(data) == 5:
                        mu_len_x = float(data[1])
                        mu_len_y = float(data[2])
                        mu_len_z = float(data[3])
                        mu_len.append([mu_len_x, mu_len_y, mu_len_z])

                        line = next(f)
                        data = line.split()

                #
                # Velocity transition electric dipole
                #
                if "excited state transition velocity dipole" in line:

                    line = u.skiplines(f, 1)
                    data = line.split()
                    while len(data) == 5:
                        mu_vel_x = float(data[1])
                        mu_vel_y = float(data[2])
                        mu_vel_z = float(data[3])
                        mu_vel.append([mu_vel_x, mu_vel_y, mu_vel_z])

                        line = next(f)
                        data = line.split()

                #
                # Transition magnetic dipole
                #
                if "excited state transition magnetic dipole" in line:

                    line = u.skiplines(f, 1)
                    data = line.split()
                    while len(data) == 4:
                        m_x = float(data[1])
                        m_y = float(data[2])
                        m_z = float(data[3])
                        mag.append([m_x, m_y, m_z])

                        data = next(f).split()

                #
                # Velocity Rotatory Strength
                #
                if "R(velocity)" in line:

                    line = u.skiplines(f)
                    data = line.split()
                    while len(data) == 5:
                        r_vel = data[-1]
                        rot_vel.append(float(r_vel))

                        data = next(f).split()

                #
                # Length Rotatory Strength
                #
                if "R(length)" in line:

                    line = u.skiplines(f)
                    data = line.split()
                    while len(data) == 5:
                        r_len = data[-1]
                        rot_len.append(float(r_len))

                        data = next(f).split()

                #
                # Transition Energies and Oscillator Strengths
                #
                if "Excited State  " in line:

                    energy = line.split()[4]
                    energies.append(float(energy))

                    f_osc = line.split()[-1].split('=')[-1]
                    oscillators.append(float(f_osc))

            Z_atoms = [ int(x[0]) for x in structure ]
            atoms = [ ELEMENTS[x[0]].symbol for x in structure ]
            coords = np.array([ x[1:] for x in structure ])
            energies = np.array(energies)
            oscillators = np.array(oscillators)
            rot_vel = np.array(rot_vel)
            rot_len = np.array(rot_len)
            mu_vel = np.array(mu_vel)
            mu_len = np.array(mu_len)
            mag = np.array(mag)

            return Z_atoms, atoms, coords, energies, oscillators, rot_vel, rot_len, mu_vel, mu_len, mag


if __name__ == '__main__':
    pass

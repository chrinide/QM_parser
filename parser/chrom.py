#!/usr/bin/env python

import os
import re
import sys
import numpy as np
from ..util import util as u
from ..util.elements import ELEMENTS

#
# Chromophore Class
#
class Chrom:
    '''A class for chromophores.'''

    def __init__(self, outfile=None):

        if outfile:
            self.name = '.'.join(os.path.split(outfile)[1].split('.')[:-1])
            self.ext = os.path.split(outfile)[1].split('.')[-1]
            self.path = os.path.split(outfile)[0]
            self.file = os.path.join(self.path, '.'.join([self.name, self.ext]))

        else:
            self.path = os.getcwd()

        # Geometrical properties
        self.Z_atoms = None
        self.atoms = None
        self.coords = None
        self.natoms = None
        self.com = None
        self.type = "Empty"

        # Electronic properties
        self.energies = None
        self.ntran = None
        self.f_osc = None
        self.r_vel = None
        self.r_len = None
        self.mu_vel = None
        self.mu_len = None
        self.mag = None
        self.trchgs = None

        pass


    def filter_props(self, idxs=None):
        '''Filters class attributes concerning electronic properties according
        to array (or list) idxs, allowing for an easy selection of transitions
        with a numbered list.'''

        #
        # Adapt to python numeration
        #
        try:
            idxs = np.array(idxs) - 1
        
        except:
            pass

        #
        # Filter Excited States properties, handling the None objects
        #
        try:
            self.energies = self.energies[idxs]
        except:
            pass

        try:
            self.ntran = len(self.energies)
        except:
            pass

        try:
            self.f_osc = self.f_osc[idxs]
        except:
            pass

        try:
            self.r_vel = self.r_vel[idxs]
        except:
            pass

        try:
            self.r_len = self.r_len[idxs]
        except:
            pass

        try:
            self.mu_vel = self.mu_vel[idxs]
        except:
            pass

        try:
            self.mu_len = self.mu_len[idxs]
        except:
            pass

        try:
            self.mag = self.mag[idxs]
        except:
            pass
        
        try:
            self.trchgs = self.trchgs[:,idxs]
        except:
            pass

        return


    def set_energies(self, new_energies, idxs=None):
        '''Sets the energies for states identified by the list (or array) idxs
        to the values specified by the user in the list (or array) new_energies.'''

        if self.energies is None:
            self.energies = np.zeros(len(idxs))
            self.mu_len = np.zeros((len(idxs), 3))
            self.mag = np.zeros((len(idxs), 3))

        if idxs is not None:
            idxs = np.array(idxs) - 1
            self.energies[idxs] = np.array(new_energies)
        
        else:
            self.energies = np.array(new_energies)

        self.ntran = len(self.energies)

        return


    def scale_dips(self, factor=2):
        '''Scales Transition Dipole Moments by a factor.'''

        #
        # Scale Transition Dipole Moments, handling the None objects
        #
        try:
            self.mu_vel *= factor
        except:
            pass

        try:
            self.mu_len *= factor
        except:
            pass

        try:
            self.mag *= factor
        except:
            pass
        
        return


    def es_analysis(self, A, B, C, mu=None, mag=None):
        '''Returns angles theta and phi with respect to the reference frame
        defined by A, B and C. A, B, C can either be integers or lists. The
        coordinates of the atoms corresponding to the A, B, C indexes (or the
        average coordinates, if they are lists) will be used to build the frame
        so that the x axis passes through A and B, the y axis is orthogonal to
        x and passes through C. The z axis is built accordingly. Theta is the
        angle of a vector in the xy plane with respect to x, phi is the angle
        with the xy plane.'''

        # Assign default Transition Dipole Moments for the analysis
        if mu is None:
            mu = self.mu_len

        if mag is None:
            mag = self.mag

        #
        # Convert idxs to python numeration
        #
        A =  np.array(A) - 1
        B =  np.array(B) - 1
        C =  np.array(C) - 1

        #
        # Average the coordinates of atoms called with idxs A, B, C
        # This is an elegant solution in my opinion: the [A,None] notation
        # adds an empty axes to the array, so that also 1d arrays can be
        # averaged along the columns. Finally, the average value is reported
        # as a 1d vector.
        #
        A1 = np.average(self.coords[A,None], axis=0).reshape(3)
        B1 = np.average(self.coords[B,None], axis=0).reshape(3)
        C1 = np.average(self.coords[C,None], axis=0).reshape(3)

        #
        # Generate reference frame, see util module
        #
        self.ref = u.refframe(A1, B1, C1)
        ref = self.ref

        #
        # Calculate theta and phi angle for each transition
        #
        theta_mu = []
        theta_mag = []
        phi_mu=[]
        phi_mag=[]
        for i in range(self.ntran):

            if mu is not None:

                # theta
                # subtract from the vector its component along the normal
                # to the plane. For the sign, check the dot prod with the
                # y plane
                t = u.v1v2_angle(mu[i] - ref[2] * np.dot(mu[i], ref[2]), ref[0])

                # dipole reorientation if needed
                if np.abs(t) > 90:
                    
                    mu[i] *= -1
                    t = u.v1v2_angle(mu[i] - ref[2] * np.dot(mu[i], ref[2]), ref[0])

                sign = np.sign(np.dot(mu[i], ref[1]))
                t = t * sign
                theta_mu.append(t)

                # phi
                # angle between the vector and the normal to xy plane
                p = 90 - u.v1v2_angle(mu[i], ref[2])
                phi_mu.append(p)

            if mag is not None:

                # theta
                # subtract from the vector its component along the normal
                # to the plane. For the sign, check the dot prod with the
                # y plane
                t = u.v1v2_angle(mag[i] - ref[2] * np.dot(mag[i], ref[2]), ref[0])

                # dipole reorientation if needed
                if np.abs(t) > 90:

                    mag[i] *= -1
                    t = u.v1v2_angle(mag[i] - ref[2] * np.dot(mag[i], ref[2]), ref[0])

                sign = np.sign(np.dot(mag[i], ref[1]))
                t = t * sign
                theta_mag.append(t)

                # phi
                # angle between the vector and the normal to xy plane
                p = 90 - u.v1v2_angle(mag[i], ref[2])
                phi_mag.append(p)

        #
        # Convert everything to np.array
        #
        theta_mu = np.array(theta_mu)
        phi_mu = np.array(phi_mu)
        theta_mag = np.array(theta_mag)
        phi_mag = np.array(phi_mag)

        return theta_mu, phi_mu, theta_mag, phi_mag


    def save_es_analysis(self, A, B, C, mu=None, mag=None, filename=None):
        '''Saves the ES_analysis.txt file, containing the orientation of
        Transition Dipole Moments in terms of Theta and Phi angles, defined
        in the es_analysis() method.'''

        # Default filename
        if filename is None:
            filename = "ES_ana_%s" % self.name

        # Assign default Transition Dipole Moments for the analysis
        if mu is None:
            mu = self.mu_len

        if mag is None:
            mag = self.mag

        # Assign to filename the same path as the logfile connected to the instance
        filename = os.path.join(self.path, filename)
        theta_mu, phi_mu, theta_mag, phi_mag = self.es_analysis(A, B, C, mu=mu, mag=mag)

        header = ("#\n"
                  "# Analysis of Excited State Transition Dipole Moments Orientation\n"
                  "#\n"
                  "# Orientation defined by atoms:\n"
                  "# A : %s\n"
                  "# B : %s\n"
                  "# C : %s\n"
                  "#\n"
                  "# X axis : A -> B\n"
                  "# Y axis : x -> C\n"
                  "#\n"
                  "# Theta : angle between Transition Dipole Moment and X axis in the XY plane\n"
                  "# Phi : angle between Transition Dipole Moment and XY plane\n"
                  "#\n#\n#\n")

        intest = tuple(["Energy (eV)", "Norm (a.u.)", "Theta (deg)", "Phi (deg)"])
        with open("%s.txt" % filename, 'w') as f:
            f.write(header % (A, B, C))

            if mu is not None:

                f.write("#\n")
                f.write("# Transition Electric Dipole Moments\n")
                f.write("#\n")
                f.write("# %12s %12s %12s %12s\n" % intest)
                f.write("#\n")

                for i in range(self.ntran):

                    data = [self.energies[i], np.linalg.norm(mu[i]), theta_mu[i], phi_mu[i]]
                    f.write("%12.6f %12.6f %12.2f %12.2f\n" % tuple(data))

            if mag is not None:

                f.write("#\n")
                f.write("#\n")
                f.write("# Transition Magnetic Dipole Moments\n")
                f.write("#\n")
                f.write("# %12s %12s %12s %12s\n" % intest)
                f.write("#\n")

                for i in range(self.ntran):

                    data = [self.energies[i], np.linalg.norm(mag[i]), theta_mag[i], phi_mag[i]]
                    f.write("%12.6f %12.6f %12.2f %12.2f\n" % tuple(data))

        pass

    
    @staticmethod
    def calc_com(atoms, coords):
        '''Calculates the Center of Mass given the atom list and the
        coordinates array.'''

        #
        # Center of Mass
        #
        masses = np.array([ ELEMENTS[atom].mass for atom in atoms ])
        com = np.dot(coords.T, masses) / np.sum(masses)

        return com


    def set_center(self, idxs=None):
        '''Calculates the Center of the molecule according to the idxs in the
        coordinates array.'''

        if idxs is not None:
            idxs = np.array(idxs) - 1
            self.com = np.average(self.coords[idxs], axis=0)
        
        return


    def save_xyz(self, filename=None):
        '''Saves an .xyz file with the structure of the chromophore.'''

        # Default filename
        if filename is None:
            filename = "geo_%s" % self.name

        # Assign to filename the same path as the logfile connected to the instance
        filename = os.path.join(self.path, filename)
        with open("%s.xyz" % filename, 'w') as f:
            f.write("%d\n" % self.natoms)
            f.write("Saved by Chrom Class\n")
            for i in range(self.natoms):
                data = [self.atoms[i], self.coords[i,0], self.coords[i,1], self.coords[i,2]]
                f.write("%-3s %12.8f %12.8f %12.8f\n" % tuple(data))

        return


    def save_visdip(self, mu=None, mag=None, center=None, filename=None):
        '''Saves a .vmd script for the visualization of Transition Dipole
        Moments centered on center. If center is not specified, the COM will
        be used. Otherwise center can be an index, which will be associated to
        an atom, or a list of indexes, which will be associated to the
        geometrical center of the atoms identified by those indexes.'''

        # Default filename
        if filename is None:
            filename = "vis_%s" % self.name

        # Assign default Transition Dipole Moments for the analysis
        if mu is None:
            mu = self.mu_len

        if mag is None:
            mag = self.mag

        # Assign to filename the same path as the logfile connected to the instance
        xyzname = filename
        filename = os.path.join(self.path, filename)
        header = ("#\n"
                  "# VMD script to draw vectors\n"
                  "#\n"
                  "menu main on\n"
                  "display projection orthographic\n"
                  "display depthcue off\n"
                  "display nearclip set 0.01\n"
                  "axes location lowerleft\n"
                  "\n"
                  "#\n"
                  "# VMD functions to draw a vector\n"
                  "#\n"
                  "proc vmd_draw_arrow {mol start end} {\n"
                  "    set length [veclength [vecsub $end $start]]\n"
                  "    set conelen [expr max(0.4,0.2*$length) ]\n"
                  "    set scale [expr max(0.5,(1.0-$conelen/$length))]\n"
                  "\n"
                  "    set middle [vecadd $start [vecscale $scale [vecsub $end $start]]]\n"
                  "    graphics $mol cylinder $start $middle radius 0.05\n"
                  "    puts [list cone $middle $end radius 0.15]\n"
                  "    graphics $mol cone $middle $end radius 0.15\n"
                  "}\n"                     
                  "\n"                                                                   
                  "proc vmd_draw_vector { mol pos val } {\n"
                  "    set end   [ vecadd $pos [ vecscale +1 $val ] ]\n"
                  "    vmd_draw_arrow $mol $pos $end\n"
                  "}\n"
                  "\n"
                  "\n")

        #
        # If center is not specified, use COM. Otherwise, get atomic coordinates
        # of the specified atom or atom list and calculate the geometrical center
        #
        if not center:
            center = self.calc_com(self.atoms, self.coords)

        else:
            center = np.array(center) - 1
            center = np.average(self.coords[center,None], axis=0).reshape(3)

        #
        # Save the geometry and write the VMD script file with transition dips
        #
        self.save_xyz(filename=xyzname)
        with open('%s.vmd' % filename, 'w') as f:
            f.write(header)
            f.write("mol new %s.xyz type xyz\n" % os.path.split(xyzname)[-1])
            f.write("\n")
            f.write("\n")

            #
            # If an Excited State Analysis has been performed, save also
            # the reference axes used for the analysis
            #
            try:
                if self.ref is not None:
                    
                    f.write("#\n")
                    f.write("# Reference Frame chosen for the Analysis\n")
                    f.write("#\n")
                    f.write("# X, Y, Z axes\n")

                    for i in range(3):
                        cmd = "graphics 0 color blue; vmd_draw_vector 0 {%8.4f %8.4f %8.4f} {%8.4f %8.4f %8.4f}\n"
                        data = [center[0], center[1], center[2], self.ref[i,0], self.ref[i,1], self.ref[i,2]]
                        f.write(cmd % tuple(data))

                    f.write("\n")
                    f.write("\n")

            except:
                pass

            for i in range(self.ntran):

                if mu is not None:

                    comment = "# Transition Electric Dipole Moment of Excited State %d\n"
                    cmd = "graphics 0 color red; vmd_draw_vector 0 {%8.4f %8.4f %8.4f} {%8.4f %8.4f %8.4f}\n"
                    if i > 0:
                        cmd = "# " + cmd
                        comment = "# " + comment

                    data = [center[0], center[1], center[2], mu[i,0], mu[i,1], mu[i,2]]
                    f.write(comment % (i + 1))
                    f.write(cmd % tuple(data))
                    f.write("\n")

                if mag is not None:

                    comment = "# Transition Magnetic Dipole Moment of Excited State %d\n"
                    cmd = "graphics 0 color green; vmd_draw_vector 0 {%8.4f %8.4f %8.4f} {%8.4f %8.4f %8.4f}\n"
                    if i > 0:
                        cmd = "# " + cmd
                        comment = "# " + comment

                    data = [center[0], center[1], center[2], mag[i,0], mag[i,1], mag[i,2]]
                    f.write(comment % (i + 1))
                    f.write(cmd % tuple(data))
                    f.write("\n")
                pass
                  

        return


    def transform(self, coords):
        '''Transform coordinates and dipoles according to the transormation
        matrix that minimizes the RMSD between the chromophore and the array
        coords.'''

        self.coords, M = u.kabsch(coords, self.coords)
        self.com = self.calc_com(self.atoms, self.coords)
        
        try:
            self.mu_vel = np.dot(self.mu_vel, M)
        except TypeError:
            pass

        try:
            self.mu_len = np.dot(self.mu_len, M)
        except TypeError:
            pass

        try:
            self.mag = np.dot(self.mag, M)
        except TypeError:
            pass

        return


if __name__ == '__main__':
    pass

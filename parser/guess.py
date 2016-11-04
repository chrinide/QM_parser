#!/usr/bin/env python

import sys

from g03 import *
from g09 import *
from gamess import *
from nwchem import *
from orca import *
from qchem import *
from xyz import *

def guess(filename):
    '''Returns the correct class needed to parse filename, if it exists.'''

    #
    # Dictionary of unique sentences in QM packages output files to guess
    # the correct parser to use
    #
    filetypes = {}
    
    filetypes["This is the Gaussian(R) 03 program."] = G03
    filetypes["This is part of the Gaussian(R) 09 program."] = G09
    filetypes["GAMESS VERSION"] = Gamess
    filetypes["Northwest Computational Chemistry Package (NWChem)"] = NWChem
    filetypes["* O   R   C   A *"] = Orca
    filetypes["A Quantum Leap Into The Future Of Chemistry"] = QChem

    filetype = None
    done = False
    with open(filename) as f:

        for line in f:
            for sentence in filetypes.keys():

                if sentence in line:
                    filetype = filetypes[sentence]
                    done = True
                    break

            # once the type has been identified, exit the cycles
            if done:
                break

    if not filetype:
        try:
            XYZ(filename)
            filetype = XYZ
        except:
            pass

    if not filetype:
        print(" %s" % filename)
        print(" File type not known")
        sys.exit()

    return filetype(filename)


if __name__ == '__main__':

    testfile = sys.argv[1]
    parser = guess(testfile)

    data = parser(testfile)

    print data.energies
    pass

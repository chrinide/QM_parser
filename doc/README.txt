Module to parse Excited State properties from various QM packages output files.

In ../data there are folder related to some QM packages, which contain example
output files used for development.

In ../parser there are the python modules containing classes to parse those
example output files. Additionally, there is a general class for chromophores
in chrom.py and a function used to guess the correct class to use in guess.py

In ../util there are a couple of modules containing utility functions, such
as a dictionary of chemical elements (elements.py) and some linear algebra and
other functions in util.py

#!/usr/bin/env python2
#
# The MIT License (MIT)
#
# Copyright (c) 2015 Anders Steen Christensen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
# Citations:
#   (1) J Rezac and P Hobza (2012) J. Chem. Theory Comput., 8:141-151
#   (2) AS Christensen (2015) https://github.com/andersx/pauli-se
#

import sys
import math

def get_coordinates(filename):

    """ Opens a standard XYZ file and returns an array with all
        hydrogen coordinates.
    """

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    hydrogens = []

    for line in lines[2:]:
        if len(line) < 3:
            break

        tokens = line.split()

        if tokens[0] == "H":

            hydrogens.append([float(tokens[1]),
                              float(tokens[2]),
                              float(tokens[3])])

    return hydrogens


def calc_interaction(a, b):

    """ Calculates the Pauli repulsion between two hydrogen atoms
        using the parameters imported from parameters.py in kcal/mol.
    """

    r = math.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)
    e = SHH * (1.0 - 1.0 / (1.0 + math.exp(-1.0 * EHH * (r/R0HH - 1.0))))

    return e


def calc_pauli(filename, SHH, EHH, R0HH):

    # Get hydrogen positions from XYZ-file.
    h = get_coordinates(filename)

    # Reset energy.
    energy = 0.0

    # Loop over all unique pairs.
    for i, a in enumerate(h):
        for b in h[i+1:]:

            # Add to sum.
            energy += calc_interaction(a, b)

    return energy


if __name__ == "__main__":

    # Import some parameters.
    from parameters import SHH, EHH, R0HH

    # Get XYZ-filename.
    filename = sys.argv[1]

    # Calculate energy, defining filename and parameters.
    energy = calc_pauli(filename, SHH, EHH, R0HH)

    # Output results in kcal/mol.
    print "Pauli repulsion: %6.4f kcal/mol" % energy

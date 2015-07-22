#!/usr/bin/env python

"""
  Format Exomol's CH4 data (Yurchenko et al. 2014) for use in Transit.

  Usage:
  ------
  To run this code, execute from the Shell:
   ./Yurchenko_CH4_format.py fileIn1 [fileIn2] ... [fileInN]

  Parameters:
  -----------
  fileIn1: String
     File name of the input CH4 opacity cross section as given by:
     http://www.exomol.com/xsecs/12C-1H4
  fileIn2 -- fileInN: String
     Optional additional file names of input CH4 opacity cross section.

  Examples:
  ---------
  For H2-H2 run from the Shell prompt:
    ./Yurchenko_CH4_format.py  /home/.../data/opacities/CH4/Yurchenko/*.sigma

  Notes:
  ------
  - The ExoMol input files have to be generated including the column with
    the wavenumber.
"""

import sys
import numpy as np

# Inout/output files:
nfiles = len(sys.argv) - 1
filein = sys.argv[1:]

# Number of temperatures:
ntemp = nfiles
temp = np.zeros(ntemp, np.double)

# Read and extract data from files:
for i in np.arange(nfiles):
  f = open(filein[i], "r")
  lines = f.readlines()
  f.close()

  if i == 0:
    # Number of wavenumber samples:
    nwave = len(lines)
    wave = np.zeros(nwave, np.double)
    # Allocate output data array:
    data = np.zeros((nwave,ntemp), np.double)

  # Extract temperature from the filename:
  temp[i] = (filein[i].split("_")[2])[:-1]

  for j in np.arange(nwave):
    val = lines[j].split()
    # Get the wavenumber only if thi is the first file:
    if i == 0:
      wave[j] = val[0]
    # Get the opacity:
    data[j,i] = val[1]

species = "CH4"

# Write to the output file:
fileout = "ExoMol_{:s}_{:.1f}-{:.1f}cm-1_{:04d}-{:04d}K.dat".format(
                   species, wave[0], wave[-1], int(temp[0]), int(temp[-1]))
fout = open(fileout, "w")

fout.write("# This file formats the tabulated ExoMol {:s} data from "
           "Yurchenko et al. (2014).\n\n".format(species))

fout.write("i {:s}\n".format(species))
fout.write("t         ")
for j in np.arange(ntemp):
    fout.write("    {:6.1f}".format(temp[j]))
fout.write("\n\n")

fout.write("# Wavenumber in cm-1, Opacity cross section in cm-1 amagat-1:\n")

# Write down the data:
for i in np.arange(nwave):
  fout.write("  {:7.1f} ".format(wave[i]))
  for j in np.arange(ntemp):
    fout.write(" {:.3e}".format(data[i,j]))
  fout.write("\n")

fout.close()

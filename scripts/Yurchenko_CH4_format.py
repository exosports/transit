#!/usr/bin/env python

import sys
import numpy as np
import scipy.constants as sc

# Amagat (in molec cm^-3):
N0 = sc.physical_constants[
        "Loschmidt constant (273.15 K, 101.325 kPa)"][0] * 1e-6

def main():
  """
  Format Exomol's CH4 data (Yurchenko et al. 2014) for use in Transit.
  http://www.exomol.com/xsecs/12C-1H4

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
  $ ./Yurchenko_CH4_format.py  /home/.../data/opacities/CH4/Yurchenko/*.sigma

  Notes:
  ------
  - The ExoMol input files have to be generated including the column with
    the wavenumber.
  """

  # Inout/output files:
  filein = sys.argv[1:]
  # Number of temperatures (one temp sample per file):
  ntemp = nfiles = len(filein)

  # Array of sampled temperatures:
  temp = np.zeros(ntemp, np.double)

  # Read and extract data from files:
  for j in np.arange(nfiles):
    f = open(filein[j], "r")
    lines = f.readlines()
    f.close()

    if j == 0:
      # Number of wavenumber samples:
      nwave = len(lines)
      wave = np.zeros(nwave, np.double)
      # Allocate output data array:
      data = np.zeros((nwave,ntemp), np.double)

    # Extract temperature from the filename:
    temp[j] = (filein[j].split("_")[2])[:-1]

    for i in np.arange(nwave):
      val = lines[i].split()
      # Get the wavenumber only if thi is the first file:
      if j == 0:
        wave[i] = val[0]
      # Get the opacity:
      data[i,j] = val[1]

  # Species name is hardcoded (undo if this script works for other datasets):
  species = "CH4"

  # Convert units from cm2 molecule-1 to cm-1 amagat-1:
  data *= N0

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


if __name__ == "__main__":
  main()

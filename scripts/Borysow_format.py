#!/usr/bin/env python

# Format Borysow CIA data files for Transit:

# Download the CIA tabulated data from Borysow's webpage:
#   http://www.astro.ku.dk/~aborysow/programs/
# to the working directory.

# Note that final_CIA_HT.dat has an incomplete column near the end that
# has to be patched before proceeding.

import sys, os
import numpy as np

def main():
  """
  This script formats a CIA data files from A. Borysow into the CIA file
  read by the Transit program.

  Usage:
  ------
  Run from the shell:
  $ ./Borysow_format.py <InputFileName> <OutputFileName> <MolName1> <MolName2>

  Examples:
  --------
  $ ./Borysow_format.py ciah2he_dh_quantmech CIA_Borysow_H2He_1000-7000K_0.5-400um.dat H2 He
  $ ./Borysow_format.py final_CIA_LT.dat     CIA_Borysow_H2H2_0060-0350K_0.6-1000um.dat H2 H2
  $ ./Borysow_format.py final_CIA_HT.dat     CIA_Borysow_H2H2_0400-1000K_0.6-1000um.dat H2 H2
  $ ./Borysow_format.py CIA.H2H2.Yi          CIA_Borysow_H2H2_1000-7000K_0.6-0500um.dat H2 H2
  """
  CIAin, CIAout, mol1, mol2 = sys.argv[1:]
  root = ""
  # Read and extract data from files:
  f = open(root + CIAin, "r")
  lines = f.readlines()
  f.close()

  # Extract temperature arrays:
  T = lines[1].split()[1:]
  # Remove the trailing "K" from the temperatures:
  for i in np.arange(len(T)):
    T[i] = T[i][:-1]
  T = np.asarray(T, np.double)

  # Number ow wavelength samples:
  size = len(lines) - 3

  # Extract wavenumber and extinction arrays:
  wn   = np.zeros(size, np.double)
  data = np.zeros((size,len(T)), np.double)

  for i in np.arange(size):
    info = lines[i+3].split()
    wn[i]   = info[0]
    data[i] = info[1:]

  # Output file name:
  fout = open(CIAout, "w")

  fout.write(
        "# This file incorporates the tabulated {:s}-{:s} CIA data from:\n"
        "#  http://www.astro.ku.dk/~aborysow/programs/{:s}\n\n".
          format(mol1, mol2, CIAin))

  fout.write("i {:s} {:s}\n".format(mol1, mol2))
  fout.write("t        ")
  for j in np.arange(len(T)):
      fout.write("      {:4.0f}".format(T[j]))
  fout.write("\n\n")

  fout.write("# Wavenumber in cm-1, CIA coefficients in cm-1 amagat-2:\n")
  for i in np.arange(len(wn)):
    fout.write(" {:7.1f} ".format(wn[i]))
    for j in np.arange(len(T)):
      fout.write(" {:.3e}".format(data[i,j]))
    fout.write("\n")

  fout.close()


if __name__ == "__main__":
  main()

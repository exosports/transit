#!/usr/bin/env python

import numpy as np
import sys, os
import scipy.constants as sc

# Amagat (in molec cm^-3):
N0 = sc.physical_constants[
        "Loschmidt constant (273.15 K, 101.325 kPa)"][0] * 1e-6

def main():
  """
  This scripts reformats a HITRAN cross-section file into the
  format read by the Transit program. Designed for use with
  Hargreaves, et al (2014) cross section files.

  Usage
  -----
  Run from the shell:
  $ ./HITRAN_csx_format.py inputfile1 [... inputfileN] outputfile

  Examples
  --------
  $ ./HITRAN_csx_format.py CH4_*C_XSC.xsc CH4_cross-section_Hargreaves.dat
  """
  infile  = sys.argv[1:-1]
  outfile = sys.argv[-1]

  # There is one temperature sample per file:
  ntemp = nfiles = len(infile)
  # Array of sampled temperatures:
  temp = np.zeros(ntemp, np.double)

  # Number of data points per line:
  linelen = 10
  # Error message:
  error = ""

  for i in np.arange(nfiles):
    # Open file for reading:
    f_in = open(infile[i], 'r')

    line = str.split(f_in.readline())
    # Set the fields from the header:
    if i == 0:
      mol     = str(  line[0])  # Molecule designation (eg CH4)
      wn_init = float(line[1])  # Initial wavenumber (cm-1)
      wn_fin  = float(line[2])  # Final wavenumber (cm-1)
      nwave   = int(  line[3])  # Number of data points
      press   = float(line[5])  # Pressure (Torr)
      res     = float(line[7])  # Resolution of measurement (cm-1)

      # Array to be filled with cross section data
      data = np.zeros((nwave, ntemp), np.double)
    else:
      # Check that the data from the other files match:
      if str(line[0]) != mol:
        error += "for different species, "
      if (float(line[1]) != wn_init or
          float(line[2]) != wn_fin  ):
        error += "with different wavelength ranges, "

      if error != "":
        print("Can't combine files {:s}".format(error))
        return

    temp[i] = line[4]  # Add temperature (in K) to list

    # For counting loop repetitions
    j = 0

    # Loop through the lines of data
    for x in f_in:
      # Loop through fields in each line
      for k in range(linelen):
        # Only read data that is there (don't read extra fields at the end)
        if k+linelen*j < nwave:
          data[k+linelen*j, i] = str.split(x)[k]
      # Increment number of loop repetitions
      j += 1

    f_in.close()

  # Convert from cm^2 molecule-1 to cm-1 amagat-1 units:
  data *= N0

  # Wavenumber array to be printed in output file
  wn = np.linspace(wn_init, wn_fin, nwave)

  # Open file for writing:
  f_out = open(outfile, 'w')

  # Write the header
  f_out.write('# This file is a reformatted version of Hargreaves et al\n')
  f_out.write('# (2015) hot methane cross section file\n\n')

  # Write the molecule and temperature lines
  f_out.write('i {:s}\n'.format(mol))
  f_out.write('t ')
  for j in np.arange(ntemp):
    f_out.write(' {:6.1f}'.format(temp[j]))
  f_out.write('\n\n')

  # Write comment on units
  f_out.write('\n # Wavenumber in cm-1, cross section in cm-1 amagat-1\n')

  # Write data in format required by Transit program
  for i in range(nwave):
    f_out.write(" {:7.2f}   ".format(wn[i]))
    for j in np.arange(ntemp):
      f_out.write(" {:.3e}".format(data[i,j]))
    f_out.write('\n')

  f_out.close()


if __name__ == "__main__":
  main()

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
  # Input and output files
  infile  = sys.argv[1:-1]
  outfile = sys.argv[-1]

  # Start and end points of header records
  mol_start     =  0
  mol_end       = 20

  wn_init_start = 20
  wn_init_end   = 30

  wn_fin_start  = 30
  wn_fin_end    = 40

  nwave_start   = 40
  nwave_end     = 47

  temp_start    = 47
  temp_end      = 54

  press_start   = 54
  press_end     = 60

  res_start     = 70
  res_end       = 75
  
  
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

    line = f_in.readline()
    # Set the fields from the header:
    if i == 0:
      mol     = str(  line[mol_start    :mol_end    ].strip())  # Molecule name (eg CH4)
      wn_init = float(line[wn_init_start:wn_init_end].strip())  # Init wavenumber (cm-1)
      wn_fin  = float(line[wn_fin_start :wn_fin_end ].strip())  # Final wavenumber (cm-1)
      nwave   = int(  line[nwave_start  :nwave_end  ].strip())  # Number of data points
      press   = float(line[press_start  :press_end  ].strip())  # Pressure (Torr)
      res     = float(line[res_start    :res_end    ].strip())  # Resolution (cm-1)
      
      # Array to be filled with cross section data
      data = np.zeros((nwave, ntemp), np.double)
    else:
      # Check that the data from the other files match:
      if str(line[mol_start:mol_end].strip()) != mol:
        error += "for different species, "
      if (float(line[wn_init_start:wn_init_end].strip()) != wn_init or
          float(line[wn_fin_start :wn_fin_end ].strip()) != wn_fin  ):
        error += "with different wavelength ranges, "

      if error != "":
        print("Can't combine files {:s}".format(error))
        return

    temp[i] = line[temp_start:temp_end].strip()  # Add temperature (in K) to list
    
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
    f_out.write(' {:6.2f}'.format(temp[j]))
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

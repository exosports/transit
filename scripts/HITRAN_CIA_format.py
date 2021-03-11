#!/usr/bin/env python

"""
  Format HITRAN CIA data (Richard et al. 2012) for use in Transit.

  Usage:
  ------
  To run this code, execute from the Shell:
   ./HITRAN_CIA_format.py fileIN fileOUT [temp_step] [wnumber_step]

  Parameters:
  -----------
  fileIN: String
     File name of the input CIA file as given by Richard et al. (2012)
  fileOUT: String
     File name of the ouput CIA file for use with Transit.
  temp_step: Float
     Optional temperature sampling step to thin the outpout array.
  wnumber_step: Float
     Optional wavenumber sampling step to thin the outpout array.

  Examples:
  ---------
  For H2-H2 run from the Shell prompt:
    ./HITRAN_CIA_format.py  H2-H2_2011.cia  CIA_HITRAN_H2H2_0200-3000K_1.0-500um.dat  50  10

  Written by:
  -----------
  - Patricio Cubillos (UCF)
  - Dylan Bruce (UCF)
"""

import sys, os
import numpy as np

# Loschmidt number (cm-3):
N0 = 2.6867774e19

# Inout/output files:
filein  = sys.argv[1]
fileout = sys.argv[2]


# User-defined sampling rates:
tstep = None           # For temperature
if len(sys.argv) > 3:
  tstep = float(sys.argv[3])

wstep = None           # For wavenumber
if len(sys.argv) > 4:
  wstep = float(sys.argv[4])

# Read and extract data from files:
f = open(filein, "r")
lines = f.readlines()
f.close()

# Calculate total number of wavenumbers for each temperature,
# and how many temperatures were calculated:
header = lines[0].split()

species = header[0].split("-")
size    = int(header[3]) + 1
number_temps = len(lines) // size
number_waves = int(header[3])

# Initialize arrays:
T    = np.zeros(number_temps, np.double)
wn   = np.zeros(number_waves, np.double)
data = np.zeros((number_waves, number_temps), np.double)

# Write temperatures, wavenumbers, and coefficients into arrays:
for i in np.arange(number_temps):
    header = lines[size*i].split()[1:]
    T[i] = header[3]
    for j in np.arange(number_waves):
        point = lines[size*i + j + 1].split()
        if(i == 0):
            wn[j] = point[0]
        data[j][i] = point[1]


# Thin the arrays if requested:
if tstep is not None  and  wstep is not None:
  # Thinned arrays:
  Tthin  = np.arange(np.amin(T),  np.amax(T) +1, tstep)
  wnthin = np.arange(np.amin(wn), np.amax(wn)+1, wstep)
  # Indices corresponding to the thinned arrays:
  itemp = np.where(np.in1d(T,  Tthin ))[0]
  iwn   = np.where(np.in1d(wn, wnthin))[0]
  # Slice the data arrays:
  T  = T [itemp]
  wn = wn[iwn]
  data = data[iwn,:]
  data = data[:, itemp]

# Scale the opacity to the correct units (cm5 molecule-1 --> cm-1 amagat-2):
data *= N0**2


# Write to the output file:
fout = open(fileout, "w")

fout.write("# This file formats the tabulated {:s}-{:s} HITRAN CIA data from:\n"
           "# {:s}\n\n".format(species[0], species[1], filein))

fout.write("i {:s} {:s}\n".format(species[0], species[1]))
fout.write("t         ")
for j in np.arange(len(T)):
    fout.write("      {:4.0f}".format(T[j]))
fout.write("\n\n")

fout.write("# Wavenumber in cm-1, CIA coefficients in cm-1 amagat-2:\n")

# Write down the data:
for i in np.arange(len(wn)):
  fout.write("  {:7.1f} ".format(wn[i]))
  for j in np.arange(len(T)):
    fout.write(" {:.3e}".format(data[i,j]))
  fout.write("\n")

fout.close()

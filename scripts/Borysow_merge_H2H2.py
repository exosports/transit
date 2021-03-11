#!/usr/bin/env python

# Merge CIA data files from 60 to 7000 K:

# Download the CIA tabulated data from Borysow to the working directory:
#  http://www.astro.ku.dk/~aborysow/programs/final_CIA_LT.dat
#  http://www.astro.ku.dk/~aborysow/programs/final_CIA_HT.dat
#  http://www.astro.ku.dk/~aborysow/programs/CIA.H2H2.Yi

# Note that final_CIA_HT.dat has an incomplete column near the end that
# has to be patched before proceeding.

# Usage:
# ------
# Run from the shell:
# $ ./Borysow_merge_H2H2.py 


import sys, os
import numpy as np

root = ""

# Read and extract data from files:
f1 = open(root + "final_CIA_LT.dat", "r")
f2 = open(root + "final_CIA_HT.dat", "r")
f3 = open(root + "CIA.H2H2.Yi",      "r")
lines1 = f1.readlines()
lines2 = f2.readlines()
lines3 = f3.readlines()
f1.close()
f2.close()
f3.close()

# Extract temperature arrays:
T1 = lines1[1].split()[1:]
T2 = lines2[1].split()[1:]
T3 = lines3[1].split()[1:]

for i in np.arange(len(T1)):
  T1[i] = T1[i][:-1]
for i in np.arange(len(T2)):
  T2[i] = T2[i][:-1]
for i in np.arange(len(T3)):
  T3[i] = T3[i][:-1]
T1 = np.asarray(T1, np.double)
T2 = np.asarray(T2, np.double)
T3 = np.asarray(T3, np.double)

# Number ow wavelength samples:
size1 = len(lines1) - 3
size2 = len(lines2) - 3
size3 = len(lines3) - 3

# Extract wavenumber and extinction arrays:
wn1 = np.zeros(size1, np.double)
wn2 = np.zeros(size2, np.double)
wn3 = np.zeros(size3, np.double)
data1 = np.zeros((size1,len(T1)), np.double)
data2 = np.zeros((size2,len(T2)), np.double)
data3 = np.zeros((size3,len(T3)), np.double)

for i in np.arange(size1):
  info = lines1[i+3].split()
  wn1[i] = info[0]
  data1[i] = info[1:]

for i in np.arange(size2):
  info = lines2[i+3].split()
  wn2[i] = info[0]
  data2[i] = info[1:]

for i in np.arange(size3):
  info = lines3[i+3].split()
  wn3[i] = info[0]
  data3[i] = info[1:]

# Intersect wavenumbes for output array:
wn = np.intersect1d(np.intersect1d(wn1, wn2), wn3)
# Concatenate temperatures for output array:
T  = np.union1d(np.union1d(T1, T2), T3)

# Index intersecting values:
indata1 = np.in1d(wn1, wn)
indata2 = np.in1d(wn2, wn)
indata3 = np.in1d(wn3, wn)

# Output extinction array:
data = np.zeros((len(wn), len(T)), np.double)

# Intersected data:
idata1 = data1[np.where(indata1)]
idata2 = data2[np.where(indata2)]
idata3 = data3[np.where(indata3)][:,1:]

# os.chdir("/home/patricio/ast/esp01/pyrat/develop/inputs/CIA")
# Output file name:
fout = open("CIA_Borysow_H2H2_0060-7000K_0.6-500um.dat", "w")

fout.write(
"# This file merges the tabulated H2-H2 CIA data from:\n"
"# http://www.astro.ku.dk/~aborysow/programs/final_CIA_LT.dat\n"
"# http://www.astro.ku.dk/~aborysow/programs/final_CIA_HT.dat\n"
"# http://www.astro.ku.dk/~aborysow/programs/CIA.H2H2.Yi\n\n")

fout.write("i H2 H2\n")
fout.write("t        ")
for j in np.arange(len(T)):
    fout.write("      {:4.0f}".format(T[j]))
fout.write("\n\n")

fout.write("# Wavenumber in cm-1, CIA coefficients in cm-1 amagat-2:\n")

for i in np.arange(len(wn)):
  fout.write(" {:7.1f} ".format(wn[i]))
  for j in np.arange(len(T1)):
    fout.write(" {:.3e}".format(idata1[i,j]))
  for j in np.arange(len(T2)):
    fout.write(" {:.3e}".format(idata2[i,j]))
  for j in np.arange(len(T3)-1):
    fout.write(" {:.3e}".format(idata3[i,j]))
  fout.write("\n")

fout.close()


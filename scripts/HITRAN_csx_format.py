#!/usr/bin/env python

import numpy as np
import sys, os

def main():
    """
    This scripts reformats a HITRAN cross section file into the
    format read by the Transit program. Designed for use with
    Hargreaves, et al (2014) cross section files.
    
    Usage
    -----
    Run from the shell:
    $ ./HITRAN_csx_format.py <inputfile> <outputfile>

    Examples
    --------
    $ ./HITRAN_csx_format.py CH4_023C_XSC.xsc methane_296K.dat
    
    Modification History
    --------------------
    2015-07-23 rchallen Initial implementation.
    """
    infile, outfile = sys.argv[1:]
    # Open file for reading
    f_in = open(infile, 'r')

    ntoa = 2.6867805e19 # Conversion factor from number density (in cm-3)
                        # to amagats
    
    # Read in the header and split into fields
    line = str.split(f_in.readline())

    # Set the fields from the header
    mol     = str(  line[0]) # Molecule designation (eg CH4)
    wn_init = float(line[1]) # Initial wavenumber (cm-1)
    wn_fin  = float(line[2]) # Final wavenumber (cm-1)
    ndata   = int(  line[3]) # Number of data points
    temp    = float(line[4]) # Temperature (K)
    press   = float(line[5]) # Pressure (Torr)
    max_cs  = float(line[6]) # Maximum cross section (cm2 molecule-1)
    res     = float(line[7]) # Resolution of measurement (cm-1)
    molname = str(  line[8]) # Full molecule name (eg Methane)

    # Number of data points per line
    linelen = 10

    # Array to be filled with cross section data
    data = np.zeros(ndata)

    # Wavenumber array to be printed in output file
    wn   = np.linspace(wn_init, wn_fin, ndata)

    # For counting loop repetitions
    j = 0

    # Loop through the lines of data
    for x in f_in:
        # Loop through fields in each line
        for i in range(linelen):
            # Only read data that is there (don't read extra fields at the end)
            if i+linelen*j < ndata:
              data[i+linelen*j] = str.split(x)[i]
        # Increment number of loop repetitions
        j +=1

    data *= ntoa # Convert to correct units
        
    f_in.close()

    # Open file for writing
    f_out = open(outfile, 'w')

    # Write the header
    f_out.write('# This file is a reformatted version of Hargreaves et al \n')
    f_out.write('# (2015) hot methane cross section file \n \n')

    # Write the molecule and temperature lines
    f_out.write('i ' + mol       + '\n')
    f_out.write('t ' + str(temp) + '\n')

    # Write comment on units
    f_out.write('\n # Wavenumber in cm-1, cross section in cm-1 amagat-1 \n')

    # Write data in format required by Transit program
    for i in range(ndata):
      f_out.write("{:.2f}".format(wn[i]) + '   ' + "{:.4e}".format(data[i]) + '\n')

    f_out.close()
    

    
if __name__ == "__main__":
    main()

# ****************************** START LICENSE ******************************
# Transit, a code to solve for the radiative-transifer equation for
# planetary atmospheres.
# 
# This project was completed with the support of the NASA Planetary
# Atmospheres Program, grant NNX12AI69G, held by Principal Investigator
# Joseph Harrington. Principal developers included graduate students
# Patricio E. Cubillos and Jasmina Blecic, programmer Madison Stemm, and
# undergraduate Andrew S. D. Foster.  The included
# 'transit' radiative transfer code is based on an earlier program of
# the same name written by Patricio Rojo (Univ. de Chile, Santiago) when
# he was a graduate student at Cornell University under Joseph
# Harrington.
# 
# Copyright (C) 2014 University of Central Florida.  All rights reserved.
# 
# This is a test version only, and may not be redistributed to any third
# party.  Please refer such requests to us.  This program is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.
# 
# Our intent is to release this software under an open-source,
# reproducible-research license, once the code is mature and the first
# research paper describing the code has been accepted for publication
# in a peer-reviewed journal.  We are committed to development in the
# open, and have posted this code on github.com so that others can test
# it and give us feedback.  However, until its first publication and
# first stable release, we do not permit others to redistribute the code
# in either original or modified form, nor to publish work based in
# whole or in part on the output of this code.  By downloading, running,
# or modifying this code, you agree to these conditions.  We do
# encourage sharing any modifications with us and discussing them
# openly.
# 
# We welcome your feedback, but do not guarantee support.  Please send
# feedback or inquiries to:
# 
# Joseph Harrington <jh@physics.ucf.edu>
# Patricio Cubillos <pcubillos@fulbrightmail.org>
# Jasmina Blecic <jasmina@physics.ucf.edu>
# 
# or alternatively,
# 
# Joseph Harrington, Patricio Cubillos, and Jasmina Blecic
# UCF PSB 441
# 4111 Libra Drive
# Orlando, FL 32816-2385
# USA
# 
# Thank you for using transit!
# ******************************* END LICENSE ******************************

import struct, time
import numpy as np
import utils as ut
import constants as c

class tioschwenke():
  """
  Notes:
  ------
  Download the linelist from:
  The binary linelist:
    http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
  and the partition function:
    http://kurucz.harvard.edu/molecules/tio/tiopart.dat

  There might be a problem with the linebreak character of the partition
  function.  One way to fix is, on vim do: :%s/\r/\r/g
  """
  def __init__(self):
    pass

  def getpf(self, dbfile, pffile, verbose):
    """
    Open, read, and store values from a partition function file pffile.
 
    Parameters:
    -----------
    dbfile: String
      Database filename.
    pffile: String
      Partition Function filename.
    verbose: Integer
      Verbosity threshold.
 
    Returns:
    --------
    PF: 2D ndarray
       A 2D array (N temperatures, N isotopes + 1) with the partition
      function values for each isotope (columns) as function of temperature
       (first column).
    DBname: String
       The database name.
    isoNames: List
       List with the P&S isotope names.
    mass: List
       List with the P&S isotope masses.

    Modification History:
    ---------------------
    2012-11-30  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2014-03-11  patricio  Adapted to work with pylineread.
    """
    # Open partition function file:
    partDB = open(pffile)
    ut.lrprint(verbose, "Parsing the partition function file.")
 
    # Read partition function file lines:
    PFlines = partDB.readlines()
    partDB.close()

    # Get isotopes names from first line:
    isoNames = PFlines[0].split()[1:]  # Skip first word
    # Skip header lines (first 5 lines):
    PFlines = PFlines[c.TS_PF_IGNORE:]
    # Number of Temperatures (number of lines):
    Ntemp = len(PFlines)
    # Allocate array for table of Temperature and PF values:
    PF = np.zeros((Ntemp, len(isoNames)+1), np.double)
 
    for i in np.arange(Ntemp):
      # Store values in array (automatic casting, wow!):
      PF[i] = PFlines[i].split()
      ut.lrprint(verbose-10, "Reading line %d: T=%4d."%(i, PF[i,0]))
 
    return PF, c.TS_NAME, isoNames, c.TS_MASS


  def dbread(self, dbfile, iwl, fwl, verbose, *args):
    """
    Read the Partridge and Schwenke H2O database (dbfile) between the
    wavelengths iwl and fwl.
 
    Parameters:
    -----------
    dblist: String
       Partridge & Schwenke database filename.
    iwl: Scalar
       Initial wavelength limit (in microns).
    fwl: Scalar
       Final wavelength limit (in microns).
    verbose: Integer
       Verbosity threshold.
    args:
       Additional arguments, not needed for pands.
 
    Returns:
    --------
    wlength: 1D ndarray (double)
      Line-transition central wavelength (microns).
    gf: 1D ndarray (double)
      gf value (unitless).
    elow: 1D ndarray (double)
      Lower-state energe (centimeter^-1).
    isoID: 2D ndarray (integer)
      Isotope index (1, 2, 3, ...).

    Modification History:
    ---------------------
    2012-11-30  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2014-03-11  patricio  Adapted to work with pylineread.
    2014-07-06  patricio  Updated return statement.
    """
 
    # Create a table of logarithms:
    ut.lrprint(verbose-1, "Creating log table.")
    tablog = 10.0**(0.001*(np.arange(c.TS_NCODIDX+1) - 16384))
 
    # Open the file:
    data = open(dbfile, "rb")
 
    # Get the number of lines in the file:
    data.seek(0, 2)                     # Set pointer at the file's end
    nlines = data.tell() / c.TS_RECSIZE # Number of lines (bytes/record_size)
 
    # Rewrite wavelength limits as given in the P&S file:
    iwav = iwl * c.MTC / c.NTC         # Microns to nanometer
    fwav = fwl * c.MTC / c.NTC
    iwav = np.log(iwav) / c.RATIOLOG
    fwav = np.log(fwav) / c.RATIOLOG
 
    # Find the positions of iwl and fwl, then jump to wl_i position:
    irec = ut.binsearch(data, iwav, 0,    nlines, 'ts', 0)
    frec = ut.binsearch(data, fwav, irec, nlines, 'ts', 1)
    nread = frec - irec + 1  # Number of records to read

    # Store data in two arrays for doubles and integers:
    # For Wavelength, Elow, and log(gf):
    wlength = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    # For Isotope index:
    isoID   = np.zeros(nread,     int)
 
    ut.lrprint(verbose, "Beginning to read Schwenke database, between "
                        "records %d and %d."%(irec, frec))
    # When the wavelength surpasses the max wavelength, stop the loop
    chk = 1  # Check-point counter
    i   = 0  # Stored record index
    wl_intvl = float((frec - irec)/20)  # Check-point interval

    iw   = np.zeros(nread, int)
    ieli = np.zeros(nread, np.short)
    ielo = np.zeros(nread, np.short)
    igf  = np.zeros(nread, np.short)

    while (i < nread):
      # Read a record:
      data.seek((irec+i)*c.TS_RECSIZE)
      iw[i], ieli[i], ielo[i], igf[i] = struct.unpack('ihhh',
                                                      data.read(c.TS_RECDATA))
 
      # Print a checkpoint statement every 1/20th interval:
      if verbose > 1:
        pos = float(data.tell()/c.TS_RECSIZE)
        if (pos/wl_intvl)%1 == 0.0:
          ut.lrprint(verbose-1, "checkpoint %d/20..."%chk)
          chk += 1
          ut.lrprint(verbose-3, "iwl: %d, ielow: %5d, igf: "
                                "%6d"%(iw[i], ielo[i], igf[i]))
          ut.lrprint(verbose-2, "Wavelength: %.3f, IsoID: %d, Elow: %.5e, "
                     "gf: %.5e"%(np.exp(iw[i] * c.RATIOLOG) * c.NTC/c.MTC,
                                 np.abs(ieli[i]) - 8950,
                                 tablog[ielo[i]], tablog[igf[i]]))
      i += 1

    # Convert wavelength to TLI format (microns):
    wlength[:] = np.exp(iw * c.RATIOLOG) * c.NTC/c.MTC
    # Get gf from log table:
    gf[:]      = tablog[igf]
    # Get lowest state energy from log table:
    elow[:]    = tablog[ielo]
    # Get isotopic index:
    isoID[:]   = np.abs(ieli) - 8950
    ut.lrprint(verbose, "Done.\n")
    data.close()
    return wlength, gf, elow, isoID

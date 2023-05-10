# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# Transit is under an open-source, reproducible-research license (see LICENSE).

import struct, time
import numpy as np
import utils as ut
import constants as c
from driver import dbdriver

class tioschwenke(dbdriver):
  """
  Notes:
  ------
  Linelist and partition function downloaded from:
    http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
    http://kurucz.harvard.edu/molecules/tio/tiopart.dat

  There might be a problem with the linebreak character of the partition
  function.  One way to fix is, on vim do: :%s/\r/\r/g
  """
  def __init__(self, dbfile, pffile):
    super(tioschwenke, self).__init__(dbfile, pffile)
    self.recsize = 16 # Record size (bytes)
    self.recdata = 10 # Useful record data size (bytes)
    self.ratiolog  = np.log(1.0 + 1.0/2000000)
    # Table of logarithms:
    self.tablog    = 10.0**(0.001*(np.arange(32769) - 16384))
    self.pf_ignore   = 1 # Lines to ignore (header) at the begining of PF file
    self.pf_isonames = 0 # PF line with isotopes names 

    # Database name:
    self.name ="Schwenke TiO (1998)"
    # Isotopes names:
    self.isotopes = self.getinfo()
    # Isotopes mass:  46Ti,  47Ti,  48Ti,  49Ti,  50Ti
    self.mass = [61.94754403, 62.94667863, 63.94286193, 64.94278573,
                 65.93970673]
    # Isotopic abundance ratio:
    self.isoratio = [0.080, 0.073, 0.738, 0.055, 0.054]
    # Molecule name:
    self.molecule = "TiO"

  def readwl(self, dbfile, irec):
    """
    Read wavelength parameter from irec record in dbfile database.

    Parameters:
    -----------
    dbfile: File object
       File where to extract the wavelength.
    irec: Integer
       Index of record.

    Returns:
    --------
    rec_wl: Float
       Wavelength value at record irec, as given in dbfile database.

    Modification History:
    ---------------------
    2014-03-05  patricio  Initial implementation, based on Madison's
                          code.          pcubillos@fulbrightmail.org
    2014-03-10  patricio  Updated dbtype to match command-line-argument
                          sytax.  Updated HITRAN data type.
    2014-03-24  patricio  Moved to db_tioschwenke.py from utils.py.
    """
    # Set pointer at required wavelength record location:
    dbfile.seek(irec*self.recsize)
    # Read and extract the wavelength:
    rec_wl = struct.unpack('ihhh', dbfile.read(self.recdata))[0]

    return rec_wl


  def getpf(self, verbose):
    """
    Open, read, and store values from the partition function file.
 
    Parameters:
    -----------
    verbose: Integer
      Verbosity threshold.
 
    Returns:
    --------
    Temp: 1D ndarray
       An array with the tabulated temperatures at which the PF is evaluated
    PF: 2D ndarray
       A 2D array (N isotopes, N temperatures) with the partition
      function values for each isotope as function of temperature

    Modification History:
    ---------------------
    2012-11-30  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2014-03-11  patricio  Adapted to work with pylineread.
    2014-03-24  patricio  Rewritten as subclass of dbdriver. 
    """
    # Open and read the partition function file:
    partDB = open(self.pffile)
    PFlines = partDB.readlines()
    partDB.close()

    # Skip header lines (first 5 lines):
    PFlines = PFlines[self.pf_ignore:]
    # Number of Temperatures (number of lines):
    Ntemp = len(PFlines)
    # Number of isotopes:
    Niso  = len(PFlines[0].split()) - 1

    # Allocate array for table of temperatures:
    Temp = np.zeros(Ntemp, np.double)
    # Allocate array for table of PF values:
    PF = np.zeros((Niso, Ntemp), np.double)
 
    for i in np.arange(Ntemp):
      values = PFlines[i].split()
      # Store values in array:
      Temp[i]  = values[0]
      PF[:, i] = values[1:]
      ut.lrprint(verbose-15, "Reading line %d: T=%4d."%(i, Temp[i]))
 
    return Temp, PF


  def getinfo(self):
    """
    Get isotopes names from partition function file

    Modification History:
    ---------------------
    2014-03-24  patricio  Initial implementation.
    """
    # Open and read the partition function file:
    partDB = open(self.pffile)
    PFlines = partDB.readlines()
    partDB.close()

    # Get isotopes names from first line:
    return PFlines[self.pf_isonames].split()[1:]  # Skip first word


  def dbread(self, iwl, fwl, verbose, *args):
    """
    Read the Partridge and Schwenke H2O database (dbfile) between the
    wavelengths iwl and fwl.
 
    Parameters:
    -----------
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
      Lower-state energy (cm^-1).
    isoID: 2D ndarray (integer)
      Isotope index (1, 2, 3, ...).

    Modification History:
    ---------------------
    2012-11-30  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2014-03-11  patricio  Adapted to work with pylineread.
    2014-03-24  patricio  Rewritten as subclass of dbdriver.
    2014-07-06  patricio  Updated return statement.
    2014-09-23  patricio  Updated for pylineread 5.0.
    """
 
    # Open the file:
    data = open(self.dbfile, "rb")
 
    # Get the number of lines in the file:
    data.seek(0, 2)                     # Set pointer at the file's end
    nlines = data.tell() // self.recsize # Number of lines (bytes/record_size)
 
    # Rewrite wavelength limits as given in the Database file:
    iwav = iwl * c.MTC / c.NTC         # Microns to nanometer
    fwav = fwl * c.MTC / c.NTC
    iwav = np.log(iwav) / self.ratiolog
    fwav = np.log(fwav) / self.ratiolog
 
    # Find the positions of iwl and fwl:
    irec = self.binsearch(data, iwav, 0,    nlines-1, 0)
    frec = self.binsearch(data, fwav, irec, nlines-1, 1)
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
    interval = (frec - irec)//20  # Check-point interval

    iw   = np.zeros(nread, int)
    ieli = np.zeros(nread, np.short)
    ielo = np.zeros(nread, np.short)
    igf  = np.zeros(nread, np.short)

    while (i < nread):
      # Read a record:
      data.seek((irec+i)*self.recsize)
      iw[i], ieli[i], ielo[i], igf[i] = struct.unpack('ihhh',
                                                      data.read(self.recdata))
 
      # Print a checkpoint statement every 1/20th interval:
      if verbose > 1:
        pos = float(data.tell()//self.recsize)
        if (pos/interval)%1 == 0.0:
          ut.lrprint(verbose-1, "checkpoint %d/20..."%chk)
          chk += 1
          ut.lrprint(verbose-3, "iwl: %d, ielow: %5d, igf: "
                                "%6d"%(iw[i], ielo[i], igf[i]))
          ut.lrprint(verbose-2, "Wavelength: %.3f, IsoID: %d, Elow: %.5e, "
                     "gf: %.5e"%(np.exp(iw[i] * self.ratiolog) * c.NTC/c.MTC,
                                 np.abs(ieli[i]) - 8950,
                                 self.tablog[ielo[i]], self.tablog[igf[i]]))
      i += 1

    # Convert wavelength to TLI format (microns):
    wlength[:] = np.exp(iw * self.ratiolog) * c.NTC/c.MTC
    # Get gf from log table:
    gf[:]      = self.tablog[igf]
    # Get lowest state energy from log table:
    elow[:]    = self.tablog[ielo]
    # Get isotopic index:
    isoID[:]   = np.abs(ieli) - 8950
    ut.lrprint(verbose, "Done.\n")
    data.close()
    return wlength, gf, elow, isoID

# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# Transit is under an open-source, reproducible-research license (see LICENSE).

import struct
import numpy as np
import utils as ut
import constants as c
from driver import dbdriver

class pands(dbdriver):
  """
  Notes:
  ------
  Linelist and partition function files downloaded from:
    http://kurucz.harvard.edu/molecules/h2o/h2ofastfix.bin
    http://kurucz.harvard.edu/molecules/h2o/h2opartfn.dat
  """
  def __init__(self, dbfile, pffile):
    """
    Initialize Basic data for the Database.

    Modification History:
    ---------------------
    2014-07-27  patricio  Added Documentation. Added self.molecule with
                          name of the molecule.
    """
    super(pands, self).__init__(dbfile, pffile)

    # Database name:
    self.name = "Partridge & Schwenke (1997)"
    # Isotopes names:
    self.isotopes  = ['1H1H16O',   '1H1H17O',   '1H1H18O',   '1H2H16O'] 
    # Isotopes masses:
    self.mass      = [18.01056468, 19.01478156, 20.01481046, 19.01684143]
    # Isotopic abundance ratio:
    self.isoratio  = [0.997000,    0.000508,    0.000508,    0.001984]

    # Molecule name:
    self.molecule  = "H2O"

    self.ratiolog  = np.log(1 + 1/2e6)
    # Table of logarithms: 
    self.tablog = 4*10.0**(0.001*(np.arange(32769) - 16384))
    self.recsize     = 8 # Record size
    self.pf_ignore   = 6 # Lines to ignore at the begining of PF file.
    self.pf_isonames = 3 # PF line with isotopes names


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
    2014-03-24  patricio  Moved to db_pands.py from utils.py.
    """
    # Set pointer at required wavelength record location:
    dbfile.seek(irec*self.recsize)
    # Read and extract the wavelength:
    rec_wl = struct.unpack('Ihh', dbfile.read(self.recsize))[0]
 
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
    2013        madison   Initial version based on P. Rojo's C code.
                                                madison.stemm@ucf.edu
    2014-03-05  patricio  Added documentation.  pcubillos@fulbrightmail.org
    2014-03-08  patricio  Moved to pands.py
    2014-03-10  patricio  Added DBname return.
    2014-03-24  patricio  Rewritten as subclass of dbdriver.
    """
    ut.lrprint(verbose, "Parsing the partition function file.")
    # Open and read partition function file:
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
      # Store values in arrays (automatic casting, so much wow!):
      Temp[i]  = values[0]
      PF[:, i] = values[1:]
      ut.lrprint(verbose-15, "Reading line %d: T=%4d."%(i, Temp[i]))
 
    return Temp, PF


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
      Lower-state energe (centimeter^-1).
    isoID: 2D ndarray (integer)
      Isotope index (1, 2, 3, ...).

    Modification History:
    ---------------------
    2013        madison   Initial implementation. madison.stemm@ucf edu
    2014-03-05  patricio  Added documentation.    pcubillos@fulbrightmail.org
    2014-03-08  patricio  Moved to pands.py
    2014-03-24  patricio  Rewritten as subclass of dbdriver.
    2014-07-06  patricio  Updated return statement.
    """
 
    # Open the binary file:
    data = open(self.dbfile, "rb")
 
    # Get the number of lines in the file:
    data.seek(0, 2)                      # Set pointer at the file's end
    nlines = data.tell()//self.recsize # Number of lines (8 bytes per line)
 
    # Rewrite wavelength limits as given in the P&S file:
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
 
    ut.lrprint(verbose, "Beginning to read P&S database, between "
                        "records %d and %d."%(irec, frec))
    # When the wavelength surpasses the max wavelength, stop the loop
    chk = 1  # Check-point counter
    i   = 0  # Stored record index
    interval = float((frec - irec)//20)  # Check-point interval

    iw   = np.zeros(nread, int)
    ielo = np.zeros(nread, np.short)
    igf  = np.zeros(nread, np.short)

    data.seek(irec*self.recsize) 
    while (i < nread):
      # Read a record:
      iw[i], ielo[i], igf[i] = struct.unpack('Ihh', data.read(self.recsize))
 
      # Print a checkpoint statement every 1/20th interval
      if verbose > 1:
        pos = float(data.tell()/self.recsize)
        if (pos/interval)%1 == 0.0:
          ut.lrprint(verbose-1, "checkpoint %d/20..."%chk)
          chk += 1
          ut.lrprint(verbose-3, "iwl: %d, ielow: %5d, gf: "
                                "%6d"%(iw[i], ielo[i], igf[i]))
          ut.lrprint(verbose-2, "Wavelength: %.3f, IsoID: %d, Elow: %.5e, "
                    "gf: %.5e"%(np.exp(iw[i] * self.ratiolog) * c.NTC/c.MTC,
                                2*(ielo[i] < 0) + 1*(igf[i] < 0),
                                np.abs(ielo[i]), self.tablog[np.abs(igf[i])]))
      i += 1

    # Convert wavelength to TLI format (microns):
    wlength[:] = np.exp(iw * self.ratiolog) * c.NTC/c.MTC
    # Get gf fom log table:
    gf[:]      = self.tablog[np.abs(igf)]
    # Set energy of lowest transition level:
    elow[:]    = np.abs(ielo)
    # Assign indices for isotopes based on Kurucz's indices-1:
    isoID[:]   = 2*(ielo < 0) + 1*(igf < 0)

    ut.lrprint(verbose, "Done.\n")
    data.close()
    return wlength, gf, elow, isoID

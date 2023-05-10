import sys, os
import six
import numpy as np
import subprocess as sp
import scipy.constants as sc
import struct

import utils as ut
import constants as c
from driver import dbdriver

#db_hitran.py dir:
DBHdir = os.path.dirname(os.path.realpath(__file__))


class repack(dbdriver):
  def __init__(self, dbfile, pffile, defn):
    """
    Initialize the basic database info.

    Parameters:
    -----------
    dbfile: String
       File with the Database info as given from HITRAN.
    pffile: String
       File with the partition function.
    defn  : String
       Path/to/file for HITRAN configuration file.
    """
    super(repack, self).__init__(dbfile, pffile)

    # Get info from file name:
    self.molecule, self.dbtype = os.path.split(dbfile)[1].split('_')[0:2]

    self.defn = defn
    self.isotopes = self.get_isotopes()
    # Get info from HITRAN configuration file:
    isotopes, mass, isoratio, gi = self.getEMinfo(self.defn)
    iso_index = []
    for iso in isotopes:
        for j,isotope in enumerate(self.isotopes):
            if iso == isotope:
                iso_index.append(j)
    iso_index = np.array(iso_index)
    self.mass = np.array(mass)[iso_index]
    self.isoratio = np.array(isoratio)[iso_index]

    # Database name: 
    self.name = 'repack ' + self.molecule
    self.iso_dict = {int(iso):i for i,iso in enumerate(self.isotopes)}

  def get_isotopes(self):
    # Open and read partition function file:
    with open(self.pffile) as partDB:
        info = partDB.readline()
        isotopes = partDB.readline().split()[1:]
    return isotopes



  def readwl(self, dbfile, irec):
    """
    Read wavelength parameter from irec record in dbfile database.

    Parameters:
    -----------
    dbfile: File object
       File where to extract the wavelength.
    irec: Integer
       Index of record.
    dbtype: String
       Database type: ['ps', 'hit']

    Returns:
    --------
    rec_wl: Float
       Wavelength value at record irec, as given in dbfile database.
    """
    # Set pointer at required wavenumber record:
    dbfile.seek(irec*28)
    # Read:
    wavenumber = struct.unpack('d', dbfile.read(8))[0]
    # Convert to float:
    return wavenumber


  def getpf(self, verbose):
    """                                    
    Calculate the partition function for a grid of temperatures for
    HITRAN molecules.

    Notes:
    ------
    This function is a wrapper of the TIPS library written by
    Patricio Cubillos, which is itself a C implementation of the
    Fortran TIPS routine by R. R. Gamache.
    The range of temperatures is limited to: 70K -- 3000K.

    Parameters:
    -----------
    verbose: Integer
       Verbosity threshold.

    Returns:
    --------
    Temp: 1D float ndarray
       The array of temeratures where the partition function was evaluated.
    PF: 2D float ndarray
       A 2D array (N temperatures, N isotopes) with the partition
       function values for each isotope (columns) as function of temperature
       (first column).
    """
    ut.lrprint(verbose, "Parsing the partition function file.")
    # Open and read partition function file:
    partDB = open(self.pffile)
    PFlines = partDB.readlines()
    partDB.close()

    # Skip header lines (first 5 lines):
    PFlines = PFlines[2:]
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


  def getMolec(self):
    """
    Get the HITRAN molecule index

    Parameters:
    -----------
    dbfile: String
       The data filename

    Modification History:
    ---------------------
    2014-03-21  patricio  Initial implementation.
    """
    # Open file and read first two characters:
    data = open(self.dbfile, "r")
    molID  = data.read(self.recmollen)
    data.close()
    # Set database name:
    return molID #self.molname[molID-1]


  def getEMinfo(self, defn):
    """
    Get HITRAN info from configuration file.

    Parameters:
    -----------
    defn : String
       Path/to/file for HITRAN configuration file

    Returns:
    --------
    isotopes: Isotopes names
    mass:     Isotopes mass
    isoratio: Isotopic abundance ratio
    gi:       State-independent statistical weight

    """
    # Read HITRAN configuration file
    hfile = open(self.defn, 'r')
    lines = hfile.readlines()
    hfile.close()
  
    isotopes = []
    mass     = []
    isoratio = []
    gi       = []
 
    # Get values for our molecule: 
    for i in np.arange(len(lines)):
      if len(lines[i].split())>2  and lines[i].split()[1] == self.molecule:
        line = lines[i].split()
        gi.      append(  int(line[4]))
        isotopes.append(      line[3] )
        isoratio.append(float(line[5]))
        mass.    append(float(line[6]))

    return isotopes, mass, isoratio, gi    

 
  def dbread(self, iwl, fwl, verbose, *args):
    """
    Read a HITRAN or HITEMP database (dbfile) between wavelengths iwl and fwl.

    Parameters:
    -----------
    dbfile: String
       A HITRAN or HITEMP database filename.
    iwl: Scalar
       Initial wavelength limit (in microns).
    fwl: Scalar
       Final wavelength limit (in microns).
    verbose: Integer
       Verbosity threshold.
    pffile: String
       Partition function filename.

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
    2012-12-10  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2014-03-10  patricio  Adapted for pylineread.
    2014-07-06  patricio  Updated to return 1D arrays.
    """
    # Open HITRAN file for reading:
    data = open(self.dbfile, "rb")
    # Get Total number of transitions in file:
    data.seek(0, 2)
    nlines   = data.tell() // 28

    # Rewrite wavelength limits as given in HITRAN file:
    fwn = 1.0 / (iwl * c.MTC)  # Microns to centimeters
    iwn = 1.0 / (fwl * c.MTC)  # Microns to centimeters

    # Find the record index for iwn and fwn:
    irec_init = self.binsearch(data, iwn, 0, nlines-1,  0)
    irec_fin  = self.binsearch(data, fwn, irec_init, nlines-1, 1)

    # Allocates arrays for values to extract:
    wnumber = np.zeros(irec_fin-irec_init+1, np.double)
    gf      = np.zeros(irec_fin-irec_init+1, np.double)
    elow    = np.zeros(irec_fin-irec_init+1, np.double)
    isoID   = np.zeros(irec_fin-irec_init+1,       int)

    i = 0  # Stored record index
    chk = 0
    interval = (irec_fin - irec_init)//20  # Check-point interval
    while irec_init + i <= irec_fin:
      data.seek( (irec_init+i) * 28 )
      wnumber[i], elow[i], gf[i], isoID[i] = struct.unpack('dddi', data.read(28))
      # Print a checkpoint statement every 1/20th interval
      if verbose > 1:
        if (i%interval) == 0.0:
          ut.lrprint(verbose-1, "checkpoint %d/20..."%chk)
          chk += 1
      i += 1

    data.close()
    isoID = np.array([self.iso_dict[iso] for iso in isoID])
    return 1e4/wnumber, gf, elow, isoID

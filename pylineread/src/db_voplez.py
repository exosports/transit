# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# Transit is under an open-source, reproducible-research license (see LICENSE).

import struct, time
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import utils as ut
import constants as c
from driver import dbdriver

class voplez(dbdriver):
  """
  Notes:
  ------
  Data obtained from private communication with B. Plez:
    http://www.pages-perso-bertrand-plez.univ-montp2.fr/
  """
  def __init__(self, dbfile, pffile):
    """
    """
    super(voplez, self).__init__(dbfile, pffile)

    # Database name:
    self.name = "Bertrand Plez VO"
    # Isotopic names:
    self.isotopes = ["16"]  # I'm Using AFGL naming convention
    # Isotopic masses:
    self.mass     = [66.941]
    # Isotopic abundance ratio:
    self.isoratio = [1.0]

    # Molecule name:
    self.molecule = "VO"

    # Partition-function polynomial coefficients:
    # (from communication with B. Pelz):
    self.PFcoeffs = np.array([ 6.62090157e+02, -4.03350494e+02,
                               9.82836218e+01, -1.18526504e+01,
                               7.08429905e-01, -1.67235124e-02])

    # Other utilities:
    self.recsize  = 53  # Record length
    self.recwnpos = 33  # Wavenumber position in record
    self.recelpos = 44  # Elow       position in record
    self.recgfpos = 21  # gf         position in record
    self.recwnlen = 10  # Record lengths
    self.recwnend = 43  # Record lengths
    self.recelend = 50
    self.recgfend = 32

  def readwl(self, dbfile, irec):
    """
    """
    # Set pointer at required wavenumber record:
    dbfile.seek(irec*self.recsize + self.recwnpos)
    # Read:
    wavenumber = dbfile.read(self.recwnlen)
    # Convert to float:
    rec_wl = 1.0 / (float(wavenumber) * c.MTC)

    return rec_wl


  def getpf(self, verbose):
    """
    Calculate the partition function for a grid of temperatures for VO.
 
    Parameters:
    -----------
    verbose: Integer
      Verbosity threshold.
 
    Returns:
    --------
    Temp: 1D float ndarray
       The array of temperatures where the partition function was evaluated.
    PF: 2D float ndarray
       A 2D array (of shape [Nisotopes,Ntemperatures]) with the partition
       function values for each isotope (columns) as function of temperature.

    Notes:
    ------
    The partition function is valid for the range of temperatures from
    1000K to 7000K. It is extrapolated down to 0K.

    Sample Return:
    Temp = np.linspace(1000., 7000., 13):
    array([ 1000.,  1500.,  2000.,  2500.,  3000.,  3500.,  4000.,  4500., 5000.,  5500.,  6000.,  6500.,  7000.])
    PF
    array([[   6696.28281847,    1000.        ],
           [  12463.29391522,    1500.        ],
           [  20096.2494866 ,    2000.        ],
           [  29606.07281774,    2500.        ],
           [  41175.25306071,    3000.        ],
           [  55080.93657951,    3500.        ],
           [  71666.40723555,    4000.        ],
           [  91330.40363125,    4500.        ],
           [ 114524.07531435,    5000.        ],
           [ 141751.43620469,    5500.        ],
           [ 173571.51471102,    6000.        ],
           [ 210601.37829024,    6500.        ],
           [ 253519.63927169,    7000.        ]])

    Modification History:
    ---------------------
    2015-06-14  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2015-06-21  sally     Calculates pf for Vanadium (II) Oxide (VO)
    2015-06-14  Jasmina	  Extrapolated pf values from 0 to 1000 K
    """
    # Temperature array:
    Temp = np.arange(1000.0, 7001.0, 50.0)
    Ntemp = len(Temp)

    # Number of isotopes:
    Niso  = 1

    # Intialize PF array:
    PF = np.zeros((Niso, Ntemp), np.double)

    # Calculate log(PF) at each Temp:
    for i in np.arange(Ntemp):
      # Formula from Irwin 1981, ApJS 45, 621 (equation #2):
      PF[0,i] = (self.PFcoeffs[0]                      +
                 self.PFcoeffs[1]* np.log(Temp[i])     +
                 self.PFcoeffs[2]*(np.log(Temp[i]))**2 +
                 self.PFcoeffs[3]*(np.log(Temp[i]))**3 +
                 self.PFcoeffs[4]*(np.log(Temp[i]))**4 +
                 self.PFcoeffs[5]*(np.log(Temp[i]))**5 )
    # Get the exponential of log(PF):
    PF = np.exp(PF)

    # Add start point for temp and PF arrays
    temp_ext = np.insert(Temp,  0, 0.0)
    PF_ext   = np.insert(PF[0], 0, 0.0)

    # Interpolate using quadratic spline
    temp_int = np.arange(50.0,1000.0, 50.0)
    s = InterpolatedUnivariateSpline(temp_ext, PF_ext, k=2)
    PF_int = s(temp_int)

    # Insert interpolated range into Temp and PF arrays
    Temp      = np.insert(temp_ext, 1, temp_int)
    PF_insert = np.insert(PF_ext  , 1, PF_int)

    # Update PF array
    PF = np.zeros((Niso, len(Temp)), np.double)
    PF[0] = PF_insert

    return Temp, PF


  def dbread(self, iwl, fwl, verbose, *args):
    """
    Read the B. Plez VO database between the wavelengths iwl and fwl.
 
    Parameters:
    -----------
    iwl: Scalar
       Initial wavelength limit (in microns).
    fwl: Scalar
       Final wavelength limit (in microns).
    verbose: Integer
       Verbosity threshold.
    args:
       Additional arguments, not needed for voplez.
 
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
    """
 
    # Open the file:
    data = open(self.dbfile, "r")
    data.seek(0, 2)
    nlines   = data.tell() // self.recsize

    # Find the record index for iwl and fwl:
    irec_init = self.binsearch(data, iwl, 0,         nlines-1, 1)
    irec_fin  = self.binsearch(data, fwl, irec_init, nlines-1, 0)

    # Number of records to read:
    nread = irec_fin - irec_init + 1

    # Store data in two arrays for doubles and integers:
    # For Wavelength, Elow, and log(gf):
    wlength = np.zeros(nread, np.double)
    gf      = np.zeros(nread, np.double)
    elow    = np.zeros(nread, np.double)
    # For Isotope index:
    isoID   = np.zeros(nread,     int)
 
    ut.lrprint(verbose, "Beginning to read Plez VO database, between "
                        "records %d and %d."%(irec_init, irec_fin))
    # When the wavelength surpasses the max wavelength, stop the loop
    chk = 1  # Check-point counter
    i   = 0  # Stored record index
    interval = (irec_fin - irec_init)//20  # Check-point interval

    wnumber = np.zeros(nread)
    gf      = np.zeros(nread)
    elow    = np.zeros(nread)

    while (i < nread):
      # Read a record:
      data.seek((irec_init+i) * self.recsize)
      line = data.read(self.recsize)
      wnumber[i] = float(line[self.recwnpos:self.recwnend])
      gf     [i] = float(line[self.recgfpos:self.recgfend])
      elow   [i] = float(line[self.recelpos:self.recelend])

      # Print a checkpoint statement every 1/20th interval:
      if verbose > 1:
        pos = float(data.tell()//self.recsize)
        if (pos/interval)%1 == 0.0:
          ut.lrprint(verbose-1, "checkpoint %d/20..."%chk)
          chk += 1
          ut.lrprint(verbose-3, "Wavenumber (cm-1): {:.2f}, Elow (eV): {:.3f}, "
                     "gf: {:.4e}".format(wnumber[i], elow[i], gf[i]))
          ut.lrprint(verbose-2, "Wavelength (um): {:.3f},  IsoID: {:d},  "
                     "Elow (cm-1): {:.5e}, gf: {:.5e}".
                     format(1.0/(wnumber[i]*c.MTC), isoID[i],
                            elow[i]*c.eV2kayser, gf[i]))
      i += 1

    # Convert wavelength to TLI format (microns):
    wlength[:] = 1.0 / (wnumber * c.MTC)
    # Convert Elow from eV to cm-1:
    elow[:] = elow * c.eV2kayser

    ut.lrprint(verbose, "Done.\n")
    data.close()
    return wlength, gf, elow, isoID

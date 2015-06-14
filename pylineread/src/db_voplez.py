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
# Copyright (C) 2015 University of Central Florida.  All rights reserved.
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
    1000K to 7000K.

    Modification History:
    ---------------------
    2015-06-14  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    """
    # FINDME: Implement me
    Temp = np.array([1000.0, 7000.0])
    PF   = np.array([[1.0, 10.0]])
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
    nlines   = data.tell() / self.recsize

    # Get database limiting wavenumbers:
    firstwl = self.readwl(data,        0)  # Highest wn in DB
    lastwl  = self.readwl(data, nlines-1)  # Lowest wn in DB

    # Find the record index for iwl:
    if iwl > firstwl:
      irec_init = self.binsearch(data, iwl, 0,         nlines, 1)
    else:
      irec_init = 0
    if fwl < lastwl:
      irec_fin  = self.binsearch(data, fwl, irec_init, nlines, 0)
    else:
      irec_fin = nlines-1

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
    interval = float((irec_fin - irec_init)/20)  # Check-point interval

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
        pos = float(data.tell()/self.recsize)
        if (pos/interval)%1 == 0.0:
          ut.lrprint(verbose-1, "checkpoint %d/20..."%chk)
          chk += 1
          ut.lrprint(verbose-3, "Wavenumber (cm-1): {:.2f}, Elow (eV): {:.3f}, "
                     "gf: {:.4e}".format(wnumber[i], elow[i], gf[i]))
          ut.lrprint(verbose-2, "Wavelength (um): {:.3f},  IsoID: {:d},  "
                     "Elow (cm-1): {:.5e}, gf: {:.5e}".format(
                       1.0/(wnumber[i]*c.MTC), 1, elow[i]*c.eV2kayser, gf[i]))
      i += 1

    # Convert wavelength to TLI format (microns):
    wlength[:] = 1.0 / (wnumber * c.MTC)
    # Convert Elow from eV to cm-1:
    elow[:] = elow * c.eV2kayser
    # Isotopic index:
    isoID[:]   = np.ones(nread, int)

    ut.lrprint(verbose, "Done.\n")
    data.close()
    return wlength, gf, elow, isoID

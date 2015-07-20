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

import sys, os
import numpy as np
import subprocess as sp
import scipy.constants as sc

import utils as ut
import constants as c
from driver import dbdriver

#db_hitran.py dir:
DBHdir = os.path.dirname(os.path.realpath(__file__))

# CTIPS
sys.path.append(DBHdir + "/ctips/lib")
import ctips as ct

class hitran(dbdriver):
  def __init__(self, dbfile, pffile):
    """
    Initialize the basic database info.

    Parameters:
    -----------
    dbfile: String
       File with the Database info as given from HITRAN.
    pffile: String
       File with the partition function.

    Modification History:
    ---------------------
    2014-08-01  patricio  Added documentation.
    """
    super(hitran, self).__init__(dbfile, pffile)

    self.recsize   =   0 # Record length (will be set in self.dbread())
    self.recwnpos  =   3 # Wavenumber     position in record
    self.recisopos =   2 # Isotope        position in record
    self.reclinpos =  15 # Line intensity position in record
    self.recApos   =  25 # Einstein coef  position in record
    self.recairpos =  35 # Air broadening position in record
    self.recelpos  =  45 # Low Energy     position in record
    self.recg2pos  = 155 # Low stat weight position in record
    self.recmollen =   2 # Molecule   record length
    self.recwnlen  =  12 # Wavenumber record length
    self.reclinend =  25 # Line intensity end position
    self.recelend  =  55 # Low Energy     end position
    self.T0       = 296.0 # K

    self.molID = self.getMolec()
    # Get info from HITRAN configuration file:
    self.molecule, self.isotopes, self.mass, self.isoratio, self.gi = \
                                                    self.getHITinfo(self.molID)
    # Database name: 
    self.name = "HITRAN " + self.molecule


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

    Modification History:
    ---------------------
    2014-03-05  patricio  Initial implementation, based on Madison's
                          code.          pcubillos@fulbrightmail.org
    2014-03-10  patricio  Updated dbtype to match command-line-argument
                          sytax.  Updated HITRAN data type.
    """
    # Set pointer at required wavenumber record:
    dbfile.seek(irec*self.recsize + self.recwnpos)
    # Read:
    wavenumber = dbfile.read(self.recwnlen)
    # Convert to float:
    rec_wl = float(wavenumber)

    return rec_wl


  def getpf(self, verbose):
    """                                    
    Calculate the partition function for a grid of temperatures for
    HITRAN molecules.

    Notes:
    ------
    This function is a wrapper of the CTIPS library written by
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

    Modification History:
    ---------------------
    2012-11-22  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2014-03-10  patricio  Adapted previous code.
    2015-06-29  aj-foster Swapped Fortran for C implementation of TIPS.
    """

    # Get molecule ID as  integer:
    molID = int(self.molID)
 
    # Get Number of isotopes:
    Niso = len(self.isotopes)

    # Array of temperatures:
    Temp = np.arange(70.0, 3000.1, 10)
    # TIPS range of temperature is: [70K, 3000K]
    Ntemp = len(Temp)

    # Output array for table of Temperature and PF values: 
    PF = np.zeros((Niso, Ntemp), np.double)
    ut.lrprint(verbose-14, "Temperature sample:\n" + str(Temp))
    ut.lrprint(verbose-5, "Ntemp: %d, Niso: %d"%(Ntemp, Niso))

    for i in np.arange(Niso):

      # Molecule ID, isotope ID, and temperature arrays:
      mol = np.repeat(molID, len(Temp))
      iso = np.repeat(int(self.isotopes[i]), len(Temp))

      # Call CTIPS with the given arrays:
      PF[i] = ct.tips(mol, iso, Temp)
      PF[i] /= self.gi[i]

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


  def getHITinfo(self, molID):
    """
    Get HITRAN info from configuration file.

    Parameters:
    -----------
    molID: Integer
       Molecule ID as given in HITRAN database

    Returns:
    --------
    molname:  Molecule's name
    isotopes: Isotopes names
    mass:     Isotopes mass
    isoratio: Isotopic abundance ratio
    gi:       State-independent statistical weight

    """
    # Read HITRAN configuration file from inputs folder:
    hfile = open(DBHdir + '/../inputs/hitran.dat', 'r')
    lines = hfile.readlines()
    hfile.close()
  
    isotopes = []
    mass     = []
    isoratio = []
    gi       = []
 
    # Get values for our molecule: 
    for i in np.arange(len(lines)):
      if lines[i][0:2] == molID:
        line = lines[i].split()
        molname  = line[1]
        gi.      append(  int(line[3]))
        isotopes.append(      line[2] )
        isoratio.append(float(line[4]))
        mass.    append(float(line[5]))

    return molname, isotopes, mass, isoratio, gi    

 
  def PFzero(self):
    """
    Calculate the partition function for the isotopes at T0 = 296K.

    Notes:
    ------
    This function is a wrapper of the CTIPS library written by
    Patricio Cubillos, which is itself a C implementation of the
    Fortran TIPS routine by R. R. Gamache.
    The range of temperatures is limited to: 70K -- 3000K.

    Parameters:
    -----------
    pffile: String
       Partition Function filename.
    molID: Integer
       Molecule ID as given in HITRAN database.

    Returns:
    --------
    PFzero: 1D ndarray
       An array (N isotopes) with the partition function evaluated at
       T0 = 296 K for each isotope.

    Modification History:
    ---------------------
    2012-03-10  patricio  Initial implementation.  pcubillos@fulbrightmail.org
    2015-06-29  aj-foster Swapped Fortran for C implementation of TIPS.
    """

    # Get molecule ID:
    molID = int(self.molID)

    # Get number of isotopes:
    Niso = len(self.isotopes)

    # Output array for table of Temperature and PF values:
    PFzero = np.zeros(Niso, np.double)

    # Molecule ID, isotope ID, and temperature arrays:
    mol = np.repeat(molID, Niso)
    iso = np.asarray(self.isotopes, dtype=int)
    temp = np.repeat(self.T0, Niso)

    PFzero = ct.tips(mol, iso, temp)
    return PFzero / self.gi


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
    data = open(self.dbfile, "r")
    # Read first line to get the record size:
    data.seek(0)
    line = data.readline()
    self.recsize = len(line)
    # Get Total number of transitions in file:
    data.seek(0, 2)
    nlines   = data.tell() / self.recsize
    # Get Molecule ID:
    molID = int(self.molID)

    # Rewrite wavelength limits as given in HITRAN file:
    fwn = 1.0 / (iwl * c.MTC)  # Microns to centimeters
    iwn = 1.0 / (fwl * c.MTC)  # Microns to centimeters

    # Find the record index for iwn and fwn:
    irec_init = self.binsearch(data, fwn, 0, nlines-1,  1)
    irec_fin  = self.binsearch(data, iwn, 0, irec_init, 0)

    # Allocates arrays for values to extract:
    wlength = np.zeros(irec_init-irec_fin+1, np.double)
    gf      = np.zeros(irec_init-irec_fin+1, np.double)
    elow    = np.zeros(irec_init-irec_fin+1, np.double)
    isoID   = np.zeros(irec_init-irec_fin+1,       int)
    # Einstein A coefficient:
    A21     = np.zeros(irec_init-irec_fin+1, np.double)
    # Lower statistical weight:
    g2      = np.zeros(irec_init-irec_fin+1, np.double)
    # Line intensity, used to calculate gf:
    S0      = np.zeros(irec_init-irec_fin+1, np.double)
    # Wavenumber:
    wnumber = np.zeros(irec_init-irec_fin+1, np.double)

    # Get the partition function at T0:
    PFzero = self.PFzero()

    i = 0  # Stored record index
    chk = 0
    interval = float((irec_fin - irec_init)/20)  # Check-point interval
    while irec_init - i >= irec_fin:
      data.seek( (irec_init-i) * self.recsize )
      # Read in wavenumber
      line = data.read(self.recsize)
      # Get the isotope index:
      isoID[i]   = float(line[self.recisopos:self.recwnpos ])
      # Get the wavenumber:
      wnumber[i] = float(line[self.recwnpos: self.reclinpos])
      # Get the line intensity:
      S0[i]      = float(line[self.reclinpos:self.reclinend])
      # Get Elow:
      elow[i]    = float(line[self.recelpos: self.recelend ])
      A21[i]     = float(line[self.recApos:  self.recairpos])
      g2[i]      = float(line[self.recg2pos: self.recsize  ])
      # Print a checkpoint statement every 1/20th interval
      if verbose > 1:
        pos = float(data.tell()/self.recsize)
        if (pos/interval)%1 == 0.0:
          ut.lrprint(verbose-1, "checkpoint %d/20..."%chk)
          chk += 1
          ut.lrprint(verbose-3, "Wavenumber: %s, S0: %s, Elow: "
                                "%s"%(wnumber[i], S0[i], elow[i]))
          ut.lrprint(verbose-2, "Wavelength: %.3f, IsoID: %d, Elow: %.3f"%(
                                1.0/(wnumber[i]*c.MTC), isoID[i], elow[i]))
      i += 1


    # Calculate the wavelength in microns:
    wlength[:] = 1.0 / (wnumber * c.MTC)
    # Set isotopic index to start counting from 0:
    isoID -= 1
    isoID[np.where(isoID < 0)] = 9 # 10th isotope had index 0 --> 10-1=9
    # Calculate gf:
    Ia = np.asarray(self.isoratio)
    gf = (S0 * c.C1 * PFzero[isoID] /  Ia[isoID] *
             np.exp( c.C2 * elow / self.T0)  /
                 (1.0-np.exp(-c.C2 * wnumber    / self.T0)) )
    # Alternative way:
    gf2 = A21 * g2 * c.C1 / (8.0 * np.pi * sc.c * 100.0) / wnumber**2.0

    # FINDME: Delete me when gf-vs-gf2 issue solved.
    # print(gf)
    # print(gf2)
    # print((gf/gf2)[:20])
    # print(isoID[:20])
    # print(np.amax(gf/gf2), np.amin(gf/gf2))
    # import matplotlib.pyplot as plt
    # plt.figure(1)
    # plt.plot(isoID, gf/gf2, "b,")
    # plt.ylim(0,4)
    # plt.xlim(-0.5, 5.5)
    # plt.savefig("gf.png")
    # plt.figure(2)
    # print(np.shape(self.gi), np.shape(PFzero))
    # ggi = np.asarray(self.gi)[isoID]
    # plt.plot(ggi, gf2/gf, "b,")
    # plt.plot([0,12], [0,12], "r")
    # plt.ylim(0,15)
    # plt.xlim(-0.5, 12.5)
    # plt.savefig("gf2.png")

    ut.lrprint(verbose-14, "GF values: " + str(gf))
    data.close()
    return wlength, gf2, elow, isoID

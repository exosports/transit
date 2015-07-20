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

#import db_pands  as ps
#import db_hitran as hit
#import db_tioschwenke as ts

class dbdriver(object):
  def __init__(self, dbfile, pffile):
    self.dbfile = dbfile
    self.pffile = pffile

  def getpf(self, verbose):
    """
      Calculate partition function for specific database type.

      Parameters:
      -----------
      dbtype: String
    """
    pass

  def dbread(self, iwl, fwl, verbose):
    """
      Read linelist values for specific database type.

      Parameters:
      -----------
      dbtype: String
    """
    pass

  def readwl(self, dbfile, irec):
    """
      Read the wavelength parameter as given in the database
    """
    pass

  def binsearch(self, dbfile, wavelength, ilo, ihi, searchup=True):
    """
    Do a binary search of record position in File dbfile that has wavelength iwl
 
    Parameters:
    -----------
    dbfile: File object
       File where to search.
    wavelength: Scalar
       Target wavelength (as given in each specific database).
    ilo: Integer
       Index of smallest record to search.
    ihi: Integer
       Index of largerst record to search.
    dbtype: String
       Database type: ['ps', 'hit']
    searchup: Boolean
       Search up (True) or down (False) the records for duplicate results
       after the binary search.
 
    Returns:
    --------
    Index of record with
 
    Modification History:
    ---------------------
    2013        madison   Initial implementation.
    2014-03-05  patricio  Added documentation, updated Madison's code, and
                          included searchup parameter
                                                    pcubillos@fulbrightmail.org
    """
    # imin and imax are the fixed boundaries where to search:
    imin, imax = ilo, ihi

    # Wavelength of record:
    rec_wl = 0
    irec   = 0
 
    # Start binary search:
    while ihi - ilo > 1:
      # Middle record index:
      irec = (ihi + ilo)/2
 
      # Read wavelength, depending on linelist format:
      rec_wl = self.readwl(dbfile, irec)
      # Update search limits:
      if rec_wl > wavelength:
        ihi = irec
      else:
        ilo = irec
 
    # Start linear search:
    if searchup:
      irec = ilo
    else:
      irec = ihi

    # Check value of contiguous entries:
    icheck  = irec  # Record index to test
    bounded = True  # Record's wavelength is still within boundaries

    while bounded:
      # Update irec:
      irec = icheck
      # Check index boundaries:
      if irec == imin or irec == imax:
        break
      # Check wavelength boundaries:
      if searchup:
        icheck += 1
        bounded = self.readwl(dbfile, icheck) < wavelength
      else:
        icheck -= 1
        bounded = self.readwl(dbfile, icheck) > wavelength


    # Move file pointer to begining:
    dbfile.seek(0)
 
    return irec


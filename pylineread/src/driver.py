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
    # Wavelength of record:
    rec_wl = 0
 
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
 
    # Search up/or down if there are multiple records with the same wavelength:
    while (self.readwl(dbfile, irec) == rec_wl):
      if searchup:
        irec += 1
      else:
        irec -= 1
 
    # Move file pointer to begining:
    dbfile.seek(0)
 
    #print("Found wavelength %d at position %d."%(rec_wl, ihi))
    return irec


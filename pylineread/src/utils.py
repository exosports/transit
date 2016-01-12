# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# Transit is under an open-source, reproducible-research license (see LICENSE).

import struct, sys
import numpy as np
import constants as c


def lrprint(threshold, text):
  """
  Lineread conditional print. Print text only if threshold is larger
  than zero.

  Parameters:
  -----------
  thereshold: Integer
     Conditional threshold to print.
  text: String
     Text to print.

  Example:
  --------
  >>>import utils as ut 
  >>>ut.lrprint(1, "Print this text")
  Print this text
  >>>ut.lrprint(0, "Don't print this one")

  Modification History:
  ---------------------
  2014-03-04  patricio  Initial version.  pcubillos@fulbrightmail.org
  """
  if threshold > 0:
    print(text)


def printexit(text):
  """
  Print message and exit the code execution.

  Parameters:
  -----------
  text: String
    Text to print.

  Modification History:
  ---------------------
  2014-03-06  patricio  Initial implementation.  pcubillos@fulbrightmail.org
  """
  print(text)
  sys.exit(0)

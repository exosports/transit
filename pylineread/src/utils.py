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

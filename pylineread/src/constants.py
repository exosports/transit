# Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
# Transit is under an open-source, reproducible-research license (see LICENSE).

import numpy as np
import scipy.constants as sc

"""
Set of constants for pylinereader package.

Modification History:
---------------------
2014-09-29  patricio  Fixed wrong equations fot C1.  Fixed wrong units
                      conversion factor in C2. Updated Version.
"""
# General Constants:
MTC  = 1e-4         # Microns to cm     (MTC = um/cm)
NTC  = 1e-7         # Nanometers to cm  (NTC = nm/cm)
e    = 4.803205e-10 # Elementary charge in statcoulombs (from Wolfram Alpha)
C1   = 4 * sc.epsilon_0 * sc.m_e * sc.c**2 / sc.e**2 * 0.01  # cm^-1
C2   = sc.h * sc.c / sc.k * 100.0        # cm / Kelvin units

# Version Constants:
TLI_VERSION = 6  # TLI version
LR_VERSION  = 0  # Lineread version
LR_REVISION = 0  # Lineread revision

# H2O Isotopic Physical Constants:
PS_CROSS = np.pi * ((3.2e-08)/2)**2.0 # Collision cross section

# Convert from eV to cm-1 (kayser):
# planck   = 6.62620e-34  # Planck constant [J * s]
# lumiere  = 2.997925e10  # speed of light  [cm / s]
# electron = 1.602192e-19 # elementary charge [Coulomb]
# kayser2eV = planck * lumiere / electron
eV2kayser = 8065.49179

import numpy as np
import scipy.constants as sc

# pylineread Constants:

# General Constants:
MTC  = 1e-4         # Microns to cm
NTC  = 1e-7         # Nanometers to cm
e    = 4.803205e-10 # Elementary charge in statcoulombs (from Wolfram Alpha)
C1   = sc.e * sc.c**2 / (sc.pi * e**2) # cm^-1
C2   = sc.h * sc.c / sc.k * 0.01       # cm / Kelvin units

# Version Constants:
TLI_VERSION = 4  # TLI version
LR_VERSION  = 4  # Lineread version
LR_REVISION = 0  # Lineread revision

# H2O Isotopic Physical Constants:
PS_CROSS = np.pi * ((3.2e-08)/2)**2.0 # Collision cross section


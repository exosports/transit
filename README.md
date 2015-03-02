### Transit
>A radiative-transfer code for planetary atmospheres  

Transit calculates the transmission or emission flux spectrum of a planetary atmosphere with application to extrasolar-planet transit and eclipse observations, respectively. Transit computes the spectra by solving for the one-dimensional line-by-line radiative-transfer equation for a plane-parallel atmospheric model.  

<a href="https://www.youtube.com/watch?v=-GHBFHyeI14" target="_blank"><img src="https://github.com/pcubillos/transit/blob/master/doc/ScreenShot.jpg" 
alt=""border="10" /></a>

### Table of Contents:
* [Team Members](#team-members)
* [Code Description](#code-description)
* [Installation and System Requirements](#installation-and-system-requirements)
* [Be Kind](#be-kind)
* [License](#license)

### Team Members:
* [Patricio Cubillos](https://github.com/pcubillos/) (author, UCF) <pcubillos@fulbrightmail.org>
* Jasmina Blecic (UCF)
* Joe Harrington (UCF)
* Patricio Rojo (U. de Chile)
* Madison Stemm (UCF)
* Andrew Foster (UCF)

### Code Description:
Get the Transit user's manual [here](doc/transitUM.pdf).

### Installation:
To obtain the Transit code download the latest stable version from the releases page (TBD). Alternatively, clone the repository to your local machine with the following terminal command:  
**git clone https://github.com/pcubillos/transit mytransit/**  

To compile the pu and transit modules execute ‘make’ on the respective folders:  
**cd mytransit/pu**  
**make**  
**cd ../mytransit/transit**  
**make**  

To compile the pylineread fortran code:  
**cd mytransit/pylineread/src/fortran**  
**make**  

To remove the program binaries, execute (from the respective directories):  
**make clean**  


### Be Kind:
Please reference these papers if you found this module useful for your research:  
  [Cubillos et al. 2015: The Bayesian Atmospheric Radiative-Transifer Code for Exoplanet Modeling I](), in preparation.  
  [Blecic et al. 2015: The Bayesian Atmospheric Radiative-Transifer Code for Exoplanet Modeling II](), in preparation.  
  [Harrington et al. 2015: The Bayesian Atmospheric Radiative-Transifer Code for Exoplanet Modeling III](), in preparation.  
Thanks!


### License:

Transit, a code to solve for the radiative-transfer equation for planetary atmospheres.  

This project was completed with the support of the NASA Planetary Atmospheres Program, grant NNX12AI69G, held by Principal Investigator Joseph Harrington. Principal developers in- cluded graduate students Patricio E. Cubillos and Jasmina Blecic, programmer Madison Stemm, and undergraduate Andrew S. D. Foster. The included ’transit’ radiative transfer code is based on an earlier program of the same name written by Patricio Rojo (Univ. de Chile, Santiago) when he was a graduate student at Cornell University under Joseph Harrington.  

Copyright (C) 2015 University of Central Florida. All rights reserved.  

This is a test version only, and may not be redistributed to any third party. Please refer such requests to us. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

Our intent is to release this software under an open-source, reproducible-research license, once the code is mature and the first research paper describing the code has been accepted for publica- tion in a peer-reviewed journal. We are committed to development in the open, and have posted this code on github.com so that others can test it and give us feedback. However, until its first pub- lication and first stable release, we do not permit others to redistribute the code in either original or modified form, nor to publish work based in whole or in part on the output of this code. By downloading, running, or modifying this code, you agree to these conditions. We do encourage sharing any modifications with us and discussing them openly.  

We welcome your feedback, but do not guarantee support. Please send feedback or inquiries to:  
Patricio Cubillos <pcubillos[at]fulbrightmail.org>  
Jasmina Blecic <jasmina[at]physics.ucf.edu>  
Joseph Harrington <jh[at]physics.ucf.edu>  

or alternatively,  
Joseph Harrington, Patricio Cubillos, and Jasmina Blecic  
UCF PSB 441  
4111 Libra Drive  
Orlando, FL 32816-2385  
USA  

Thank you for using transit!

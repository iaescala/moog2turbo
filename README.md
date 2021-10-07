# moog2turbo
Code to convert MOOG linelists into a format compatible with Turbospectrum without sourcing additional information from e.g. VALD

Pulls total angular momentum quantum numbers from Barklem et al. 2000, 2005
(same source as Barklem.dat in MOOG according to Gammabark.f) and damping parameters from Barklem.dat

Note that the formatting in the convert_moog_linelist() function is 
heavily based on Alex Ji's turbopy/linelists.py

It is recommended that you merge your standard MOOG linelist and strong linelist prior to conversion

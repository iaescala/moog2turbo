# moog2turbo
Code to convert MOOG linelists into a format compatible with Turbospectrum without sourcing additional information from e.g. VALD

Pulls total angular momentum quantum numbers from Barklem et al. 2000, 2005
(same source as Barklem.dat in MOOG according to Gammabark.f) and damping parameters from Barklem.dat

Uses Barklem VdW parameter where available, or VdW parameter
if specified in MOOG linelist following dampingopt = 1 (according to Damping.f)

Includes optional keyword fdamp_flag = False (default), where
fdamp_flag = True uses Unsold correction factors from vald3line-BPz-freeformat.f
in Turbospectrum2019 where VdW is otherwise zero

Note that the formatting in the convert_moog_linelist() function is 
heavily based on Alex P. Ji's turbopy/linelists.py

It is recommended that you merge your standard MOOG linelist and strong linelist prior to conversion

For now, if you use this code in your research, please cite Escala et al. 2023, in prep, and link to this repository.

"""
Author: Ivanna Escala (Carnegie Observatories)

Program to convert MOOG linelists to a Turbospectrum compatible format
without sourcing additional information from e.g. VALD

Pulls total angular momentum quantum numbers from Barklem et al. 2000, 2005
(same source as Barklem.dat in MOOG according to Gammabark.f)

Uses Barklem VdW parameter where available, or VdW parameter
if specified in MOOG linelist following dampingopt = 1 (according to Damping.f).
Includes optional keyword fdamp_flag = False (default), where
fdamp_flag = True uses Unsold correction factors from vald3line-BPz-freeformat.f
in Turbospectrum2019 where VdW is otherwise zero

Note that the formatting in the convert_moog_linelist() function is
heavily based on Alex P. Ji's turbopy/linelists.py
"""

import numpy as np
import os
from astropy.table import Table
from astropy.io import fits

def get_elem(species, ion):
    _all_elems = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
                  "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
                  "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
                  "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
                  "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
                  "Fr","Ra","Ac","Th","Pa","U"]
    _ion_state = ["I", "II"]
    _swap = ["NH", "NC"]

    intspecies = int(species)
    if intspecies < 100: #if an atom
        return f"{_all_elems[intspecies-1]} {_ion_state[ion-1]}"
    else: #if a diatomic molecule
        strspecies = str(intspecies)
        species_lo, species_hi = strspecies[0], strspecies[1:]
        Z_lo, Z_hi = _all_elems[int(species_lo)-1], _all_elems[int(species_hi)-1]
        if Z_lo != Z_hi:
            if f"{Z_hi}{Z_lo}" in _swap:
                Z_hi, Z_lo = Z_lo, Z_hi
            return f"{Z_hi}{Z_lo} {_ion_state[ion-1]}"
        else:
            return f"{Z_lo}2 {_ion_state[ion-1]}"


def load_barklem(root=os.getcwd()):
    """"
    Read in the Barklem data
    """

    data = []
    with open(os.path.join(root,"Barklem.dat"), "r") as f:
        lines = f.readlines()

        for line in lines:
            l_strip = line.strip().split()

            if len(l_strip) == 5:
                wavebk, idbk, gammabk, alphabk, gammarad = l_strip

            if len(l_strip) == 4:
                wave, idbk, gammabk, alphabk = l_strip
                gammarad = -99.99

            data.append((float(wavebk), float(idbk), float(gammabk), float(alphabk), float(gammarad)))

    #In Barklem.dat, gammabk is the base 10 log of the FWHM per perturber at 10,000 K
    #e.g., log(2w/N) following Barklem & O'Mara 1998
    #alphabk is the velocity parameter
    cols = ["wavebk", "idbk", "gammabk", "alphabk", "gammarad"]
    tabbk = Table(rows=data, names=cols)

    return tabbk

def barklem(tab, root=os.getcwd()):
    """
    Load the data from Barklem.dat in MOOG and match the data
    to the linelist following Gammabark.f
    """

    tabbk = load_barklem(root=root)
    tabbk = get_Jnumber(tabbk, root=root)

    #Identify the Barklem list positions of the wavelength limits of the input linelist
    wavemin = np.min(tab["wave"]); wavemax = np.max(tab["wave"])

    nummin = np.where((wavemin - tabbk["wavebk"]) < 1.)[0][0]
    nummax = np.where((tabbk["wavebk"] - wavemax > 1.))[0][0]

    #Search for Barklem data

    #Initalize arrays -- moog uses -1 but uses FWHM
    #If using log(FWHM) zero seems fine for Turbospec

    gambark = np.zeros(tab["wave"].size)
    #alpbark = np.full(tab["wave"].size, -1.)
    gamrad = np.zeros(tab["wave"].size)
    jhi = np.zeros(tab["wave"].size)

    watom = tab["species_moog"].astype(float) < 100.

    #species_moog = tab["species_moog"][watom].astype(float)
    #species_bk = tabbk["idbk"][nummin:nummax+1]

    #This weirdness is to account for isotopes in the species designation
    species_moog = 10. * tab["species_moog"][watom].astype(float) + 0.0001
    species_moog = species_moog.astype(int)
    species_bk = 10. * tabbk["idbk"][nummin:nummax+1] + 0.001
    species_bk = species_bk.astype(int)

    #Do just for the atoms
    for j in range(tab["wave"][watom].size):

        waverror = (tab["wave"][watom][j] - tabbk["wavebk"][nummin:nummax+1])/tabbk["wavebk"][nummin:nummax+1]
        #wavediff = tab["wave"][watom][j] - tabbk["wavebk"][nummin:nummax+1]

        ww = np.where( (np.abs(waverror) < 5.e-6) & (species_bk == species_moog[j]) )[0]
        #ww = np.where( (wavediff == 0.) & (species_bk == species_moog[watom][j]) )[0]

        if len(ww) == 0:

            #Use the approximation from line 166 in Damping.f in MOOG in the case where
            #GAMRAD = 0, to circumvent the approximation used in line 999 of bsyn.f
            #in Turbospectrum, which depends on the GU parameter
            #P.S. we don't know where this constant comes from
            gamrad[j] = 2.223e15/tab["wave"][watom][j]**2

            continue
        elif len(ww) > 1:
            ww = ww[np.argmin(waverror[ww])]
        else:
            ww = ww[0]

        #X < 0 where X = a Van der Waals damping parameter
        #implies the parameter is log(FWHM) -- as is used by VALD2
        #in constrast MOOG uses FWHM (not the log)
        gambark[j] = tabbk["gammabk"][nummin:nummax+1][ww]

        #moog uses the temperature dependence exponent here (1-alpha)/2.
        #alpbark[j] = (1. - tabbk["alphabk"][k])/2

        gamrad[j] = tabbk["gammarad"][nummin:nummax+1][ww]

        jhi[j] = tabbk["jhi"][nummin:nummax+1][ww]

    #Molecules
    wzero = gamrad == 0.
    gamrad[wzero] = 2.223e15/tab["wave"][wzero]**2.


    tab["fdamp"] = gambark
    tab["raddamp"] = gamrad
    tab["jhi"] = jhi
    tab["gu"] = 2. * jhi + 1.

    return tab


def get_Jnumber(tab_bk, root=os.getcwd()):

    #Read in the two Barklem data sets
    tab00 = Table.read(os.path.join(root, "Barklem_2000.fits"))
    tab05 = Table.read(os.path.join(root, "Barklem_2005.fits"))

    #Rename the relevant columns and join the tables
    tab05['Ion'].name = 'spec'
    tab05['Jupp'].name = 'uppJ'

    tab_j = Table()
    tab_j['lambda'] = np.concatenate((tab00['lambda'], tab05['lambda']))
    tab_j['spec'] = np.concatenate((tab00['spec'], tab05['spec']))
    tab_j['uppJ'] = np.concatenate((tab00['uppJ'], tab05['uppJ']))

    #Some re-formatting the species ID in the Barklem data from Vizier
    spec_str_len = np.array([len(str(x)) for x in tab_j['spec']])
    wion = spec_str_len == 5
    tab_j['spec'][wion] = np.round(tab_j['spec'][wion]+0.09,decimals=1)

    #Perform the matching for the J number
    jhi = np.zeros(tab_bk["wavebk"].size)

    for i in range(tab_bk['wavebk'].size):

        waverror = (tab_j["lambda"] - tab_bk["wavebk"][i])/tab_bk["wavebk"][i]

        ww = np.where( (np.abs(waverror) < 5.e-6) &\
                      (tab_j['spec'] == tab_bk['idbk'][i] ) )[0]

        if len(ww) == 0:
            continue
        elif len(ww) > 1:
            ww = ww[np.argmin(waverror[ww])]
        else:
            ww = ww[0]

        jhi[i] = tab_j['uppJ'][ww]

    tab_bk["jhi"] = jhi

    return tab_bk


def load_moog_list(filename, skipheader=0, root=os.getcwd()):
    """"
    Load in the MOOG format linelist
    """

    with open(os.path.join(root, filename), "r") as f:

        lines = f.readlines()
        lines = lines[skipheader:]

        data = []
        for line in lines:

            l_strip = line.strip().split()

            if len(l_strip) == 4: #an atom

                wave, species, expot, loggf = l_strip
                d0 = np.nan
                vdw = np.nan

            if len(l_strip) == 5: #a molecule, OR an atom with an VdW factor specified (dampingopt = 1)

                wave, species, expot, loggf, vdw_or_d0 = l_strip

                if float(species) > 100.: # a molecule
                    d0 = vdw_or_d0
                    vdw = np.nan

                else: #an atom with unsold factor specified
                    d0 = np.nan
                    vdw = vdw_or_d0

            ion_iso = species.split('.')[-1]
            if len(ion_iso) > 1:
                ion = int(ion_iso[0])+1
            else:
                ion = int(ion_iso)+1

            data.append((float(wave), species, ion, float(expot), float(loggf), float(vdw), float(d0)))

    cols = ["wave", "species_moog", "ion", "expot", "loggf", "vdw", "d0"]
    tab = Table(rows=data, names=cols)

    return tab


def convert_species_format(tab):

    #Get the species in a format for Turbospectrum
    watom = tab["species_moog"].astype(float) < 100.
    species = np.zeros(tab["wave"].size)

    #if no isotope specified, floor the atomic species (e.g. 25.1 --> 25.0)
    species_atom = np.floor(tab["species_moog"][watom].astype(float))
    #if isotope specified, then keep the same
    wiso = np.array([len(x.split('.')[-1]) > 1 for x in tab["species_moog"][watom]])
    species_atom[wiso] = tab["species_moog"][watom][wiso].astype(float)

    #Assign the atomic species to the new TS species array
    species[watom] = species_atom

    #Molecular species
    species_molec = tab["species_moog"][~watom]

    species_molec_ts = []
    for spec in species_molec:
        ion_iso = spec.split('.')[-1]
        if len(ion_iso) > 1:
            els = spec.split('.')[0]
            iso1, iso2 = ion_iso[1:3], ion_iso[3:5]
            if iso1 == '01': iso1 = '00'
            species_molec_ts.append( f"{els}.0{iso1}0{iso2}" )
        else:
            species_molec_ts.append( spec )

    #Assign the molecular species to the new TS species array
    species[~watom] = np.array(species_molec_ts).astype(float)
    tab["species"] = species

    #Formatting
    tab["species"][watom] = [f"{x:6.3f}" for x in tab["species"][watom]]
    tab["species"][~watom] = [f"{x:10.6f}" for x in tab["species"][~watom]]

    return tab


def convert_moog_linelist(filename, skipheader=0, outfilename=None, root=os.getcwd(),
                          fdamp_flag=False):

    """
    Note: It is recommended that you merge your standard MOOG linelist and strong MOOG linelist
    into a single list prior to performing this conversion
    """

    tab = load_moog_list(filename, skipheader=skipheader, root=root)

    tab = convert_species_format(tab)

    #Sort the table according to species, and then wavelength within a given species
    tab["sortspecies"] = tab["species"].astype(float) + 0.0000001 * tab["ion"]
    tab.sort(["sortspecies", "wave"])

    ## Parameters that require sourcing from elsewhere: fdamp, gu, raddamp
    tabbk = barklem(tab, root=root)

    tab["gu"] = tabbk["gu"]
    tab["fdamp"] = tabbk["fdamp"]
    tab["raddamp"] = tabbk["raddamp"]

    ##If a corrective factor is provided in the MOOG linelist, use this instead
    ## (assuming dampingopt = 1)
    ## Corrective factors for Unsold approx in TS for 0 < VdW < 20
    wvdw = ~np.isnan(tab["vdw"]) & (tab["vdw"] > 0.) & (tab["vdw"] < 20.)
    tab["fdamp"][wvdw] = tab["vdw"][wvdw]

    #If this flags is True, then use the Unsold correction factors
    #specified in vald3line-BPz-freeformat.f in Turbospectrum2019 package
    #Factors from BDP (A&A 275,101) Edvarsson et al. 1993
    if fdamp_flag:

        fdampdict1 = {11: 2.0, 14: 1.3, 20: 1.8, 26: 1.4} # neutral damping, the rest are 2.5
        fdampdict2 = {20: 1.4, 38: 1.8, 56: 3.0} # ionized damping, the rest are 2.5
        unsold_ts_factor = 2.5 #everything else

        intspecies = tab["species"].astype(int)

        for key in fdampdict1.keys():
            wkey = (tab["fdamp"] == 0.) & (intspecies == key) & (tab["ion"] == 1)
            tab["fdamp"][wkey] = fdampdict1[key]

        for key in fdampdict2.keys():
            wkey = (tab["fdamp"] == 0.) & (intspecies == key) & (tab["ion"] == 2)
            tab["fdamp"][wkey] = fdampdict2[key]

        w_still_zero = tab["fdamp"] == 0.
        tab["fdamp"][w_still_zero] = unsold_ts_factor

    #Remove bad entries that TS doesn't like
    iibad = (tab["sortspecies"] < 3) | (tab["ion"] > 2) | (tab["ion"] < 1) | \
            (tab["expot"] > 15.) | (tab["loggf"] < -10.) | (tab["loggf"] > 100.)
    tab = tab[~iibad]

    if outfilename is not None:

        if os.path.exists(os.path.join(root, outfilename)):
            os.remove(os.path.join(root, outfilename))

        with open(os.path.join(root, outfilename), "w") as f:

            def write(x):
                f.write(f"{x}\n")

            for sortspecies in np.unique(tab["sortspecies"]):

                t = tab[tab["sortspecies"]==sortspecies]
                N = len(t)
                write(f"' {t[0]['species']}'      {t[0]['ion']}       {N}")

                elem_str = get_elem(t[0]['species'], t[0]['ion'])
                write(f"' {elem_str}'")

                for row in t:
                    #fmt = "{wave:10.3f} {expot:6.3f} {loggf:6.3f} {fdamp:8.3f} {gu:6.1f} {raddmp:9.2e} '{levlo}' '{levup}'"
                    write(f"{row['wave']:10.3f} {row['expot']:6.3f} {row['loggf']:6.3f} {row['fdamp']:8.3f} "+\
                          f"{row['gu']:6.1f} {row['raddamp']:9.2e} 'X' 'X'")

    return tab

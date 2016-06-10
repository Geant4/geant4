#$Id: emcalculator.py 66892 2013-01-17 10:57:59Z gunter $
"""
# ==================================================================
#   Python module (Python3)
#
#   Calculation of photon cross section and stopping power for
#   chared particles
#
#                                              Q, 2005
# ==================================================================
"""
from Geant4 import *

# ==================================================================
# public symbols
# ==================================================================
__all__ = [ 'CalculatePhotonCrossSection', 'CalculateDEDX' ]

# ==================================================================
# Photon Cross Section
# ==================================================================
def CalculatePhotonCrossSection(mat, elist, verbose=0,
                                plist=["compt", "", "phot", "conv"]):
  """
  Calculate photon cross section for a given material and
  a list of energy, returing a list of cross sections for
  the components of "Copmton scattering", "rayleigh scattering",
  "photoelectric effect", "pair creation" and total one.

  Arguments:
    mat:      material name
    elist:    list of energy
    verbose:  verbose level [0]
    plist:    list of process name
              (compton/rayleigh/photoelectic/conversion) [StandardEM set]
    
  Keys of index:
    "compt":     Compton Scattering
    "rayleigh":  Rayleigh Scattering
    "phot" :     photoelectric effect
    "conv" :     pair Creation
    "tot"  :     total

  Example:
    xsec_list= CalculatePhotonCrossSection(...)
    value= xsec_list[energy_index]["compt"]
  """
  if(verbose>0):
    print("-------------------------------------------------------------------")
    print("                  Photon Cross Section (", mat, ")")
    print("Energy      Compton     Raleigh     Photo-      Pair        Total")
    print("            Scattering  Scattering  electric    Creation")
    print("(MeV)       (cm2/g)     (cm2/g)     (cm2/g)     (cm2/g)     (cm2/g)")
    print("-------------------------------------------------------------------")

  xsection_list= []
  for ekin in elist:
    xsec= {}
    xsec["compt"] \
      = gEmCalculator.ComputeCrossSectionPerVolume(ekin, "gamma", plist[0],
                                                   mat) * cm2/g
    xsec["rayleigh"] \
      = gEmCalculator.ComputeCrossSectionPerVolume(ekin, "gamma", plist[1],
                                                   mat) * cm2/g
    
    xsec["phot"] \
      = gEmCalculator.ComputeCrossSectionPerVolume(ekin, "gamma", plist[2],
                                                   mat) * cm2/g
    xsec["conv"] \
      = gEmCalculator.ComputeCrossSectionPerVolume(ekin, "gamma", plist[3],
                                                   mat) * cm2/g

    xsec["tot"]= xsec["compt"] + xsec["rayleigh"] + xsec["phot"] + xsec["conv"]
    
    xsection_list.append((ekin, xsec))

    if(verbose>0):    
      print(" %8.3e   %8.3e   %8.3e   %8.3e   %8.3e   %8.3e" \
            % (ekin/MeV, xsec["compt"]/(cm2/g), xsec["rayleigh"]/(cm2/g), 
               xsec["phot"]/(cm2/g), xsec["conv"]/(cm2/g), xsec["tot"]/(cm2/g)))
      
  return xsection_list


# ==================================================================
# Stopping Power
# ==================================================================
def CalculateDEDX(part, mat, elist, verbose=0, 
                  plist=["eIoni", "eBrem", "muIoni", "muBrems", "hIoni"]):
  """
  Calculate stopping powers for a give particle, material and
  a list of energy, returing stopping power for the components of
  "Ionization", "Radiation" and total one.
  
  Arguments:
    part:     particle name
    mat:      material name
    elist:    list of energy
    verbose:  verbose level [0]
    plist:    list of process name
              (electron ionization/electron brems/
               muon ionization/muon brems/hadron ionization) [StandardEM set]
               
  Keys of index:
    "ioni":   ionization
    "brems":  Bremsstrahlung
    "tot":    total

  Example:
    dedx_list= CalculateDEDX(...)
    value= dedx_list[energy_index]["ioni"]    
  """
  if(verbose>0):
    print("------------------------------------------------------")
    print("       Stopping Power (", part, ",", mat, ")")
    print("  Energy       Ionization    Radiation     Total")
    print("  (MeV)        (MeVcm2/g)    (MeVcm2/g)    (MeVcm2/g)")
    print("------------------------------------------------------")  

  procname_brems= ""
  procname_ioni= ""
  if ( part=="e+" or part=="e-" ):
    procname_ioni= plist[0]
    procname_brems= plist[1]
  elif ( part=="mu+" or part=="mu-"):
    procname_ioni= plist[2]
    procname_brems= plist[3]
  else:
    procname_ioni= plist[4]
    procname_brems= ""

  dedx_list= []
  for ekin in elist:
    dedx= {}
    dedx["ioni"] \
    = gEmCalculator.ComputeDEDX(ekin, part, procname_ioni, mat) * MeV*cm2/g
    dedx["brems"] \
    = gEmCalculator.ComputeDEDX(ekin, part, procname_brems, mat) * MeV*cm2/g
    dedx["tot"]= dedx["ioni"]+ dedx["brems"]

    if(verbose>0):
      print(" %8.3e     %8.3e     %8.3e     %8.3e" \
            % (ekin/MeV, dedx["ioni"]/(MeV*cm2/g), 
               dedx["brems"]/(MeV*cm2/g), dedx["tot"]/(MeV*cm2/g) ))
        

    dedx_list.append((ekin, dedx))    

  return dedx_list


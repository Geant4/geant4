//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// ----------------------
// Class description:
//
// The table of elements.
// In the constructor an array of strings with element name,
// symbol and A are created. 
// The G4Element instances of given Z are created when
// GetEle(G4double Z) is called using the string array.
// For each Z the G4Element is created only once.

// ----------------------

#ifndef G4ELETABLE_HH
#define G4ELETABLE_HH 1

#include "G3toG4Defs.hh"
#include "globals.hh"
#include "G4Element.hh"

class G3EleTable
{

public:  // with description

  G3EleTable();
  virtual ~G3EleTable();
  G4Element* GetEle(G4double Z);

private:

  void LoadUp();
  G4int parse(G4double& Z, char* name, char* sym, G4double& A); 

private:

  G4Element** _Ele;
  static const G4int _MaxEle = 109;
  G4String _EleNames[_MaxEle] = { "Hydrogen H 1.00794", "Helium He 4.0026",
                            "Lithium Li 6.941", "Beryllium Be 9.012182",
                            "Boron B 10.811", "Carbon C 12.011",
                            "Nitrogen N 14.00674", "Oxygen O 15.9994",
                            "Fluorine F 18.9984032", "Neon Ne 20.1797",

                            "Sodium Na 22.989768", "Magnesium Mg 24.3050",
                            "Aluminum Al 26.981539", "Silicon Si 28.0855",
                            "Phosphorus P 30.973762", "Sulfur S 32.066",
                            "Chlorine Cl 35.4527", "Argon Ar 39.948",
                            "Potassium K 39.0983", "Calcium Ca 40.078",

                            "Scandium Sc 44.955910", "Titanium Ti 47.867",
                            "Vanadium V 50.9415", "Chromium Cr 51.9961",
                            "Manganese Mn 54.93805", "Iron Fe 55.845",
                            "Cobalt Co 58.93320", "Nickel Ni 58.6934",
                            "Copper Cu 63.546", "Zinc Zn 65.39",

                            "Gallium Ga 69.723", "Germanium Ge 72.61",
                            "Arsenic As 74.92159", "Selenium Se 78.96",
                            "Bromine Br 79.904", "Krypton Kr 83.80",
                            "Rubidium Rb 85.4678", "Strontium Sr 87.62",
                            "Yttrium Y 88.90585", "Zirconium Zr 91.224",

                            "Niobium Nb 92.90638", "Molybdenum Mo 95.94",
                            "Technetium Tc 97.907215", "Ruthenium Ru 101.07",
                            "Rhodium Rh 102.90550", "Palladium Pd 106.42",
                            "Silver Ag 107.8682", "Cadmium Cd 112.41",
                            "Indium In 114.818", "Tin Sn 118.710",

                            "Antimony Sb 121.760", "Tellurium Te 127.60",
                            "Iodine I 126.90447", "Xenon Xe 131.29",
                            "Cesium Cs 132.90543", "Barium Ba 137.27",
                            "Lanthanum La 138.9055", "Cerium Ce 140.115",
                            "Praeseodymium Pr 140.90765", "NeoDymium Nd 144.24",

                            "Promethium Pm 144.912745", "Samarium Sm 150.36",
                            "Europium Eu 151.965", "Gadolinium Gd 157.25",
                            "Terbium Tb 158.92534", "Dysprosium Dy 162.50",
                            "Holmium Ho 164.93032", "Erbium Er 167.26",
                            "Thulium Tm 168.93421", "Ytterbium Yb 173.04",

                            "Lutetium Lu 174.967", "Hafnium Hf 178.49",
                            "Tantalum Ta 180.9479", "Tungsten W 183.84",
                            "Rhenium Re 186.207", "Osmium Os 190.23",
                            "Iridium Ir 192.217", "Platinum Pt 195.08",
                            "Gold Au 196.96654", "Mercury Hg 200.59",

                            "Thallium Tl 204.3833", "Lead Pb 207.2",
                            "Bismuth Bi 208.98037", "Polonium Po 208.982415",
                            "Astatine At 209.987131", "Radon Rn 222.017570",
                            "Francium Fr 223.019731", "Radium Ra 226.025402",
                            "Actinium Ac 227.027747", "Thorium Th 232.0381",

                            "Protactinium Pa 231.03588", "Uranium U 238.0289",
                            "Neptunium Np 237.048166", "Plutonium Pu 244.064197",
                            "Americium Am 243.061372", "Curium Cm 247.070346",
                            "Berkelium Bk 247.070298", "Californium Cf 251.079579",
                            "Einsteinium Es 252.08297", "Fermium Fm 257.095096",

                            "Mendelevium Md 258.098427", "Nobelium No 259.1011",
                            "Lawrencium Lr 262.1098", "Rutherfordium Rf 261.1089",
                            "Hahnium Ha 262.1144", "Seaborgium Sg 263.1186",
                            "Nielsborium Ns 262.1231", "Hassium Hs 265.1306",
                            "Meitnerium Mt 266.1378" };
};

extern G3G4DLL_API G3EleTable G3Ele;
#endif

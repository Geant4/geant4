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
// $Id: G4PSTARStopping.hh 82967 2014-07-21 15:54:41Z vnivanch $

#ifndef G4NISTStoppingData_h
#define G4NISTStoppingData_h 1

//---------------------------------------------------------------------------
//
// ClassName:   G4NISTStoppingData
//
// Description: Material names for data on stopping power
//
// Author:      V. Ivantchenko 22.07.2013
//
//----------------------------------------------------------------------------
//
// Class Description:
//
// Data on Stopping Powers from the NIST Data Base  
// http://physics.nist.gov/PhysRefData/STAR
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "CLHEP/Units/SystemOfUnits.h"

static const G4String nameNIST[74] = {
  "G4_A-150_TISSUE", "G4_ACETYLENE","G4_ADIPOSE_TISSUE_ICRP","G4_Ag","G4_AIR",
  "G4_Al","G4_ALUMINUM_OXIDE","G4_Ar","G4_Au","G4_B-100_BONE", //0 - 9
  "G4_Be","G4_BONE_COMPACT_ICRU","G4_C","G4_GRAPHITE_POROUS","G4_ETHYLENE",
  "G4_C-552","G4_CARBON_DIOXIDE","G4_CALCIUM_FLUORIDE","G4_CERIC_SULFATE",
  "G4_CELLULOSE_NITRATE", // 10-19
  "G4_BONE_CORTICAL_ICRP","G4_CESIUM_IODIDE","G4_Cu","G4_Fe",
  "G4_FERROUS_SULFATE",
  "G4_Gd","G4_Ge","G4_Pyrex_Glass","G4_H","G4_He", //20-29
  "G4_KAPTON","G4_Kr","G4_LITHIUM_TETRABORATE","G4_LITHIUM_FLUORIDE",
  "G4_M3_WAX",
  "G4_MS20_TISSUE","G4_METHANE","G4_Mo","G4_MUSCLE_WITH_SUCROSE",
  "G4_MUSCLE_WITHOUT_SUCROSE", // 30 -39
  "G4_MUSCLE_SKELETAL_ICRP","G4_MUSCLE_STRIATED_ICRU","G4_N",
  "G4_SODIUM_IODIDE","G4_Ne",
  "G4_NYLON-6-6","G4_O","G4_PARAFFIN","G4_Pb","G4_PHOTO_EMULSION", // 40-49
  "G4_PLASTIC_SC_VINYLTOLUENE","G4_POLYCARBONATE","G4_POLYETHYLENE",
  "G4_MYLAR","G4_PLEXIGLASS",
  "G4_POLYPROPYLENE","G4_POLYSTYRENE","G4_TEFLON","G4_POLYVINYL_CHLORIDE",
  "G4_PROPANE", // 50-59
  "G4_Pt","G4_Si","G4_SILICON_DIOXIDE","G4_STILBENE","G4_Ti",
  "G4_Sn","G4_TISSUE-METHANE","G4_TISSUE-PROPANE","G4_TOLUENE","G4_U",//60-69 
  "G4_W","G4_WATER","G4_WATER_VAPOR","G4_Xe"};

static const G4int numberOfMolecula = 12;

static const G4String molecularName[numberOfMolecula] = {
  "Al_2O_3",                 "CO_2",                      "CH_4",
  "(C_2H_4)_N-Polyethylene", "(C_2H_4)_N-Polypropylene",  "(C_8H_8)_N",
  "C_3H_8",                  "SiO_2",                     "CsI",
  "H_2O",                    "H_2O-Gas",                  "Graphite" };

static const G4int molecularIndex[numberOfMolecula] = {
  6, 16, 36, 52, 55, 54, 56, 62, 21, 71, 72, 13};

static const G4double fac = CLHEP::MeV*CLHEP::cm2/CLHEP::g;

static const G4double mvx = CLHEP::MeV;

static const G4double T0[78] = { 
    0.001*mvx, 0.0015*mvx, 0.002*mvx, 0.0025*mvx, 0.003*mvx, 
    0.004*mvx, 0.005*mvx,  0.006*mvx, 0.007*mvx,  0.008*mvx,
    0.009*mvx, 0.01*mvx,   0.0125*mvx,0.015*mvx,  0.0175*mvx, 
    0.02*mvx,  0.0225*mvx, 0.025*mvx, 0.0275*mvx, 0.03*mvx, 
    0.035*mvx, 0.04*mvx,   0.045*mvx, 0.05*mvx,   0.055*mvx, 
    0.06*mvx,  0.065*mvx,  0.07*mvx,  0.075*mvx,  0.08*mvx, 
    0.085*mvx, 0.09*mvx,   0.095*mvx, 0.1*mvx,    0.125*mvx, 
    0.15*mvx,  0.175*mvx,  0.2*mvx,   0.225*mvx,  0.25*mvx, 
    0.275*mvx, 0.3*mvx,    0.35*mvx,  0.4*mvx,    0.45*mvx,  
    0.5*mvx,   0.55*mvx,   0.6*mvx,   0.65*mvx,   0.7*mvx, 
    0.75*mvx,  0.8*mvx,    0.85*mvx,  0.9*mvx,    0.95*mvx,  
    1.*mvx,    1.25*mvx,   1.5*mvx,   1.75*mvx,   2.*mvx,   
    2.25*mvx,  2.5*mvx,    2.75*mvx,  3.*mvx,     3.5*mvx, 
    4.*mvx,    4.5*mvx,    5.*mvx,    5.5*mvx,    6.*mvx, 
   6.5*mvx,    7.*mvx,     7.5*mvx,   8.*mvx,     8.5*mvx, 
    9.*mvx,     9.5*mvx,   10.*mvx  }; 


#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

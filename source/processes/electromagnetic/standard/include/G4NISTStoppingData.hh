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

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
//   Author:            Mathieu Fontaine           Rachid Mazini
//                      fontaine@lps.umontreal.ca  Rachid.Mazini@cern.ch
//   Language:          C++
//   Tested on:         g++
//   Prerequisites:     None
//   Purpose:           Header file for FCALHadModule.cc, which defines
//                      the  geometry of the FCAL HadModule 0.
//   Developped:        10-March-2000   M.F.
//
//----------------------------------------------------------------------------

#ifndef FCALHadModule_h
#define FCALHadModule_h 1

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include "FCALHadModuleSD.hh"

class FCALHadModule
{
public:

  FCALHadModule();
  ~FCALHadModule();

public:

  G4LogicalVolume * Construct();

  void InitializeGeometry(); 
  G4int GetF2TileID(G4int); 

private:
  G4int NF2LarGap;
  G4int* F2LArGapID;
  G4int* F2LArIX;
  G4int* F2LArJY;
  G4int* F2LArITile;
  G4double* F2LArGapPosX;
  G4double* F2LArGapPosY;

private:

  G4double HadModuleRMin, HadModuleRMax, HadModuleLenght;
  G4double HadModuleStartPhi, HadModuleDPhi;
  G4double FCAL2HadSmart;

  G4double WAbsorberRMin, WAbsorberRMax, WAbsorberLenght;
  G4double WAbsorberStartPhi, WAbsorberDPhi;

  G4double CuPlateLenght, CuPlateAPosZ, CuPlateBPosZ;

  G4int NCableTroff;
  G4double F2TroffRmin, F2TroffRmax, F2TroffMainLenght;
  G4double F2TroffABLenght, F2TroffStartPhi,  F2TroffDphi; 
  G4double F2TroffRotZ, F2TroffABPosZ; 

  G4double F2LArGapRmin, F2LArGapRmax, F2LArGapLenght;
  G4double F2LArGapStartPhi, F2LArGapDphi;
  G4double F2RodRmin, F2RodRmax, F2RodLenght, F2RodStartPhi, F2RodDphi;

  FCALHadModuleSD* FcalHadModuleSD;

};

#endif  /* FCALHadModule.hh */


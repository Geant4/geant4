//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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

  G4double HadModuleRMin, HadModuleRMax, HadModuleLength;
  G4double HadModuleStartPhi, HadModuleDPhi;
  G4double FCAL2HadSmart;

  G4double WAbsorberRMin, WAbsorberRMax, WAbsorberLength;
  G4double WAbsorberStartPhi, WAbsorberDPhi;

  G4double CuPlateLength, CuPlateAPosZ, CuPlateBPosZ;

  G4int NCableTroff;
  G4double F2TroffRmin, F2TroffRmax, F2TroffMainLength;
  G4double F2TroffABLength, F2TroffStartPhi,  F2TroffDphi; 
  G4double F2TroffRotZ, F2TroffABPosZ; 

  G4double F2LArGapRmin, F2LArGapRmax, F2LArGapLength;
  G4double F2LArGapStartPhi, F2LArGapDphi;
  G4double F2RodRmin, F2RodRmax, F2RodLength, F2RodStartPhi, F2RodDphi;

  FCALHadModuleSD* FcalHadModuleSD;

};

#endif  /* FCALHadModule.hh */

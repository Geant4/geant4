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
//   Purpose:           Header file for FCALEMModule.cc, which defines
//                      the  geometry of the FCAL EMModule 0.
//   Developped:        10-March-2000   M.F.
//
//----------------------------------------------------------------------------

#ifndef FCALEMModule_h
#define FCALEMModule_h 1

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "FCALEMModuleSD.hh"


class FCALEMModule
{
public:

  FCALEMModule();
  ~FCALEMModule();

public:

  G4LogicalVolume * Construct();
  
  void InitializeGeometry();
  G4int GetF1TileID(G4int);
  G4double GetF1LArGapPosX(G4int);

private:
  G4int NF1LarGap;
  G4int* F1LArGapID;
  G4int* F1LArIX;
  G4int* F1LArJY;
  G4int* F1LArITile;
  G4double* F1LArGapPosX;
  G4double* F1LArGapPosY;

private:

  G4double EmModuleRMin, EmModuleRMax, EmModuleLength, CableTroffLength,
           EmModuleStartPhi, EmModuleDPhi;
  G4double FCALEmSmart;

  G4int NCableTroff;
  G4double F1CableTroffRMin, F1CableTroffRMax, F1CableTroffLength, F1CableTroffStartPhi, 
           F1CableTroffDPhi;
  G4double F1CableTroffRotZ;

  G4double F1LArGapRmin, F1LArGapRmax, F1LArGapLength, F1LArGapStartPhi, F1LArGapDPhi;

  FCALEMModuleSD* FcalEmModuleSD;

  //  G4double GapPosX, GapPosY;
  //G4int GapID, IX, IY, ITile;

};

#endif  /* FCALEMModule.hh */

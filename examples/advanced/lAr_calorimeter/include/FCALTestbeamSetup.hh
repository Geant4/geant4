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
//   Purpose:           Header file for FCALFrontVolume.cc, which defines
//                      the  volumes in the testbeam front.
//   Developped:        10-March-2000   M.F.
//
//----------------------------------------------------------------------------


#ifndef FCALTestbeamSetup_h
#define FCALTestbeamSetup_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

//class FCALFrontVolumes;
//class FCALTailVolumes;
//class FCALCryostatVolumes;

class G4LogicalVolume;
class G4VPhysicalVolume;

class FCALTestbeamSetup : public G4VUserDetectorConstruction
{

public:

  FCALTestbeamSetup();
  ~FCALTestbeamSetup();

public:

  G4VPhysicalVolume* Construct();
  void ConstructSDandField();
    
private:

  G4double MotherSizeX, MotherSizeY, MotherSizeZ;

  G4double MWPCSizeX, MWPCSizeY, MWPCSizeZ;
  G4double MWPCPosX, MWPCPosY, MWPCPosZ[5];

  G4double ScintS1andS3SizeX, ScintS1andS3SizeY, ScintS1andS3SizeZ;
  G4double ScintS2SizeX, ScintS2SizeY, ScintS2SizeZ;
  G4double ScintS1_S3PosX, ScintS1_S3PosY, ScintS1PosZ,ScintS2PosZ, ScintS3PosZ;

  G4double HoleCntrSizeX, HoleCntrSizeY, HoleCntrScintSizeZ, HoleCntrAbsrbSizeZ;
  G4double HoleCntrScintPosX, HoleCntrScintPosY, HoleCntrScintPosZ;
  G4double HoleCntrPbPosX, HoleCntrPbPosY, HoleCntrPbPosZ;
  G4double HoleCntrAlPosX, HoleCntrAlPosY, HoleCntrAlPosZ;
  G4double ScintHoleRmin, ScintHoleRmax, ScintHoleLenght;
  G4double AbsrbHoleRmin, AbsrbHoleRmax, AbsrbHoleLenght;
  G4double HoleStartPhi, HoleDPhi;
  G4double HolePosX, HolePosY, HolePosZ;

  G4double LeadWallSizeX, LeadWallSizeY, LeadWallSizeZ;
  G4double LeadWallSlitSizeX, LeadWallSlitSizeY, LeadWallSlitSizeZ;
  G4double LeadWallPosX,LeadWallPosY, LeadWallPosZ;

  G4double IronWallSizeX, IronWallSizeY, IronWallSizeZ;
  G4double IronWallSlitSizeX, IronWallSlitSizeY, IronWallSlitSizeZ;
  G4double IronWallPosX,IronWallPosY, IronWallPosZ;

  G4int NBigScint, NSmallScint, NBigIron, NSmallIron;
  G4double BigScintSizeX, BigScintSizeY, SmallScintSizeX, SmallScintSizeY, ScintSizeZ;
  G4double ScintPosX, ScintPosY, ScintPosZ[7];
  G4double BigIronSizeX, BigIronSizeY, SmallIronSizeX, SmallIronSizeY, IronSizeZ;
  G4double IronPosX, IronPosY, IronPosZ[6];

  G4double ConcWallSizeX, ConcWallSizeY, ConcWallSizeZ; 
  G4double ConcWallPosX, ConcWallPosY, ConcWallAPosZ, ConcWallBPosZ;
  G4double ConcWallInsSizeX, ConcWallInsSizeY, ConcWallInsSizeZ;
  G4double ConcWallInsPosZ;

  G4double MuCntrSIzeX, MuCntrSIzeY, MuCntrSIzeZ;
  G4double MuCntrPosX, MuCntrPosY, MuCntrPosZ;
  

  G4double CryostatPosX, CryostatPosY, CryostatPosZ;


  //FCALTestbeamSetupSD* FCALTBSetupSD;  // Senstive detector

  /* 
  G4double TailPosX;
  G4double TailPosY;
  G4double TailPosZ;

  G4double FrontPosX;
  G4double FrontPosY;
  G4double FrontPosZ;
  */
  /*
  G4double EMModulePosX;
  G4double EMModulePosY;
  G4double EMModulePosZ;

  G4double HadModulePosX;
  G4double HadModulePosY;
  G4double HadModulePosZ;
  */

};

#endif   /* FCALTestbeamSetup.hh */


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
//   Purpose:           Header file for FCALCryostatVolume.cc, which defines
//                      the  volumes in the cryostat.
//   Developped:        10-March-2000   M.F.
//
//----------------------------------------------------------------------------

#ifndef FCALCryostatVolumes_h
#define FCALCryostatVolumes_h 1

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "FCALEMModule.hh"
#include "FCALHadModule.hh"

class FCALCryostatVolumes
{
public:

  FCALCryostatVolumes();
  ~FCALCryostatVolumes();

public:

  G4LogicalVolume * Construct();

private:

  G4double CryostatRMin, CryostatRMax, CryostatLength, StartingPhi, DPhi;

  G4double InsulationRMin, InsulationRMax, InsulationLength;
  G4double InsulationPosX, InsulationPosY, InsulationPosZ;

  G4double LArgRMin, LArgRMax, LArgLength;
  G4double LArgPosX, LArgPosY, LArgPosZ;

  G4double FrontExcluderSizeX, FrontExcluderSizeY, FrontExcluderSizeZ;
  G4double FrontExcluderPosX, FrontExcluderPosY, FrontExcluderPosZ;

  G4double BackExcluderSize1X, BackExcluderSize1Y, 
           BackExcluderSize2X, BackExcluderSize2Y,
           BackExcluderSizeZ;
  G4double BackExcluderPosX, BackExcluderPosY, BackExcluderPosZ,
           BackExcluderRotX;

  G4double FCALEnvelopeRMin, FCALEnvelopeRMax, FCALEnvelopeLength, 
           FCALEnvelopeStartPhi, FCALEnvelopeDPhi;
  G4double FCALEnvelopePosX, FCALEnvelopePosY, FCALEnvelopePosZ,
           FCALEnvelopeRotX, FCALEnvelopeRotY;

  G4double FCALEmModulePosX, FCALEmModulePosY, FCALEmModulePosZ;
  
  G4double ModuleRotZ;

  G4double FCALHadModulePosX, FCALHadModulePosY, FCALHadModulePosZ;

  FCALEMModule * EmModule;
  FCALHadModule * HadModule;


};

#endif /* FCALcryostatVolumes.hh */


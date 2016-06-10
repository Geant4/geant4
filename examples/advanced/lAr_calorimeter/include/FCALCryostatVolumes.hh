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

  G4double CryostatRMin, CryostatRMax, CryostatLenght, StartingPhi, DPhi;

  G4double InsulationRMin, InsulationRMax, InsulationLenght;

  G4double LArgRMin, LArgRMax, LArgLenght;
  G4double LArgPosX, LArgPosY, LArgPosZ;

  G4double FrontExcluderSizeX, FrontExcluderSizeY, FrontExcluderSizeZ;
  G4double FrontExcluderPosX, FrontExcluderPosY, FrontExcluderPosZ;

  G4double BackExcluderSize1X, BackExcluderSize1Y, 
           BackExcluderSize2X, BackExcluderSize2Y,
           BackExcluderSizeZ;
  G4double BackExcluderPosX, BackExcluderPosY, BackExcluderPosZ,
           BackExcluderRotX;

  G4double FCALEnvelopeRMin, FCALEnvelopeRMax, FCALEnvelopeLenght, 
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


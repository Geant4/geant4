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
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ---------------- GFlashHitMaker ----------------
//
// Authors: E.Barberio & Joanna Weng 
// ------------------------------------------------------------

#include "G4ios.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHandle.hh"
#include "G4VGFlashSensitiveDetector.hh"

#include "GFlashHitMaker.hh"
#include "G4GFlashSpot.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"

GFlashHitMaker::GFlashHitMaker()
{
  fTouchableHandle   = new G4TouchableHistory(); // talk to ?@@@
  fpNavigator        = new G4Navigator();
  fNaviSetup         = false;
  fWorldWithSdName   = "";
  fpSpotS = new G4Step();
  fpSpotP = new G4StepPoint();
  // N.B. Pre and Post step points are common.
  fpSpotS->SetPreStepPoint(fpSpotP);
  fpSpotS->SetPostStepPoint(fpSpotP);
}

GFlashHitMaker::~GFlashHitMaker()
{
  delete fpNavigator;
  delete fpSpotP;
  fpSpotS->ResetPreStepPoint();
  fpSpotS->ResetPostStepPoint();
  delete fpSpotS;
}

void GFlashHitMaker::make(GFlashEnergySpot * aSpot, const G4FastTrack * aT)
{
  // Locate the spot
  if (!fNaviSetup)
  {
    // Choose the world volume that contains the sensitive detector based on its name (empty name for mass geometry)
    G4VPhysicalVolume* worldWithSD = nullptr;
    if(fWorldWithSdName.empty()) {
      worldWithSD = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
      } else {
      worldWithSD = G4TransportationManager::GetTransportationManager()->GetParallelWorld(fWorldWithSdName);
    }
    fpNavigator->SetWorldVolume(worldWithSD);
    fpNavigator->
      LocateGlobalPointAndUpdateTouchable(aSpot->GetPosition(),
                                          fTouchableHandle(), false);
    fNaviSetup = true;
  }
  else
  {
    fpNavigator->
      LocateGlobalPointAndUpdateTouchable(aSpot->GetPosition(),
                                          fTouchableHandle());
  }
  
  //--------------------------------------
  // Produce Hits
  // call sensitive part: taken/adapted from the stepping:
  // Send G4Step information to Hit/Dig if the volume is sensitive
  //--------------G4TouchableHistory----------------------------------------
  
  G4VPhysicalVolume* pCurrentVolume = fTouchableHandle()->GetVolume();    
  G4VSensitiveDetector* pSensitive;
  if( pCurrentVolume != 0 )
  {
    pSensitive = pCurrentVolume->GetLogicalVolume()->GetSensitiveDetector();
    G4VGFlashSensitiveDetector * gflashSensitive = 
                   dynamic_cast<G4VGFlashSensitiveDetector * > (pSensitive);
    if( gflashSensitive )
    {
      // set spot information:
      G4GFlashSpot theSpot(aSpot, aT, fTouchableHandle);
      gflashSensitive->Hit(&theSpot);
    }
    else if( pSensitive )
    {
      fpSpotS->SetTotalEnergyDeposit(aSpot->GetEnergy());
      fpSpotS->SetTrack(const_cast<G4Track*>(aT->GetPrimaryTrack()));
      fpSpotP->SetWeight(aT->GetPrimaryTrack()->GetWeight());
      fpSpotP->SetPosition(aSpot->GetPosition());
      fpSpotP->SetGlobalTime(aT->GetPrimaryTrack()->GetGlobalTime());
      fpSpotP->SetLocalTime(aT->GetPrimaryTrack()->GetLocalTime());
      fpSpotP->SetProperTime(aT->GetPrimaryTrack()->GetProperTime());
      fpSpotP->SetTouchableHandle(fTouchableHandle);
      fpSpotP->SetProcessDefinedStep(fpProcess);
      fpSpotP->SetStepStatus(fUserDefinedLimit);
      pSensitive->Hit(fpSpotS);
    }
  }
  else
  {     
    #ifdef GFLASH_DEBUG
    G4cout << "GFlashHitMaker::Out of volume  "<< G4endl;
    #endif
  }
}

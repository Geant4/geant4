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
// $Id: GFlashHitMaker.cc 98338 2016-07-07 08:20:39Z gcosmo $
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

GFlashHitMaker::GFlashHitMaker()
{
  fTouchableHandle   = new G4TouchableHistory(); // talk to ?@@@
  fpNavigator        = new G4Navigator();
  fNaviSetup         = false;
}

GFlashHitMaker::~GFlashHitMaker()
{
  delete fpNavigator;
}

void GFlashHitMaker::make(GFlashEnergySpot * aSpot, const G4FastTrack * aT)
{
  // Locate the spot
  if (!fNaviSetup)
  {
    fpNavigator->
      SetWorldVolume(G4TransportationManager::GetTransportationManager()->
                     GetNavigatorForTracking()->GetWorldVolume() );
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
  // Fills attribute of the G4Step needed
  // by our sensitive detector:
  //-------------------------------------
  // set spot information:
  G4GFlashSpot theSpot(aSpot, aT, fTouchableHandle);
  ///Navigator
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
      gflashSensitive->Hit(&theSpot);
    }
    else if (( pSensitive ) && 
             ( pCurrentVolume->GetLogicalVolume()->GetFastSimulationManager() )
            ) // Using gflash without implementing the 
              // gflashSensitive detector interface -> not allowed!
    
    {    
      G4cerr << "ERROR - GFlashHitMaker::make()" << G4endl
             << "        It is required to implement the "<< G4endl
             << "        G4VGFlashSensitiveDetector interface in "<< G4endl
             << "        addition to the usual SensitiveDetector class."
             << G4endl;
      G4Exception("GFlashHitMaker::make()", "InvalidSetup", FatalException, 
                  "G4VGFlashSensitiveDetector interface not implemented.");
    }
  }
  else
  {     
    #ifdef GFLASH_DEBUG
    G4cout << "GFlashHitMaker::Out of volume  "<< G4endl;
    #endif
  }
}

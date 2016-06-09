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
//
// $Id: GFlashHitMaker.cc,v 1.6 2005/10/04 09:08:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
    else if ( (!gflashSensitive ) && 
             ( pSensitive ) && 
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

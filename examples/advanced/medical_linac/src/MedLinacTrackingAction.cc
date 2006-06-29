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


#include "MedLinacTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
  
#include "MedLinacDetectorConstruction.hh" 
#include "MedLinacTargetAndFilterDecorator.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

void MedLinacTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
 {

  
    // global geometry navigator
    gNavigator = G4TransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking();


   // Create trajectory only for primaries

      // check if particle is gamma
   G4ParticleDefinition* particleType = aTrack->GetDefinition();
   G4ThreeVector direction = aTrack->GetMomentumDirection();
   G4double theta  = std::acos(std::abs(direction.z()));

   //if(aTrack->GetParentID()==0)
   //fpTrackingManager->SetStoreTrajectory(true); 

       if(particleType==G4Gamma::GammaDefinition()) {
	if(theta>=20.*deg){
        // check if particle is in target
        G4ThreeVector pos = aTrack->GetPosition();
        G4ThreeVector *ptr = NULL;
        G4VPhysicalVolume *theVolume;
	theVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,false);
        if(theVolume->GetName() == "targetA"||theVolume->GetName() == "targetB") {
    
	  fpTrackingManager->SetStoreTrajectory(false); }
	}
}
       else
	 {fpTrackingManager->SetStoreTrajectory(true);} 
}

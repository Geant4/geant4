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

#include "G4Track.hh"

#include "F04SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "F04SteppingActionMessenger.hh"

#include "G4ParticleTypes.hh"

#include "F04UserTrackInformation.hh"

F04SteppingAction::F04SteppingAction()
{
  steppingMessenger = new F04SteppingActionMessenger(this);
}

F04SteppingAction::~F04SteppingAction()
{
  delete steppingMessenger ;
}

void F04SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  G4Track* theTrack = theStep->GetTrack();

  if (theTrack->GetParentID()==0) {
    //This is a primary track
    G4String theVolumeName = theStep->GetPreStepPoint()->
                                  GetPhysicalVolume()->GetName();
    if (theVolumeName != "Target") {
       theTrack->SetTrackStatus(fStopAndKill);
       return;
    }
  }

//  G4ParticleDefinition* particleType = theTrack->GetDefinition();

  // check if it is entering the test volume
  G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
  G4String thePrePVname = thePrePV->GetName();

  G4String thePostPVname = " ";
  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();

  if (thePostPoint) {
     G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
     if (thePostPV) thePostPVname = thePostPV->GetName();
  }

  if (thePrePVname  != "TestPlane" &&
      thePostPVname == "TestPlane") {

//     G4double x = theTrack->GetPosition().x();
//     G4double y = theTrack->GetPosition().y();

     G4ThreeVector theMomentumDirection = theTrack->
                          GetDynamicParticle()->GetMomentumDirection();
//     G4double theTotalMomentum = theTrack->GetDynamicParticle()->
//                                                    GetTotalMomentum();

     // then kill the track
     theTrack->SetTrackStatus(fStopAndKill);
     return;

  }

//  G4double z = theTrack->GetPosition().z();

  F04UserTrackInformation* trackInformation =
                        (F04UserTrackInformation*)theTrack->GetUserInformation();

  if (trackInformation->GetTrackStatusFlag() != reverse) {
     if (thePrePVname  != "Target" ) {
        if ( theTrack->      GetMomentumDirection().z()>0.0 &&
             theTrack->GetVertexMomentumDirection().z()<0.0     )
        {
             trackInformation->SetTrackStatusFlag(reverse);
        }
     }
  }

  // check if it is alive
  if (theTrack->GetTrackStatus() == fAlive) { return; }

  if (thePostPoint->GetProcessDefinedStep() != 0) {
     if (thePostPoint->GetProcessDefinedStep()->GetProcessName() != "Decay") return;
  }

}

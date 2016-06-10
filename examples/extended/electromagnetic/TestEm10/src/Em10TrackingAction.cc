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
// $Id: Em10TrackingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
/// \file electromagnetic/TestEm10/src/Em10TrackingAction.cc
/// \brief Implementation of the Em10TrackingAction class
//
 #include "Em10TrackingAction.hh"

 #include "G4TrackingManager.hh"
 #include "G4Track.hh"
 #include "G4RunManager.hh"

 #include "G4UImanager.hh"
 #include "G4SystemOfUnits.hh"

 Em10TrackingAction::Em10TrackingAction() 
 : G4UserTrackingAction()
 { }

 void Em10TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
 {

  if( aTrack->GetParentID() == 1 && aTrack->GetKineticEnergy() > 100.*GeV ){
    G4cout << "[Em10TrackingAction::DEBUG]" << G4endl;
    G4cout << " Track ID:          " << aTrack->GetTrackID() << G4endl;
    G4cout << " particle:          " << aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() << G4endl;
    G4cout << " Parent ID:         " << aTrack->GetParentID() << G4endl;
    G4cout << " created by:        " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
    G4cout << " kin. energy (TeV): " << aTrack->GetKineticEnergy() / TeV << G4endl;
    G4cout << " volume:            " << aTrack->GetVolume()->GetName() << G4endl;
    G4cout << " global time:       " << aTrack->GetGlobalTime() << G4endl;

    G4cout << " Killing event..." << G4endl;
    if( aTrack->GetTrackID() != 1 )
      const_cast<G4Track*>(aTrack)->SetTrackStatus( fKillTrackAndSecondaries );
  }

  return;
 }

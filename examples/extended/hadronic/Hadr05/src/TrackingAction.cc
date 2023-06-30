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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"

#include "G4RunManager.hh"
#include "G4StepStatus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, EventAction* evt)
:G4UserTrackingAction(),fDetector(det),fEventAct(evt)
{ }
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track )
{
  //get Run
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
             
  // Energy flow initialisation for primary particle
  //
  if (track->GetTrackID() == 1) {
    G4int Idnow = 1;
    if (track->GetVolume() != fDetector->GetphysiWorld()) {
      // unique identificator of layer+absorber
      const G4VTouchable* touchable = track->GetTouchable();
      G4int absorNum = touchable->GetCopyNumber();
      G4int layerNum = touchable->GetReplicaNumber(1);
      Idnow = (fDetector->GetNbOfAbsor())*layerNum + absorNum;
    }
    
    G4double Eflow = track->GetKineticEnergy();
    ///if (track->GetDefinition() == G4Positron::Positron()) {
    ///  Eflow += 2*electron_mass_c2;
    ///}
         
    //flux artefact, if primary vertex is inside the calorimeter   
    for (G4int pl=1; pl<=Idnow; ++pl) {run->SumEnergyFlow(pl, Eflow);}
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{             
 // energy leakage
 G4StepStatus status = aTrack->GetStep()->GetPostStepPoint()->GetStepStatus();
 if (status == fWorldBoundary) { 
    G4int parentID = aTrack->GetParentID();
    G4int index = 0; if (parentID > 0) index = 1;    //primary=0, secondaries=1
    G4double eleak = aTrack->GetKineticEnergy();
    fEventAct->SumEnergyLeak(eleak,index); 
 }               
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


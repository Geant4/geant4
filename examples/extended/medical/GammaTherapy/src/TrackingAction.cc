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
// $Id: TrackingAction.cc 103469 2017-04-11 07:29:36Z gcosmo $
//
/// \file medical/GammaTherapy/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//

//---------------------------------------------------------------------------
//
// ClassName:   TrackingAction
//
// Description: Implementation file for MC truth.
//
// Author:      V.Ivanchenko 17/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TrackingAction.hh"
#include "G4RunManager.hh"
#include "Run.hh"

#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackingAction::TrackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackingAction::~TrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  run->ScoreNewTrack(aTrack);

  if(run->GetVerbose()) {
    G4cout << "New track #"
           << aTrack->GetTrackID() << " of " 
           << aTrack->GetParticleDefinition()->GetParticleName()
           << " Ekin(MeV)= " << aTrack->GetKineticEnergy()/MeV
           << " parent # " << aTrack->GetParentID()
           << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackingAction::PostUserTrackingAction(const G4Track*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
/// \file hadronic/Hadr01/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
//
/////////////////////////////////////////////////////////////////////////
//
// StackingAction
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#include "StackingAction.hh"
#include "HistoManager.hh"
#include "StackingMessenger.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
 : G4UserStackingAction(),
   fHistoManager(0), fStackMessenger(0), fParticle(0) 
{
  fStackMessenger = new StackingMessenger(this);
  fHistoManager   = HistoManager::GetPointer();
  fKillSecondary  = false;
  fParticle       = 0; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete fStackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4ClassificationOfNewTrack status = fUrgent;

  if (aTrack->GetTrackStatus() == fAlive) {
    fHistoManager->ScoreNewTrack(aTrack);
  }

  const G4ParticleDefinition* part = aTrack->GetDefinition();

  if(fHistoManager->GetVerbose() > 1 ) {
    G4cout << "Track #"
           << aTrack->GetTrackID() << " of " << part->GetParticleName()
           << " E(MeV)= " << aTrack->GetKineticEnergy()/MeV
           << " produced by Track ID= " << aTrack->GetParentID()
           << G4endl;
  }

  //stack or delete secondaries
  if(aTrack->GetTrackID() > 1) {  
    if (fKillSecondary || fParticle == part) { status = fKill; }
  }
  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::SetKillStatus(G4bool value)    
{ 
  fKillSecondary = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::SetKill(const G4String& name)  
{ 
  fParticle = G4ParticleTable::GetParticleTable()->FindParticle(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

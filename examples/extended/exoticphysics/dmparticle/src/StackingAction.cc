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
/// \file exoticphysics/dmparticle/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
// $Id: StackingAction.cc,v 1.8 2009-03-06 18:04:23 maire Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "Run.hh"
#include "StackingMessenger.hh"
#include "G4LDMPhoton.hh"
#include "G4LDMHi.hh"
#include "G4LDMHiBar.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
  : G4UserStackingAction(), fKillSecondary(false)
{
  fStackMessenger = new StackingMessenger(this);
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
  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;

  //keep primary particle
  if(aTrack->GetParentID() == 0) { return status; }
  
  const G4ParticleDefinition* part = aTrack->GetDefinition();
  if(part == G4LDMPhoton::LDMPhoton() ||
     part == G4LDMHi::LDMHi() ||
     part == G4LDMHiBar::LDMHiBar()) { 
    G4cout << "### New exotic particle produced: " 
           << part->GetParticleName()
           << " Ekin(GeV)= " << aTrack->GetKineticEnergy()/CLHEP::GeV
           << " Mass(GeV)= " << part->GetPDGMass()/CLHEP::GeV
           << " TrackId= " << aTrack->GetTrackID() 
           << G4endl;
    return status; 
  }

  if(fKillSecondary) { status = fKill; }
  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
// $Id: StackingAction.cc,v 1.1 2008-11-20 08:55:42 antoni Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// StackingAction
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
// 

#include "StackingAction.hh"

#include "HistoManager.hh"
#include "StackingMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
{
  stackMessenger = new StackingMessenger(this);
  histoManager   = HistoManager::GetPointer();
  killSecondary  = false;
  pname          = ""; 
  nBiased        = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete stackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4ClassificationOfNewTrack status = fUrgent;

  if (aTrack->GetTrackStatus() == fAlive) {
    histoManager->ScoreNewTrack(aTrack);
  }

  if(histoManager->GetVerbose() > 1 ) {
    const G4String name = aTrack->GetDefinition()->GetParticleName();
    G4cout << "Track #"
	   << aTrack->GetTrackID() << " of " << name
	   << " E(MeV)= " << aTrack->GetKineticEnergy()/MeV
	   << " produced by " 
	   << " ID= " << aTrack->GetParentID()
	   << G4endl;
  }
  if(aTrack->GetTrackID() == 1) { 
    return status; 

  } else if(0 < nBiased) {
    for(G4int i=0; i<nBiased; ++i) {
      if(aTrack->GetDefinition() == biasedParticle[i]) {
        if(aTrack->GetKineticEnergy() < biasedEnergy[i]) {
	  if(G4UniformRand() < rrProbability[i]) {
	    const_cast<G4Track*>(aTrack)->SetWeight(aTrack->GetWeight()/rrProbability[i]);
	  } else {
	    status = fKill;
	  }
	}
	break;
      }
    }
  }

  //stack or delete secondaries
  //if (killSecondary)      { status = fKill; }
  //else if(pname == name)  { status = fKill; }

  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
StackingAction::ActivateSecondaryBiasing(const G4String& pname, 
					 G4double f, G4double e)
{
  G4ParticleDefinition* part = 
    G4ParticleTable::GetParticleTable()->FindParticle(pname);
  if(part) {

    biasedParticle.push_back(part);
    rrProbability.push_back(f);
    biasedEnergy.push_back(e);
    ++nBiased;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
// $Id: StackingAction.cc,v 1.1 2006-06-02 19:00:02 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "HistoManager.hh"
#include "StackingMessenger.hh"

#include "G4Track.hh"
#include "G4Neutron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
{
  verbose        = 1;
  killSecondary  = false;
  stackMessenger = new StackingMessenger(this);
  histoManager   = HistoManager::GetPointer();
  tmin           = 90.*keV; 
  tmax           = 110.*keV; 
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
  if (aTrack->GetTrackStatus() == fAlive) 
    histoManager->ScoreNewTrack(aTrack);
  //keep primary particle
  if (aTrack->GetParentID() == 0) return fUrgent;

  //
  //energy spectrum of secondaries
  //
  G4double e = aTrack->GetKineticEnergy();
  const G4ParticleDefinition* p = aTrack->GetDefinition();
  //  G4double m = p->GetPDGMass();

  if(verbose > 1 && p == G4Neutron::Neutron()
     && ( e < 2.*eV ||( e < tmax && e> tmin)) ) {

    G4int pid = aTrack->GetParentID();
    const G4String name = aTrack->GetDefinition()->GetParticleName();
    G4cout << "Track #"
	   << aTrack->GetTrackID() << " of " << name
	   << " E(keV)= " << e/keV
	   << " produced by " 
	   << histoManager->CurrentDefinition()->GetParticleName()
	   << " ID= " << pid
	   << " with E(MeV)= " << histoManager->CurrentKinEnergy()/MeV
	   << G4endl;
  }

  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;
  if (killSecondary)         status = fKill;
  
  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

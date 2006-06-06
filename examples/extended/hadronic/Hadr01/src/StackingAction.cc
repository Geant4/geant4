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
// $Id: StackingAction.cc,v 1.2 2006-06-06 19:48:38 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// StackingAction
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#include "StackingAction.hh"

#include "HistoManager.hh"
#include "StackingMessenger.hh"

#include "G4Track.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
{
  stackMessenger = new StackingMessenger(this);
  histoManager   = HistoManager::GetPointer();
  killSecondary  = false;
  pname          = ""; 
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

  if (aTrack->GetTrackStatus() == fAlive) 
    histoManager->ScoreNewTrack(aTrack);

  const G4String name = aTrack->GetDefinition()->GetParticleName();

  if(histoManager->GetVerbose() > 1 ) {
    G4cout << "Track #"
	   << aTrack->GetTrackID() << " of " << name
	   << " E(MeV)= " << aTrack->GetKineticEnergy()/MeV
	   << " produced by " 
	   << histoManager->CurrentDefinition()->GetParticleName()
	   << " ID= " << aTrack->GetParentID()
	   << " with E(MeV)= " << histoManager->CurrentKinEnergy()/MeV
	   << G4endl;
  }

  //stack or delete secondaries
  if (killSecondary)      status = fKill;
  else if(pname == name)  status = fKill; 

  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

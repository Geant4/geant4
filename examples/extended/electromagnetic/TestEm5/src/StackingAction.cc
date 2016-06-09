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
//
// $Id: StackingAction.cc,v 1.1 2003/08/11 10:21:33 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "StackingMessenger.hh"

#include "G4Track.hh"

#ifdef G4ANALYSIS_USE
 #include "AIDA/IHistogram1D.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(RunAction* RA, EventAction* EA, HistoManager* HM)
:runaction(RA), eventaction(EA), histoManager(HM)
{
 killSecondary = false;
 stackMessenger = new StackingMessenger(this);
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
  //keep primary particle 
  if (aTrack->GetParentID() == 0) return fUrgent;
  
  //count secondary particles
  runaction->CountParticles(aTrack->GetDefinition());
  
#ifdef G4ANALYSIS_USE
  //  
  //energy spectrum of secondaries
  //
  G4double energy = aTrack->GetKineticEnergy();
  G4double charge = aTrack->GetDefinition()->GetPDGCharge();
  
  if (charge != 0.) {
    if (histoManager->GetHisto(2)) {
      G4double unit = histoManager->GetHistoUnit(2);
      histoManager->GetHisto(2)->fill(energy/unit);
    }
  }
      
  if (aTrack->GetDefinition() == G4Gamma::Gamma()) {
    if (histoManager->GetHisto(3)) {
      histoManager->GetHisto(3)->fill(log10(energy/MeV));
    }
  }  
#endif
      
  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;
  if (killSecondary)         status = fKill;
  
  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

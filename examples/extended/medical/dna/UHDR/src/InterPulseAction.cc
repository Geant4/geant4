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
// author: Le Tuan Anh
/// \file InterPulseAction.cc
/// \brief Implementation of the InterPulseAction class

#include "InterPulseAction.hh"

#include "CLHEP/Random/RandGeneral.h"
#include "Scorer.hh"

#include "G4EventManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

InterPulseAction::InterPulseAction(const G4String& pulse, G4bool useHisto, G4double pulsePeriod,
                                   G4int npulses)
  : PulseAction(pulse, useHisto), fNumberOfPulse(npulses), fPulsePeriod(pulsePeriod)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void InterPulseAction::PreUserTrackingAction(const G4Track* atrk)
{
  if (IsActivedPulse() && GetPulseFileName().empty()) {
    return;
  }
  if (atrk->GetParentID() == 0) {
    G4double DelayedTimeInPulse = RandomizeInPulse();

    G4int activePulse = WhichPulse();
    fDelayedTime = DelayedTimeInPulse + G4double(activePulse - 1) * fPulsePeriod;

    if (GetVerbose() > 1) {
      G4cout << "Particle comes at : " << G4BestUnit(fDelayedTime, "Time")
             << " in Pulse : " << activePulse << G4endl;
    }
    if (GetLonggestDelayedTime() < fDelayedTime) {
      SetLonggestDelayedTime(fDelayedTime);
    }
  }
  auto pPulseInfo = new PulseInfo(fDelayedTime);
  ((G4Track*)atrk)->SetUserInformation(pPulseInfo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int InterPulseAction::WhichPulse() const
{
  G4int random = std::floor(1 + fNumberOfPulse * G4UniformRand());
  G4int output = fNumberOfPulse == 1 ? 1 : random;
  return output;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

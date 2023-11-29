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
//
#include "TimeStepAction.hh"
#include "G4DNAEventScheduler.hh"
#include "G4DNAGillespieDirectMethod.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4ITLeadingTracks.hh"
#include "G4ITTransportationManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4Scheduler.hh"
#include "PulseAction.hh"
#include "G4MoleculeCounter.hh"
#include "G4VChemistryWorld.hh"

TimeStepAction::TimeStepAction(const G4VChemistryWorld *pChemWorld, PulseAction *pPulse)
    : G4UserTimeStepAction(), fpPulse(pPulse), fpChemWorld(pChemWorld) {
  fpEventScheduler = std::make_unique<G4DNAEventScheduler>();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TimeStepAction::UserPreTimeStepAction() {}

void TimeStepAction::UserPostTimeStepAction() {
  G4double T1 = 5 * CLHEP::ns;
  if (fpPulse != nullptr && fpPulse->IsActivedPulse()) {
    G4MoleculeCounter::Instance()->Use(false);
    //we don't count molecules during the pulse
    T1 = fpPulse->GetLonggestDelayedTime() + 5 * CLHEP::ns;
  }
  // T1: time to start mesoscopic model
  if (G4Scheduler::Instance()->GetGlobalTime() >= T1) {
    CompartmentBased();
  }
}

void TimeStepAction::UserReactionAction(
    const G4Track & /*a*/, const G4Track & /*b*/,
    const std::vector<G4Track *> * /*products*/) {}

void TimeStepAction::EndProcessing() {}

void TimeStepAction::CompartmentBased() {
  SetInitialPixel();
  fpEventScheduler->SetVerbose(G4Scheduler::Instance()->GetVerbose());
  fpEventScheduler->SetStartTime(G4Scheduler::Instance()->GetGlobalTime()); // continue from globalTime
  fpEventScheduler->SetEndTime(G4Scheduler::Instance()->GetEndTime() - 1 * ps);
  fpEventScheduler->SetChangeMesh(true);
  fpEventScheduler->Initialize(*fpChemWorld->GetChemistryBoundary(), fPixel);
  fpEventScheduler->Run();
  G4Scheduler::Instance()->Stop();
}

G4DNAEventScheduler *TimeStepAction::GetEventScheduler() const {
  return fpEventScheduler.get();
}

void TimeStepAction::SetInitialPixel() {
  auto pBoundingBox = fpChemWorld->GetChemistryBoundary();
  G4double Box = pBoundingBox->halfSideLengthInX();
  if (Box == 1.6 * um) {
    fPixel = 512;// for CONV
  } else if (Box == 0.8 * um) {
    fPixel = 256;// for FLASH
  } else {
    G4cout << "Box  : " << *pBoundingBox << "  Pixel : " << fPixel << G4endl;
    G4Exception("This chem volume is not optimized and the result may be incorrect.",
                "TimeStepAction::TimeStepAction", FatalException, "");
  }
  // 512 : for conventional dose rate
  // 256 : for higher dose rate
}

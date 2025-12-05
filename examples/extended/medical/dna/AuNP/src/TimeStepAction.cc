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
/// \file TimeStepAction.cc
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"

#include "DetectorConstruction.hh"
#include "G4AnalysisManager.hh"
#include "G4IT.hh"
#include "G4ITTrackHolder.hh"
#include "G4RunManager.hh"
#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"

TimeStepAction::TimeStepAction()
{
  fpDetector = dynamic_cast<const DetectorConstruction *>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::Save(const G4MolecularConfiguration *molconf)
{
  const G4int moleculeID = molconf->GetMoleculeID();
  const G4String &moleculeName = molconf->GetFormatedName();
  G4TrackList *trackList = G4ITTrackHolder::Instance()->GetMainList(moleculeID);

  if (trackList == nullptr) return;

  for(const G4Track* track : *trackList)
  {
    SaveMoleculeInfo(track, moleculeID, moleculeName);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::
SaveMoleculeInfo(const G4Track *track, const G4int molID,
                 const G4String & /*moleculeName*/)
{
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  if (!analysisManager->IsActive()) {
    return;
  }

  const G4ThreeVector &position = track->GetPosition();
  // H_{3}O^{1} ID = 0
  // OH^{-1}    ID = 1
  // OH^{0}     ID = 2
  // e_{aq}^{1} ID = 3
  // H0         ID = 4
  // H_{2}^{0}  ID = 5
  // H2O2       ID = 6
  // H_{2}O^{0} or H_{2}O^{1} ID=7-17

  const G4double xp = position.x();
  const G4double yp = position.y();
  const G4double zp = position.z();
  const G4double R = std::sqrt(xp * xp + yp * yp + zp * zp) / nm;
  if (molID < 7) {
    constexpr G4int offset = 10;
    analysisManager->FillH1(molID + offset, R);
  } else {
    analysisManager->FillH1(17, R);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserPreTimeStepAction()
{
  // Loop over defined molecules
  G4ConfigurationIterator it =
    G4MoleculeTable::Instance()->GetConfigurationIterator();

  const G4double time = G4Scheduler::Instance()->GetGlobalTime();

  if (time == 1. * us) {
    while (it()) {
      const G4MolecularConfiguration *molconf = it.value();
      Save(molconf);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

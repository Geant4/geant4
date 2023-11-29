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
/// \file scavenger/src/Run.cc
/// \brief Implementation of the scavenger::Run class

#include "Run.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "ScoreSpecies.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"

namespace scavenger
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

Run::Run()
  : G4Run(){
  G4MultiFunctionalDetector *mfdet = dynamic_cast<G4MultiFunctionalDetector *>
  (G4SDManager::GetSDMpointer()->FindSensitiveDetector("mfDetector"));
  G4int CollectionIDspecies =
    G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/Species");
  fScorerRun = mfdet->GetPrimitive(CollectionIDspecies);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void Run::RecordEvent(const G4Event *event) {
  if (event->IsAborted()) {
    return;
  }
  G4int CollectionID =
    G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/Species");
  //Hits collections
  G4HCofThisEvent *pHCE = event->GetHCofThisEvent();
  if (!pHCE) {
    return;
  }
  auto evtMap = dynamic_cast<G4THitsMap<G4double> *>(pHCE->GetHC(CollectionID));
  auto map = evtMap->GetMap();
  for (const auto &it : *map) {
    G4double edep = *(it.second);
    fSumEne += edep;
  }
  G4Run::RecordEvent(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void Run::Merge(const G4Run *aRun) {
  if (aRun == this) {
    return;
  }
  const auto localRun = dynamic_cast<const Run *>(aRun);
  fSumEne += localRun->fSumEne;
  auto masterScorer = dynamic_cast<ScoreSpecies *>(this->fScorerRun);
  auto localScorer = dynamic_cast<ScoreSpecies *>(localRun->fScorerRun);
  masterScorer->AbsorbResultsFromWorkerScorer(localScorer);
  G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

}
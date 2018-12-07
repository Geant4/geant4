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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file Run.cc
// \brief Implementation of the Run class

#include "Run.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "ScoreSpecies.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"

#include <map>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
 : G4Run(),
   fSumEne(0),
   fScorerRun(0)
{
  G4MultiFunctionalDetector* mfdet =  dynamic_cast<G4MultiFunctionalDetector*>
    (G4SDManager::GetSDMpointer()->FindSensitiveDetector("mfDetector"));
  G4int CollectionIDspecies =
    G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/Species");

  fScorerRun = mfdet->GetPrimitive(CollectionIDspecies);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* event)
{
  if(event->IsAborted()) return;

  G4int CollectionID =
    G4SDManager::GetSDMpointer()->GetCollectionID("mfDetector/Species");
  // G4int evtNb = event->GetEventID();
  
  //  G4cout << G4endl << "---> end of event: " << evtNb << G4endl;
  //Hits collections
  //
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if(!HCE) return;

  G4THitsMap<G4double>* evtMap =
  static_cast<G4THitsMap<G4double>*>(HCE->GetHC(CollectionID));

  std::map<G4int,G4double*>::iterator itr;
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++){
    G4double edep = *(itr->second);
    fSumEne+=edep;
    // G4cout<<"Energy for this event: "<<edep/eV<<" eV"<<G4endl;
  }

  G4Run::RecordEvent(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* aRun)
{
  if(aRun == this){
    return;
  }

  const Run* localRun = static_cast<const Run*>(aRun);
  fSumEne += localRun->fSumEne;

  ScoreSpecies* masterScorer =
  dynamic_cast<ScoreSpecies*>(this->fScorerRun);

  ScoreSpecies* localScorer =
  dynamic_cast<ScoreSpecies*>(localRun->fScorerRun);

  masterScorer->AbsorbResultsFromWorkerScorer(localScorer);
  G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

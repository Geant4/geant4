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
/// @file RunAction.cc
/// @brief Describe run actions

#include <G4VUserMPIrunMerger.hh>
#include "G4MPImanager.hh"
#include <stdio.h>
#include "G4Threading.hh"
#include "Analysis.hh"
#include "RunAction.hh"
#include "Run.hh"

#include "G4MPIscorerMerger.hh"
#include "toolx/mpi/wrmpi"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::RunAction(G4bool useNtuple, G4bool mergeNtuple)
 : G4UserRunAction()
{
  // Book analysis in ctor
  Analysis* myana = Analysis::GetAnalysis();
  myana->SetUseNtuple(useNtuple);
  myana->SetMergeNtuple(mergeNtuple);
  G4cout<<"Book analysis on rank: " << G4MPImanager::GetManager()-> GetRank() << G4endl;;
  myana->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Run* RunAction::GenerateRun()
{
  return new Run;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::BeginOfRunAction(const G4Run*)
{
  G4cout << "RunAction::BeginOfRunAction" << G4endl;

  Analysis* myana = Analysis::GetAnalysis();

  // ntuples will be written on each rank
  G4int rank = G4MPImanager::GetManager()->GetRank();
  std::ostringstream fname;
  fname<<"dose-rank"<<rank;
  myana->OpenFile(fname.str());
     // OpenFile triggeres creating collecting/sending ntuples objects;
     // must be called at BeginOfRunAction
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::EndOfRunAction(const G4Run*)
{
  G4cout << "RunAction::EndOfRunAction" << G4endl;

  Analysis* myana = Analysis::GetAnalysis();
  myana-> Save();
}

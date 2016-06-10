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
// $Id: RunAction.cc 82385 2014-06-18 09:25:07Z gcosmo $
//
/// @file RunAction.cc
/// @brief Describe run actions

#include "G4MPImanager.hh"
#include <stdio.h>
#include "G4Threading.hh"
#include "Analysis.hh"
#include "RunAction.hh"

#include "G4MPIRunMerger.hh"
#include "G4MPIScorerMerger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::BeginOfRunAction(const G4Run*)
{
  Analysis* myana = Analysis::GetAnalysis();
  myana-> Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::EndOfRunAction(const G4Run* arun)
{
  G4bool is_worker = G4Threading::IsWorkerThread();
  if ( is_worker ) return;

  G4int rank = G4MPImanager::GetManager()-> GetRank();
  
  G4cout<<"================================================"<<G4endl;
  G4cout<<"My rank is: "<<rank<<G4endl;
  G4cout<<"================================================"<<G4endl;
  G4MPIRunMerger rm( arun );
  //rm.SetVerbosity( 3 );
  rm.Merge();
  //G4cout<<"Num events: "<<arun->GetNumberOfEvent()<<G4endl;
  if (  G4ScoringManager::GetScoringManagerIfExist() ) {
    G4MPIScorerMerger sm( G4ScoringManager::GetScoringManagerIfExist() );
    //sm.SetVerbosity(4);
    sm.Merge();
  }
  G4cout<<"================================================"<<G4endl;
  G4cout<<"================================================"<<G4endl;

  char str[64];
  sprintf(str, "dose-%03d.root", rank);
  G4String fname(str);

  Analysis* myana = Analysis::GetAnalysis();
  myana-> Save(fname);
}

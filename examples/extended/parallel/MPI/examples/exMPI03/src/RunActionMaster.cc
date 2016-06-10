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
// $Id: RunActionMaster.cc 82310 2014-06-14 00:27:39Z adotti $
//
/// @file RunActionMaster.cc
/// @brief Describe run actions

#include "G4MPImanager.hh"
#include <stdio.h>
#include "G4Threading.hh"
#include "Analysis.hh"
#include "RunActionMaster.hh"
#include "g4root.hh" //ROOT Format for output
#include "RunMerger.hh"
#include "G4MPIscorerMerger.hh"
#include "G4MPIhistoMerger.hh"    // New code with use of g4analysis
#include "Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunActionMaster::RunActionMaster()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Run* RunActionMaster::GenerateRun()
{
  return new Run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunActionMaster::~RunActionMaster()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
RunActionMaster::BeginOfRunAction(const G4Run*)
{
  Analysis* myana = Analysis::GetAnalysis();
  myana-> Clear();
  myana->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void
RunActionMaster::EndOfRunAction(const G4Run* arun)
{
  //This is executed by master thread. Workers have already
  //merged their histograms into this master threads
  G4int rank = G4MPImanager::GetManager()-> GetRank();

  //Save histograms before MPI merging
  if (rank == 0)
    {
      G4String fname("dose-rank000");
      Analysis* myana = Analysis::GetAnalysis();
      myana-> Save(fname);
    }

  G4cout << "======================================";
  G4cout <<  "===========" << G4endl;
  G4cout << "Start EndOfRunAction for master thread in rank: " << rank<<G4endl;
  G4cout << "======================================";
  G4cout << "===========" << G4endl;

  //Merging of G4Run object:
  //All ranks > 0 merge to rank #0
  RunMerger rm(static_cast<const Run*>(arun));
  G4int ver = 0 ; //Use 4 for lots of info
  if ( rank == 0 && ver == 0) ver = 1;

  if ( ver > 1 ) {
  G4cout<<"Before merge the Run object has test counter value: "
      <<static_cast<const Run*>(arun)->GetCounter()<<G4endl;
  }
  rm.SetVerbosity( ver );
  rm.Merge();
  if ( ver > 1 ) {
  G4cout<<"After merge the Run object has test counter value: "
      <<static_cast<const Run*>(arun)->GetCounter()
      <<" (with 1 thread== number of ranks)"<<G4endl;
  }

  //Merge of scorers
  //ver = 0;
  if (G4ScoringManager::GetScoringManagerIfExist())
    {
      G4MPIscorerMerger sm(G4ScoringManager::GetScoringManagerIfExist());
      sm.SetVerbosity(ver);
      sm.Merge();
    }

  //Merge of g4analysis objects
  ver=0;
  G4MPIhistoMerger hm(G4AnalysisManager::Instance());
  hm.SetVerbosity(ver);
  hm.Merge();

  G4cout << "======================================";
  G4cout <<  "==========" << G4endl;
  G4cout << "End EndOfRunAction for master thread in rank: " << rank << G4endl;
  G4cout << "======================================";
  G4cout <<  "==========" << G4endl;

  //Save g4analysis objects to a file
  //NB: It is important that the save is done *after* MPI-merging of histograms

  //One can save all ranks or just rank0: remember in case of all ranks,
  //the file of rank0 contains the sum of everything
  if (true /*rank == 0*/)
    {
      char str[64];
      sprintf(str, "dose-rank%03d", rank);
      G4String fname(str);
      if (rank == 0)
        fname = "dose-merged";
      Analysis* myana = Analysis::GetAnalysis();
      myana-> Save(fname);
    }
    Analysis* myana = Analysis::GetAnalysis();
    myana-> Close();
}

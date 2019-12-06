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
  //This is executed by master thread only. Worker threads have already
  //merged their results into this master threads.

  //We are going to merge the results via MPI for:
  // 1. User-defined "Run" object
  // 2. Command line scorers, if exists
  // 3. G4Analysis objects

  // For debugging purposes in the following we write out results twice:
  // BEFORE and AFTER the merging, so that rank 0 actually
  // writes two files, the one called "*rank000*" will contain the partial
  // results only from rank #0,
  // the file *merged* contains the reduction from all ranks.
  // It should be very easy to adapt this code

  const G4int rank = G4MPImanager::GetManager()-> GetRank();

  G4cout << "=====================================================" << G4endl;
  G4cout << "Start EndOfRunAction for master thread in rank: " << rank<<G4endl;
  G4cout << "=====================================================" << G4endl;

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
  ver = 0;
  if (G4ScoringManager::GetScoringManagerIfExist()) {
      const auto scor = G4ScoringManager::GetScoringManager();
      G4MPIscorerMerger sm(scor);
      sm.SetVerbosity(ver);
      //Debug!
      auto debugme = [&scor](){
        for ( size_t idx = 0 ; idx < scor->GetNumberOfMesh() ; ++idx) {
            const auto m = scor->GetMesh(idx);
            const auto map = m->GetScoreMap();
            std::for_each(map.begin(),map.end(),
              [](const G4VScoringMesh::MeshScoreMap::value_type& e) {
                G4cout<<e.first<<"("<<e.second<<"):"<<G4endl;
                const auto data = e.second->GetMap();
                for( auto it = data->begin() ; it != data->end() ; ++it ) {
                    G4cout<<it->first<<" => G4StatDouble(n,sum_w,sum_w2,sum_wx,sum_wx2): "
                    <<it->second->n()<<" "<<it->second->sum_w()<<" "
                    <<it->second->sum_w2()<<" "<<it->second->sum_wx()<<" "
                    <<it->second->sum_wx2()<<G4endl;
                }
            });
        }
      };
      //Debug!
      if ( ver > 4 ) {
          G4cout<<"Before merging: Meshes dump"<<G4endl;
          debugme();
      }
      //Write partial scorers from single ranks *before* merging
      //Do not rely on UI command to write out scorers, because rank-specific
      //files will have same file name: need to add rank # to file name
      if ( true ) {
       for ( size_t idx = 0 ; idx < scor->GetNumberOfMesh() ; ++idx) {
              const auto m = scor->GetMesh(idx);
              const auto& mn = m->GetWorldName();
              std::ostringstream fname;
              fname<<"scorer-"<<mn<<"-rank"<<rank<<".csv";
              scor->DumpAllQuantitiesToFile(mn,fname.str());
       }
      }

      //Now reduce all scorers to rank #0
      sm.Merge();

      //Debug!
      if ( ver > 4 ) {
          G4cout<<"After merging: Meshes dump"<<G4endl;
          debugme();
      }
      //For rank #0 write out the merged files
      if ( rank == 0 ) {
        for ( size_t idx = 0 ; idx < scor->GetNumberOfMesh() ; ++idx) {
            const auto m = scor->GetMesh(idx);
            const auto& mn = m->GetWorldName();
            std::ostringstream fname;
            fname<<"scorer-"<<mn<<"-merged.csv";
            scor->DumpAllQuantitiesToFile(mn,fname.str());
        }
    }
  }

  //Save histograms *before* MPI merging for rank #0
  if (rank == 0)
    {
      G4String fname("dose-rank0");
      Analysis* myana = Analysis::GetAnalysis();
      myana-> Save(fname);
      myana-> Close(false); // close file withour resetting data
    }
  //Merge of g4analysis objects
  ver=0;
  G4MPIhistoMerger hm(G4AnalysisManager::Instance());
  hm.SetVerbosity(ver);
  hm.Merge();

  //Save g4analysis objects to a file
  //NB: It is important that the save is done *after* MPI-merging of histograms

  //One can save all ranks or just rank0, chane the if
  if (true /*rank == 0*/)
    {
      std::ostringstream fname;
      fname<<"dose-rank"<<rank;
      if (rank == 0) {
        fname.str("dose-merged");
      }
      Analysis* myana = Analysis::GetAnalysis();
      myana-> Save(fname.str());
    }
    Analysis* myana = Analysis::GetAnalysis();
    myana-> Close();

    G4cout << "===================================================" << G4endl;
    G4cout << "End EndOfRunAction for master thread in rank: " << rank << G4endl;
    G4cout << "===================================================" << G4endl;

}

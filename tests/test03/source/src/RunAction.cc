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
// $Id$
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "DetectorParameters.hh"
#include "ApplicationParameters.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

using namespace ApplicationParameters;
using namespace DetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction(),
   fFileName("testAnalysis")
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Test reading
  if ( TestWrite ) TestWriting();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  if ( TestWrite ) delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TestWriting() const
{
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);

  // Book histograms, ntuple
  //
  
  // H1 properties
  
  if ( TestH1 ) {
    // ID = 0
    analysisManager->CreateH1("Edep","Edep in absorber", 
                               100, 0., 100*MeV);  
    // ID = 1
    analysisManager->CreateH1("1_TrackL",
                              "TrackL",
                              100, 0., NofZLayers*LayerZSize);
    // ID = 2  unit
    analysisManager->CreateH1("2_TrackL_Unit",
                              "TrackL unit=cm",
                              100, 0., NofZLayers*LayerZSize, "cm");
    // ID = 3  fcn (log10)
    analysisManager->CreateH1("3_TrackL_Fcn",
                              "TrackL fcn=log10", 
                              100, 1., NofZLayers*LayerZSize, "none", "log10");
    // ID = 4  binScheme
    analysisManager->CreateH1("4_TrackL_BinScheme",
                              "TrackL binScheme=log",
                              100, 1., NofZLayers*LayerZSize, "none", "none", "log");
    // ID = 5  unit & fcn (log10)
    analysisManager->CreateH1("5_TrackL_Unit_Fcn",
                              "TrackL unit=cm & fcn=log10", 
                              100, 1., NofZLayers*LayerZSize, "cm", "log10");
    // ID = 6  unit & binScheme
    analysisManager->CreateH1("6_TrackLength_Unit_BinScheme",
                              "TrackL unit=cm & binScheme=log", 
                              100, 1., NofZLayers*LayerZSize, "cm", "none", "log");

    // ID = 7  as ID =1 
    analysisManager->CreateH1("7_TrackL",
                              "Not in test",
                              100, 0., NofZLayers*LayerZSize);
    // ID = 8  as ID =1 
    analysisManager->CreateH1("8_TrackL",
                              "Not in test",
                              100, 0., NofZLayers*LayerZSize);
    // ID = 9  as ID =1 
    analysisManager->CreateH1("9_TrackL",
                              "Not in test",
                              100, 0., NofZLayers*LayerZSize);
    // Histogram dimensions
  
    // 1D histograms
    //
    analysisManager->CreateH1("PX", "PX in absorber layer #2",              // ID = 10
                              100, -10.*MeV, 10*MeV);
    analysisManager->CreateH1("PX_set", "PX in absorber layer #2",          // ID = 11
                              100, -10.*MeV, 10*MeV);
  }
  
  if ( TestH2 ) {
    // 2D histograms
    //
    analysisManager->CreateH2("PX_PY", "PX PY in absorber layer #2",         // ID = 0
                              100, -10.*MeV, 10*MeV, 100, -10.*MeV, 10*MeV);

    analysisManager->CreateH2("PX_PY_set", "PX PY in absorber layer #2",     // ID = 1
                              100, -10.*MeV, 10*MeV, 100, -10.*MeV, 10*MeV);
  }                              

  if ( TestH3 ) {
    // 3D histograms
    //
    analysisManager->CreateH3("PX_PY_PZ", "PX PY PZ in absorber layer #2",     // ID = 0
                              100, -10.*MeV, 10*MeV, 100, -10.*MeV, 10*MeV,  
                              100, -10.*MeV, 10*MeV);
    analysisManager->CreateH3("PX_PY_PZ_set", "PX PY PZ in absorber layer #2", // ID = 1 
                              100, -10.*MeV, 10*MeV, 100, -10.*MeV, 10*MeV, 
                              100, -10.*MeV, 10*MeV);
  }                              
                            
  if ( TestP1 ) {
    // 1D profiles
    //
    analysisManager->CreateP1("Edep_Vs_Z", "Longit Edep pseudo-profile",      // ID = 0
                              NofZLayers, 0, NofZLayers*LayerZSize, 0., 10.*GeV);
    analysisManager->CreateP1("Edep_Vs_Z_set", "Longit Edep pseudo-profile",  // ID = 1
                            NofZLayers, 0, NofZLayers*LayerZSize, 0., 10.*GeV);
  }                            
  
  if ( TestP2 ) {
    // 2D profiles
    analysisManager->CreateP2("Edep_Vs_XY", "Longit Edep pseudo-profile",      // ID = 0
                              NofXYLayers, - NofXYLayers*LayerXYSize/2, NofXYLayers*LayerXYSize/2,
                              NofXYLayers, - NofXYLayers*LayerXYSize/2, NofXYLayers*LayerXYSize/2,
                              0, 10*GeV);
    analysisManager->CreateP2("Edep_Vs_XY_set", "Longit Edep pseudo-profile",   // ID = 1
                              NofXYLayers, - NofXYLayers*LayerXYSize/2, NofXYLayers*LayerXYSize/2,
                              NofXYLayers, - NofXYLayers*LayerXYSize/2, NofXYLayers*LayerXYSize/2,
                              0, 10*GeV);
  }                              
 
  if ( TestH1 ) {
    // Ilegal definitions 
    // (objects are not created)
    //
    // log10 and xmin = 0
    analysisManager->CreateH1("7_TrackL_Fcn",
                              "TrackL fcn=log10", 
                              100, 0., NofZLayers*LayerZSize, "none", "log10");

    // fcn (log10) & binScheme
    analysisManager->CreateH1("8_TrackLength_Unit_Fcn_BinScheme",
                              "TrackL fcn=log10 & binScheme=log", 
                              100, 1., NofZLayers*LayerZSize, "none", "log10", "log");
  }                              

  if ( TestNtuple ) {
    // Creating ntuples
    // ntupleId = 0
    analysisManager->CreateNtuple("EDep", "Edep and Track Length");
    analysisManager->CreateNtupleDColumn("Eabs");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->CreateNtupleSColumn("Label");
    analysisManager->FinishNtuple();
    // ntupleId = 1
    analysisManager->CreateNtuple("TrackL", "Track Length");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->FinishNtuple();
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TestReading() const
{
  // Read something from another analysis file
  G4AnalysisReader* analysisReader = G4AnalysisReader::Instance();
  analysisReader->SetVerboseLevel(1);

  // Define base file name
  analysisReader->SetFileName(fFileName);
    
  if ( TestH1 ) {
     G4int h1Id = analysisReader->ReadH1(G4String("Edep"));
     if ( h1Id >= 0 ) {
       G4H1* h1 = analysisReader->GetH1(h1Id);
       if ( h1 ) {
         G4cout << "   H1: " 
                << "   mean: " << h1->mean() << " rms: " << h1->rms() << G4endl;
       }  
     }  
  }  
    
  if ( TestH2 ) {
     G4int h2Id = analysisReader->ReadH2("PX_PY");
     if ( h2Id >= 0 ) {
       G4H2* h2 =  analysisReader->GetH2(h2Id);
       if ( h2 ) {
         G4cout << "   H2: " 
                << "   mean_x: " << h2->mean_x() << " rms_x: " << h2->rms_x() 
                << G4endl << "       "
                << "   mean_y: " << h2->mean_y() << " rms_y: " << h2->rms_y() 
                << G4endl;
       }
     }           
  }    

  if ( TestH3 ) {
     G4int h3Id = analysisReader->ReadH3("PX_PY_PZ");
     if ( h3Id >= 0 ) {
       G4H3* h3 =  analysisReader->GetH3(h3Id);
       if ( h3 ) {
         G4cout << "   H3: " 
                << "   mean_x: " << h3->mean_x() << " rms_x: " << h3->rms_x() 
                << G4endl << "       "
                << "   mean_y: " << h3->mean_y() << " rms_y: " << h3->rms_y() 
                << G4endl << "       "
                << "   mean_z: " << h3->mean_z() << " rms_z: " << h3->rms_z() 
                << G4endl;
       }
     }
  }  

  if ( TestP1 ) {
     G4int p1Id = analysisReader->ReadP1("Edep_Vs_Z");
     if ( p1Id >= 0 ) {
       G4P1* p1 =  analysisReader->GetP1(p1Id);
       if ( p1 ) {
         G4cout << "   P1: " 
                << "   mean: " << p1->mean() << " rms: " << p1->rms() 
                << G4endl;
       }
     }
  }   

  if ( TestP2 ) {
     G4int p2Id = analysisReader->ReadP2("Edep_Vs_XY");
     if ( p2Id >= 0 ) {
       G4P2* p2 =  analysisReader->GetP2(p2Id);
       if ( p2 ) {
         G4cout << "   P2: " 
                << "   mean_x: " << p2->mean_x() << " rms_x: " << p2->rms_x() 
                << G4endl << "       "
                << "   mean_y: " << p2->mean_y() << " rms_y: " << p2->rms_y() 
                << G4endl;
       }
     }           
   }  

  if ( TestNtuple ) {
    // TrackL
    G4int ntupleId; 
    if ( ! G4Threading::IsMultithreadedApplication() ||
         ( G4Threading::IsMultithreadedApplication() && ! isMaster ) ) {
         // In MT application only ntuples are written only on workers
      ntupleId 
        = analysisReader->GetNtuple("TrackL");
            // file name can be omitted then the file name will be
            // define from the base file name like in writing phase
    }
    else {        
      // In MT application on master thread 
      // test reading from a file with explicitly given name
      // (specific to the output format)
#ifdef TEST_ANALYSIS_ROOT 
      ntupleId 
        = analysisReader->GetNtuple("TrackL", "testAnalysis_t1");
#endif
#ifdef TEST_ANALYSIS_XML 
      ntupleId 
        = analysisReader->GetNtuple("TrackL", "testAnalysis_nt_TrackL_t1");
#endif
#ifdef TEST_ANALYSIS_CSV 
      ntupleId 
        = analysisReader->GetNtuple("TrackL", "testAnalysis_nt_TrackL_t1");
#endif
    }
    if ( ntupleId >= 0 ) {
      G4double trackL;
      analysisReader->SetNtupleDColumn("Labs", trackL);
      G4int counter = 0;
      G4cout << "Ntuple TrackL, reading selected column Labs" << G4endl;
      while ( analysisReader->GetNtupleRow() && counter < 10 ) {
        G4cout << counter++ << "th entry: "
               << "  TrackL: " << trackL << std::endl;
      }
    }

    // EDep
    // Read column of string type
    // Note when handlind more than one ntuple, ntupleId has to be passed
    // to SetNtupleSColumn(..) and GetNtupleRow(..) calls
    if ( ! G4Threading::IsMultithreadedApplication() ||
         ( G4Threading::IsMultithreadedApplication() && ! isMaster ) ) {
         // In MT application only ntuples are written only on workers
      ntupleId 
        = analysisReader->GetNtuple("EDep");
            // file name can be omitted then the file name will be
            // define from the base file name like in writing phase
    }
    else {        
      // In MT application on master thread 
      // test reading from a file with explicitly given name
      // (specific to the output format)
#ifdef TEST_ANALYSIS_ROOT 
      ntupleId 
        = analysisReader->GetNtuple("EDep", "testAnalysis_t1");
#endif
#ifdef TEST_ANALYSIS_XML 
      ntupleId 
        = analysisReader->GetNtuple("EDep", "testAnalysis_nt_EDep_t1");
#endif
#ifdef TEST_ANALYSIS_CSV
      ntupleId 
        = analysisReader->GetNtuple("EDep", "testAnalysis_nt_EDep_t1");
#endif
    }
    if ( ntupleId >= 0 ) {
      G4String label;
      analysisReader->SetNtupleSColumn(ntupleId, "Label", label);
      G4int counter = 0;
      G4cout << "Ntuple EDep, reading selected column Label" << G4endl;
      while ( analysisReader->GetNtupleRow(ntupleId) && counter < 10 ) {
        G4cout << counter++ << "th entry: "
               << "  Label: " << label << std::endl;
      }
    }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintStatistics() const
{
  // print histogram statistics

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if(isMaster) {
    if ( TestH1 ) {
      G4H1* h1 = analysisManager->GetH1(0);
      if ( h1 ) {
        G4cout << "   H1: " 
               << "   mean: " << h1->mean() << " rms: " << h1->rms() << G4endl;
      }
    }
    if ( TestH2 ) {
      G4H2* h2 =  analysisManager->GetH2(0);
      if ( h2 ) {
         G4cout << "   H2: " 
                << "   mean_x: " << h2->mean_x() << " rms_x: " << h2->rms_x() 
                << G4endl << "       "
                << "   mean_y: " << h2->mean_y() << " rms_y: " << h2->rms_y() 
                << G4endl;
      }
    }    
    
    if ( TestH3 ) {
      G4H3* h3 =  analysisManager->GetH3(0);
      if ( h3 ) {
        G4cout << "   H3: " 
               << "   mean_x: " << h3->mean_x() << " rms_x: " << h3->rms_x() 
               << G4endl << "       "
               << "   mean_y: " << h3->mean_y() << " rms_y: " << h3->rms_y() 
               << G4endl << "       "
               << "   mean_z: " << h3->mean_z() << " rms_z: " << h3->rms_z() 
               << G4endl;
      }         
    }  

    if ( TestP1 ) {
      G4P1* p1 =  analysisManager->GetP1(0);
      if ( p1 ) {
        G4cout << "   P1: "
               << "   mean: " << p1->mean() << " rms: " << p1->rms() 
               << G4endl;
      }
    }   

    if ( TestP2 ) {
      G4P2* p2 =  analysisManager->GetP2(0);
      if ( p2 ) {
        G4cout << "   P2: "
               << "   mean_x: " << p2->mean_x() << " rms_x: " << p2->rms_x() 
               << G4endl << "       "
               << "   mean_y: " << p2->mean_y() << " rms_y: " << p2->rms_y() 
               << G4endl;
      }
    }  
  }
}    


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Test reading
  if ( TestRead ) TestReading();
  
  if ( TestWrite ) {
    // Define file name and
    // append it with "W" if both Read and Write tests are selected
    G4String fileName = fFileName;
    if ( TestRead ) fileName.append("W");
    
    // Open the output file
    G4AnalysisManager::Instance()->OpenFile(fileName);
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{

  if ( TestWrite ) {
    // print some information
    PrintStatistics();

    // save histograms & ntuple
    G4AnalysisManager::Instance()->Write();
    G4AnalysisManager::Instance()->CloseFile();
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

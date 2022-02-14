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
///////////////////////////////////////////////////////////////////////////////
// File: CCalRunAction.cc
// Description: A class for providing user actions at begin and end of run
///////////////////////////////////////////////////////////////////////////////
#include "CCalRunAction.hh"

#include "globals.hh"
#include "G4Run.hh"

#include "G4ios.hh"

#include "G4AnalysisManager.hh"
#include "G4Threading.hh"


CCalRunAction::CCalRunAction()
{
  numberOfTimeSlices = 200;
  Book();
}


CCalRunAction::~CCalRunAction()
{
}


<<<<<<< HEAD
void CCalRunAction::BeginOfRunAction(const G4Run* aRun) {

=======
void CCalRunAction::BeginOfRunAction(const G4Run* aRun) 
{
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->SetDefaultFileType("root");

  //cleanup
  G4int timeHist = analysis->GetH1Id("h300");
  for (G4int i=0; i<numberOfTimeSlices; ++i) {
    analysis->GetH1(timeHist+i)->reset();
  }

  // Open an output file
  analysis->OpenFile("ccal");
  G4cout << "********************************************" << G4endl
         << "* o/p file ccal"  << G4endl
         << "********************************************" << G4endl 
         << G4endl;
}


void CCalRunAction::EndOfRunAction(const G4Run* aRun) {

  G4cout << "### Run " << aRun->GetRunID() << " end." << G4endl;

  // Close-out analysis: save histograms and ntuple
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}   


void CCalRunAction::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);
  analysisManager->SetFirstNtupleId(1);

  // Note: merging ntuples is available only with Root output
  if ( G4Threading::IsMultithreadedApplication() ) analysisManager->SetNtupleMerging(true);
 
  // Create ntuple
  analysisManager->CreateNtuple("ntuple1", "Event info");
  for (G4int i=0;i<28;++i) {
      G4String tupleidString = "hcal" + std::to_string( i );
      analysisManager->CreateNtupleFColumn( tupleidString.c_str() );
  }
  for (G4int i=0; i<49; ++i) {
      G4String tupleidString = "ecal" + std::to_string( i );
      analysisManager->CreateNtupleFColumn( tupleidString.c_str() );
  }
  analysisManager->CreateNtupleFColumn("ELAB");
  analysisManager->CreateNtupleFColumn("XPOS");
  analysisManager->CreateNtupleFColumn("YPOS");
  analysisManager->CreateNtupleFColumn("ZPOS");
  analysisManager->CreateNtupleFColumn("EDEP");
  analysisManager->CreateNtupleFColumn("EDEC");
  analysisManager->CreateNtupleFColumn("EHDC");
  analysisManager->FinishNtuple();

  //
  //Create histograms. Save the ID of the first histogram in each block
  //
  //Energy deposit in Hcal layers 
  for (G4int i = 0; i<28; ++i) {
    G4String idString = "h" + std::to_string( i+100 );
    G4String ntupletagString = "Energy Deposit in Hcal Layer" + std::to_string( i ) + "  in GeV";
    analysisManager->CreateH1( idString.c_str(), ntupletagString.c_str(), 100, 0., 1.0 );    
  }
  // Energy deposits in Ecal towers
  for (G4int i = 0; i<49; ++i) {
    G4String idString = "h" + std::to_string( i+200 );
    G4String ntupletagString = "Energy Deposit in Ecal Tower" + std::to_string( i ) + "  in GeV";
    analysisManager->CreateH1( idString.c_str(), ntupletagString.c_str(), 100, 0., 1.0 );      
  }
  // Total energy deposit
  analysisManager->CreateH1( "h4000", "Total energy deposited in GeV", 100, 0., 100.0 );

  // Time slices          
  for (G4int i=0; i<numberOfTimeSlices; ++i){
    G4String idString = "h" + std::to_string( i+300 );
    G4String ntupletagString = "Time slice " + std::to_string( i ) + " nsec energy profile in GeV";
    analysisManager->CreateH1( idString.c_str(), ntupletagString.c_str(), 100, 0., 100.0 );   
  }

  // Profile of lateral energy deposit in Hcal
  for (G4int i = 0; i<70; ++i) {
    G4String idString = "h" + std::to_string( i+500 );
    G4String ntupletagString = "Lateral energy profile at " + std::to_string( i ) + " cm  in GeV";
    analysisManager->CreateH1( idString.c_str(), ntupletagString.c_str(), 100, 0., 10.0 );    
  }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

}   


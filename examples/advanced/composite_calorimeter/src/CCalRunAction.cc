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

#include "CCalAnalysis.hh"


void CCalRunAction::BeginOfRunAction(const G4Run* aRun) 
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  numberOfTimeSlices = 200;

  if (aRun->GetRunID() == 0) //first run
    Book();
      
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //cleanup
  G4int timeHist = analysis->GetH1Id("h300");
  for (G4int i=0; i<numberOfTimeSlices; ++i) {
    analysis->GetH1(timeHist+i)->reset();
  }
}


void CCalRunAction::EndOfRunAction(const G4Run* aRun) 
{
  G4cout << "### Run " << aRun->GetRunID() << " end." << G4endl;
}   

void CCalRunAction::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);
  analysisManager->SetFirstNtupleId(1);
  
  // Open an output file
  analysisManager->OpenFile("ccal");
  G4cout << "********************************************" << G4endl
         << "* o/p file ccal"  << G4endl
         << "********************************************" << G4endl 
         << G4endl;
  
  
  // Create a tuple :
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

  // Time profile 
  analysisManager->CreateH1( "h901", "Time Profile in Sensitive Detector", 200, 0., 200. );
  analysisManager->CreateH1( "h902", "Time Profile in Sensitive+Passive",  200, 0., 200. );

  return;
}

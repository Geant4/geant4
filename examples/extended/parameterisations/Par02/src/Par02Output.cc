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
/// \file Par02Output.cc
/// \brief Implementation of the Par02Output class

#include "Par02Output.hh"
#include "Par02EventInformation.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02Output* Par02Output::fPar02Output = 0;
G4ThreadLocal G4int Par02Output::fCurrentNtupleId = 0;
G4ThreadLocal G4int Par02Output::fCurrentID = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02Output::Par02Output() : fFileNameWithRunNo( false ) {
  fFileName = "DefaultOutput.root";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02Output::~Par02Output() {
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02Output* Par02Output::Instance() {
  if ( ! fPar02Output ) {
    fPar02Output = new Par02Output();
  }
  return fPar02Output;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02Output::SetFileName( G4String aName ) {
  fFileName = aName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02Output::AppendName( G4bool aApp ) {
  fFileNameWithRunNo = aApp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Par02Output::GetFileName() {
  return fFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02Output::StartAnalysis( G4int aRunID ) {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( fFileNameWithRunNo ) {
    fFileName +=  "_run";
    fFileName += G4UIcommand::ConvertToString( aRunID );
  }
  analysisManager->SetVerboseLevel( 1 );
  analysisManager->SetFileName( fFileName );
  analysisManager->OpenFile( fFileName );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02Output::EndAnalysis() {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02Output::CreateNtuples() {
  const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();
  G4String evName = "Event_";
  evName += G4UIcommand::ConvertToString( event->GetEventID() );
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  fCurrentNtupleId = analysisManager->CreateNtuple( evName, evName );

  analysisManager->CreateNtupleIColumn( "particleID" );  // column Id = 0
  analysisManager->CreateNtupleIColumn( "PID" );  // column Id = 1
  analysisManager->CreateNtupleDColumn( "MC_pX" );  // column Id = 2
  analysisManager->CreateNtupleDColumn( "MC_pY" ); // column Id = 3
  analysisManager->CreateNtupleDColumn( "MC_pZ" ); // column Id = 4

  analysisManager->CreateNtupleDColumn( "tracker_res" ); // column Id = 5
  analysisManager->CreateNtupleDColumn( "tracker_eff" ); // column Id = 6
  analysisManager->CreateNtupleDColumn( "tracker_pX" );  // column Id = 7
  analysisManager->CreateNtupleDColumn( "tracker_pY" ); // column Id = 8
  analysisManager->CreateNtupleDColumn( "tracker_pZ" ); // column Id = 9

  analysisManager->CreateNtupleDColumn( "emcal_res" ); // column Id = 10
  analysisManager->CreateNtupleDColumn( "emcal_eff" );  // column Id = 11
  analysisManager->CreateNtupleDColumn( "emcal_X" );  // column Id = 12
  analysisManager->CreateNtupleDColumn( "emcal_Y" ); // column Id = 13
  analysisManager->CreateNtupleDColumn( "emcal_Z" ); // column Id = 14
  analysisManager->CreateNtupleDColumn( "emcal_E" ); // column Id = 15

  analysisManager->CreateNtupleDColumn( "hcal_res" ); // column Id = 16
  analysisManager->CreateNtupleDColumn( "hcal_eff" ); // column Id = 17
  analysisManager->CreateNtupleDColumn( "hcal_X" );  // column Id = 18
  analysisManager->CreateNtupleDColumn( "hcal_Y" ); // column Id = 19
  analysisManager->CreateNtupleDColumn( "hcal_Z" ); // column Id = 20
  analysisManager->CreateNtupleDColumn( "hcal_E" );  // column Id = 21

  analysisManager->FinishNtuple( fCurrentNtupleId );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02Output::CreateHistograms()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->CreateH1( "Pdiff", "momentum smeared in tracker", 100, 0.8, 1.2 );
  analysisManager->SetH1XAxisTitle( 0, "p_{smeared}/p_{true}" );
  analysisManager->SetH1YAxisTitle( 0, "Entries" );
  analysisManager->CreateH1( "EMCalEdiff", "energy smeared in EMCal", 100, 0.8, 1.2 );
  analysisManager->SetH1XAxisTitle( 1, "E_{smeared}/E_{true}" );
  analysisManager->SetH1YAxisTitle( 1, "Entries" );
  analysisManager->CreateH1( "HCalEdiff", "energy smeared in HCal", 100, 0.0, 2.0 );
  analysisManager->SetH1XAxisTitle( 2, "E_{smeared}/E_{true}" );
  analysisManager->SetH1YAxisTitle( 2, "Entries" );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02Output::SaveTrack( SaveType aWhatToSave, G4int aPartID,  G4int aPDG,
                             G4ThreeVector aVector, G4double aResolution, 
                             G4double aEfficiency, G4double aEnergy ) {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  switch( aWhatToSave ) {
    case Par02Output::eNoSave :
      break;
    case Par02Output::eSaveMC : {
      analysisManager->FillNtupleIColumn( fCurrentNtupleId, 0, aPartID );
      analysisManager->FillNtupleIColumn( fCurrentNtupleId, 1, aPDG );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 2, aVector.x() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 3, aVector.y() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 4, aVector.z() );
      fCurrentID = aPartID;
      break;
    }
    case Par02Output::eSaveTracker : {
      if ( aPartID != fCurrentID ) G4cout << 
        " Wrong particle - trying to save Tracker information of different particle"
        << G4endl;
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 5, aResolution );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 6, aEfficiency );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 7, aVector.x() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 8, aVector.y() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 9, aVector.z() );
      break;
    }
    case Par02Output::eSaveEMCal : {
      if ( aPartID != fCurrentID ) G4cout <<
        " Wrong particle - trying to save EMCal information of different particle"
        << G4endl;
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 10, aResolution );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 11, aEfficiency );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 12, aVector.x() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 13, aVector.y() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 14, aVector.z() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 15, aEnergy );
      break;
    }
    case Par02Output::eSaveHCal : {
      if ( aPartID != fCurrentID ) G4cout <<
        " Wrong particle - trying to save HCal information of different particle"
        << G4endl;
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 16, aResolution );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 17, aEfficiency );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 18, aVector.x() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 19, aVector.y() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 20, aVector.z() );
      analysisManager->FillNtupleDColumn( fCurrentNtupleId, 21, aEnergy );
      analysisManager->AddNtupleRow( fCurrentNtupleId );
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02Output::FillHistogram( G4int aHistNo, G4double aValue ) const {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1( aHistNo, aValue );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


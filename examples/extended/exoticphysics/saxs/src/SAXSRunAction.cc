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
/// \file SAXSRunAction.cc
/// \brief Implementation of the SAXSRunAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "SAXSRunAction.hh"
#include "SAXSRunActionMessenger.hh"
#include "SAXSDetectorConstruction.hh"
#include "SAXSAnalysis.hh"

#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSRunAction::SAXSRunAction():
  G4UserRunAction(),
  fIsFileOpened(false)
{
  //set the number of events to print progress during the simulation
  G4RunManager::GetRunManager()->SetPrintProgress(10000);
  
  //define the messenger  
  fMessenger = new SAXSRunActionMessenger(this);
  
  //default output filename (can be set through macro)
  fFileName = "output";
  
  //Create the analysis manager  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetNtupleMerging(true); //only for root
  analysisManager->SetVerboseLevel(1);
  
  //Creating the SD scoring ntuple
  analysisManager->CreateNtuple("part","Particle");
  analysisManager->CreateNtupleDColumn("e");
  analysisManager->CreateNtupleDColumn("posx");
  analysisManager->CreateNtupleDColumn("posy");
  analysisManager->CreateNtupleDColumn("posz");
  analysisManager->CreateNtupleDColumn("momx");
  analysisManager->CreateNtupleDColumn("momy");
  analysisManager->CreateNtupleDColumn("momz");
  analysisManager->CreateNtupleDColumn("t");
  analysisManager->CreateNtupleIColumn("type");
  analysisManager->CreateNtupleIColumn("trackID");
  analysisManager->CreateNtupleIColumn("NRi");
  analysisManager->CreateNtupleIColumn("NCi");
  analysisManager->CreateNtupleIColumn("NDi");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleDColumn("weight");
  analysisManager->FinishNtuple();
    
  //Creating ntuple for scattering
  analysisManager->CreateNtuple("scatt","Scattering");
  analysisManager->CreateNtupleIColumn("processID");
  analysisManager->CreateNtupleDColumn("e");
  analysisManager->CreateNtupleDColumn("theta");
  analysisManager->CreateNtupleDColumn("weight");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSRunAction::~SAXSRunAction() {
   //write results
  if (fIsFileOpened)
    {
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
      analysisManager->Write();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSRunAction::BeginOfRunAction(const G4Run*)
{
  //open the output file         
  if (!fIsFileOpened)
    {
      G4AnalysisManager::Instance()->OpenFile(fFileName);
      fIsFileOpened = true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSRunAction::EndOfRunAction(const G4Run* run)
{               
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
      
  //print
  if (IsMaster()) {
    G4cout
         << G4endl
         << "--------------------End of Global Run-----------------------"
         << G4endl
         << " The run had " << nofEvents << " events";
   } else {
     G4cout
          << G4endl
          << "--------------------End of Local Run------------------------"
          << G4endl
          << " The run had " << nofEvents << " events";
  }      
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSRunAction::SetFileName(G4String filename)
{
  //method to set the output filename
  if (filename != "") fFileName = filename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


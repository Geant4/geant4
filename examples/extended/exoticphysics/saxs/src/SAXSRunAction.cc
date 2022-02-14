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
#include "G4AnalysisManager.hh"

#include "SAXSRunAction.hh"
#include "SAXSRunActionMessenger.hh"
#include "SAXSDetectorConstruction.hh"

#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSRunAction::SAXSRunAction(): G4UserRunAction()
{
  //define the messenger  
  fMessenger = new SAXSRunActionMessenger(this);
  
  //default output filename (can be set through macro)
  fFileName = "output";
  
  //Create the analysis manager  
  fAnalysisManager = G4AnalysisManager::Instance();

  fAnalysisManager->SetDefaultFileType("root");
  fAnalysisManager->SetFileName(fFileName);
  fAnalysisManager->SetNtupleMerging(true); //only for root
  fAnalysisManager->SetVerboseLevel(1);
  
  //Creating the SD scoring ntuple
  fAnalysisManager->CreateNtuple("part","Particle");
  fAnalysisManager->CreateNtupleDColumn("e");
  fAnalysisManager->CreateNtupleDColumn("posx");
  fAnalysisManager->CreateNtupleDColumn("posy");
  fAnalysisManager->CreateNtupleDColumn("posz");
  fAnalysisManager->CreateNtupleDColumn("momx");
  fAnalysisManager->CreateNtupleDColumn("momy");
  fAnalysisManager->CreateNtupleDColumn("momz");
  fAnalysisManager->CreateNtupleDColumn("t");
  fAnalysisManager->CreateNtupleIColumn("type");
  fAnalysisManager->CreateNtupleIColumn("trackID");
  fAnalysisManager->CreateNtupleIColumn("NRi");
  fAnalysisManager->CreateNtupleIColumn("NCi");
  fAnalysisManager->CreateNtupleIColumn("NDi");
  fAnalysisManager->CreateNtupleIColumn("eventID");
  fAnalysisManager->CreateNtupleDColumn("weight");
  fAnalysisManager->FinishNtuple();
    
  //Creating ntuple for scattering
  fAnalysisManager->CreateNtuple("scatt","Scattering");
  fAnalysisManager->CreateNtupleIColumn("processID");
  fAnalysisManager->CreateNtupleDColumn("e");
  fAnalysisManager->CreateNtupleDColumn("theta");
  fAnalysisManager->CreateNtupleDColumn("weight");
  fAnalysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSRunAction::~SAXSRunAction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSRunAction::BeginOfRunAction(const G4Run*)
{
  //open the output file         
  if (!fIsFileOpened)
    {
      fAnalysisManager->OpenFile(fFileName);
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
  if (fIsFileOpened) {
    fAnalysisManager->Write();
    fAnalysisManager->CloseFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSRunAction::SetFileName(const G4String& filename)
{
  //method to set the output filename
  if (filename != "") fFileName = filename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


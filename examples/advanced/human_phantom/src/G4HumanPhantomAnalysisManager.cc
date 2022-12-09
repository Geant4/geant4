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
/*
Author: Susanna Guatelli, University of Wollongong, Australia
*/
// The class G4HumanPhantomAnalysisManager creates and manages ntuples

// The analysis was included in this application following the extended Geant4
// example analysis/AnaEx01

#include "G4HumanPhantomAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

G4HumanPhantomAnalysisManager::G4HumanPhantomAnalysisManager()
{
fFactoryOn = false;

// Initialization ntuple
  for (G4int k=0; k<MaxNtCol; k++) {
    fNtColId[k] = 0;
  }  
}

void G4HumanPhantomAnalysisManager::book() 
{  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  analysisManager->SetVerboseLevel(2);
 
  // Create a root file
  G4String fileName = "human_phantom.root";

  // Create directories  
  analysisManager->SetNtupleDirectoryName("human_phantom_ntuple");
  
  G4bool fileOpen = analysisManager->OpenFile(fileName);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " 
           << fileName
           << G4endl;
    return;
  }

  //creating a ntuple, containg 3D energy deposition in the phantom
  analysisManager->SetFirstNtupleId(1);
  analysisManager->CreateNtuple("1", "3Dedep");
  fNtColId[0] = analysisManager->CreateNtupleDColumn("organID");
  fNtColId[1] = analysisManager->CreateNtupleDColumn("edep");

  analysisManager->FinishNtuple();
  
  fFactoryOn = true;    
}

void G4HumanPhantomAnalysisManager::FillNtupleWithEnergyDeposition(G4int organ,G4double energyDep)
{
  if (energyDep !=0)
 {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(1, fNtColId[0], organ);
  analysisManager->FillNtupleDColumn(1, fNtColId[1], energyDep);
  analysisManager->AddNtupleRow(1);  
  G4cout << "Analysis: organ " << organ << " edep: "<< energyDep << G4endl;  
}
 }

void G4HumanPhantomAnalysisManager::save() 
{  
 if (fFactoryOn) 
   {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();  
    analysisManager->Clear();  
    fFactoryOn = false;
   }
}













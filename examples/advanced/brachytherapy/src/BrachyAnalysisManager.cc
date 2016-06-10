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
Author: Susanna Guatelli
*/
// The class BrachyAnalysisManager creates and manages histograms and ntuples

// The analysis was included in this application following the extended Geant4
// example analysis/AnaEx01

#include <stdlib.h>
#include "BrachyAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

BrachyAnalysisManager::BrachyAnalysisManager()
{

factoryOn = false;

// Initialization
// histograms
for (G4int k=0; k<MaxHisto; k++) fHistId[k] = 0;

// Initialization ntuple
  for (G4int k=0; k<MaxNtCol; k++) {
    fNtColId[k] = 0;
  }  

primaryParticleSpectrum = 0;

}

BrachyAnalysisManager::~BrachyAnalysisManager() 
{ 
}

void BrachyAnalysisManager::book() 
{  
  G4AnalysisManager* AnalysisManager = G4AnalysisManager::Instance();
  
  AnalysisManager->SetVerboseLevel(2);
 
  // Create a root file
  G4String fileName = "brachytherapy.root";

  // Create directories  
  AnalysisManager->SetHistoDirectoryName("brachy_histo");
  AnalysisManager->SetNtupleDirectoryName("brachy_ntuple");
  

  G4bool fileOpen = AnalysisManager->OpenFile(fileName);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " 
           << fileName[1] 
           << G4endl;
    return;
  }

//creating a 1D histograms
  AnalysisManager->SetFirstHistoId(1);
  
  // Histogram containing the primary particle energy (MeV) 
  fHistId[1]  = AnalysisManager -> CreateH1("1", 
                                            "Initial energy", 
                                            1000, 0., 1000.);

  //Parameters of CreateH1: histoID, histo name, bins' number, xmin, xmax 
  
  primaryParticleSpectrum = AnalysisManager-> GetH1(fHistId[1]);
  
  //creating a ntuple, containg 3D energy deposition in the phantom
  AnalysisManager -> CreateNtuple("1", "3Dedep");
  fNtColId[0] = AnalysisManager->CreateNtupleDColumn("xx");
  fNtColId[1] = AnalysisManager->CreateNtupleDColumn("yy");
  fNtColId[2] = AnalysisManager->CreateNtupleDColumn("zz");
  fNtColId[3] = AnalysisManager->CreateNtupleDColumn("edep");
  AnalysisManager->FinishNtuple();
  
 factoryOn = true;    
}

void BrachyAnalysisManager::FillPrimaryParticleHistogram(G4double primaryParticleEnergy)
{
 // 1DHistogram: energy spectrum of primary particles  
  primaryParticleSpectrum -> fill(primaryParticleEnergy);
}
 
void BrachyAnalysisManager::FillNtupleWithEnergyDeposition(G4double xx,
                                                     G4double yy, 
                                                     G4double zz,
                                                     G4double energyDep)
{
  G4AnalysisManager* AnalysisManager = G4AnalysisManager::Instance();
  AnalysisManager->FillNtupleDColumn(fNtColId[0], xx);
  AnalysisManager->FillNtupleDColumn(fNtColId[1], yy);
  AnalysisManager->FillNtupleDColumn(fNtColId[2], zz);
  AnalysisManager->FillNtupleDColumn(fNtColId[3], energyDep);
  AnalysisManager->AddNtupleRow();  
}

void BrachyAnalysisManager::save() 
{  
 if (factoryOn) 
   {
    G4AnalysisManager* AnalysisManager = G4AnalysisManager::Instance();    
    AnalysisManager->Write();
    AnalysisManager->CloseFile();  
      
    delete G4AnalysisManager::Instance();
    factoryOn = false;
   }
}













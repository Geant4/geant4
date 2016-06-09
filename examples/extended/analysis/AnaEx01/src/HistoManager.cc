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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fileName[0] = "AnaEx01";
  factoryOn = false;
  
  // histograms
  for (G4int k=0; k<MaxHisto; k++) {
    fHistId[k] = 0;
    fHistPt[k] = 0;    
  }
  // ntuple
  for (G4int k=0; k<MaxNtCol; k++) {
    fNtColId[k] = 0;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(2);
  G4String extension = analysisManager->GetFileType();
  fileName[1] = fileName[0] + "." + extension;
      
  // Create directories 
  analysisManager->SetHistoDirectoryName("histo");
  analysisManager->SetNtupleDirectoryName("ntuple");
    
  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile(fileName[0]);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " << fileName[1] 
           << G4endl;
    return;
  }
  
  // create selected histograms
  //
  analysisManager->SetFirstHistoId(1);

  fHistId[1] = analysisManager->CreateH1("1","Edep in absorber (MeV)",
                                              100, 0., 800*MeV);
  fHistPt[1] = analysisManager->GetH1(fHistId[1]);
                                           
  fHistId[2] = analysisManager->CreateH1("2","Edep in gap (MeV)",
                                              100, 0., 100*MeV);
  fHistPt[2] = analysisManager->GetH1(fHistId[2]);
                                           
  fHistId[3] = analysisManager->CreateH1("3","trackL in absorber (mm)",
                                              100, 0., 1*m);
  fHistPt[3] = analysisManager->GetH1(fHistId[3]);
                                           
  fHistId[4] = analysisManager->CreateH1("4","trackL in gap (mm)",
                                              100, 0., 50*cm);
  fHistPt[4] = analysisManager->GetH1(fHistId[4]);
                                  
  // Create 1 ntuple
  //    
  analysisManager->CreateNtuple("101", "Edep and TrackL");
  fNtColId[0] = analysisManager->CreateNtupleDColumn("Eabs");
  fNtColId[1] = analysisManager->CreateNtupleDColumn("Egap");
  fNtColId[2] = analysisManager->CreateNtupleDColumn("Labs");
  fNtColId[3] = analysisManager->CreateNtupleDColumn("Lgap");
  analysisManager->FinishNtuple();
  
  factoryOn = true;       
  G4cout << "\n----> Histogram Tree is opened in " << fileName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{
  if (factoryOn) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();  
    G4cout << "\n----> Histogram Tree is saved in " << fileName[1] << G4endl;
      
    delete G4AnalysisManager::Instance();
    factoryOn = false;
  }                    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << "does note xist; xbin= " << xbin << " w= " << weight << G4endl;
    return;
  }

  if (fHistPt[ih]) fHistPt[ih]->fill(xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << "  fac= " << fac << G4endl;
    return;
  }

  if (fHistPt[ih]) fHistPt[ih]->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap,
                              G4double trackLAbs, G4double trackLGap)
{                
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(fNtColId[0], energyAbs);
  analysisManager->FillNtupleDColumn(fNtColId[1], energyGap);
  analysisManager->FillNtupleDColumn(fNtColId[2], trackLAbs);
  analysisManager->FillNtupleDColumn(fNtColId[2], trackLGap);
  analysisManager->AddNtupleRow();  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
  if(factoryOn) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
       << " EAbs : mean = " << G4BestUnit(fHistPt[1]->mean(), "Energy") 
               << " rms = " << G4BestUnit(fHistPt[1]->rms(),  "Energy") 
               << G4endl;
    G4cout                
       << " EGap : mean = " << G4BestUnit(fHistPt[2]->mean(), "Energy") 
               << " rms = " << G4BestUnit(fHistPt[2]->rms(),  "Energy") 
               << G4endl;
    G4cout 
       << " LAbs : mean = " << G4BestUnit(fHistPt[3]->mean(), "Length") 
               << " rms = " << G4BestUnit(fHistPt[3]->rms(),  "Length") 
               << G4endl;
    G4cout 
       << " LGap : mean = " << G4BestUnit(fHistPt[4]->mean(), "Length") 
               << " rms = " << G4BestUnit(fHistPt[4]->rms(),  "Length") 
               << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



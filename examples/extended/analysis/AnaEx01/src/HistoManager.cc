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
// $Id: HistoManager.cc,v 1.1 2010-09-16 16:26:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fileName = "AnaEx01";

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
  
  // Open an output file
  //
  analysisManager->OpenFile(fileName);

  // create selected histograms
  //
  analysisManager->SetFirstHistoId(1);
  
  fHistId[1] = analysisManager->CreateH1("1","Edep in absorber",
                                              100, 0., 800*MeV);
  fHistPt[1] = analysisManager->GetH1(fHistId[1]);
  					 
  fHistId[2] = analysisManager->CreateH1("2","Edep in gap",
                                              100, 0., 100*MeV);
  fHistPt[2] = analysisManager->GetH1(fHistId[2]);
  					 
  fHistId[3] = analysisManager->CreateH1("3","trackL in absorber",
                                              100, 0., 1*m);
  fHistPt[3] = analysisManager->GetH1(fHistId[3]);
  					 
  fHistId[4] = analysisManager->CreateH1("4","trackL in gap",
                                              100, 0., 50*cm);
  fHistPt[4] = analysisManager->GetH1(fHistId[4]);
  				
  // Create 1 ntuple
  //
  analysisManager->SetFirstNtupleId(0);
    
  analysisManager->CreateNtuple("101", "Edep and TrackL");
  fNtColId[0] = analysisManager->CreateNtupleDColumn("Eabs");
  fNtColId[1] = analysisManager->CreateNtupleDColumn("Egap");
  fNtColId[2] = analysisManager->CreateNtupleDColumn("Labs");
  fNtColId[3] = analysisManager->CreateNtupleDColumn("Lgap");
  analysisManager->FinishNtuple();
          
  G4cout << "\n----> Histogram Tree is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->Write();
  analysisManager->CloseFile();
  
  G4cout << "\n----> Histogram Tree is saved in " << fileName << G4endl;
      
  // complete cleanup
  //
  delete G4AnalysisManager::Instance();         
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

void HistoManager::FillNtuple(G4int column, G4double value)
{
  if (column > 3) {
    G4cout << "---> warning from HistoManager::FillNtuple : " 
     << "column=" << column << " value=" << value << G4endl;
    return;
  }
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(fNtColId[column], value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddRowNtuple()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->AddNtupleRow();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
  if(fHistPt[1]) {
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



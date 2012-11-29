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
/// \file electromagnetic/TestEm6/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Open an output file
  //
  G4String fileName = "testem6";
  analysisManager->OpenFile(fileName);    
  analysisManager->SetVerboseLevel(1);
  G4String extension = analysisManager->GetFileType();
  fileName = fileName + "." + extension;
    
  // Creating histograms
  //
  analysisManager->SetFirstHistoId(1);  
  analysisManager->CreateH1("1","1/(1+(theta+[g]+)**2)",100, 0 ,1.);
  analysisManager->CreateH1("2","log10(theta+ [g]+)",   100,-3.,1.);
  analysisManager->CreateH1("3","log10(theta- [g]-)",   100,-3.,1.);
  analysisManager->CreateH1("4","log10(theta+ [g]+ -theta- [g]-)", 100,-3.,1.);
  analysisManager->CreateH1("5","xPlus" ,100,0.,1.);
  analysisManager->CreateH1("6","xMinus",100,0.,1.);
    
  G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  analysisManager->Write();
  analysisManager->CloseFile();

  delete G4AnalysisManager::Instance();    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();
  fProcCounter = new ProcessesCount;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName)
{
  //does the process  already encounted ?
  size_t nbProc = fProcCounter->size();
  size_t i = 0;
  while ((i<nbProc)&&((*fProcCounter)[i]->GetName()!=procName)) i++;
  if (i == nbProc) fProcCounter->push_back( new OneProcessCount(procName));
  
  (*fProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
  //total number of process calls
  G4double countTot = 0.;
  G4cout << "\n Number of process calls --->";
  for (size_t i=0; i< fProcCounter->size();i++) {
        G4String procName = (*fProcCounter)[i]->GetName();
        if (procName != "Transportation") {
          G4int count    = (*fProcCounter)[i]->GetCounter(); 
          G4cout << "\t" << procName << " : " << count;
          countTot += count;
        }
  }
  G4cout << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

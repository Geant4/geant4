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
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "FCALRunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

#include "FCALAnalysisManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALRunAction::FCALRunAction()
{
    // Get/create analysis manager
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    
    // Open an output file
    G4cout << "Opening output file " << man->GetFileName() << " ... ";
    man->SetFileName("fcal");
    man->SetFirstHistoId(1);
    G4cout << " done" << G4endl;
    
    // Create histogram(s)
    man->CreateH1("1","Number of Out Of World", 100,0.,10.);
    man->CreateH1("2","Number of Secondaries", 100,0.,100.);
    man->CreateH1("3","Electromagnetic Energy/MeV", 100,0.,100.);
    man->CreateH1("4","hadronic Energy/MeV", 100,10.,60.);
    

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALRunAction::~FCALRunAction()
{
 delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALRunAction::BeginOfRunAction(const G4Run* aRun)
{
 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

//  if (G4VVisManager::GetConcreteInstance())
//    {
//      G4UImanager* UI = G4UImanager::GetUIpointer(); 
//      UI->ApplyCommand("/vis/scene/notifyHandlers");
//    }
    
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALRunAction::EndOfRunAction(const G4Run* )
{

//  if (G4VVisManager::GetConcreteInstance()) {
//     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
//  }

  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  // Complete clean-up
  delete G4AnalysisManager::Instance();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







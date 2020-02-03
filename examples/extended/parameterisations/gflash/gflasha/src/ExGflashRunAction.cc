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
/// \file ExGflashRunAction.cc
/// \brief Implementation of the ExGflashRunAction class
//
#include "ExGflashRunAction.hh"
#include "G4Run.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashRunAction::ExGflashRunAction(ExGflashDetectorConstruction* det)
  //  : G4UserRunAction(), fRunID(0), fDetector(det)
  : G4UserRunAction(), fDetector(det)
{
  fHistoManager = new ExGflashHistoManager(fDetector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashRunAction::~ExGflashRunAction()
{
  delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashRunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //histograms file
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  if (analysis->IsActive()) analysis->OpenFile();
  
  if (IsMaster()) fRunTimer.Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashRunAction::EndOfRunAction(const G4Run* aRun)
{
  if (IsMaster()) {
    // For MT we need Timer Merge
    fRunTimer.Stop();
    G4cout << G4endl;
    G4cout << "******************************************";
    G4cout << G4endl;
    G4cout << "Run Real Elapsed Time is: "<< fRunTimer.GetRealElapsed();
    G4cout << G4endl;
    G4cout << "Run System Elapsed Time: " << fRunTimer.GetSystemElapsed();
    G4cout << G4endl;
    G4cout << "Run GetUserElapsed Time: " << fRunTimer.GetUserElapsed();
    G4cout << G4endl;
    G4cout << "******************************************"<< G4endl;
    
    G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  }
    //save histograms
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();   
    if (analysis->IsActive()) {   
      analysis->Write();
    analysis->CloseFile();
  }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
/// \file hadronic/Hadr01/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
// -------------------------------------------------------------
//  
//    GEANT4 class file
//    RunAction
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "HistoManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4Timer.hh"
#include "G4NuclearLevelData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{
  timer = new G4Timer();
  timer->Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int id = aRun->GetRunID();
  G4cout << "### Run " << id << " start" << G4endl;
  timer->Stop();
  G4cout << "RunAction::BeginOfRunAction:  "  << *timer << G4endl;

  delete timer;
  timer = new G4Timer();
  timer->Start();
  G4NuclearLevelData::GetInstance();
  timer->Stop();
  G4cout << "NuclearData:  "  << *timer << G4endl;

  (HistoManager::GetPointer())->BeginOfRun();

#ifdef G4VIS_USE
  G4UImanager* UI = G4UImanager::GetUIpointer();

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  }
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{

  G4cout << "RunAction: End of run actions are started" << G4endl;
  (HistoManager::GetPointer())->EndOfRun();

#ifdef G4VIS_USE
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
/// \file runAndEvent/RE03/src/RE03ActionInitialization.cc
/// \brief Implementation of the RE03ActionInitialization class
//
//
//

#include "RE03ActionInitialization.hh"
#include "G4RunManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "RE03PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE03ActionInitialization::RE03ActionInitialization()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE03ActionInitialization::~RE03ActionInitialization()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE03ActionInitialization::BuildForMaster() const
{
  auto* runManager = G4RunManager::GetRunManager();

  // check if in sub-event parallel mode
  if(runManager->GetRunManagerType()==G4RunManager::subEventMasterRM)
  {
    // defining size of sub-event
    runManager->RegisterSubEventType(0,100);
    // sending electron, positron and gamma to worker thread
    runManager->SetDefaultClassification(G4Electron::Definition(),fSubEvent_0);
    runManager->SetDefaultClassification(G4Positron::Definition(),fSubEvent_0);
    runManager->SetDefaultClassification(G4Gamma::Definition(),fSubEvent_0);

    // primary generator action is defined to the master thread 
    SetUserAction(new RE03PrimaryGeneratorAction);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE03ActionInitialization::Build() const
{
  auto* runManager = G4RunManager::GetRunManager();

  // check if NOT in the sub-event parallel mode
  if(runManager->GetRunManagerType()!=G4RunManager::subEventWorkerRM)
  {
    // primary generator action is defined to the worker thread
    SetUserAction(new RE03PrimaryGeneratorAction);
  }
}

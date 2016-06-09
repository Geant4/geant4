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
// Code developed by:
//  S.Larsson
//
//    ******************************
//    *                            *
//    *    PurgMagRunAction.cc     *
//    *                            *
//    ******************************
//
// $Id$
//

#include "PurgMagRunAction.hh"

#include "G4SteppingManager.hh"

#include "G4Run.hh"
#include "G4Material.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#ifdef G4ANALYSIS_USE
#include "PurgMagAnalysisManager.hh"
#endif

#include <assert.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagRunAction::PurgMagRunAction(PurgMagDetectorConstruction* det)
:Detector(det)
{   
  saveRndm = 1;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagRunAction::~PurgMagRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagRunAction::BeginOfRunAction(const G4Run* aRun)
{  

#ifdef G4ANALYSIS_USE
PurgMagAnalysisManager* analysis = PurgMagAnalysisManager::getInstance();
   analysis->book();
#endif  

  G4cout << "---> Run " << aRun->GetRunID() << " start." << G4endl;
  
  
  // save Rndm status
  if (saveRndm > 0)
    { CLHEP::HepRandom::showEngineStatus();
      CLHEP::HepRandom::saveEngineStatus("beginOfRun.rndm");
    }
       
 
  //Drawing
  // 
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagRunAction::EndOfRunAction(const G4Run* aRun)
{     
  
#ifdef G4ANALYSIS_USE
  PurgMagAnalysisManager* analysis = PurgMagAnalysisManager::getInstance();
#endif
  
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  
#ifdef G4ANALYSIS_USE      
      analysis->finish();
#endif

  //drawing
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
       
  // save Rndm status
  if (saveRndm == 1)
    { CLHEP::HepRandom::showEngineStatus();
      CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
    }                         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....








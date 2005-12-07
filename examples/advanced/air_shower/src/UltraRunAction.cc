//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraRunAction.cc
//    ****************************************************
//
//    RunAction class for Ultra; it has also a Messenger Class
//
#include "UltraRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "UltraRunActionMessenger.hh"
#include "Randomize.hh"
#ifdef G4ANALYSIS_USE
#include "UltraAnalysisManager.hh"
#endif

#include "Randomize.hh"
#include "G4ios.hh"

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraRunAction::UltraRunAction()
{

   G4int seed = 1;                                         // RANLUX seed
   G4int luxury = 3;                                       // RANLUX luxury level (3 is default)
   CLHEP::HepRandom::setTheSeed(seed,luxury);
   CLHEP::HepRandom::showEngineStatus();

  theRunActMessenger = new UltraRunActionMessenger(this);
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraRunAction::~UltraRunAction()
{ 
  delete theRunActMessenger;
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraRunAction::BeginOfRunAction(const G4Run* aRun)
{

#ifdef G4ANALYSIS_USE
  UltraAnalysisManager* analysis = UltraAnalysisManager::getInstance();
  analysis->book();
#endif

  G4cout << "### Run " << aRun->GetRunID() << " start..." << G4endl;

  // save Rndm status
  if (saveRndm > 0)
    { CLHEP::HepRandom::showEngineStatus();
      CLHEP::HepRandom::saveEngineStatus("beginOfRun.rndm");
    }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraRunAction::EndOfRunAction(const G4Run* aRun)
{

  // save Rndm status
  if (saveRndm == 1)
    { CLHEP::HepRandom::showEngineStatus();
      CLHEP::HepRandom::saveEngineStatus("endOfRun.rndm");
    }

  // Write histograms to file
  G4cout << "Close and Save Histograms" << G4endl;
  G4cout << "### Run " << aRun->GetRunID() << " ended." << G4endl;

  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis/drawVolume");
  } 

#ifdef G4ANALYSIS_USE
  UltraAnalysisManager* analysis = UltraAnalysisManager::getInstance();
  analysis->finish();
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

    void UltraRunAction::MySetRunID(G4int myrun)
{
     runID = myrun ;
}

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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISARunAction class                                              *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************


#include "LISARunAction.hh"

#include "LISARunActionMessenger.hh"
#include "LISAEventAction.hh"

#ifdef G4ANALYSIS_USE
#include "LISAAnalysisManager.hh"
#endif

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "globals.hh"
#include <strstream>

#include "Randomize.hh"
#include <time.h>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISARunAction::LISARunAction() {

  // defaults
  autoSeed = false;

  // create messenger
  runMessenger = new LISARunActionMessenger(this);

#ifdef G4ANALYSIS_USE
  // Book histograms and ntuples
  LISAAnalysisManager* analysis = LISAAnalysisManager::getInstance();
  analysis->Init();
#endif


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISARunAction::~LISARunAction() {

#ifdef G4ANALYSIS_USE
  // delete analysis
  LISAAnalysisManager* analysis = LISAAnalysisManager::getInstance();
  analysis->Dispose();
#endif

  delete runMessenger;


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISARunAction::BeginOfRunAction(const G4Run* aRun) {

  // run id
  G4int run_id = aRun->GetRunID();
  G4cout << "### Run " << run_id << " start." << G4endl;


  char filename[100];
  std::ostrstream os(filename,100);

  if(autoSeed) {
    // automatic (time-based) random seeds and filenames for each run
    G4cout << "*******************" << G4endl;
    G4cout << "*** AUTOSEED ON ***" << G4endl;
    G4cout << "*******************" << G4endl;
    long seeds[2];
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    // G4cout << "seed1: " << seeds[0] << "; seed2: " << seeds[1] << G4endl;
    HepRandom::setTheSeeds(seeds);
    HepRandom::showEngineStatus();

    // form filename (eg run00_7324329387_3284798343.out)
    os << "run" << std::setw(2) << std::setfill('0') << run_id
       << "_" << seeds[0] << "_" << seeds[1] << std::ends;
  }
  else {
    // default filename, seeds set by /random/reset
    os << "charge" << std::ends;
  }
  //  G4cout << "Filename: " << G4String(filename) << G4endl;


  // send filename to eventAction
  LISAEventAction* eventAction = (LISAEventAction*)
    G4RunManager::GetRunManager()->GetUserEventAction();
  eventAction->SetFilename( G4String(filename)+G4String(".out") );


#ifdef G4ANALYSIS_USE
  // Book histograms and ntuples
  LISAAnalysisManager* analysis = LISAAnalysisManager::getInstance();
  analysis->bookRun( G4String(filename)+G4String(".hbook") );
#endif


  //   if (G4VVisManager::GetConcreteInstance()) {
  //     G4UImanager* UI = G4UImanager::GetUIpointer(); 
  //     UI->ApplyCommand("/vis/scene/notifyHandlers");
  //   } 

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISARunAction::EndOfRunAction(const G4Run* aRun) {

#ifdef G4ANALYSIS_USE
  LISAAnalysisManager* analysis = LISAAnalysisManager::getInstance();
  analysis->FinishRun();
#endif

  //    if (G4VVisManager::GetConcreteInstance()) {
  //      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  //    }
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

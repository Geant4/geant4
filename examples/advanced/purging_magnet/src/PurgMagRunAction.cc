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
// Code developed by:
//  S.Larsson
//
//    ******************************
//    *                            *
//    *    PurgMagRunAction.cc     *
//    *                            *
//    ******************************
//
// $Id: PurgMagRunAction.cc,v 1.2 2004/06/18 09:18:00 gunter Exp $
// GEANT4 tag $Name: geant4-06-02 $
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
    { HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("beginOfRun.rndm");
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
    { HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("endOfRun.rndm");
    }                         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....








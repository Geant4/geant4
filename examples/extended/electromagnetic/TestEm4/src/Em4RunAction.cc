// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em4RunAction.cc,v 1.3 1999-12-15 14:49:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em4RunAction.hh"
#include "Em4RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "Randomize.hh"
#include "CLHEP/Hist/HBookFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em4RunAction::Em4RunAction()
{
  bookHisto();
  
  runMessenger = new Em4RunActionMessenger(this);
  saveFlag = "off";
  saveRndm = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em4RunAction::~Em4RunAction()
{ 
 // Delete HBOOK stuff
  delete [] histo;
  delete hbookManager;
  
  delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em4RunAction::bookHisto()
{
  // init hbook
  hbookManager = new HBookFile("TestEm4.histo", 68);

  // book histograms
  histo[0] = hbookManager->histogram("total energy deposit in C6F6 (MeV)"
                                   , 100,0.,10.);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em4RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  if (saveRndm > 0)
    { HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("beginOfRun.rndm");
    }  

  G4UImanager* UI = G4UImanager::GetUIpointer();
   
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em4RunAction::EndOfRunAction(const G4Run* aRun)
{  
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
    
  // Write histogram file 
  if (saveFlag == "on") hbookManager->write();
  
  // save Rndm status
  if (saveRndm == 1)
    { HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("endOfRun.rndm");
    }      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelRunAction.cc                             *
// * -------                                                            *
// *                                                                    *
// * Version:           0.5                                             *
// * Date:              08/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R. Nartallo
// - First implementation of RunAction
// - Based on Chandra and XMM models
//
// 08.11.2000 R. Nartallo
// - Modified "/vis/****" commands to the G4UIManager in 
//   BeginOfRunAction and EndOfRunAction
//
// **********************************************************************

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "g4std/fstream"
#include "g4std/vector"

#include "XrayTelRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelRunAction::XrayTelRunAction(G4std::vector<G4double*> *enEnergy,
				   G4std::vector<G4ThreeVector*> *enDirect,
				   G4bool* dEvent)
  :EnteringEnergy(enEnergy),
   EnteringDirection(enDirect),drawEvent(dEvent)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelRunAction::~XrayTelRunAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int RunN = aRun->GetRunID();
  if ( RunN % 1000 == 0 ) 
    G4cout << "### Run : " << RunN << G4endl;

  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager* UI = G4UImanager::GetUIpointer(); 
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 

  EnteringEnergy->clear();
  EnteringDirection->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelRunAction::EndOfRunAction(const G4Run* )
{
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/show/view");
  }

  G4std::ofstream outscat("detector.hist", ios::app);

  G4cout << "End of Run summary" << G4endl << G4endl;

  G4double TotEnteringEnergy = 0.0;

  for (G4int i=0;i< EnteringEnergy->size();i++)
    TotEnteringEnergy += *(*EnteringEnergy)[i];
  G4cout << "Total Entering Detector : " << EnteringEnergy->size()  << G4endl;
  G4cout << "Total Entering Detector Energy : " << TotEnteringEnergy  << G4endl;

  for (i=0;i<EnteringEnergy->size();i++) {
    outscat << "  "
	    << *(*EnteringEnergy)[i]
            << "  "
            << (*EnteringDirection)[i]->x()
	    << "  "
            << (*EnteringDirection)[i]->y()
	    << "  "
            << (*EnteringDirection)[i]->z()
	    << G4endl;
  }
  outscat.close();
}




























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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelRunAction.cc                             *
// * -------                                                            *
// *                                                                    *
// * Version:           0.5                                             *
// * Date:              16/10/01                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 16.10.2001 R. Nartallo
// - Updated "/vis" commands to new versions
// - Clean up code to avoid 'pedantic' and 'ANSI' compiler warnings 
//
// 30.11.2000 R. Nartallo
// - Add pre-processor directives to compile without analysis option
//
// 16.11.2000 A. Pfeiffer
// - Implementation of analysis manager call
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope Physics list
// - Based on Chandra and XMM models
// 
//
// **********************************************************************

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "g4std/fstream"
#include "g4std/vector"

#include "XrayTelRunAction.hh"
#include "XrayTelAnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelRunAction::XrayTelRunAction(G4std::vector<G4double*> *enEnergy,
				   G4std::vector<G4ThreeVector*> *enDirect,
				   G4bool* dEvent,
				   XrayTelAnalysisManager* aAnalysisManager)
  :enteringEnergy(enEnergy),
   enteringDirection(enDirect),drawEvent(dEvent),
   fAnalysisManager(aAnalysisManager)
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

  enteringEnergy->clear();
  enteringDirection->clear();

#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->BeginOfRun();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelRunAction::EndOfRunAction(const G4Run* )
{

  size_t i;

  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  G4std::ofstream outscat("detector.hist", ios::app);

  G4cout << "End of Run summary" << G4endl << G4endl;

  G4double totEnteringEnergy = 0.0;

  for (i=0;i< enteringEnergy->size();i++)
    totEnteringEnergy += *(*enteringEnergy)[i];
  G4cout << "Total Entering Detector : " << enteringEnergy->size()  << G4endl;
  G4cout << "Total Entering Detector Energy : " << totEnteringEnergy  << G4endl;

  for (i=0;i<enteringEnergy->size();i++) {
    outscat << "  "
	    << *(*enteringEnergy)[i]
            << "  "
            << (*enteringDirection)[i]->x()
	    << "  "
            << (*enteringDirection)[i]->y()
	    << "  "
            << (*enteringDirection)[i]->z()
	    << G4endl;
  }
  outscat.close();

#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->EndOfRun();
#endif
}




























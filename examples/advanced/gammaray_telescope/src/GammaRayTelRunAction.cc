// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelRunAction.cc,v 1.3 2000-12-06 16:53:14 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelRunAction  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelRunAction.hh"

#include <stdlib.h>
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

extern ofstream outFile;

#ifdef G4ANALYSIS_USE
GammaRayTelRunAction::GammaRayTelRunAction(GammaRayTelAnalysisManager* aMgr)
  :analysisManager(aMgr)
{
}
#else
GammaRayTelRunAction::GammaRayTelRunAction()
{
}
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelRunAction::~GammaRayTelRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelRunAction::BeginOfRunAction(const G4Run* aRun)
{  
  char name[15];

  // Open the file for the tracks of this run
  sprintf(name,"Tracks_%d.dat", aRun->GetRunID());
  outFile.open(name);

  // Prepare the visualization
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 

  // If analysis is used reset the histograms
#ifdef G4ANALYSIS_USE
  analysisManager->BeginOfRun();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelRunAction::EndOfRunAction(const G4Run* aRun)
{
  // Run ended, update the visualization
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  // Close the file with the hits information
  outFile.close();

  // If analysis is used, print out the histograms
#ifdef G4ANALYSIS_USE
  analysisManager->EndOfRun(aRun->GetRunID());
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









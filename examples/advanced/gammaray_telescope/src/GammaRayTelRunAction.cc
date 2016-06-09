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
// $Id: GammaRayTelRunAction.cc,v 1.10 2003/06/16 16:46:31 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelRunAction  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
// 18.11.2001 G.Santin
// - Modified the analysis management according to the new design
//
// ************************************************************


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelRunAction.hh"

#ifdef  G4ANALYSIS_USE
#include "GammaRayTelAnalysis.hh"
#endif

#include <stdlib.h>
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

extern std::ofstream outFile;

GammaRayTelRunAction::GammaRayTelRunAction()
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelRunAction::~GammaRayTelRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelRunAction::BeginOfRunAction(const G4Run* aRun)
{  

  // Open the file for the tracks of this run

  char name[15];
  sprintf(name,"Tracks_%d.dat", aRun->GetRunID());

#ifdef G4STORE_DATA
  outFile.open(name);
#endif

  // Prepare the visualization
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 

  // If analysis is used reset the histograms
#ifdef G4ANALYSIS_USE
  GammaRayTelAnalysis* analysis = GammaRayTelAnalysis::getInstance();
  //  analysis->BeginOfRun(aRun->GetRunID());
  analysis->BeginOfRun();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelRunAction::EndOfRunAction(const G4Run* aRun)
{
  char name[15];
  sprintf(name,"Tracks_%d.dat", aRun->GetRunID());
  G4cout << "End of Run " << G4endl;
  G4cout << "File " << name << G4endl;

  // Run ended, update the visualization
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  // Close the file with the hits information
#ifdef G4STORE_DATA
  outFile.close();
#endif

  // If analysis is used, print out the histograms
#ifdef G4ANALYSIS_USE
  GammaRayTelAnalysis* analysis = GammaRayTelAnalysis::getInstance();
  analysis->EndOfRun(aRun->GetRunID());
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









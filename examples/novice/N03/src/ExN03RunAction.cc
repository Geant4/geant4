// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03RunAction.cc,v 1.1 1999-01-07 16:05:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "ExN03RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExN03RunAction::ExN03RunAction()
{
  runIDcounter = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExN03RunAction::~ExN03RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);
   
  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis~/clear/view");
      UI->ApplyCommand("/vis~/draw/current");
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03RunAction::EndOfRunAction(G4Run* aRun)
{
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

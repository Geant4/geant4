// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02RunAction.cc,v 1.1 1999-01-07 16:05:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "ExN02RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

ExN02RunAction::ExN02RunAction()
{
  runIDcounter = 0;
}

ExN02RunAction::~ExN02RunAction()
{
}

void ExN02RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);
   
  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/tracking/storeTrajectory 1"); 
   
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    UI->ApplyCommand("/vis~/clear/view");
    UI->ApplyCommand("/vis~/draw/current");
  } 
}



void ExN02RunAction::EndOfRunAction(G4Run* aRun)
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  
  if(pVVisManager)
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
  }
}


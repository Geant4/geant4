// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0RunAction.cc,v 1.4 1999-12-15 14:51:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em0RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0RunAction::Em0RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0RunAction::~Em0RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  NbOfTraks0 = 0; NbOfTraks1 = 0; NbOfSteps0 = 0; NbOfSteps1 = 0;   

  G4UImanager* UI = G4UImanager::GetUIpointer();  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if (pVVisManager)
    {
      UI->ApplyCommand("/vis~/clear/view");
      UI->ApplyCommand("/vis~/draw/current");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0RunAction::EndOfRunAction(const G4Run* aRun)
{
  //statistic of nb of tracks and steps
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents)
    { G4double dNbOfEvents = double(NbOfEvents);
    
      G4long oldform = G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
      G4int  oldprec = G4cout.precision(2);
      
      G4cout << "\n nb tracks/event"
                      << "   neutral: " << G4std::setw(7) << NbOfTraks0/dNbOfEvents
                      << "   charged: " << G4std::setw(7) << NbOfTraks1/dNbOfEvents     
             << "\n nb  steps/event"
                      << "   neutral: " << G4std::setw(7) << NbOfSteps0/dNbOfEvents
                      << "   charged: " << G4std::setw(7) << NbOfSteps1/dNbOfEvents
             << G4endl;
             
      G4cout.setf(oldform,G4std::ios::floatfield);
      G4cout.precision(oldprec);       
    }         
                             
  //draw the events
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  
  if (pVVisManager) 
      G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

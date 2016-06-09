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
#include "ExGflashRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"




ExGflashRunAction::ExGflashRunAction(): runID(0)
{
}


ExGflashRunAction::~ExGflashRunAction()
{
}


void ExGflashRunAction::BeginOfRunAction(const G4Run* aRun)
{  
   ((G4Run *)(aRun))->SetRunID(runID++);
   
   std::cout << "### Run " << aRun->GetRunID() << " start." << std::endl;
   
   if (G4VVisManager::GetConcreteInstance())
     {
       G4UImanager* UI = G4UImanager::GetUIpointer();
       UI->ApplyCommand("/vis/scene/notifyHandlers");
     } 
}


void ExGflashRunAction::EndOfRunAction(const G4Run* aRun)
{ 
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
    }

  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;

}












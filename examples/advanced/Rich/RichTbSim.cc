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
// Rich Test Beam Simulation   Main program
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
// ----------------------------------------------------------------
#include <iostream>
#include "RichTbRunAction.hh"
#include "RichTbEventAction.hh"
#include "RichTbDetectorConstruction.hh"
#include "RichTbPrimaryGeneratorAction.hh"
#include "RichTbStackingAction.hh"
#include "RichTbSteppingAction.hh"
#include "RichTbTrackingAction.hh"
#include "RichTbPhysicsList.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "RichTbRunConfig.hh"
#include "RichTbIOData.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif
#include "Randomize.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4ios.hh"
#include <stdlib.h>


int main(int argc,char** argv) {

  // Seed the random number generator manually
  // ------------------------------------------

  G4long myseed = 345354;

   HepRandom::setTheSeed(myseed);

  // Run manager

   G4RunManager * runManager = new G4RunManager;
   //Job and Run  options.
   RichTbRunConfig* rConfig= new RichTbRunConfig();
   // Datafile streams for input and output
   RichTbIOData* rIOData = new  RichTbIOData( rConfig );
   RichTbDetectorConstruction* RichTbDet 
     = new  RichTbDetectorConstruction(rConfig);
   
   runManager->SetUserInitialization(RichTbDet);
   
   RichTbPhysicsList* RichTbPhy
     =new RichTbPhysicsList(rConfig);

   runManager-> SetUserInitialization(RichTbPhy);

   

// UserAction classes - optional

  G4VisManager* visManager = new G4VisExecutive();
  visManager->SetVerboseLevel(0);
  visManager->Initialize();

   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   G4cout<<" PVVisManager "<<pVVisManager<<G4endl;
   G4cout<<"VisManager "<<visManager<<G4endl;

   runManager->SetUserAction(new RichTbRunAction);

   RichTbPrimaryGeneratorAction* PrimaryGenAction =
     new RichTbPrimaryGeneratorAction(RichTbDet);

   runManager->SetUserAction( PrimaryGenAction );

   RichTbEventAction* eventAction=
     new RichTbEventAction(rConfig,visManager,
			   rIOData);

   runManager->SetUserAction(eventAction);

   runManager->SetUserAction(new RichTbStackingAction);

   RichTbSteppingAction* StepAction= new RichTbSteppingAction(rConfig, PrimaryGenAction );
   runManager->SetUserAction(StepAction);
   
   runManager->SetUserAction(new RichTbTrackingAction);

   G4UImanager* UI = G4UImanager::GetUIpointer();  

   G4UIsession* session=0;

   //Initialize G4 kernel
   runManager -> Initialize();

   if(argc==1){

       session = new G4UIterminal;
   }

  if (session){    // Interactive mode  

    UI->ApplyCommand("/run/verbose 0");
    UI->ApplyCommand("/event/verbose 0");
    UI->ApplyCommand("/tracking/verbose 0");
    UI->ApplyCommand("/particle/process/verbose 0");

    session->SessionStart();
    delete session;

  }
  else    // Batch mode
  {
   G4UImanager * UI = G4UImanager::GetUIpointer();
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UI->ApplyCommand(command+fileName);
  }
#ifdef G4VIS_USE
  delete visManager;
#endif

  G4cout << "\nVisManager deleted..\n" <<G4endl;

   delete runManager;

  G4cout << "\nRunManager deleted..\n" <<G4endl;

  return 0;
}







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
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4String.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"

#include "LXeDetectorConstruction.hh"
#include "LXePhysicsList.hh"
#include "LXePrimaryGeneratorAction.hh"
#include "LXeEventAction.hh"
#include "LXeStackingAction.hh"
#include "LXeSteppingAction.hh"
#include "LXeTrackingAction.hh"
#include "LXeRunAction.hh"
#include "LXeSteppingVerbose.hh"

#include "RecorderBase.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
int main(int argc, char** argv)
{
  G4VSteppingVerbose::SetInstance(new LXeSteppingVerbose);

  G4RunManager* runManager = new G4RunManager;

  runManager->SetUserInitialization(new LXeDetectorConstruction);
  runManager->SetUserInitialization(new LXePhysicsList);

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  RecorderBase* recorder = NULL;//No recording is done in this example

  runManager->SetUserAction(new LXePrimaryGeneratorAction);
  runManager->SetUserAction(new LXeStackingAction);
  
  runManager->SetUserAction(new LXeRunAction(recorder));
  runManager->SetUserAction(new LXeEventAction(recorder));
  runManager->SetUserAction(new LXeTrackingAction(recorder));
  runManager->SetUserAction(new LXeSteppingAction(recorder));

  runManager->Initialize();
 
  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  
  if(argc==1){
    G4UIsession* session=0;
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);
#else
    session = new G4UIterminal();
#endif

    //execute vis.mac
    UI->ApplyCommand("/control/execute vis.mac");

    session->SessionStart();
    delete session;

  }
  else{
    G4String command = "/control/execute ";
    G4String filename = argv[1];
    UI->ApplyCommand(command+filename);
  }

  if(recorder)delete recorder;

#ifdef G4VIS_USE
  delete visManager;
#endif

  // job termination
  delete runManager;
  return 0;
}



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
// $Id: G4VisManagerRegisterMessengers.cc,v 1.42 2001-11-14 14:27:32 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4VisManager::RegisterMessengers - John Allison 30/July/1998.

#include "G4VisManager.hh"
#include "G4VVisCommand.hh"
#include "G4VisCommands.hh"
#include "G4VisCommandsCompound.hh"
#include "G4VisCommandsScene.hh"
#include "G4VisCommandsSceneAdd.hh"
#include "G4VisCommandsSceneHandler.hh"
#include "G4VisCommandsViewer.hh"
#include "G4VisCommandsViewerSet.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"

void G4VisManager::RegisterMessengers () {

  // Instantiate individual messengers/commands (often one command per
  // messenger).

  G4VVisCommand::SetVisManager (this);

  G4UIcommand* directory;

  directory = new G4UIdirectory ("/vis/");
  directory -> SetGuidance ("Visualization commands.");
  fDirectoryList.push_back (directory);
  fMessengerList.push_back (new G4VisCommandEnable);
  fMessengerList.push_back (new G4VisCommandVerbose);

  directory = new G4UIdirectory ("/vis/scene/");
  directory -> SetGuidance ("Operations on Geant4 scenes.");
  fDirectoryList.push_back (directory);
  fMessengerList.push_back (new G4VisCommandSceneCreate);
  fMessengerList.push_back (new G4VisCommandSceneEndOfEventAction);
  fMessengerList.push_back (new G4VisCommandSceneList);
  fMessengerList.push_back (new G4VisCommandSceneNotifyHandlers);
  fMessengerList.push_back (new G4VisCommandSceneRemove);
  fMessengerList.push_back (new G4VisCommandSceneSelect);

  directory = new G4UIdirectory ("/vis/scene/add/");
  directory -> SetGuidance ("Add model to current scene.");
  fDirectoryList.push_back (directory);
  fMessengerList.push_back (new G4VisCommandSceneAddAxes);
  fMessengerList.push_back (new G4VisCommandSceneAddGhosts);
  fMessengerList.push_back (new G4VisCommandSceneAddHits);
  fMessengerList.push_back (new G4VisCommandSceneAddLogicalVolume);
  fMessengerList.push_back (new G4VisCommandSceneAddScale);
  fMessengerList.push_back (new G4VisCommandSceneAddText);
  fMessengerList.push_back (new G4VisCommandSceneAddTrajectories);
  fMessengerList.push_back (new G4VisCommandSceneAddVolume);

  directory = new G4UIdirectory ("/vis/sceneHandler/");
  directory -> SetGuidance ("Operations on Geant4 scene handlers.");
  fDirectoryList.push_back (directory);
  fMessengerList.push_back (new G4VisCommandSceneHandlerAttach);
  fMessengerList.push_back (new G4VisCommandSceneHandlerCreate);
  fMessengerList.push_back (new G4VisCommandSceneHandlerList);
  fMessengerList.push_back (new G4VisCommandSceneHandlerRemove);
  fMessengerList.push_back (new G4VisCommandSceneHandlerSelect);

  directory = new G4UIdirectory ("/vis/viewer/");
  directory -> SetGuidance ("Operations on Geant4 viewers.");
  fDirectoryList.push_back (directory);
  fMessengerList.push_back (new G4VisCommandViewerClear);
  fMessengerList.push_back (new G4VisCommandViewerCreate);
  fMessengerList.push_back (new G4VisCommandViewerDolly);
  fMessengerList.push_back (new G4VisCommandViewerFlush);
  fMessengerList.push_back (new G4VisCommandViewerLights);
  // DEPRECATED - moved to /vis/viewer/set/.
  fMessengerList.push_back (new G4VisCommandViewerList);
  fMessengerList.push_back (new G4VisCommandViewerPan);
  fMessengerList.push_back (new G4VisCommandViewerRefresh);
  fMessengerList.push_back (new G4VisCommandViewerRemove);
  fMessengerList.push_back (new G4VisCommandViewerReset);
  fMessengerList.push_back (new G4VisCommandViewerSelect);
  fMessengerList.push_back (new G4VisCommandViewerUpdate);
  fMessengerList.push_back (new G4VisCommandViewerViewpoint);
  // DEPRECATED - moved to /vis/viewer/set/.
  fMessengerList.push_back (new G4VisCommandViewerZoom);

  directory = new G4UIdirectory ("/vis/viewer/set/");
  directory -> SetGuidance ("Set view parameters of current viewer.");
  fDirectoryList.push_back (directory);
  fMessengerList.push_back (new G4VisCommandsViewerSet);

  // Compound commands...
  fMessengerList.push_back (new G4VisCommandDrawTree);
  fMessengerList.push_back (new G4VisCommandDrawVolume);
  fMessengerList.push_back (new G4VisCommandDrawView);
  fMessengerList.push_back (new G4VisCommandOpen);
  fMessengerList.push_back (new G4VisCommandSpecify);

}

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManagerRegisterMessengers.cc,v 1.8 1999-11-05 16:07:32 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4VisManager::RegisterMessengers - John Allison 30/July/1998.

#include "G4VisManager.hh"
#include "G4VVisCommand.hh"
#include "G4VisCommandsScene.hh"
#include "G4VisCommandsSceneAdd.hh"
#include "G4VisCommandsSceneInclude.hh"
#include "G4VisCommandsSceneHandler.hh"
#include "G4VisCommandsViewer.hh"

// OLD STYLE!!
#include "G4VisManMessenger.hh"
#include "G4VisCommandTemplates.hh"
#include "G4VisCommandsCamera.hh"
#include "G4VisCommandsClear.hh"
#include "G4VisCommandsCopy.hh"
#include "G4VisCommandsCreateScene.hh"
#include "G4VisCommandsCreateView.hh"
#include "G4VisCommandsDelete.hh"
#include "G4VisCommandsDraw.hh"
#include "G4VisCommandsLights.hh"
#include "G4VisCommandsPrint.hh"
#include "G4VisCommandsRefresh.hh"
#include "G4VisCommandsSet.hh"
#include "G4VisCommandsShow.hh"

#include "G4ios.hh"
#include <assert.h>

void G4VisManager::RegisterMessengers () {

  /******************************************************************

(Extract for working note.)

Introduction
============

This note defines the concepts of scene, scene handler and viewers.

Scenes
======

We introduce the concept of a graphics-system-independent scene, which
is a list of visualizable objects, such as detector components
(G4PhysicalVolumeModels), hits, trajectories.  The G4VisManager
manages a list of scenes.

* means "not implemented yet".

/vis/scene/create [<scene-name>]
/vis/scene/list [-l] [<scene-name>]
/vis/scene/select [<scene-name>]
* /vis/scene/edit
/vis/scene/remove <scene-name>
/vis/scene/add/ghosts [<particle>]
* /vis/scene/add/ghost [<particle>] [<physical-volume-name>]
*                      [<copy-no>] [<depth>]
/vis/scene/add/logicalVolume <logical-volume-name> [<depth>] 
/vis/scene/add/volume [<physical-volume-name>] [<copy-no>] [<depth>] 
/vis/scene/include/hits [<sensitive-volume-name>] (argument not impl'd yet.)
/vis/scene/include/trajectories [<sensitive-volume-name>] (do.)
* /vis/scene/include/transientObjects
* /vis/scene/set/hitOption accumulate|byEvent
* /vis/scene/set/notifyOption immediate|delayed
* /vis/scene/set/modelingStyle [<modeling-style>]
/vis/scene/notifyHandlers


Scene Handlers
==============

A scene handler is a graphics-system-specific thing that turns a scene
into something that the graphics system can understand.

The commands would look something like:

/vis/sceneHandler/create
       <graphics-system> [<scene-handler-name>] [<scene-name>]
/vis/sceneHandler/attach [<scene-name>]
/vis/sceneHandler/list
/vis/sceneHandler/select [<scene-handler-name>]
/vis/sceneHandler/remove <scene-handler-name>
* /vis/sceneHandler/processScene
* /vis/sceneHandler/notifyEndOfProcessing


Viewers
=======

A viewer opens windows and draws to the screen, hardcopy, etc.  It can
be dumb (a non-interactive window) or intelligent (respond to mouse
clicks, spawn other windows, change viewpoint, etc.).

/vis/viewer/create [<scene-handler>] [<viewer-name>]
/vis/viewer/list
/vis/viewer/select [<viewer-name>]
/vis/viewer/remove <viewer-name>
* /vis/viewer/set/style wireframe|surface|logical
* /vis/viewer/set/viewpoint
* /vis/viewer/set/notifyOption immediate|delayed
* /vis/viewer/notifyHandler
* /vis/viewer/clone
/vis/viewer/update [<viewer-name>]


Global Commands
===============

* /vis/copy/sceneAndView <from-viewer-name> <to-viewer-name>


Compound Commands
=================

* /vis/create/view OGLSXM

would be

/vis/scene/create
/vis/sceneHandler/create $1
/vis/sceneHandler/attach
/vis/viewer/create

and

* /vis/draw Calorimiter

would be

/vis/scene/add/volume $1
/vis/scene/notifyHandlers
/vis/viewer/update

or some such.

  ******************************************************************/

  // Instantiate individual messengers/commands (usually one command
  // per messenger).

  //A The GEANT4 commands available with the GEANT4 Visualization
  //A System.

  // Introduction

  //I The GEANT4 Visualization System consists of a Visualization
  //I Manager and a Visualization Interface.  The problem is to design
  //I a set of commands which offer the GEANT4 user good functionality
  //I and, at the same time, maximally exploit the power of modern
  //I graphics systems.  The chief concept is the scene, a
  //I graphics-system-independent list of visualizable GEANT4 objects,
  //I which can be rendered to any graphics system in a variety of
  //I ways.  If the graphics system has its own graphical database,
  //I this is used.

  G4VVisCommand::SetVisManager (this);

  G4UIcommand* command;
  command = new G4UIdirectory ("/vis/");
  command -> SetGuidance ("Visualization commands.");

  //S \subsection {Scenes}
  //S 
  //S A GEANT4 scene is a set of visualizable objects which can be created,
  //S edited and copied independent of any graphics system.

  command = new G4UIdirectory ("/vis/scene/");
  command -> SetGuidance ("Operations on Geant4 scenes.");
  fMessengerList.append (new G4VisCommandSceneCreate);
  // fMessengerList.append (new G4VisCommandSceneEdit);
  fMessengerList.append (new G4VisCommandSceneList);
  fMessengerList.append (new G4VisCommandSceneNotifyHandlers);
  fMessengerList.append (new G4VisCommandSceneSelect);
  fMessengerList.append (new G4VisCommandSceneRemove);
  command = new G4UIdirectory ("/vis/scene/add/");
  command -> SetGuidance ("Add model to current scene.");
  fMessengerList.append (new G4VisCommandSceneAddVolume);
  fMessengerList.append (new G4VisCommandSceneAddGhosts);
  fMessengerList.append (new G4VisCommandSceneAddLogicalVolume);
  command = new G4UIdirectory ("/vis/scene/include/");
  command -> SetGuidance ("Include drawing option in current scene.");
  fMessengerList.append (new G4VisCommandSceneIncludeHits);
  fMessengerList.append (new G4VisCommandSceneIncludeTrajectories);
  command = new G4UIdirectory ("/vis/sceneHandler/");
  command -> SetGuidance ("Operations on Geant4 scene handlers.");
  fMessengerList.append (new G4VisCommandSceneHandlerAttach);
  fMessengerList.append (new G4VisCommandSceneHandlerCreate);
  fMessengerList.append (new G4VisCommandSceneHandlerList);
  fMessengerList.append (new G4VisCommandSceneHandlerSelect);
  fMessengerList.append (new G4VisCommandSceneHandlerRemove);
  command = new G4UIdirectory ("/vis/viewer/");
  command -> SetGuidance ("Operations on Geant4 viewers.");
  fMessengerList.append (new G4VisCommandViewerCreate);
  fMessengerList.append (new G4VisCommandViewerList);
  fMessengerList.append (new G4VisCommandViewerSelect);
  fMessengerList.append (new G4VisCommandViewerRemove);
  fMessengerList.append (new G4VisCommandViewerUpdate);

  // Camera - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandCamera>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandCameraReset>);

  // Clear - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandClear>);
//  fMessengerList.append
//    (new G4VisSimpleCommandMessenger <G4VisCommandClearScene>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandClearView>);
//  fMessengerList.append
//    (new G4VisSimpleCommandMessenger <G4VisCommandClearViewAndScene>);

  // Copy - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandCopy>);
//  fMessengerList.append
//    (new G4VisSimpleCommandMessenger <G4VisCommandCopyAll>);
//  fMessengerList.append
//    (new G4VisSimpleCommandMessenger <G4VisCommandCopyScene>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandCopyView>);

  // Create Scene - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandCreateScene>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandCreateSceneClear>);

  // Create View - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandCreateView>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandCreateViewNewScene>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandCreateViewNewView>);

  // Delete - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandDelete>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandDeleteScene>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandDeleteView>);

  // Draw - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandDraw>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandDrawCurrent>);

  // Lights - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandLights>);
  fMessengerList.append
    (new G4VisButtonCommandMessenger <G4VisCommandLightsMoveWithCamera>);

  // Print - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandPrint>);

  // Refresh - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandRefresh>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandRefreshView>);

  // Set - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandSet>);
  fMessengerList.append
    (new G4VisButtonCommandMessenger <G4VisCommandSetCulling>);
  fMessengerList.append
    (new G4VisButtonCommandMessenger <G4VisCommandSetCullCoveredDaughters>);
  fMessengerList.append
    (new G4VisButtonCommandMessenger <G4VisCommandSetCullInvisible>);

  // Show - OLD STYLE!!
  fMessengerList.append
    (new G4VisCommandDirectoryMessenger <G4VisCommandShow>);
  fMessengerList.append
    (new G4VisSimpleCommandMessenger <G4VisCommandShowView>);

}

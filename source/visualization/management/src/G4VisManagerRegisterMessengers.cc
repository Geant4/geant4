// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManagerRegisterMessengers.cc,v 1.26 2001-02-06 23:37:03 johna Exp $
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
#include "G4VisCommandsSceneInclude.hh"  // Deprecated
#include "G4VisCommandsSceneHandler.hh"
#include "G4VisCommandsViewer.hh"
#include "G4VisCommandsViewerSet.hh"

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
// END OLD STYLE!!

#include "G4ios.hh"
#include <assert.h>

void G4VisManager::RegisterMessengers () {

  /******************************************************************

(Extract for working note.)

NI means "not implemented yet".


Introduction
============

This note defines the concepts of scene, scene handler and viewer.

A scene handler can handle only one scene at any one time but can have
any number of viewers.

A scene can be attached to any number of scene handlers.

It is important to realise that there is also the concept of a "current"
scene, "current" scene handler and "current" viewer which are maintained
by the G4VisManager.  It is usually the last object created or selected
or operated upon.

Note that there is also a concept of a "standard view" which is that which
comfortably includes all components of the scene.  From this follow the
concept of standard target point, etc.  It is the responsibility of
each viewer to apply its view parameters relatively, taking into account
the scene it represents.  This is a dynamic operation and the scene handler
and its viewers must be smart enough to know when to recalculate
graphics-system-dependent quantities when the scene or the view parameters
change.


Scenes
======

A scene is a list of visualizable objects, such as detector components,
hits, trajectories, axes, etc.

Scenes are graphics-system-independent.

The G4VisManager has a list of scenes.

/vis/scene/create [<scene-name>]
  default       auto-generated name
  This scene becomes current.

/vis/scene/list [<scene-name>] [<verbosity>]
  default:           all             0
  Current scene remains current.

/vis/scene/select [<scene-name>]
  default:      current scene name
  This scene becomes current.

NI /vis/scene/edit

/vis/scene/remove <scene-name>
  Current scene can change or become invalid.

/vis/scene/add/ghosts [<particle>]
  default:                 all
  Adds to current scene.

NI /vis/scene/add/ghost [<particle>] [<physical-volume-name>]
NI                      [<copy-no>] [<depth>]

/vis/scene/add/logicalVolume <logical-volume-name> [<depth>] 
  default:                                             1
  Adds logical volume to current scene.

/vis/scene/add/volume [<physical-volume-name>] [<copy-no>] [<depth>]
  default:                     world                0         -1
  Adds physical volume to current scene.

/vis/scene/add/hits [<sensitive-volume-name>]
  default:              (argument not impl'd yet.)
  Causes hits, if any, to be drawn at the end of processiing an event.
  Adds to current scene.

/vis/scene/add/trajectories [<sensitive-volume-name>]
  default:                      (argument not impl'd yet.)
  Causes trajectories, if any, to be drawn at the end of processiing an event.
  Adds to current scene.

NI /vis/scene/add/transientObjects
  Adds in current scene.

NI /vis/scene/set/hitOption accumulate|byEvent
  Affects current scene.

NI /vis/scene/set/notifyOption immediate|delayed
  Affects current scene.

NI /vis/scene/set/modelingStyle [<modeling-style>]
  Affects current scene.

/vis/scene/notifyHandlers
  Refreshes all viewers of current scene.
  Does not issue "update" (see /vis/viewer/update).

NI /vis/scene/add/axes
  Adds to current scene.

NI /vis/scene/add/text
  Adds to current scene.


Scene Handlers
==============

A scene handler is an object which knows how to interpret a scene for a
specific graphics system.

Each scene handler handles one scene and has, in general, any number of
viewers.

The G4VisManager has a list of scene handlers.

/vis/sceneHandler/create [<graphics-system-name>] [<scene-handler-name>]
  default:                     error               auto-generated name
  (The first default simply triggers a list of possibilities.)
  This scene handler becomes current.
  The current scene, if any, is attached.

/vis/sceneHandler/attach [<scene-name>]
  default:             current scene name
  Attaches scene to current scene handler.

/vis/sceneHandler/list [<scene-handler-name>] [<verbosity>]
  default:                     all                  0
  Current scene handler remains current.

/vis/sceneHandler/select [<scene-handler-name>]
  default:            current scene handler name
  This scene handler becomes current.

/vis/sceneHandler/remove <scene-handler-name>
  Current scene handler can change or become invalid.

NI /vis/sceneHandler/processScene
  Refreshes all viewers of current scene handler.
  Does not issue "update" (see /vis/viewer/update).

NI /vis/sceneHandler/notifyEndOfProcessing
  Issues "update" for each viewer of current scene handler.


Viewers
=======

A viewer opens windows and draws to the screen or writes to file for
off-line viewing or hardcopy, etc.  It can be dumb (a non-interactive
window) or intelligent (respond to mouse clicks, spawn other windows,
change viewpoint, etc.).

Most viewer commands respond to the viewer "short name", which is the
name up to the first space character, if any.  Thus, a viewer name can
contain spaces but must be unique up to the first space.

/vis/viewer/clear
  Affects current viewer.

NI /vis/viewer/clone
  Creates a clone.  Clone becomes current viewer.

/vis/viewer/create [<scene-handler>]           [<viewer-name>]     [<pixels>]
  default:     current scene handler name    auto-generated name      600
  Pixel size of square window (hint only).
  This viewer becomes current.

/vis/viewer/dolly [<increment>] [<unit>]
  default:        current-value    m
  Moves the camera incrementally in by this distance.
  Affects current viewer.

/vis/viewer/dollyTo [<distance>] [<unit>]
  default:         current-value    m
  Moves the camera in this distance relative to standard target point.
  Affects current viewer.

/vis/viewer/lightsThetaPhi [<theta>] [<phi>] [deg|rad]
  default: current as default.
/vis/viewer/lightsVector [<x>] [<y>] [<z>]
  default: current as default.
  Affects current viewer.

/vis/viewer/list   [<viewer-name>]    [<verbosity>]
  default:             all                  0
  Current viewer remains current.

/vis/viewer/pan [<right-increment>] [<up-increment>] [<unit>]
  default:                 0           0
  Moves the camera incrementally right and up by these amounts.
  Affects current viewer.

/vis/viewer/panTo [<right>] [<up>] [<unit>]
  default:                0     0
  Moves the camera to this position right and up relative to standard target
    point.
  Affects current viewer.

/vis/viewer/select <viewer-name>
  default:          no default
  This viewer becomes current.

/vis/viewer/refresh [<viewer-name>]
  default:       current viewer name
  Re-traverses graphical data, which is enough in some cases to refresh the
  view.  In some cases post-processing is required to complete the view -
  see /vis/viewer/update.  This viewer becomes current.

/vis/viewer/remove <viewer-name>
  default:          no default
  Current viewer can change or become invalid.

/vis/viewer/reset [<viewer-name>]
  default:      current viewer name
  Resets view parameters to defaults.
  This viewer becomes current.

/vis/viewer/viewpointThetaPhi [<theta>] [<phi>] [deg|rad]
  default: current as default.
/vis/viewer/viewpointVector [<x>] [<y>] [<z>]
  default: current as default.
  Set direction from target to camera.  Also changes lightpoint direction if
  lights are set to move with camera.
  Affects current viewer.

/vis/viewer/set/zoom [<factor>]
  default:                 1
  Multiplies magnification by this factor.
  Affects current viewer.

/vis/viewer/set/zoomTo [<factor>]
  default:                   1
  Magnifies by this factor relative to standard view.
  Affects current viewer.

(Now follow some viewer/set commands.  A possible syntax is:

NI /vis/viewer/set <parameter-name> <parameter-value> [<parameter-value>]...

but this would mean re-inventing parameter parsing which is already
contained in many G4UIcommand subclasses.  To re-use this code we have
to know the type of the parameter, so we need separate commands.)

/vis/viewer/set/all <from-viewer-name>
  Copies view parameters from from-viewer to current viewer.
  Affects current viewer.

/vis/viewer/set/autoRefresh [true|false]
  default:                    false
  View is automatically refreshed after a change of view parameters.
  Affects current viewer.

/vis/viewer/set/culling
  g[lobal]|c[overedDaughters]|i[nvisible]|d[ensity] [true|false]
  [density] [unit]
  default: none true 0.01 g/cm3
  Affects current viewer.

NI /vis/viewer/set/cutawayPlane ...
  Set plane(s) for cutaway views.
  Affects current viewer.

/vis/viewer/set/edge [true|false]
  default:              true
  Affects current viewer.

/vis/viewer/set/hiddenEdge  [true|false]
  default:                       true
  Affects current viewer.

/vis/viewer/set/hiddenMarker  [true|false]
  default:                       true
  Affects current viewer.

/vis/viewer/set/lightsMove with-camera|with-object
  Affects current viewer.

/vis/viewer/set/lineSegmentsPerCircle        [<number-of-sides-per-circle>]
  default:                                         24
  Number of sides per circle in polygon/polyhedron graphical representation
    of objects with curved lines/surfaces.
  Affects current viewer.

/vis/viewer/set/projection
  o[rthogonal]|p[erspective] [<field-half-angle>] [deg|rad]
  default: none                       30             deg
  Affects current viewer.

#specialNI /vis/viewer/set/orbit (parameters)
  Orbits the scene about the up-vector, lights fixed to the scene.  Draws N
    frames, the camera rotated Delta-beta about the up-vector each frame.
  Affects current viewer.

/vis/viewer/set/sectionPlane ...
  Set plane for drawing section (DCUT).  Specify plane by x y z units nx ny nz,
  e.g., for a y-z plane at x = 1 cm:
  /vis/viewer/set/section_plane on 1 0 0 cm 1 0 0
  Affects current viewer.

#specialNI /vis/viewer/set/spin (parameters)
#orNI      /vis/viewer/advanced?
  Spins the scene about the up-vector, lights fixed to the camera.  Draws N
  frames, the scene rotated Delta-beta about the up-vector each frame.
  Affects current viewer.

/vis/viewer/set/style w[ireframe]|s[urface]
  Affects current viewer.

NI /vis/viewer/notifyOption immediate|delayed ?Issue of "update" after "set"?

NI /vis/viewer/notifyHandler ??

/vis/viewer/update [<viewer-name>]
  default:     current viewer name
  Initiates post-processing if required.  This viewer becomes current.


Attributes (nothing implemented yet)
==========

The G4VisManager also keeps a list of visualization attributes which can
be created and changed and attributed to visualizable objects.

NI /vis/attributes/create [<vis-attributes-name>]
  default:                 auto-generated name
  These attributes become current.

NI /vis/attributes/set/colour [<vis-attributes-name>] [<r>]  [<g>] [<b>] [<o>]
NI /vis/attributes/set/color  [<vis-attributes-name>] [<r>]  [<g>] [<b>] [<o>]
  default:                   current attributes name   1      1     1     1
  Sets colour (red, green, blue, opacity).
  Affects current vis attributes.

NI /vis/scene/set/attributes <logical-volume-name>
  Associates current vis attributes with logical volume.  (Do we need to
    provide possibility of resetting to original attributes?)
  (Move to scene when implemented.)


General Commands
================

/vis/enable [true|false]
  default:       true
/vis/disable
  Enables/disables visualization system.

/vis/verbose [<verbosity-integer>]
  default:               0


Compound Commands
=================

/vis/open [<graphics-system-name>] [<[pixels>]
default:          error                600
  /vis/sceneHandler/create $1
  /vis/viewer/create ! ! $2

NI /vis/draw <physical-volume-name> clashes with old /vis~/draw/, so...
/vis/drawVolume [<physical-volume-name>]
default:             world
  /vis/scene/create
  /vis/scene/add/volume $1
  /vis/sceneHandler/attach
  /vis/viewer/refresh

/vis/drawView [<theta-deg>] [<phi-deg>]
              [<pan-right>] [<pan-up>] [<pan-unit>]
              [<zoom-factor>]
              [<dolly>] [<dolly-unit>]
default: 0 0 0 0 cm 1 0 cm
  /vis/viewer/viewpointThetaPhi $1 $2 deg
  /vis/viewer/panTo $3 $4 $5
  /vis/viewer/zoomTo $6
  /vis/viewer/dollyTo $7 $8

/vis/specify <logical-volume-name>
  /geometry/print $1
  /vis/scene/create
  /vis/scene/add/logicalVolume $1
  /vis/sceneHandler/attach
  /vis/viewer/refresh

  ******************************************************************/

  // Instantiate individual messengers/commands (usually one command
  // per messenger).

  G4VVisCommand::SetVisManager (this);

  G4UIcommand* command;
  command = new G4UIdirectory ("/vis/");
  command -> SetGuidance ("Visualization commands.");
  fMessengerList.append (new G4VisCommandEnable);
  fMessengerList.append (new G4VisCommandVerbose);

  command = new G4UIdirectory ("/vis/scene/");
  command -> SetGuidance ("Operations on Geant4 scenes.");
  fMessengerList.append (new G4VisCommandSceneCreate);
  fMessengerList.append (new G4VisCommandSceneList);
  fMessengerList.append (new G4VisCommandSceneNotifyHandlers);
  fMessengerList.append (new G4VisCommandSceneSelect);
  fMessengerList.append (new G4VisCommandSceneRemove);

  command = new G4UIdirectory ("/vis/scene/add/");
  command -> SetGuidance ("Add model to current scene.");
  fMessengerList.append (new G4VisCommandSceneAddGhosts);
  fMessengerList.append (new G4VisCommandSceneAddHits);
  fMessengerList.append (new G4VisCommandSceneAddLogicalVolume);
  fMessengerList.append (new G4VisCommandSceneAddTrajectories);
  fMessengerList.append (new G4VisCommandSceneAddVolume);

  command = new G4UIdirectory ("/vis/scene/include/");
  command -> SetGuidance ("Deprecated commands; now in /vis/scene/add/.");
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
  fMessengerList.append (new G4VisCommandViewerClear);
  fMessengerList.append (new G4VisCommandViewerCreate);
  fMessengerList.append (new G4VisCommandViewerDolly);
  fMessengerList.append (new G4VisCommandViewerLights);
  fMessengerList.append (new G4VisCommandViewerList);
  fMessengerList.append (new G4VisCommandViewerPan);
  fMessengerList.append (new G4VisCommandViewerRefresh);
  fMessengerList.append (new G4VisCommandViewerRemove);
  fMessengerList.append (new G4VisCommandViewerReset);
  fMessengerList.append (new G4VisCommandViewerSelect);
  fMessengerList.append (new G4VisCommandViewerUpdate);
  fMessengerList.append (new G4VisCommandViewerViewpoint);
  fMessengerList.append (new G4VisCommandViewerZoom);

  command = new G4UIdirectory ("/vis/viewer/set/");
  command -> SetGuidance ("Set view parameters of current viewer.");
  fMessengerList.append (new G4VisCommandsViewerSet);

  // Compound commands...
  fMessengerList.append (new G4VisCommandDrawVolume);
  fMessengerList.append (new G4VisCommandDrawView);
  fMessengerList.append (new G4VisCommandOpen);
  fMessengerList.append (new G4VisCommandSpecify);

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

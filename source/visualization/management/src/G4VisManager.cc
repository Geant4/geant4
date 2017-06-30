//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4VisManager.cc 103926 2017-05-03 13:43:27Z gcosmo $
//
// 
// GEANT4 Visualization Manager - John Allison 02/Jan/1996.

#include "G4VisManager.hh"

#include "G4VisCommands.hh"
#include "G4VisCommandsCompound.hh"
#include "G4VisCommandsGeometry.hh"
#include "G4VisCommandsGeometrySet.hh"
#include "G4VisCommandsMultithreading.hh"
#include "G4VisCommandsSet.hh"
#include "G4VisCommandsScene.hh"
#include "G4VisCommandsSceneAdd.hh"
#include "G4VisCommandsSceneHandler.hh"
#include "G4VisCommandsTouchable.hh"
#include "G4VisCommandsTouchableSet.hh"
#include "G4VisCommandsViewer.hh"
#include "G4VisCommandsViewerDefault.hh"
#include "G4VisCommandsViewerSet.hh"
#include "G4UImanager.hh"
#include "G4VisStateDependent.hh"
#include "G4UIdirectory.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VSceneHandler.hh"
#include "G4VViewer.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Vector3D.hh"
#include "G4Point3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Polyline.hh"
#include "G4Polyhedron.hh"
#include "G4NullModel.hh"
#include "G4ModelingParameters.hh"
#include "G4TransportationManager.hh"
#include "G4VisCommandModelCreate.hh"
#include "G4VisCommandsListManager.hh"
#include "G4VisModelManager.hh"
#include "G4VModelFactory.hh"
#include "G4VisFilterManager.hh"
#include "G4VTrajectoryModel.hh"
#include "G4TrajectoryDrawByCharge.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include <map>
#include <set>
#include <vector>
#include <sstream>

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4GeometryWorkspacePool.hh"
#include "G4SolidsWorkspacePool.hh"
#include <deque>
#include <typeinfo>
#ifdef G4VIS_USE_STD11
#include <chrono>
#include <thread>
#endif
#endif

G4VisManager* G4VisManager::fpInstance = 0;

G4VisManager::Verbosity G4VisManager::fVerbosity = G4VisManager::warnings;

G4VisManager::G4VisManager (const G4String& verbosityString):
  fVerbose         (1),
  fInitialised     (false),
  fpGraphicsSystem (0),
  fpScene          (0),
  fpSceneHandler   (0),
  fpViewer         (0),
  fpStateDependent (0),
  fEventRefreshing          (false),
  fTransientsDrawnThisRun   (false),
  fTransientsDrawnThisEvent (false),
  fNoOfEventsDrawnThisRun   (0),
  fNKeepRequests            (0),
  fEventKeepingSuspended    (false),
  fKeptLastEvent            (false),
  fpRequestedEvent          (0),
  fAbortReviewKeptEvents    (false),
  fIsDrawGroup              (false),
  fDrawGroupNestingDepth    (0),
  fIgnoreStateChanges       (false)
#ifdef G4MULTITHREADED
, fMaxEventQueueSize        (100)
, fWaitOnEventQueueFull     (true)
#endif
  // All other objects use default constructors.
{
  fpTrajDrawModelMgr = new G4VisModelManager<G4VTrajectoryModel>("/vis/modeling/trajectories");
  fpTrajFilterMgr = new G4VisFilterManager<G4VTrajectory>("/vis/filtering/trajectories");
  fpHitFilterMgr = new G4VisFilterManager<G4VHit>("/vis/filtering/hits");
  fpDigiFilterMgr = new G4VisFilterManager<G4VDigi>("/vis/filtering/digi");

  VerbosityGuidanceStrings.push_back
    ("Simple graded message scheme - digit or string (1st character defines):");
  VerbosityGuidanceStrings.push_back
    ("  0) quiet,         // Nothing is printed.");
  VerbosityGuidanceStrings.push_back
    ("  1) startup,       // Startup and endup messages are printed...");
  VerbosityGuidanceStrings.push_back
    ("  2) errors,        // ...and errors...");
  VerbosityGuidanceStrings.push_back
    ("  3) warnings,      // ...and warnings...");
  VerbosityGuidanceStrings.push_back
    ("  4) confirmations, // ...and confirming messages...");
  VerbosityGuidanceStrings.push_back
    ("  5) parameters,    // ...and parameters of scenes and views...");
  VerbosityGuidanceStrings.push_back
    ("  6) all            // ...and everything available.");

  if (fpInstance) {
    G4Exception
      ("G4VisManager::G4VisManager",
       "visman0001", FatalException,
       "Attempt to Construct more than one VisManager");
  }

  fpInstance = this;
  SetConcreteInstance(this);

  fpStateDependent = new G4VisStateDependent (this);
  // No need to delete this; G4StateManager does this.

  fVerbosity = GetVerbosityValue(verbosityString);
  if (fVerbosity >= startup) {
      G4cout
	<< "Visualization Manager instantiating with verbosity \""
	<< VerbosityString(fVerbosity)
	<< "\"..." << G4endl;
  }

  // Note: The specific graphics systems must be instantiated in a
  // higher level library to avoid circular dependencies.  Also,
  // some specifically need additional external libararies that the
  // user must supply.  Therefore we ask the user to implement
  // RegisterGraphicsSystems() and RegisterModelFactories()
  // in a subclass.  We have to wait for the subclass to instantiate
  // so RegisterGraphicsSystems() cannot be called from this
  // constructor; it is called from Initialise().  So we ask the
  // user:
  //   (a) to write a subclass and implement  RegisterGraphicsSystems()
  //       and RegisterModelFactories().  See
  //       visualization/include/G4VisExecutive.hh/icc as an example.
  //   (b) instantiate the subclass.
  //   (c) invoke the Initialise() method of the subclass.
  // For example:
  //   ...
  // #ifdef G4VIS_USE
  //   // Instantiate and initialise Visualization Manager.
  //   G4VisManager* visManager = new G4VisExecutive;
  //   visManager -> SetVerboseLevel (Verbose);
  //   visManager -> Initialise ();
  // #endif
  //   // (Don't forget to delete visManager;)
  //   ...

  // Make top level command directory...
  // Vis commands should *not* be broadcast to threads (2nd argument).
  G4UIcommand* directory = new G4UIdirectory ("/vis/",false);
  directory -> SetGuidance ("Visualization commands.");
  fDirectoryList.push_back (directory);

  // Instantiate *basic* top level commands so that they can be used
  // immediately after instantiation of the vis manager.  Other top
  // level and lower level commands are instantiated later in
  // RegisterMessengers.
  G4VVisCommand::SetVisManager (this);  // Sets shared pointer
  RegisterMessenger(new G4VisCommandVerbose);
  RegisterMessenger(new G4VisCommandInitialize);
}

G4VisManager::~G4VisManager () {
  fpInstance = 0;
  size_t i;
  for (i = 0; i < fSceneList.size (); ++i) {
    delete fSceneList[i];
  }
  for (i = 0; i < fAvailableSceneHandlers.size (); ++i) {
  if (fAvailableSceneHandlers[i] != NULL) {
    delete fAvailableSceneHandlers[i];
  }
  }
  for (i = 0; i < fAvailableGraphicsSystems.size (); ++i) {
    if (fAvailableGraphicsSystems[i]) {
      delete fAvailableGraphicsSystems[i];
    }
  }
  if (fVerbosity >= startup) {
    G4cout << "Graphics systems deleted." << G4endl;
    G4cout << "Visualization Manager deleting..." << G4endl;
  }
  for (i = 0; i < fMessengerList.size (); ++i) {
    delete fMessengerList[i];
  }
  for (i = 0; i < fDirectoryList.size (); ++i) {
    delete fDirectoryList[i];
  }

  delete fpDigiFilterMgr;
  delete fpHitFilterMgr;
  delete fpTrajFilterMgr;
  delete fpTrajDrawModelMgr;
}

G4VisManager* G4VisManager::GetInstance () {
  if (!fpInstance) {
    G4Exception
      ("G4VisManager::GetInstance",
       "visman0002", FatalException, "VisManager not yet instantiated");
  }
  return fpInstance;
}

void G4VisManager::Initialise () {

  if (fInitialised && fVerbosity >= warnings) {
    G4cout << "WARNING: G4VisManager::Initialise: already initialised."
	   << G4endl;
    return;
  }

  if (fVerbosity >= startup) {
    G4cout << "Visualization Manager initialising..." << G4endl;
  }

  if (fVerbosity >= parameters) {
    G4cout <<
      "\nYou have instantiated your own Visualization Manager, inheriting"
      "\n  G4VisManager and implementing RegisterGraphicsSystems(), in which"
      "\n  you should, normally, instantiate drivers which do not need"
      "\n  external packages or libraries, and, optionally, drivers under"
      "\n  control of environment variables."
      "\n  Also you should implement RegisterModelFactories()."
      "\n  See visualization/management/include/G4VisExecutive.hh/icc, for example."
      "\n  In your main() you will have something like:"
      "\n  #ifdef G4VIS_USE"
      "\n    G4VisManager* visManager = new G4VisExecutive;"
      "\n    visManager -> SetVerboseLevel (Verbose);"
      "\n    visManager -> Initialize ();"
      "\n  #endif"
      "\n  (Don't forget to delete visManager;)"
      "\n"
	 << G4endl;
  }

  if (fVerbosity >= startup) {
    G4cout << "Registering graphics systems..." << G4endl;
  }

  RegisterGraphicsSystems ();

  if (fVerbosity >= startup) {
    G4cout <<
      "\nYou have successfully registered the following graphics systems."
	 << G4endl;
    PrintAvailableGraphicsSystems (fVerbosity);
    G4cout << G4endl;
  }

  // Make command directories for commands instantiated in the
  // modeling subcategory...
  G4UIcommand* directory;
  directory = new G4UIdirectory ("/vis/modeling/");
  directory -> SetGuidance ("Modeling commands.");
  fDirectoryList.push_back (directory);
  directory = new G4UIdirectory ("/vis/modeling/trajectories/");
  directory -> SetGuidance ("Trajectory model commands.");
  fDirectoryList.push_back (directory);
  directory = new G4UIdirectory ("/vis/modeling/trajectories/create/");
  directory -> SetGuidance ("Create trajectory models and messengers.");
  fDirectoryList.push_back (directory);

  // Filtering command directory
  directory = new G4UIdirectory ("/vis/filtering/");
  directory -> SetGuidance ("Filtering commands.");
  fDirectoryList.push_back (directory);
  directory = new G4UIdirectory ("/vis/filtering/trajectories/");
  directory -> SetGuidance ("Trajectory filtering commands.");
  fDirectoryList.push_back (directory);
  directory = new G4UIdirectory ("/vis/filtering/trajectories/create/");
  directory -> SetGuidance ("Create trajectory filters and messengers.");
  fDirectoryList.push_back (directory);
  directory = new G4UIdirectory ("/vis/filtering/hits/");
  directory -> SetGuidance ("Hit filtering commands.");
  fDirectoryList.push_back (directory);
  directory = new G4UIdirectory ("/vis/filtering/hits/create/");
  directory -> SetGuidance ("Create hit filters and messengers.");
  fDirectoryList.push_back (directory);
  directory = new G4UIdirectory ("/vis/filtering/digi/");
  directory -> SetGuidance ("Digi filtering commands.");
  fDirectoryList.push_back (directory);
  directory = new G4UIdirectory ("/vis/filtering/digi/create/");
  directory -> SetGuidance ("Create digi filters and messengers.");
  fDirectoryList.push_back (directory);

  RegisterMessengers ();

  if (fVerbosity >= startup) {
    G4cout << "Registering model factories..." << G4endl;
  }

  RegisterModelFactories();

  if (fVerbosity >= startup) {
    G4cout <<
      "\nYou have successfully registered the following model factories."
	   << G4endl;
    PrintAvailableModels (fVerbosity);
    G4cout << G4endl;
  }

  if (fVerbosity >= startup) {
    PrintAvailableUserVisActions (fVerbosity);
    G4cout << G4endl;
  }

  if (fVerbosity >= startup) {
    PrintAvailableColours (fVerbosity);
    G4cout << G4endl;
  }

  fInitialised = true;
}

void G4VisManager::RegisterMessengers () {
  
  // Instantiate individual messengers/commands (often - but not
  // always - one command per messenger).
  
  G4UIcommand* directory;
  
  // *Basic* top level commands were instantiated in the constructor
  // so that they can be used immediately after instantiation of the
  // vis manager.  Other top level and lower level commands are
  // instantiated here.
  
  // Other top level commands...
  RegisterMessenger(new G4VisCommandAbortReviewKeptEvents);
  RegisterMessenger(new G4VisCommandEnable);
  RegisterMessenger(new G4VisCommandList);
  RegisterMessenger(new G4VisCommandReviewKeptEvents);
  
  // Compound commands...
  RegisterMessenger(new G4VisCommandDrawTree);
  RegisterMessenger(new G4VisCommandDrawView);
  RegisterMessenger(new G4VisCommandDrawVolume);
  RegisterMessenger(new G4VisCommandOpen);
  RegisterMessenger(new G4VisCommandSpecify);
  
  directory = new G4UIdirectory ("/vis/geometry/");
  directory -> SetGuidance("Operations on vis attributes of Geant4 geometry.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandGeometryList);
  RegisterMessenger(new G4VisCommandGeometryRestore);
  
  directory = new G4UIdirectory ("/vis/geometry/set/");
  directory -> SetGuidance("Set vis attributes of Geant4 geometry.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandGeometrySetColour);
  RegisterMessenger(new G4VisCommandGeometrySetDaughtersInvisible);
  RegisterMessenger(new G4VisCommandGeometrySetLineStyle);
  RegisterMessenger(new G4VisCommandGeometrySetLineWidth);
  RegisterMessenger(new G4VisCommandGeometrySetForceAuxEdgeVisible);
  RegisterMessenger(new G4VisCommandGeometrySetForceLineSegmentsPerCircle);
  RegisterMessenger(new G4VisCommandGeometrySetForceSolid);
  RegisterMessenger(new G4VisCommandGeometrySetForceWireframe);
  RegisterMessenger(new G4VisCommandGeometrySetVisibility);
  
#ifdef G4MULTITHREADED
  directory = new G4UIdirectory ("/vis/multithreading/");
  directory -> SetGuidance("Commands unique to multithreading mode.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandMultithreadingActionOnEventQueueFull);
  RegisterMessenger(new G4VisCommandMultithreadingMaxEventQueueSize);
#endif

  directory = new G4UIdirectory ("/vis/set/");
  directory -> SetGuidance
    ("Set quantities for use in future commands where appropriate.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandSetArrow3DLineSegmentsPerCircle);
  RegisterMessenger(new G4VisCommandSetColour);
  RegisterMessenger(new G4VisCommandSetLineWidth);
  RegisterMessenger(new G4VisCommandSetTextColour);
  RegisterMessenger(new G4VisCommandSetTextLayout);
  RegisterMessenger(new G4VisCommandSetTextSize);
  RegisterMessenger(new G4VisCommandSetTouchable);
  
  directory = new G4UIdirectory ("/vis/scene/");
  directory -> SetGuidance ("Operations on Geant4 scenes.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandSceneActivateModel);
  RegisterMessenger(new G4VisCommandSceneCreate);
  RegisterMessenger(new G4VisCommandSceneEndOfEventAction);
  RegisterMessenger(new G4VisCommandSceneEndOfRunAction);
  RegisterMessenger(new G4VisCommandSceneList);
  RegisterMessenger(new G4VisCommandSceneNotifyHandlers);
  RegisterMessenger(new G4VisCommandSceneSelect);
  
  directory = new G4UIdirectory ("/vis/scene/add/");
  directory -> SetGuidance ("Add model to current scene.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandSceneAddArrow);
  RegisterMessenger(new G4VisCommandSceneAddArrow2D);
  RegisterMessenger(new G4VisCommandSceneAddAxes);
  RegisterMessenger(new G4VisCommandSceneAddDate);
  RegisterMessenger(new G4VisCommandSceneAddDigis);
  RegisterMessenger(new G4VisCommandSceneAddEventID);
  RegisterMessenger(new G4VisCommandSceneAddExtent);
  RegisterMessenger(new G4VisCommandSceneAddFrame);
  RegisterMessenger(new G4VisCommandSceneAddGPS);
  RegisterMessenger(new G4VisCommandSceneAddHits);
  RegisterMessenger(new G4VisCommandSceneAddLine);
  RegisterMessenger(new G4VisCommandSceneAddLine2D);
  RegisterMessenger(new G4VisCommandSceneAddLogicalVolume);
  RegisterMessenger(new G4VisCommandSceneAddLogo);
  RegisterMessenger(new G4VisCommandSceneAddLogo2D);
  RegisterMessenger(new G4VisCommandSceneAddMagneticField);
  RegisterMessenger(new G4VisCommandSceneAddPSHits);
  RegisterMessenger(new G4VisCommandSceneAddScale);
  RegisterMessenger(new G4VisCommandSceneAddText);
  RegisterMessenger(new G4VisCommandSceneAddText2D);
  RegisterMessenger(new G4VisCommandSceneAddTrajectories);
  RegisterMessenger(new G4VisCommandSceneAddUserAction);
  RegisterMessenger(new G4VisCommandSceneAddVolume);
  
  directory = new G4UIdirectory ("/vis/sceneHandler/");
  directory -> SetGuidance ("Operations on Geant4 scene handlers.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandSceneHandlerAttach);
  RegisterMessenger(new G4VisCommandSceneHandlerCreate);
  RegisterMessenger(new G4VisCommandSceneHandlerList);
  RegisterMessenger(new G4VisCommandSceneHandlerSelect);
  
  directory = new G4UIdirectory ("/vis/touchable/");
  directory -> SetGuidance ("Operations on touchables.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandsTouchable);
  
  directory = new G4UIdirectory ("/vis/touchable/set/");
  directory -> SetGuidance ("Set vis attributes of current touchable.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandsTouchableSet);

  directory = new G4UIdirectory ("/vis/viewer/");
  directory -> SetGuidance ("Operations on Geant4 viewers.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandViewerAddCutawayPlane);
  RegisterMessenger(new G4VisCommandViewerChangeCutawayPlane);
  RegisterMessenger(new G4VisCommandViewerClear);
  RegisterMessenger(new G4VisCommandViewerClearCutawayPlanes);
  RegisterMessenger(new G4VisCommandViewerClearTransients);
  RegisterMessenger(new G4VisCommandViewerClearVisAttributesModifiers);
  RegisterMessenger(new G4VisCommandViewerClone);
  RegisterMessenger(new G4VisCommandViewerCopyViewFrom);
  RegisterMessenger(new G4VisCommandViewerCreate);
  RegisterMessenger(new G4VisCommandViewerDolly);
  RegisterMessenger(new G4VisCommandViewerFlush);
  RegisterMessenger(new G4VisCommandViewerInterpolate);
  RegisterMessenger(new G4VisCommandViewerList);
  RegisterMessenger(new G4VisCommandViewerPan);
  RegisterMessenger(new G4VisCommandViewerRebuild);
  RegisterMessenger(new G4VisCommandViewerRefresh);
  RegisterMessenger(new G4VisCommandViewerReset);
  RegisterMessenger(new G4VisCommandViewerSave);
  RegisterMessenger(new G4VisCommandViewerScale);
  RegisterMessenger(new G4VisCommandViewerSelect);
  RegisterMessenger(new G4VisCommandViewerUpdate);
  RegisterMessenger(new G4VisCommandViewerZoom);
  
  directory = new G4UIdirectory ("/vis/viewer/default/");
  directory -> SetGuidance("Set default values for future viewers.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandViewerDefaultHiddenEdge);
  RegisterMessenger(new G4VisCommandViewerDefaultStyle);
  
  directory = new G4UIdirectory ("/vis/viewer/set/");
  directory -> SetGuidance ("Set view parameters of current viewer.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandsViewerSet);
  
  // List manager commands
  RegisterMessenger(new G4VisCommandListManagerList< G4VisModelManager<G4VTrajectoryModel> >
		    (fpTrajDrawModelMgr, fpTrajDrawModelMgr->Placement()));
  RegisterMessenger(new G4VisCommandListManagerSelect< G4VisModelManager<G4VTrajectoryModel> >
		    (fpTrajDrawModelMgr, fpTrajDrawModelMgr->Placement()));
  
  // Trajectory filter manager commands
  RegisterMessenger(new G4VisCommandListManagerList< G4VisFilterManager<G4VTrajectory> >
                    (fpTrajFilterMgr, fpTrajFilterMgr->Placement()));
  RegisterMessenger(new G4VisCommandManagerMode< G4VisFilterManager<G4VTrajectory> >
                    (fpTrajFilterMgr, fpTrajFilterMgr->Placement()));
  
  // Hit filter manager commands
  RegisterMessenger(new G4VisCommandListManagerList< G4VisFilterManager<G4VHit> >
                    (fpHitFilterMgr, fpHitFilterMgr->Placement()));
  RegisterMessenger(new G4VisCommandManagerMode< G4VisFilterManager<G4VHit> >
                    (fpHitFilterMgr, fpHitFilterMgr->Placement()));
  
  // Digi filter manager commands
  RegisterMessenger(new G4VisCommandListManagerList< G4VisFilterManager<G4VDigi> >
                    (fpDigiFilterMgr, fpDigiFilterMgr->Placement()));
  RegisterMessenger(new G4VisCommandManagerMode< G4VisFilterManager<G4VDigi> >
                    (fpDigiFilterMgr, fpDigiFilterMgr->Placement()));
}

void G4VisManager::Enable() {
  if (IsValidView ()) {
    SetConcreteInstance(this);
    if (fVerbosity >= confirmations) {
      G4cout << "G4VisManager::Enable: visualization enabled." << G4endl;
    }
    if (fVerbosity >= warnings) {
      auto nKeptEvents =
      G4RunManager::GetRunManager()->GetCurrentRun()->GetEventVector()->size();
      G4cout <<
  "There are " << nKeptEvents << " kept events."
  "\n  \"/vis/reviewKeptEvents\" to review them one by one."
  "\n  \"/vis/viewer/flush\" or \"/vis/viewer/rebuild\" to see them accumulated."
      << G4endl;
    }
  }
  else {
    if (fVerbosity >= warnings) {
      G4cout <<
	"G4VisManager::Enable: WARNING: visualization remains disabled for"
	"\n  above reasons.  Rectifying with valid vis commands will"
	"\n  automatically enable."
	     << G4endl;
    }
  }
}

void G4VisManager::Disable() {
  SetConcreteInstance(0);
  if (fVerbosity >= confirmations) {
    G4cout <<
    "G4VisManager::Disable: visualization disabled."
    "\n  The pointer returned by GetConcreteInstance will be zero."
    "\n  Note that it will become enabled after some valid vis commands."
	   << G4endl;
  }
  if (fVerbosity >= warnings) {
    G4int currentTrajectoryType =
    G4RunManagerKernel::GetRunManagerKernel()->GetTrackingManager()->GetStoreTrajectory();
    if (currentTrajectoryType > 0) {
    G4cout <<
      "You may wish to disable trajectory production too:"
      "\n  \"/tracking/storeTrajectory 0\""
      "\nbut don't forget to re-enable with"
      "\n  \"/vis/enable\""
      "\n  \"/tracking/storeTrajectory " << currentTrajectoryType << "\" (for your case)."
      << G4endl;
    }
  }
}

const G4GraphicsSystemList& G4VisManager::GetAvailableGraphicsSystems () {
  G4int nSystems = fAvailableGraphicsSystems.size ();
  if (nSystems == 0) {
    if (fVerbosity >= warnings) {
      G4cout << "G4VisManager::GetAvailableGraphicsSystems: WARNING: no"
	"\n graphics system available!"
	"\n  1) Did you have environment variables G4VIS_BUILD_xxxx_DRIVER set"
	"\n     when you compiled/built the visualization code?"
	"\n  2) Did you instantiate your own Visualization Manager and forget"
	"\n     to implement RegisterGraphicsSystems correctly?"
	"\n  3) You can register your own graphics system, e.g.,"
	"\n     visManager->RegisterGraphicsSystem(new MyGraphicsSystem);)"
	"\n     after instantiating your vis manager and before"
	"\n     visManager->Initialize()."
	     << G4endl;
    }
  }
  return fAvailableGraphicsSystems;
}

G4bool G4VisManager::RegisterGraphicsSystem (G4VGraphicsSystem* pSystem) {
  G4bool happy = true;
  if (pSystem) {
    fAvailableGraphicsSystems.push_back (pSystem);
    if (fVerbosity >= confirmations) {
      G4cout << "G4VisManager::RegisterGraphicsSystem: "
	     << pSystem -> GetName ();
      if (pSystem -> GetNickname () != "") {
	G4cout << " (" << pSystem -> GetNickname () << ")";
      }
      G4cout << " registered." << G4endl;
    }
  }
  else {
    if (fVerbosity >= errors) {
      G4cout << "G4VisManager::RegisterGraphicsSystem: null pointer!"
	     << G4endl;
    }
    happy=false;
  }
  return happy;
}

const G4VTrajectoryModel*
G4VisManager::CurrentTrajDrawModel() const
{
  assert (0 != fpTrajDrawModelMgr);

  const G4VTrajectoryModel* model = fpTrajDrawModelMgr->Current();

  if (0 == model) {
    // No model was registered with the trajectory model manager.
    // Use G4TrajectoryDrawByCharge as a fallback.
    fpTrajDrawModelMgr->Register(new G4TrajectoryDrawByCharge("DefaultModel"));
    if (fVerbosity >= warnings) {
      G4cout<<"G4VisManager: Using G4TrajectoryDrawByCharge as fallback trajectory model."<<G4endl;
      G4cout<<"See commands in /vis/modeling/trajectories/ for other options."<<G4endl;
    }
  }

  model = fpTrajDrawModelMgr->Current();
  assert (0 != model); // Should definitely exist now

  return model;
}

void G4VisManager::RegisterModel(G4VTrajectoryModel* model)
{
  fpTrajDrawModelMgr->Register(model);
}

void
G4VisManager::RegisterModelFactory(G4TrajDrawModelFactory* factory) 
{
  fpTrajDrawModelMgr->Register(factory);
}

void G4VisManager::RegisterModel(G4VFilter<G4VTrajectory>* model)
{
  fpTrajFilterMgr->Register(model);
}

void
G4VisManager::RegisterModelFactory(G4TrajFilterFactory* factory)
{
  fpTrajFilterMgr->Register(factory);
}

void G4VisManager::RegisterModel(G4VFilter<G4VHit>* model)
{
  fpHitFilterMgr->Register(model);
}

void
G4VisManager::RegisterModelFactory(G4HitFilterFactory* factory)
{
  fpHitFilterMgr->Register(factory);
}

void G4VisManager::RegisterModel(G4VFilter<G4VDigi>* model)
{
  fpDigiFilterMgr->Register(model);
}

void
G4VisManager::RegisterModelFactory(G4DigiFilterFactory* factory)
{
  fpDigiFilterMgr->Register(factory);
}

void G4VisManager::SelectTrajectoryModel(const G4String& model) 
{
   fpTrajDrawModelMgr->SetCurrent(model);
}

void G4VisManager::BeginDraw (const G4Transform3D& objectTransform)
{
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  fDrawGroupNestingDepth++;
  if (fDrawGroupNestingDepth > 1) {
    G4Exception
      ("G4VisManager::BeginDraw",
       "visman0008", JustWarning,
       "Nesting detected. It is illegal to nest Begin/EndDraw."
       "\n Ignored");
    return;
  }
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fIsDrawGroup = true;
  }
}

void G4VisManager::EndDraw ()
{
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  fDrawGroupNestingDepth--;
  if (fDrawGroupNestingDepth != 0) {
    if (fDrawGroupNestingDepth < 0) fDrawGroupNestingDepth = 0;
    return;
  }
  if (IsValidView ()) {
    fpSceneHandler -> EndPrimitives ();
  }
  fIsDrawGroup = false;
}

void G4VisManager::BeginDraw2D (const G4Transform3D& objectTransform)
{
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  fDrawGroupNestingDepth++;
  if (fDrawGroupNestingDepth > 1) {
    G4Exception
      ("G4VisManager::BeginDraw2D",
       "visman0009", JustWarning,
       "Nesting detected. It is illegal to nest Begin/EndDraw2D."
       "\n Ignored");
    return;
  }
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives2D (objectTransform);
    fIsDrawGroup = true;
  }
}

void G4VisManager::EndDraw2D ()
{
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  fDrawGroupNestingDepth--;
  if (fDrawGroupNestingDepth != 0) {
    if (fDrawGroupNestingDepth < 0) fDrawGroupNestingDepth = 0;
    return;
  }
  if (IsValidView ()) {
    fpSceneHandler -> EndPrimitives2D ();
  }
  fIsDrawGroup = false;
}

template <class T> void G4VisManager::DrawT
(const T& graphics_primitive, const G4Transform3D& objectTransform) {
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  if (fIsDrawGroup) {
    if (objectTransform != fpSceneHandler->GetObjectTransformation()) {
      G4Exception
	("G4VSceneHandler::DrawT",
	 "visman0010", FatalException,
	 "Different transform detected in Begin/EndDraw group.");
    }
    fpSceneHandler -> AddPrimitive (graphics_primitive);
  } else {
    if (IsValidView ()) {
      ClearTransientStoreIfMarked();
      fpSceneHandler -> BeginPrimitives (objectTransform);
      fpSceneHandler -> AddPrimitive (graphics_primitive);
      fpSceneHandler -> EndPrimitives ();
    }
  }
}

template <class T> void G4VisManager::DrawT2D
(const T& graphics_primitive, const G4Transform3D& objectTransform) {
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  if (fIsDrawGroup) {
    if (objectTransform != fpSceneHandler->GetObjectTransformation()) {
      G4Exception
	("G4VSceneHandler::DrawT",
	 "visman0011", FatalException,
	 "Different transform detected in Begin/EndDraw2D group.");
    }
    fpSceneHandler -> AddPrimitive (graphics_primitive);
  } else {
    if (IsValidView ()) {
      ClearTransientStoreIfMarked();
      fpSceneHandler -> BeginPrimitives2D (objectTransform);
      fpSceneHandler -> AddPrimitive (graphics_primitive);
      fpSceneHandler -> EndPrimitives2D ();
    }
  }
}

void G4VisManager::Draw (const G4Circle& circle,
			 const G4Transform3D& objectTransform)
{
  DrawT (circle, objectTransform);
}

void G4VisManager::Draw (const G4Polyhedron& polyhedron,
			 const G4Transform3D& objectTransform)
{
  DrawT (polyhedron, objectTransform);
}

void G4VisManager::Draw (const G4Polyline& line,
			 const G4Transform3D& objectTransform)
{
  DrawT (line, objectTransform);
}

void G4VisManager::Draw (const G4Polymarker& polymarker,
			 const G4Transform3D& objectTransform)
{
  DrawT (polymarker, objectTransform);
}

void G4VisManager::Draw (const G4Scale& scale,
			 const G4Transform3D& objectTransform)
{
  DrawT (scale, objectTransform);
}

void G4VisManager::Draw (const G4Square& square,
			 const G4Transform3D& objectTransform)
{
  DrawT (square, objectTransform);
}

void G4VisManager::Draw (const G4Text& text,
			 const G4Transform3D& objectTransform)
{
  DrawT (text, objectTransform);
}

void G4VisManager::Draw2D (const G4Circle& circle,
			   const G4Transform3D& objectTransform)
{
  DrawT2D (circle, objectTransform);
}

void G4VisManager::Draw2D (const G4Polyhedron& polyhedron,
			   const G4Transform3D& objectTransform)
{
  DrawT2D (polyhedron, objectTransform);
}

void G4VisManager::Draw2D (const G4Polyline& line,
			   const G4Transform3D& objectTransform)
{
  DrawT2D (line, objectTransform);
}

void G4VisManager::Draw2D (const G4Polymarker& polymarker,
			   const G4Transform3D& objectTransform)
{
  DrawT2D (polymarker, objectTransform);
}

void G4VisManager::Draw2D (const G4Square& square,
			   const G4Transform3D& objectTransform)
{
  DrawT2D (square, objectTransform);
}

void G4VisManager::Draw2D (const G4Text& text,
			   const G4Transform3D& objectTransform)
{
  DrawT2D (text, objectTransform);
}

void G4VisManager::Draw (const G4VHit& hit) {
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  if (fIsDrawGroup) {
    fpSceneHandler -> AddCompound (hit);
  } else {
    if (IsValidView ()) {
      ClearTransientStoreIfMarked();
      fpSceneHandler -> AddCompound (hit);
    }
  }
}

void G4VisManager::Draw (const G4VDigi& digi) {
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  if (fIsDrawGroup) {
    fpSceneHandler -> AddCompound (digi);
  } else {
    if (IsValidView ()) {
      ClearTransientStoreIfMarked();
      fpSceneHandler -> AddCompound (digi);
    }
  }
}

void G4VisManager::Draw (const G4VTrajectory& traj) {
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  // A trajectory needs a trajectories model to provide G4Atts, etc.
  static G4TrajectoriesModel trajectoriesModel;
  trajectoriesModel.SetCurrentTrajectory(&traj);
  G4RunManager* runManager = G4RunManager::GetRunManager();
#ifdef G4MULTITHREADED
  if (G4Threading::IsMultithreadedApplication()) {
    runManager = G4MTRunManager::GetMasterRunManager();
  }
#endif
  const G4Run* currentRun = runManager->GetCurrentRun();
  if (currentRun) {
    trajectoriesModel.SetRunID(currentRun->GetRunID());
  }
  const G4Event* currentEvent =
  G4EventManager::GetEventManager()->GetConstCurrentEvent();
  if (currentEvent) {
    trajectoriesModel.SetEventID(currentEvent->GetEventID());
  }
  if (fIsDrawGroup) {
    fpSceneHandler -> SetModel (&trajectoriesModel);
    fpSceneHandler -> AddCompound (traj);
    fpSceneHandler -> SetModel (0);
  } else {
    if (IsValidView ()) {
      ClearTransientStoreIfMarked();
      fpSceneHandler -> SetModel (&trajectoriesModel);
      fpSceneHandler -> AddCompound (traj);
      fpSceneHandler -> SetModel (0);
    }
  }
}

void G4VisManager::Draw (const G4LogicalVolume& logicalVol,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  // Find corresponding solid.
  G4VSolid* pSol = logicalVol.GetSolid ();
  Draw (*pSol, attribs, objectTransform);
}

void G4VisManager::Draw (const G4VSolid& solid,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  if (fIsDrawGroup) {
    fpSceneHandler -> PreAddSolid (objectTransform, attribs);
    solid.DescribeYourselfTo (*fpSceneHandler);
    fpSceneHandler -> PostAddSolid ();
  } else {
    if (IsValidView ()) {
      ClearTransientStoreIfMarked();
      fpSceneHandler -> PreAddSolid (objectTransform, attribs);
      solid.DescribeYourselfTo (*fpSceneHandler);
      fpSceneHandler -> PostAddSolid ();
    }
  }
}

void G4VisManager::Draw (const G4VPhysicalVolume& physicalVol,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
  // Note: It is tempting to use a temporary model here, as for
  // trajectories, in order to get at the G4Atts of the physical
  // volume.  I tried it (JA).  But it's not easy to pass the
  // vis attributes.  Also other aspects of the model seem not to
  // be properly set up.  So, the idea has been abandoned for the time
  // being.  The model pointer will be null.  So when picking there
  // will be no G4Atts from this physical volume.
  //
  // If this is called from DrawHit, for example, the user may G4Atts to the
  // hit and these will be available with "/vis/scene/add/hits".
  //
  // Find corresponding logical volume and solid.
  G4LogicalVolume* pLV  = physicalVol.GetLogicalVolume ();
  G4VSolid*        pSol = pLV -> GetSolid ();
  Draw (*pSol, attribs, objectTransform);
}

void G4VisManager::CreateSceneHandler (const G4String& name) {
  if (!fInitialised) Initialise ();
  if (fpGraphicsSystem) {
    G4VSceneHandler* pSceneHandler =
      fpGraphicsSystem -> CreateSceneHandler (name);
    if (pSceneHandler) {
      fAvailableSceneHandlers.push_back (pSceneHandler);
      fpSceneHandler = pSceneHandler;                         // Make current.
    }
    else {
      if (fVerbosity >= errors) {
	G4cout << "ERROR in G4VisManager::CreateSceneHandler during "
	       << fpGraphicsSystem -> GetName ()
	       << " scene handler creation.\n  No action taken."
	       << G4endl;
      }
    }
  }
  else PrintInvalidPointers ();
}

void G4VisManager::CreateViewer
(const G4String& name, const G4String& XGeometry)
{

  if (!fInitialised) Initialise ();

  if (!fpSceneHandler) {
    PrintInvalidPointers ();
    return;
  }

  G4VViewer* p = fpGraphicsSystem -> CreateViewer (*fpSceneHandler, name);

  if (!p) {
    if (fVerbosity >= errors) {
      G4cerr << "ERROR in G4VisManager::CreateViewer: null pointer during "
	     << fpGraphicsSystem -> GetName ()
	     << " viewer creation.\n  No action taken."
	     << G4endl;
    }
    return;
  }

  if (p -> GetViewId() < 0) {
    if (fVerbosity >= errors) {
      G4cerr << "ERROR in G4VisManager::CreateViewer during "
	     << fpGraphicsSystem -> GetName ()
	     << " viewer instantiation.\n  No action taken."
	     << G4endl;
    }
    return;
  }

  // Viewer is created, now we can set geometry parameters
  // Before 12/2008, it was done in G4VViewer.cc but it did not have to be there!
    
  G4ViewParameters initialvp = p -> GetViewParameters();
  initialvp.SetXGeometryString(XGeometry); //parse string and store parameters
  p -> SetViewParameters(initialvp);
  p -> Initialise ();  // (Viewer itself may change view parameters further.)
  if (p -> GetViewId() < 0) {
    if (fVerbosity >= errors) {
      G4cerr << "ERROR in G4VisManager::CreateViewer during "
	     << fpGraphicsSystem -> GetName ()
	     << " viewer initialisation.\n  No action taken."
	     << G4endl;
    }
    return;
  }

  fpViewer = p;                             // Make current.
  fpSceneHandler -> AddViewerToList (fpViewer);
  fpSceneHandler -> SetCurrentViewer (fpViewer);
  if (fVerbosity >= confirmations) {
    G4cout << "G4VisManager::CreateViewer: new viewer created."
	   << G4endl;
  }

  const G4ViewParameters& vp = fpViewer->GetViewParameters();
  if (fVerbosity >= parameters) {
    G4cout << " view parameters are:\n  " << vp << G4endl;
  }

  if (vp.IsCulling () && vp.IsCullingInvisible ()) {
    static G4bool warned = false;
    if (fVerbosity >= confirmations) {
      if (!warned) {
	G4cout <<
  "NOTE: objects with visibility flag set to \"false\""
  " will not be drawn!"
  "\n  \"/vis/viewer/set/culling global false\" to Draw such objects."
  "\n  Also see other \"/vis/viewer/set\" commands."
	       << G4endl;
	warned = true;
      }
    }
  }
  if (vp.IsCullingCovered ()) {
    static G4bool warned = false;
    if (fVerbosity >= warnings) {
      if (!warned) {
	G4cout <<
  "WARNING: covered objects in solid mode will not be rendered!"
  "\n  \"/vis/viewer/set/culling coveredDaughters false\" to reverse this."
  "\n  Also see other \"/vis/viewer/set\" commands."
	       << G4endl;
	warned = true;
      }
    }
  }
}

void G4VisManager::GeometryHasChanged () {
  if (fVerbosity >= confirmations) {
    G4cout << "G4VisManager::GeometryHasChanged() called." << G4endl;
  }

  // Change the world...
  G4VPhysicalVolume* pWorld =
    G4TransportationManager::GetTransportationManager ()
    -> GetNavigatorForTracking () -> GetWorldVolume ();
  if (!pWorld) {
    if (fVerbosity >= warnings) {
      G4cout << "WARNING: There is no world volume!" << G4endl;
    }
  }

  // Check scenes.
  G4SceneList& sceneList = fSceneList;
  G4int iScene, nScenes = sceneList.size ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    G4Scene* pScene = sceneList [iScene];
    std::vector<G4Scene::Model>& modelList = pScene -> SetRunDurationModelList ();
    if (modelList.size ()) {
      G4bool modelInvalid;
      do {  // Remove, if required, one at a time.
	modelInvalid = false;
	std::vector<G4Scene::Model>::iterator iterModel;
	for (iterModel = modelList.begin();
	     iterModel != modelList.end();
	     ++iterModel) {
	  modelInvalid = !(iterModel->fpModel->Validate(fVerbosity>=warnings));
	  if (modelInvalid) {
	    // Model invalid - remove and break.
	    if (fVerbosity >= warnings) {
	      G4cout << "WARNING: Model \""
		     << iterModel->fpModel->GetGlobalDescription ()
		     <<
		"\" is no longer valid - being removed\n  from scene \""
		     << pScene -> GetName () << "\""
		     << G4endl;
	    }
	    modelList.erase (iterModel);
	    break;
	  }
	}
      } while (modelInvalid);

      if (modelList.size () == 0) {
	if (fVerbosity >= warnings) {
	  G4cout << "WARNING: No models left in this scene \""
		 << pScene -> GetName ()
		 << "\"."
		 << G4endl;
	}
      }
      else {
	pScene->CalculateExtent();
	G4UImanager::GetUIpointer () ->
	  ApplyCommand (G4String("/vis/scene/notifyHandlers " + pScene->GetName()));
      }
    }
  }

  // Check the manager's current scene...
  if (fpScene && fpScene -> GetRunDurationModelList ().size () == 0) {
    if (fVerbosity >= warnings) {
      G4cout << "WARNING: The current scene \""
	     << fpScene -> GetName ()
	     << "\" has no models."
	     << G4endl;
    }
  }
}

void G4VisManager::NotifyHandlers () {

  if (fVerbosity >= confirmations) {
    G4cout << "G4VisManager::NotifyHandler() called." << G4endl;
  }

  // Check scenes.
  G4SceneList& sceneList = fSceneList;
  G4int iScene, nScenes = sceneList.size ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    G4Scene* pScene = sceneList [iScene];
    std::vector<G4Scene::Model>& modelList = pScene -> SetRunDurationModelList ();
    
    if (modelList.size ()) {
      pScene->CalculateExtent();
      G4UImanager::GetUIpointer () ->
        ApplyCommand (G4String("/vis/scene/notifyHandlers " + pScene->GetName()));
    }
  }

  // Check the manager's current scene...
  if (fpScene && fpScene -> GetRunDurationModelList ().size () == 0) {
    if (fVerbosity >= warnings) {
      G4cout << "WARNING: The current scene \""
	     << fpScene -> GetName ()
	     << "\" has no models."
	     << G4endl;
    }
  }
}

G4bool G4VisManager::FilterTrajectory(const G4VTrajectory& trajectory)
{
  return fpTrajFilterMgr->Accept(trajectory);
}   

G4bool G4VisManager::FilterHit(const G4VHit& hit)
{
  return fpHitFilterMgr->Accept(hit);
}   

G4bool G4VisManager::FilterDigi(const G4VDigi& digi)
{
  return fpDigiFilterMgr->Accept(digi);
}   

void G4VisManager::DispatchToModel(const G4VTrajectory& trajectory)
{
  G4bool visible(true);

  // See if trajectory passes filter
  G4bool passed = FilterTrajectory(trajectory);

  if (!passed) {
    // Draw invisible trajectory if trajectory failed filter and
    // are filtering in soft mode
    if (fpTrajFilterMgr->GetMode() == FilterMode::Soft) visible = false;
    else {return;}
  }

  // Go on to draw trajectory
  assert (0 != fpTrajDrawModelMgr);

  const G4VTrajectoryModel* trajectoryModel = CurrentTrajDrawModel();

  assert (0 != trajectoryModel); // Should exist

  if (IsValidView()) {
      trajectoryModel->Draw(trajectory, visible);
  }
}

void G4VisManager::SetUserAction
(G4VUserVisAction* pVisAction,
 const G4VisExtent& extent) {
  if (fVerbosity >= warnings) {
    G4cout <<
  "WARNING: SetUserAction is deprecated. Use RegisterRunDurationUserVisAction."
    << G4endl;
  }
  RegisterRunDurationUserVisAction("SetUserAction",pVisAction,extent);
}

void G4VisManager::RegisterRunDurationUserVisAction
(const G4String& name,
 G4VUserVisAction* pVisAction,
 const G4VisExtent& extent) {
  fRunDurationUserVisActions.push_back(UserVisAction(name,pVisAction));
  if (extent.GetExtentRadius() > 0.) {
    fUserVisActionExtents[pVisAction] = extent;
  } else {
    if (fVerbosity >= warnings) {
      G4cout << 
	"WARNING: No extent set for user vis action \"" << name << "\"."
	     << G4endl;
    }
  }
}

void G4VisManager::RegisterEndOfEventUserVisAction
(const G4String& name,
 G4VUserVisAction* pVisAction,
 const G4VisExtent& extent) {
  fEndOfEventUserVisActions.push_back(UserVisAction(name,pVisAction));
  if (extent.GetExtentRadius() > 0.) {
    fUserVisActionExtents[pVisAction] = extent;
  } else {
    if (fVerbosity >= warnings) {
      G4cout << 
	"WARNING: No extent set for user vis action \"" << name << "\"."
	     << G4endl;
    }
  }
}

void G4VisManager::RegisterEndOfRunUserVisAction
(const G4String& name,
 G4VUserVisAction* pVisAction,
 const G4VisExtent& extent) {
  fEndOfRunUserVisActions.push_back(UserVisAction(name,pVisAction));
  if (extent.GetExtentRadius() > 0.) {
    fUserVisActionExtents[pVisAction] = extent;
  } else {
    if (fVerbosity >= warnings) {
      G4cout << 
	"WARNING: No extent set for user vis action \"" << name << "\"."
	     << G4endl;
    }
  }
}

void G4VisManager::SetCurrentScene (G4Scene* pScene) {
  if (pScene != fpScene) {
    // A change of scene.  Therefore reset transients drawn flags.  All
    // memory of previous transient proceessing thereby erased...
    ResetTransientsDrawnFlags();
  }
  fpScene = pScene;
}

void G4VisManager::SetCurrentGraphicsSystem (G4VGraphicsSystem* pSystem) {
  fpGraphicsSystem = pSystem;
  if (fVerbosity >= confirmations) {
    G4cout << "G4VisManager::SetCurrentGraphicsSystem: system now "
	   << pSystem -> GetName () << G4endl;
  }
  // If current scene handler is of same graphics system, leave unchanged.
  // Else find the most recent scene handler of same graphics system.
  // Or clear pointers.
  if (!(fpSceneHandler && fpSceneHandler -> GetGraphicsSystem () == pSystem)) {
    const G4SceneHandlerList& sceneHandlerList = fAvailableSceneHandlers;
    G4int nSH = sceneHandlerList.size ();  // No. of scene handlers.
    G4int iSH;
    for (iSH = nSH - 1; iSH >= 0; iSH--) {
      if (sceneHandlerList [iSH] -> GetGraphicsSystem () == pSystem) break;
    }
    if (iSH >= 0) {
      fpSceneHandler = sceneHandlerList [iSH];
      if (fVerbosity >= confirmations) {
	G4cout << "  Scene Handler now "
	       << fpSceneHandler -> GetName () << G4endl;
      }
      if (fpScene != fpSceneHandler -> GetScene ()) {
	fpScene = fpSceneHandler -> GetScene ();
	if (fVerbosity >= confirmations) {
	  G4cout << "  Scene now \""
		 << fpScene -> GetName () << "\"" << G4endl;
	}
      }
      const G4ViewerList& viewerList = fpSceneHandler -> GetViewerList ();
      if (viewerList.size ()) {
	fpViewer = viewerList [0];
	if (fVerbosity >= confirmations) {
	  G4cout << "  Viewer now " << fpViewer -> GetName () << G4endl;
	}
      }
      else {
	fpViewer = 0;
      }
    }
    else {
      fpSceneHandler = 0;
      fpViewer = 0;
    }
  }
}

void G4VisManager::SetCurrentSceneHandler (G4VSceneHandler* pSceneHandler) {
  fpSceneHandler = pSceneHandler;
  if (fVerbosity >= confirmations) {
    G4cout << "G4VisManager::SetCurrentSceneHandler: scene handler now \""
	   << pSceneHandler -> GetName () << "\"" << G4endl;
  }
  if (fpScene != fpSceneHandler -> GetScene ()) {
    fpScene = fpSceneHandler -> GetScene ();
    if (fVerbosity >= confirmations) {
      G4cout << "  Scene now \""
	     << fpScene -> GetName () << "\"" << G4endl;
    }
  }
  if (fpGraphicsSystem != pSceneHandler -> GetGraphicsSystem ()) {
    fpGraphicsSystem = pSceneHandler -> GetGraphicsSystem ();
    if (fVerbosity >= confirmations) {
      G4cout << "  Graphics system now \""
	     << fpGraphicsSystem -> GetName () << "\"" << G4endl;
    }
  }
  const G4ViewerList& viewerList = fpSceneHandler -> GetViewerList ();
  G4int nViewers = viewerList.size ();
  if (nViewers) {
    G4int iViewer;
    for (iViewer = 0; iViewer < nViewers; iViewer++) {
      if (fpViewer == viewerList [iViewer]) break;
    }
    if (iViewer >= nViewers) {
      fpViewer = viewerList [0];
      if (fVerbosity >= confirmations) {
	G4cout << "  Viewer now \"" << fpViewer -> GetName () << "\""
	       << G4endl;
      }
    }
    if (!IsValidView ()) {
      if (fVerbosity >= warnings) {
	G4cout <<
  "WARNING: Problem setting scene handler - please report circumstances."
	       << G4endl;
      }
    }
  }
  else {
    fpViewer = 0;
    if (fVerbosity >= warnings) {
      G4cout <<
	"WARNING: No viewers for this scene handler - please create one."
	     << G4endl;
    }
  }
}

void G4VisManager::SetCurrentViewer (G4VViewer* pViewer) {
  fpViewer  = pViewer;
  if (fVerbosity >= confirmations) {
    G4cout << "G4VisManager::SetCurrentViewer: viewer now "
	   << pViewer -> GetName ()
	   << G4endl;
  }
  fpSceneHandler = fpViewer -> GetSceneHandler ();
  if (!fpSceneHandler) {
    if (fVerbosity >= warnings) {
      G4cout <<
      "WARNING: No scene handler for this viewer - please create one."
      << G4endl;
    }
    return;
  }
  fpViewer->SetView();
  fpSceneHandler -> SetCurrentViewer (pViewer);
  fpScene = fpSceneHandler -> GetScene ();
  fpGraphicsSystem = fpSceneHandler -> GetGraphicsSystem ();
  if (!IsValidView ()) {
    if (fVerbosity >= warnings) {
      G4cout <<
	"WARNING: Problem setting viewer - please report circumstances."
	     << G4endl;
    }
  }
}

void G4VisManager::PrintAvailableGraphicsSystems (Verbosity verbosity) const
{
  G4cout << "Current available graphics systems are:\n";
  if (fAvailableGraphicsSystems.size ()) {
    for (const auto& gs: fAvailableGraphicsSystems) {
      const G4String& name = gs->GetName();
      const std::vector<G4String>& nicknames = gs->GetNicknames();
      if (verbosity <= warnings) {
        // Brief output
        G4cout << name << " (";
        for (size_t i = 0; i < nicknames.size(); ++i) {
          if (i != 0) {
            G4cout << ", ";
          }
          G4cout << nicknames[i];
        }
        G4cout << ')';
      } else {
        // Full output
        G4cout << *gs;
      }
      G4cout << G4endl;
    }
  } else {
    G4cout << "\n  NONE!!!  None registered - yet!  Mmmmm!" << G4endl;
  }
}

void G4VisManager::PrintAvailableModels (Verbosity verbosity) const
{
  {
    //fpTrajDrawModelMgr->Print(G4cout);
    G4cout << "Registered model factories:" << G4endl;
    const std::vector<G4VModelFactory<G4VTrajectoryModel>*>& factoryList =
      fpTrajDrawModelMgr->FactoryList();
    if (factoryList.empty()) G4cout << "  None" << G4endl;
    else {
      std::vector<G4VModelFactory<G4VTrajectoryModel>*>::const_iterator i;
      for (i = factoryList.begin(); i != factoryList.end(); ++i)
	(*i)->Print(G4cout);
    }
    const G4VisListManager<G4VTrajectoryModel>* listManager =
      fpTrajDrawModelMgr->ListManager();
    const std::map<G4String, G4VTrajectoryModel*>& modelMap =
      listManager->Map();
    if (!modelMap.empty()) {
      G4cout << "\nRegistered models:" << G4endl;
      std::map<G4String, G4VTrajectoryModel*>::const_iterator i;
      for (i = modelMap.begin(); i != modelMap.end(); ++i) {
	G4cout << "  " << i->second->Name();
	if (i->second == listManager->Current()) G4cout << " (Current)";
	G4cout << G4endl;
	if (verbosity >= parameters) i->second->Print(G4cout);
      }
    }
  }

  G4cout << G4endl;

  {
    //fpTrajFilterMgr->Print(G4cout);
    G4cout << "Registered filter factories:" << G4endl;
    const std::vector<G4VModelFactory<G4VFilter<G4VTrajectory> >*>&
      factoryList = fpTrajFilterMgr->FactoryList();
    if (factoryList.empty()) G4cout << "  None" << G4endl;
    else {
      std::vector<G4VModelFactory<G4VFilter<G4VTrajectory> >*>::const_iterator i;
      for (i = factoryList.begin(); i != factoryList.end(); ++i)
	(*i)->Print(G4cout);
    }
    const std::vector<G4VFilter<G4VTrajectory>*>&
      filterList = fpTrajFilterMgr->FilterList();
    if (!filterList.empty()) {
      G4cout << "\nRegistered filters:" << G4endl;
      std::vector<G4VFilter<G4VTrajectory>*>::const_iterator i;
      for (i = filterList.begin(); i != filterList.end(); ++i) {
	G4cout << "  " << (*i)->GetName() << G4endl;
	if (verbosity >= parameters) (*i)->PrintAll(G4cout);
      }
    }
  }
}

void G4VisManager::PrintAvailableUserVisActions (Verbosity) const
{
  G4cout <<
    "You have successfully registered the following user vis actions."
	 << G4endl;
  G4cout << "Run Duration User Vis Actions:";
  if (fRunDurationUserVisActions.empty()) G4cout << " none" << G4endl;
  else {
    G4cout << G4endl;
    for (size_t i = 0; i < fRunDurationUserVisActions.size(); i++) {
      const G4String& name = fRunDurationUserVisActions[i].fName;
      G4cout << "  " << name << G4endl;
    }
  }

  G4cout << "End of Event User Vis Actions:";
  if (fEndOfEventUserVisActions.empty()) G4cout << " none" << G4endl;
  else {
    G4cout << G4endl;
    for (size_t i = 0; i < fEndOfEventUserVisActions.size(); i++) {
      const G4String& name = fEndOfEventUserVisActions[i].fName;
      G4cout << "  " << name << G4endl;
    }
  }

  G4cout << "End of Run User Vis Actions:";
  if (fEndOfRunUserVisActions.empty()) G4cout << " none" << G4endl;
  else {
    G4cout << G4endl;
    for (size_t i = 0; i < fEndOfRunUserVisActions.size(); i++) {
      const G4String& name = fEndOfRunUserVisActions[i].fName;
      G4cout << "  " << name << G4endl;
    }
  }
}

void G4VisManager::PrintAvailableColours (Verbosity) const {
  G4cout <<
    "Some /vis commands (optionally) take a string to specify colour."
    "\nAvailable colours:\n  ";
  const std::map<G4String, G4Colour>& map = G4Colour::GetMap();
  for (std::map<G4String, G4Colour>::const_iterator i = map.begin();
       i != map.end();) {
    G4cout << i->first;
    if (++i != map.end()) G4cout << ", ";
  }
  G4cout << G4endl;
}

void G4VisManager::PrintInvalidPointers () const {
  if (fVerbosity >= errors) {
    G4cerr << "ERROR: G4VisManager::PrintInvalidPointers:";
    if (!fpGraphicsSystem) {
      G4cerr << "\n null graphics system pointer.";
    }
    else {
      G4cerr << "\n  Graphics system is " << fpGraphicsSystem -> GetName ()
	     << " but:";
      if (!fpScene)
	G4cerr <<
	  "\n  Null scene pointer. Use \"/vis/drawVolume\" or"
	  " \"/vis/scene/create\".";
      if (!fpSceneHandler)
	G4cerr <<
	  "\n  Null scene handler pointer. Use \"/vis/open\" or"
	  " \"/vis/sceneHandler/create\".";
      if (!fpViewer )
	G4cerr <<
	  "\n  Null viewer pointer. Use \"/vis/viewer/create\".";
    }
    G4cerr << G4endl;
  }
}

#ifdef G4MULTITHREADED

namespace {
  G4bool mtRunInProgress = false;
  std::deque<const G4Event*> mtVisEventQueue;
  G4Thread* mtVisSubThread = 0;
  G4Mutex mtVisSubThreadMutex = G4MUTEX_INITIALIZER;
}

G4ThreadFunReturnType G4VisManager::G4VisSubThread(G4ThreadFunArgType p)
{
  G4VisManager* pVisManager = (G4VisManager*)p;
  G4VSceneHandler* pSceneHandler = pVisManager->GetCurrentSceneHandler();
  if (!pSceneHandler) return 0;
  G4Scene* pScene = pSceneHandler->GetScene();
  if (!pScene) return 0;
  G4VViewer* pViewer = pVisManager->GetCurrentViewer();
  if (!pViewer) return 0;

  G4UImanager::GetUIpointer()->SetUpForSpecialThread("G4VIS");

//  G4cout << "G4VisManager::G4VisSubThread: thread: "
//  << G4Threading::G4GetThreadId() << std::endl;

  // Set up geometry and navigation for a thread
  G4GeometryWorkspacePool::GetInstance()->CreateAndUseWorkspace();
  G4SolidsWorkspacePool::GetInstance()->CreateAndUseWorkspace();
  G4Navigator* navigator =
  G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  navigator->SetWorldVolume
  (G4MTRunManager::GetMasterRunManagerKernel()->GetCurrentWorld());

  pViewer->SwitchToVisSubThread();

  while (true) {

    G4MUTEXLOCK(&mtVisSubThreadMutex);
    G4int eventQueueSize = mtVisEventQueue.size();
    G4MUTEXUNLOCK(&mtVisSubThreadMutex);
//    G4cout << "Event queue size (A): " << eventQueueSize << G4endl;

    while (eventQueueSize) {

      G4MUTEXLOCK(&mtVisSubThreadMutex);
      const G4Event* event = mtVisEventQueue.front();
      G4MUTEXUNLOCK(&mtVisSubThreadMutex);
//      G4int eventID = event->GetEventID();
//      G4cout <<
//      "G4VisManager::G4VisSubThread: Vis sub-thread: Dealing with event: "
//      << eventID << G4endl;

      // Here comes the event drawing
      pVisManager->SetTransientsDrawnThisEvent(false);
      pSceneHandler->SetTransientsDrawnThisEvent(false);

      // We are about to draw the event (trajectories, etc.), but first we
      // have to clear the previous event(s) if necessary.  If this event
      // needs to be drawn afresh, e.g., the first event or any event when
      // "accumulate" is not requested, the old event has to be cleared.
      // We have postponed this so that, for normal viewers like OGL, the
      // previous event(s) stay on screen until this new event comes
      // along.  For a file-writing viewer the geometry has to be drawn.
      // See, for example, G4HepRepFileSceneHandler::ClearTransientStore.
      pVisManager->ClearTransientStoreIfMarked();

      // Now draw the event...
      pSceneHandler->DrawEvent(event);
      ++pVisManager->fNoOfEventsDrawnThisRun;

      if (pScene->GetRefreshAtEndOfEvent()) {

        // ShowView guarantees the view is flushed to the screen.  It also
        // triggers other features such picking (if enabled) and allows
        // file-writing viewers to close the file.
        pViewer->ShowView();
        pSceneHandler->SetMarkForClearingTransientStore(true);

      }

      // Testing.
//      std::this_thread::sleep_for(std::chrono::seconds(5));

      // Then pop and release event
      G4MUTEXLOCK(&mtVisSubThreadMutex);
      mtVisEventQueue.pop_front();
      event->PostProcessingFinished();
      eventQueueSize = mtVisEventQueue.size();
      G4MUTEXUNLOCK(&mtVisSubThreadMutex);
//      G4cout << "Event queue size (B): " << eventQueueSize << G4endl;
    }

    G4MUTEXLOCK(&mtVisSubThreadMutex);
    G4int runInProgress = mtRunInProgress;
    G4MUTEXUNLOCK(&mtVisSubThreadMutex);
    if (!runInProgress) {
      // EndOfRun on master thread has signalled end of run.  There is
      // nothing to draw so...
      break;
    }

    // Run still in progress but nothing to draw, so wait a while.
#ifdef G4VIS_USE_STD11
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
#else
    G4THREADSLEEP(1);
#endif
  }

  // Inform viewer that we have finished all sub-thread drawing
  pViewer->DoneWithVisSubThread();
  pViewer->MovingToMasterThread();
//  G4cout << "G4VisManager::G4VisSubThread: Vis sub-thread: ending" << G4endl;
  return /*(G4ThreadFunReturnType)*/0;
}

namespace {
  //  G4Mutex visBeginOfRunMutex = G4MUTEX_INITIALIZER;
  //  G4Mutex visBeginOfEventMutex = G4MUTEX_INITIALIZER;
  G4Mutex visEndOfEventMutex = G4MUTEX_INITIALIZER;
  //  G4Mutex visEndOfRunMutex = G4MUTEX_INITIALIZER;
}

#endif

void G4VisManager::BeginOfRun ()
{
  if (fIgnoreStateChanges) return;

  if (!GetConcreteInstance()) return;

#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif
//  G4cout << "G4VisManager::BeginOfRun: thread: "
//  << G4Threading::G4GetThreadId() << G4endl;

  G4RunManager* runManager = G4RunManager::GetRunManager();
#ifdef G4MULTITHREADED
  if (G4Threading::IsMultithreadedApplication()) {
    runManager = G4MTRunManager::GetMasterRunManager();
  }
#endif

  // For a fake run...
  G4int nEventsToBeProcessed = runManager->GetNumberOfEventsToBeProcessed();
  if (nEventsToBeProcessed == 0) return;

  fNKeepRequests = 0;
  fKeptLastEvent = false;
  fEventKeepingSuspended = false;
  fTransientsDrawnThisRun = false;
  if (fpSceneHandler) fpSceneHandler->SetTransientsDrawnThisRun(false);
  fNoOfEventsDrawnThisRun = 0;

#ifdef G4MULTITHREADED
//   There is a static method G4Threading::IsMultithreadedApplication()
//   that returns true only if G4MTRunManager is instantiated with MT
//   installation. Thus method returns false if G4RunManager base class is
//   instantiated even with the MT installation, or of course with sequential
//   installation.
  if (G4Threading::IsMultithreadedApplication()) {

    // Inform viewer that we have finished all master thread drawing for now...
    if (fpViewer) fpViewer->DoneWithMasterThread();

    // Start vis sub-thread
//    G4cout << "G4VisManager::BeginOfRun: Starting vis sub-thread" << G4endl;
    G4MUTEXLOCK(&mtVisSubThreadMutex);
    mtRunInProgress = true;
    G4MUTEXUNLOCK(&mtVisSubThreadMutex);
    mtVisSubThread = new G4Thread;
    // Launch vis thread
    G4THREADCREATE(mtVisSubThread,G4VisSubThread,this);

    // Tricky things for some viewers (e.g., Qt):
    // - Launch the vis thread
    // - Wait for the vis thread to set its QThread
    // - Then move current QOpenGL context (if Qt) to this Qthread
    // - Go ahead
    if (fpViewer) fpViewer->MovingToVisSubThread();
  }
#endif
}

void G4VisManager::BeginOfEvent ()
{
  if (fIgnoreStateChanges) return;

  if (!GetConcreteInstance()) return;

//  G4cout << "G4VisManager::BeginOfEvent: thread: "
//  << G4Threading::G4GetThreadId() << G4endl;

  // Some instructions that should NOT be in multithreaded version.
#ifndef G4MULTITHREADED
  // These instructions are in G4VisSubThread for multithreading.
  fTransientsDrawnThisEvent = false;
  if (fpSceneHandler) fpSceneHandler->SetTransientsDrawnThisEvent(false);
#endif
}

void G4VisManager::EndOfEvent ()
{
  if (fIgnoreStateChanges) return;

  if (!GetConcreteInstance()) return;

  // Don't call IsValidView unless there is a scene handler.  This
  // avoids WARNING message at end of event and run when the user has
  // not instantiated a scene handler, e.g., in batch mode.
  G4bool valid = fpSceneHandler && IsValidView();
  if (!valid) return;

//  G4cout << "G4VisManager::EndOfEvent: thread: "
//  << G4Threading::G4GetThreadId() << G4endl;

#ifdef G4MULTITHREADED
  G4AutoLock al(&visEndOfEventMutex);
  // Testing.
//  std::this_thread::sleep_for(std::chrono::seconds(5));
#endif

  G4RunManager* runManager = G4RunManager::GetRunManager();
#ifdef G4MULTITHREADED
  if (G4Threading::IsMultithreadedApplication()) {
    runManager = G4MTRunManager::GetMasterRunManager();
  }
#endif

  const G4Run* currentRun = runManager->GetCurrentRun();
  if (!currentRun) return;

  // This gets the thread-local event manager
  G4EventManager* eventManager = G4EventManager::GetEventManager();
  const G4Event* currentEvent = eventManager->GetConstCurrentEvent();
  if (!currentEvent) return;

  if (G4Threading::IsMultithreadedApplication()) {

#ifdef G4MULTITHREADED

    // Wait if too many events in the queue.
    G4MUTEXLOCK(&mtVisSubThreadMutex);
    G4int eventQueueSize = mtVisEventQueue.size();
    G4MUTEXUNLOCK(&mtVisSubThreadMutex);
//    G4cout << "Event queue size (1): " << eventQueueSize << G4endl;

    G4bool eventQueueFull = false;
    while (fMaxEventQueueSize > 0 && eventQueueSize >= fMaxEventQueueSize) {

//      G4cout << "Event queue size (2): " << eventQueueSize << G4endl;
      if (fWaitOnEventQueueFull) {
        static G4bool warned = false;
        if (!warned) {
          G4cout <<
          "WARNING: The number of events in the visualisation queue has exceeded"
          "\n  the maximum, "
          << fMaxEventQueueSize <<
          ".\n  If, during a multithreaded run, the simulation gets ahead of the"
          "\n  visualisation by more than this maximum, the simulation is delayed"
          "\n  until the vis sub-thread has drawn a few more events and removed them"
          "\n  from the queue.  You may change this maximum number of events with"
          "\n  \"/vis/multithreading/maxEventQueueSize <N>\", where N is the maximum"
          "\n  number you wish to allow.  N <= 0 means \"unlimited\"."
          "\n  Alternatively you may choose to discard events for drawing by setting"
          "\n  \"/vis/multithreading/actionOnEventQueueFull discard\"."
          "\n  To avoid visualisation altogether: \"/vis/disable\"."
          "\n  And maybe \"/tracking/storeTrajectories 0\"."
          << G4endl;
          warned = true;
        }
        //      G4cout << "Event queue size (3): " << eventQueueSize << G4endl;
        // Wait a while to give event drawing time to reduce the queue...
#ifdef G4VIS_USE_STD11
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
#else
        G4THREADSLEEP(1);
#endif
        //      G4cout << "Event queue size (4): " << eventQueueSize << G4endl;
      } else {
        static G4bool warned = false;
        if (!warned) {
          G4cout <<
          "WARNING: The number of events in the visualisation queue has exceeded"
          "\n  the maximum, "
          << fMaxEventQueueSize <<
          ".\n  Some events have been discarded for drawing.  You may change this"
          "\n  behaviour with \"/vis/multithreading/actionOnEventQueueFull wait\"."
          "\n  To avoid visualisation altogether: \"/vis/disable\"."
          "\n  And maybe \"/tracking/storeTrajectories 0\"."
          << G4endl;
          warned = true;
        }
        eventQueueFull = true;  // Causes event to be discarded for drawing.
        break;
      }

      G4MUTEXLOCK(&mtVisSubThreadMutex);
      eventQueueSize = mtVisEventQueue.size();
      G4MUTEXUNLOCK(&mtVisSubThreadMutex);
    }

    if (!eventQueueFull) {
      G4MUTEXLOCK(&mtVisSubThreadMutex);
      // Keep event for processing and put event on vis event queue
      currentEvent->KeepForPostProcessing();
      mtVisEventQueue.push_back(currentEvent);
      G4MUTEXUNLOCK(&mtVisSubThreadMutex);
    }

//    G4MUTEXLOCK(&mtVisSubThreadMutex);
//    G4int eQS = mtVisEventQueue.size();
//    G4MUTEXUNLOCK(&mtVisSubThreadMutex);
//    G4cout << "Event queue size (5): " << eQS << G4endl;

#endif

  } else {

    // Sequential mode

    G4int nEventsToBeProcessed = 0;
    G4int nKeptEvents = 0;
    G4int eventID = -2;  // (If no run manager, triggers ShowView as normal.)
    if (currentRun) {
      nEventsToBeProcessed = runManager->GetNumberOfEventsToBeProcessed();
      eventID = currentEvent->GetEventID();
      const std::vector<const G4Event*>* events = currentRun->GetEventVector();
      if (events) nKeptEvents = events->size();
    }

    // We are about to draw the event (trajectories, etc.), but first we
    // have to clear the previous event(s) if necessary.  If this event
    // needs to be drawn afresh, e.g., the first event or any event when
    // "accumulate" is not requested, the old event has to be cleared.
    // We have postponed this so that, for normal viewers like OGL, the
    // previous event(s) stay on screen until this new event comes
    // along.  For a file-writing viewer the geometry has to be drawn.
    // See, for example, G4HepRepFileSceneHandler::ClearTransientStore.
    ClearTransientStoreIfMarked();

    // Now draw the event...
    fpSceneHandler->DrawEvent(currentEvent);
    ++fNoOfEventsDrawnThisRun;

    if (fpScene->GetRefreshAtEndOfEvent()) {

      // Unless last event (in which case wait end of run)...
      if (eventID < nEventsToBeProcessed - 1) {
        // ShowView guarantees the view is flushed to the screen.  It also
        // triggers other features such picking (if enabled) and allows
        // file-writing viewers to close the file.
        fpViewer->ShowView();
      } else {  // Last event...
                // Keep, but only if user has not kept any...
        if (nKeptEvents == 0) {
          eventManager->KeepTheCurrentEvent();
          fNKeepRequests++;
          fKeptLastEvent = true;
        }
      }
      fpSceneHandler->SetMarkForClearingTransientStore(true);

    }
  }

  // Both modes - sequential and MT

  if (!(fpScene->GetRefreshAtEndOfEvent())) {

    //  Accumulating events...

    G4int maxNumberOfKeptEvents = fpScene->GetMaxNumberOfKeptEvents();

    if (maxNumberOfKeptEvents > 0 &&
        fNKeepRequests >= maxNumberOfKeptEvents) {

      fEventKeepingSuspended = true;
      static G4bool warned = false;
      if (!warned) {
        if (fVerbosity >= warnings) {
          G4cout <<
          "WARNING: G4VisManager::EndOfEvent: Automatic event keeping suspended."
          "\n  The number of events exceeds the maximum, "
          << maxNumberOfKeptEvents <<
          ", that may be kept by\n  the vis manager."
          << G4endl;
        }
        warned = true;
      }

    } else if (maxNumberOfKeptEvents != 0) {

      // If not disabled nor suspended.
      if (GetConcreteInstance() && !fEventKeepingSuspended) {
//        G4cout <<
//        "Requesting keeping event " << currentEvent->GetEventID()
//        << G4endl;
        eventManager->KeepTheCurrentEvent();
        fNKeepRequests++;
      }
      
    }
  }
}

void G4VisManager::EndOfRun ()
{
  if (fIgnoreStateChanges) return;

  if (!GetConcreteInstance()) return;

#ifdef G4MULTITHREADED
  if (G4Threading::IsWorkerThread()) return;
#endif

//  G4cout << "G4VisManager::EndOfRun: thread: "
//  << G4Threading::G4GetThreadId() << G4endl;

  G4RunManager* runManager = G4RunManager::GetRunManager();
#ifdef G4MULTITHREADED
  if (G4Threading::IsMultithreadedApplication()) {
    runManager = G4MTRunManager::GetMasterRunManager();
  }
#endif

  // For a fake run...
  G4int nEventsToBeProcessed = runManager->GetNumberOfEventsToBeProcessed();
  if (nEventsToBeProcessed == 0) return;

  const G4Run* currentRun = runManager->GetCurrentRun();
  if (!currentRun) return;

#ifdef G4MULTITHREADED
  //  G4AutoLock al(&visEndOfRunMutex);  ???
  if (G4Threading::IsMultithreadedApplication()) {
    // Reset flag so that sub-thread exits when it has finished processing.
    G4MUTEXLOCK(&mtVisSubThreadMutex);
    mtRunInProgress = false;
    G4MUTEXUNLOCK(&mtVisSubThreadMutex);
    // Wait for sub-thread to finish.
    G4THREADJOIN(*mtVisSubThread);
    delete mtVisSubThread;
    if (fpViewer) fpViewer->SwitchToMasterThread();
  }
#endif

#ifdef G4MULTITHREADED
  // Print warning about discarded events, if any.
  // Don't call IsValidView unless there is a scene handler.  This
  // avoids WARNING message from IsValidView() when the user has
  // not instantiated a scene handler, e.g., in batch mode.
  if (fpSceneHandler && IsValidView()) {  // Events should have been drawn
    G4int noOfEventsRequested = runManager->GetNumberOfEventsToBeProcessed();
    if (fNoOfEventsDrawnThisRun != noOfEventsRequested) {
      if (!fWaitOnEventQueueFull && fVerbosity >= warnings) {
        G4cout
        << "WARNING: Number of events drawn this run, "
        << fNoOfEventsDrawnThisRun << ", is different to number requested, "
        << noOfEventsRequested <<
        ".\n  (This is because you requested \"/vis/multithreading/actionOnEventQueueFull discard\".)"
        << G4endl;
      }
    }
  }
#endif

  G4int nKeptEvents = 0;
  const std::vector<const G4Event*>* events = currentRun->GetEventVector();
  if (events) nKeptEvents = events->size();
  if (nKeptEvents && !fKeptLastEvent) {
    if (fVerbosity >= warnings) {
      G4cout << "WARNING: " << nKeptEvents;
      if (nKeptEvents == 1) G4cout << " event has";
      else G4cout << " events have";
      G4cout << " been kept for refreshing and/or reviewing." << G4endl;
      if (nKeptEvents != fNKeepRequests) {
        G4cout << "  (Note: ";
        if (fNKeepRequests == 0) {
          G4cout << "No keep requests were";
        } else if (fNKeepRequests == 1) {
          G4cout << "Only 1 keep request was";
        } else {
          G4cout << "Only " << fNKeepRequests << " keep requests were";
        }
        G4cout << " made by the vis manager.)" << G4endl;
      }
      G4cout <<
  "\n  \"/vis/reviewKeptEvents\" to review them one by one."
  "\n  \"/vis/viewer/flush\" or \"/vis/viewer/rebuild\" to see them accumulated."
      << G4endl;
    }
    //    static G4bool warned = false;
    //    if (!valid && fVerbosity >= warnings && !warned) {
    //      G4cout <<
    //	"  Only useful if before starting the run:"
    //	"\n    a) trajectories are stored (\"/vis/scene/add/trajectories [smooth|rich]\"), or"
    //	"\n    b) the Draw method of any hits or digis is implemented."
    //	"\n  To view trajectories, hits or digis:"
    //	"\n    open a viewer, draw a volume, \"/vis/scene/add/trajectories\""
    //	"\n    \"/vis/scene/add/hits\" or \"/vis/scene/add/digitisations\""
    //	"\n    and, possibly, \"/vis/viewer/flush\"."
    //	"\n  To see all events: \"/vis/scene/endOfEventAction accumulate\"."
    //  "\n  (You may need \"/vis/viewer/flush\" or even \"/vis/viewer/rebuild\".)"
    //	"\n  To see events individually: \"/vis/reviewKeptEvents\"."
    //	     << G4endl;
    //      warned = true;
    //    }
  }

  if (fEventKeepingSuspended && fVerbosity >= warnings) {
    G4cout <<
    "WARNING: G4VisManager::EndOfRun: Automatic event keeping was suspended."
    "\n  The number of events in the run exceeded the maximum, "
    << fpScene->GetMaxNumberOfKeptEvents() <<
    ", that may be\n  kept by the vis manager." <<
    "\n  The number of events kept by the vis manager can be changed with"
    "\n  \"/vis/scene/endOfEventAction accumulate <N>\", where N is the"
    "\n  maximum number you wish to allow.  N < 0 means \"unlimited\"."
    << G4endl;
  }

  // Don't call IsValidView unless there is a scene handler.  This
  // avoids WARNING message at end of event and run when the user has
  // not instantiated a scene handler, e.g., in batch mode.
  G4bool valid = fpSceneHandler && IsValidView();
  if (GetConcreteInstance() && valid) {
//    // ???? I can't remember why
//    // if (!fpSceneHandler->GetMarkForClearingTransientStore()) {
//    // is here.  It prevents ShowView at end of run, which seems to be OK
//    // for sequential mode, but MT mode seems to need it (I have not
//    // figured out why). ???? JA ????
//    if (!fpSceneHandler->GetMarkForClearingTransientStore()) {
      if (fpScene->GetRefreshAtEndOfRun()) {
	fpSceneHandler->DrawEndOfRunModels();
        // An extra refresh for auto-refresh viewers.
        // ???? I DON'T WHY THIS IS NECESSARY ???? JA ????
        if (fpViewer->GetViewParameters().IsAutoRefresh()) {
          fpViewer->RefreshView();
        }
        // ShowView guarantees the view is flushed to the screen.  It also
        // triggers other features such picking (if enabled) and allows
        // file-writing viewers to close the file.
        fpViewer->ShowView();
	fpSceneHandler->SetMarkForClearingTransientStore(true);
      } else {
        if (fpGraphicsSystem->GetFunctionality() ==
            G4VGraphicsSystem::fileWriter) {
          if (fVerbosity >= warnings) {
            G4cout << "\"/vis/viewer/update\" to close file." << G4endl;
          }
        }
      }
//    }
  }
  fEventRefreshing = false;
}

void G4VisManager::ClearTransientStoreIfMarked(){
  // Assumes valid view.
  if (fpSceneHandler->GetMarkForClearingTransientStore()) {
    fpSceneHandler->SetMarkForClearingTransientStore(false);
    fpSceneHandler->ClearTransientStore();
  }
  // Record if transients drawn.  These local flags are only set
  // *after* ClearTransientStore.  In the code in G4VSceneHandler
  // triggered by ClearTransientStore, use these flags so that
  // event refreshing is not done too early.
  fTransientsDrawnThisEvent = fpSceneHandler->GetTransientsDrawnThisEvent();
  fTransientsDrawnThisRun = fpSceneHandler->GetTransientsDrawnThisRun();
}

void G4VisManager::ResetTransientsDrawnFlags()
{
  fTransientsDrawnThisRun = false;
  fTransientsDrawnThisEvent = false;
  G4SceneHandlerListConstIterator i;
  for (i = fAvailableSceneHandlers.begin();
       i != fAvailableSceneHandlers.end(); ++i) {
    (*i)->SetTransientsDrawnThisEvent(false);
    (*i)->SetTransientsDrawnThisRun(false);
  }
}

G4String G4VisManager::ViewerShortName (const G4String& viewerName) const {
  G4String viewerShortName (viewerName);
  viewerShortName = viewerShortName (0, viewerShortName.find (' '));
  return viewerShortName.strip ();
}

G4VViewer* G4VisManager::GetViewer (const G4String& viewerName) const {
  G4String viewerShortName = ViewerShortName (viewerName);
  size_t nHandlers = fAvailableSceneHandlers.size ();
  size_t iHandler, iViewer;
  G4VViewer* viewer = 0;
  G4bool found = false;
  for (iHandler = 0; iHandler < nHandlers; iHandler++) {
    G4VSceneHandler* sceneHandler = fAvailableSceneHandlers [iHandler];
    const G4ViewerList& viewerList = sceneHandler -> GetViewerList ();
    for (iViewer = 0; iViewer < viewerList.size (); iViewer++) {
      viewer = viewerList [iViewer];
      if (viewerShortName == viewer -> GetShortName ()) {
	found = true;
	break;
      }
    }
    if (found) break;
  }
  if (found) return viewer;
  else return 0;
}

std::vector<G4String> G4VisManager::VerbosityGuidanceStrings;

G4String G4VisManager::VerbosityString(Verbosity verbosity) {
  G4String rs;
  switch (verbosity) {
  case         quiet: rs = "quiet (0)"; break;
  case       startup: rs = "startup (1)"; break;
  case        errors: rs = "errors (2)"; break;
  case      warnings: rs = "warnings (3)"; break;
  case confirmations: rs = "confirmations (4)"; break;
  case    parameters: rs = "parameters (5)"; break;
  case           all: rs = "all (6)"; break;
  }
  return rs;
}

G4VisManager::Verbosity
G4VisManager::GetVerbosityValue(const G4String& verbosityString) {
  G4String ss(verbosityString); ss.toLower();
  Verbosity verbosity;
  if      (ss(0) == 'q') verbosity = quiet;
  else if (ss(0) == 's') verbosity = startup;
  else if (ss(0) == 'e') verbosity = errors;
  else if (ss(0) == 'w') verbosity = warnings;
  else if (ss(0) == 'c') verbosity = confirmations;
  else if (ss(0) == 'p') verbosity = parameters;
  else if (ss(0) == 'a') verbosity = all;
  else {
    G4int intVerbosity;
    std::istringstream is(ss);
    is >> intVerbosity;
    if (!is) {
      G4cerr << "ERROR: G4VisManager::GetVerbosityValue: invalid verbosity \""
	     << verbosityString << "\"";
      for (size_t i = 0; i < VerbosityGuidanceStrings.size(); ++i) {
	G4cerr << '\n' << VerbosityGuidanceStrings[i];
      }
      verbosity = warnings;
      G4cerr << "\n  Returning " << VerbosityString(verbosity)
	     << G4endl;
    }
    else {
      verbosity = GetVerbosityValue(intVerbosity);
    }
  }
  return verbosity;
}

G4VisManager::Verbosity G4VisManager::GetVerbosityValue(G4int intVerbosity) {
  Verbosity verbosity;
  if      (intVerbosity < quiet) verbosity = quiet;
  else if (intVerbosity > all)   verbosity = all;
  else                           verbosity = Verbosity(intVerbosity);
  return verbosity;
}

G4VisManager::Verbosity G4VisManager::GetVerbosity () {
  return fVerbosity;
}

void G4VisManager::SetVerboseLevel (G4int intVerbosity) {
  fVerbosity = GetVerbosityValue(intVerbosity);
}

void G4VisManager::SetVerboseLevel (const G4String& verbosityString) {
  fVerbosity = GetVerbosityValue(verbosityString);
}

G4bool G4VisManager::IsValidView () {

  if (!fInitialised) Initialise ();

  static G4bool noGSPrinting = true;
  if (!fpGraphicsSystem) {
    // Limit printing - we do not want printing if the user simply does
    // not want to use graphics, e.g., in batch mode.
    if (noGSPrinting) {
      noGSPrinting = false;
      if (fVerbosity >= warnings) {
	G4cout <<
  "WARNING: G4VisManager::IsValidView(): Attempt to draw when no graphics system"
  "\n  has been instantiated.  Use \"/vis/open\" or \"/vis/sceneHandler/create\"."
  "\n  Alternatively, to avoid this message, suppress instantiation of vis"
  "\n  manager (G4VisExecutive), possibly by setting G4VIS_NONE, and ensure"
  "\n  drawing code is executed only if G4VVisManager::GetConcreteInstance()"
  "\n  is non-zero."
	       << G4endl;
      }
    }
    return false;
  }

  if ((!fpScene) || (!fpSceneHandler) || (!fpViewer)) {
    if (fVerbosity >= errors) {
      G4cerr <<
	"ERROR: G4VisManager::IsValidView(): Current view is not valid."
	     << G4endl;
      PrintInvalidPointers ();
    }
    return false;
  }

  if (fpScene != fpSceneHandler -> GetScene ()) {
    if (fVerbosity >= errors) {
      G4cerr << "ERROR: G4VisManager::IsValidView ():";
      if (fpSceneHandler -> GetScene ()) {
	G4cout <<
	  "\n  The current scene \""
	       << fpScene -> GetName ()
	       << "\" is not handled by"
	  "\n  the current scene handler \""
	       << fpSceneHandler -> GetName ()
	       << "\""
	  "\n  (it currently handles scene \""
	       << fpSceneHandler -> GetScene () -> GetName ()
	       << "\")."
	  "\n  Either:"
	  "\n  (a) attach it to the scene handler with"
	  "\n      /vis/sceneHandler/attach "
	       << fpScene -> GetName ()
	       <<	", or"
	  "\n  (b) create a new scene handler with "
	  "\n      /vis/sceneHandler/create <graphics-system>,"
	  "\n      in which case it should pick up the the new scene."
	       << G4endl;
      }
      else {
	G4cout << "\n  Scene handler \""
	       << fpSceneHandler -> GetName ()
	       << "\" has null scene pointer."
	  "\n  Attach a scene with /vis/sceneHandler/attach [<scene-name>]"
	       << G4endl;
      }
    }
    return false;
  }

  const G4ViewerList& viewerList = fpSceneHandler -> GetViewerList ();
  if (viewerList.size () == 0) {
    if (fVerbosity >= errors) {
      G4cerr <<
	"ERROR: G4VisManager::IsValidView (): the current scene handler\n  \""
	     << fpSceneHandler -> GetName ()
	     << "\" has no viewers.  Do /vis/viewer/create."
	     << G4endl;
    }
    return false;
  }

  G4bool isValid = true;
  if (fpScene -> IsEmpty ()) {  // Add world by default if possible...
    G4bool warn(fVerbosity >= warnings);
    G4bool successful = fpScene -> AddWorldIfEmpty (warn);
    if (!successful || fpScene -> IsEmpty ()) {        // If still empty...
      if (fVerbosity >= errors) {
	G4cerr << "ERROR: G4VisManager::IsValidView ():";
	G4cerr <<
	  "\n  Attempt at some drawing operation when scene is empty."
	  "\n  Maybe the geometry has not yet been defined."
	  "  Try /run/initialize."
          "\n  Or use \"/vis/scene/add/extent\"."
	       << G4endl;
      }
      isValid = false;
    }
    else {
      G4UImanager::GetUIpointer()->ApplyCommand ("/vis/scene/notifyHandlers");
      if (fVerbosity >= warnings) {
	G4cout <<
	  "WARNING: G4VisManager: the scene was empty, \"world\" has been"
	  "\n  added and the scene handlers notified.";
	G4cout << G4endl;
      }
    }
  }
  return isValid;
}

void
G4VisManager::RegisterModelFactories() 
{
  if (fVerbosity >= warnings) {
    G4cout<<"G4VisManager: No model factories registered with G4VisManager."<<G4endl;
    G4cout<<"G4VisManager::RegisterModelFactories() should be overridden in derived"<<G4endl;
    G4cout<<"class. See G4VisExecutive for an example."<<G4endl;
  }
}

#ifdef G4MULTITHREADED
void G4VisManager::SetUpForAThread()
{
  new G4VisStateDependent(this); 
}
#endif

void G4VisManager::IgnoreStateChanges(G4bool val)
{
  fIgnoreStateChanges = val; 
}

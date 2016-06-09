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
// $Id: G4VisManager.cc,v 1.98 2006/06/29 21:30:00 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// GEANT4 Visualization Manager - John Allison 02/Jan/1996.

#include "G4VisManager.hh"

#include "G4VisCommands.hh"
#include "G4VisCommandsCompound.hh"
#include "G4VisCommandsGeometry.hh"
#include "G4VisCommandsGeometrySet.hh"
#include "G4VisCommandsScene.hh"
#include "G4VisCommandsSceneAdd.hh"
#include "G4VisCommandsSceneHandler.hh"
#include "G4VisCommandsViewer.hh"
#include "G4VisCommandsViewerSet.hh"
#include "G4UImanager.hh"
#include "G4VisStateDependent.hh"
#include "G4UIdirectory.hh"
#include "G4VisFeaturesOfFukuiRenderer.hh"
#include "G4VisFeaturesOfDAWNFILE.hh"
#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VisFeaturesOfOpenInventor.hh"
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
#include "G4NURBS.hh"
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
#include "G4Run.hh"
#include "G4Event.hh"

#include <sstream>

G4VisManager* G4VisManager::fpInstance = 0;

G4VisManager::G4VisManager ():
  fInitialised     (false),
  fpUserVisAction  (0),
  fpGraphicsSystem (0),
  fpScene          (0),
  fpSceneHandler   (0),
  fpViewer         (0),
  fVerbosity       (warnings),
  fVerbose         (1),
  fpStateDependent (0),
  fEventCount      (0),
  fReprocessing (false),
  fReprocessingLastEvent (false),
  fLastRunID       (0),
  fLastEventID     (0),
  fTransientsDrawnThisRun   (false),
  fTransientsDrawnThisEvent (false)
  // All other objects use default constructors.
{
  fpTrajDrawModelMgr = new G4VisModelManager<G4VTrajectoryModel>("/vis/modeling/trajectories");
  fpTrajFilterMgr = new G4VisFilterManager<G4VTrajectory>("/vis/filtering/trajectories");

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
      ("G4VisManager: attempt to Construct more than one VisManager.");
  }
  else {

    fpInstance = this;
    SetConcreteInstance(this);

    fpStateDependent = new G4VisStateDependent (this);
    // No need to delete this; G4StateManager does this.

    if (fVerbosity >= startup) {
      G4cout << "Visualization Manager instantiating..." << G4endl;
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
  }
}

G4VisManager::~G4VisManager () {
  fpInstance = 0;
  size_t i;
  for (i = 0; i < fSceneList.size (); ++i) {
    delete fSceneList[i];
  }
  for (i = 0; i < fAvailableSceneHandlers.size (); ++i) {
    delete fAvailableSceneHandlers[i];
  }
  for (i = 0; i < fAvailableGraphicsSystems.size (); ++i) {
    delete fAvailableGraphicsSystems[i];
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

  delete fpTrajDrawModelMgr;
  delete fpTrajFilterMgr;
}

G4VisManager* G4VisManager::GetInstance () {
  if (!fpInstance) {
    G4Exception
      ("G4VisManager::GetInstance: VisManager not yet instantiated!");
  }
  return fpInstance;
}

void G4VisManager::Initialise () {

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
      "\n  See visualization/include/G4VisExecutive.hh/icc, for example."
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
    PrintAvailableGraphicsSystems ();
    G4cout << G4endl;
  }

  // Make top level directory...
  G4UIcommand* directory;
  directory = new G4UIdirectory ("/vis/");
  directory -> SetGuidance ("Visualization commands.");
  fDirectoryList.push_back (directory);

  // ... and make command directory for commands instantiated in the
  // modeling subcategory...
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

  RegisterMessengers ();

  if (fVerbosity >= startup) {
    G4cout << "Registering model factories..." << G4endl;
  }

  RegisterModelFactories();

  if (fVerbosity >= startup) {
    G4cout <<
      "\nYou have successfully registered the following model factories."
	 << G4endl;
    PrintAvailableModels ();
    G4cout << G4endl;
  }

  fInitialised = true;
}

void G4VisManager::Enable() {
  if (IsValidView ()) {
    if (fVerbosity >= confirmations) {
      G4cout << "G4VisManager::Enable: visualization enabled." << G4endl;
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
    // Use G4TrajectoryDrawByCharge as a default.
    fpTrajDrawModelMgr->Register(new G4TrajectoryDrawByCharge("AutoGenerated"));

    if (fVerbosity >= warnings) {
      G4cout<<"G4VisManager: Using G4TrajectoryDrawByCharge as default trajectory model."<<G4endl;
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

void G4VisManager::SelectTrajectoryModel(const G4String& model) 
{
   fpTrajDrawModelMgr->SetCurrent(model);
}

void G4VisManager::Draw (const G4Circle& circle,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (circle);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4NURBS& nurbs,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (nurbs);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Polyhedron& polyhedron,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (polyhedron);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Polyline& line,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (line);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Polymarker& polymarker,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (polymarker);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Scale& scale,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (scale);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Square& square,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (square);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Text& text,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (text);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw2D (const G4Text& text)
{
  if (IsValidView()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> BeginPrimitives2D();
    fpSceneHandler -> AddPrimitive(text);
    fpSceneHandler -> EndPrimitives2D();
  }
}

void G4VisManager::Draw (const G4VHit& hit) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> AddCompound (hit);
  }
}

void G4VisManager::Draw (const G4VTrajectory& traj,
			 G4int i_mode) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> SetModel (&dummyTrajectoriesModel);
    dummyTrajectoriesModel.SetDrawingMode(i_mode);
    fpSceneHandler -> AddCompound (traj);
  }
}

void G4VisManager::Draw (const G4LogicalVolume& logicalVol,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
  // Find corresponding solid.
  G4VSolid* pSol = logicalVol.GetSolid ();
  Draw (*pSol, attribs, objectTransform);
}

void G4VisManager::Draw (const G4VSolid& solid,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    fpSceneHandler -> PreAddSolid (objectTransform, attribs);
    solid.DescribeYourselfTo (*fpSceneHandler);
    fpSceneHandler -> PostAddSolid ();
  }
}

void G4VisManager::Draw (const G4VPhysicalVolume& physicalVol,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
  // Find corresponding logical volume and solid.
  G4LogicalVolume* pLV  = physicalVol.GetLogicalVolume ();
  G4VSolid*        pSol = pLV -> GetSolid ();
  Draw (*pSol, attribs, objectTransform);
}

void G4VisManager::CreateSceneHandler (G4String name) {
  if (!fInitialised) Initialise ();
  if (fpGraphicsSystem) {
    G4VSceneHandler* pSceneHandler =
      fpGraphicsSystem -> CreateSceneHandler (name);
    if (pSceneHandler) {
      fAvailableSceneHandlers.push_back (pSceneHandler);
      fpSceneHandler = pSceneHandler;                         // Make current.
    }
    else {
      if(fVerbosity >= errors) {
	G4cout << "ERROR in G4VisManager::CreateSceneHandler during "
	       << fpGraphicsSystem -> GetName ()
	       << " scene creation.\n  No action taken."
	       << G4endl;
      }
    }
  }
  else PrintInvalidPointers ();
}

void G4VisManager::CreateViewer (G4String name) {

  if (!fInitialised) Initialise ();

  if (fpSceneHandler) {
    G4VViewer* p = fpGraphicsSystem -> CreateViewer (*fpSceneHandler, name);
    if (p) {
      fpViewer = p;                             // Make current.
      fpViewer -> Initialise ();
      fpSceneHandler -> AddViewerToList (fpViewer);
      fpSceneHandler -> SetCurrentViewer (fpViewer);

      const G4ViewParameters& vp = fpViewer->GetViewParameters();
      G4bool warn = false;
      if (vp.IsCulling () && vp.IsCullingInvisible ()) {
	warn = true;
	if (fVerbosity >= confirmations) {
	  G4cout << "G4VisManager::CreateViewer: new viewer created:"
		 << G4endl;
	}
	if (fVerbosity >= parameters) {
	  G4cout << " view parameters are:\n  " << vp << G4endl;
	}
	if (fVerbosity >= warnings) {
	  G4cout <<
	    "WARNING: objects with visibility flag set to \"false\""
	    " will not be drawn!"
	    "\n  \"/vis/viewer/set/culling global false\" to Draw such objects."
		 << G4endl;
	}
      }
      if (vp.IsCullingCovered ()) {
	if (!warn) {
	  if (fVerbosity >= confirmations) {
	    G4cout << "G4VisManager::CreateViewer: new viewer created:"
		   << G4endl;
	  }
	}
	warn = true;
	if (fVerbosity >= warnings) {
	  G4cout <<
	    "WARNING: covered objects in solid mode will not be rendered!"
	    "\n  \"/vis/viewer/set/culling coveredDaughters false\" to reverse this."
		 << G4endl;
	}
      }
      if (warn) {
	if (fVerbosity >= warnings) {
	  G4cout << "  Also see other \"/vis/viewer/set\" commands."
		 << G4endl;
	}
      }
    }
    else {
      if (fVerbosity >= errors) {
	G4cout << "ERROR in G4VisManager::CreateViewer during "
	       << fpGraphicsSystem -> GetName ()
	       <<	" viewer creation.\n  No action taken."
	       << G4endl;
      }
    }
  }
  else PrintInvalidPointers ();
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
    std::vector<G4VModel*>& modelList = pScene -> SetRunDurationModelList ();

    if (modelList.size ()) {
      G4bool modelInvalid;
      do {  // Remove, if required, one at a time.
	modelInvalid = false;
	std::vector<G4VModel*>::iterator iterModel;
	for (iterModel = modelList.begin();
	     iterModel != modelList.end();
	     ++iterModel) {
	  modelInvalid = !((*iterModel) -> Validate (fVerbosity >= warnings));
	  if (modelInvalid) {
	    // Model invalid - remove and break.
	    if (fVerbosity >= warnings) {
	      G4cout << "WARNING: Model \""
		     << (*iterModel) -> GetGlobalDescription ()
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

G4bool G4VisManager::FilterTrajectory(const G4VTrajectory& trajectory)
{
  return fpTrajFilterMgr->Accept(trajectory);
}   

void G4VisManager::DispatchToModel(const G4VTrajectory& trajectory, G4int i_mode)
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

  const G4VTrajectoryModel* model = CurrentTrajDrawModel();

  assert (0 != model); // Should exist

  model->Draw(trajectory, i_mode, visible);
} 

void G4VisManager::SetUserAction
(G4VUserVisAction* pVisAction,
 const G4VisExtent& extent) {
  fpUserVisAction = pVisAction;
  fUserVisActionExtent = extent;
  if (extent.GetExtentRadius() <= 0.) {
    if (fVerbosity >= warnings) {
      G4cout << 
	"WARNING: No extent set for user vis action.  (You may"
	"\n  set it later when adding with /vis/scene/add/userAction.)"
	     << G4endl;
    }
  }
}

void G4VisManager::SetCurrentScene (G4Scene* pScene) {
  if (pScene != fpScene) {
    // A change of scene.  Therefore reset transients drawn flags.  All
    // memory of previous transient proceessing thereby erased...
    fTransientsDrawnThisRun = false;
    if (fpSceneHandler) fpSceneHandler->SetTransientsDrawnThisRun(false);
    fTransientsDrawnThisEvent = false;
    if (fpSceneHandler) fpSceneHandler->SetTransientsDrawnThisEvent(false);
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
    IsValidView ();  // Checks.
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
  fpSceneHandler -> SetCurrentViewer (pViewer);
  fpScene = fpSceneHandler -> GetScene ();
  fpGraphicsSystem = fpSceneHandler -> GetGraphicsSystem ();
  IsValidView ();  // Checks.
}

void G4VisManager::RegisterMessengers () {

  // Instantiate individual messengers/commands (often - but not
  // always - one command per messenger).

  G4VVisCommand::SetVisManager (this);  // Sets shared pointer to vis manager.

  G4UIcommand* directory;

  // Top level commands...
  RegisterMessenger(new G4VisCommandEnable);
  RegisterMessenger(new G4VisCommandList);
  RegisterMessenger(new G4VisCommandVerbose);

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
  RegisterMessenger(new G4VisCommandGeometrySetForceSolid);
  RegisterMessenger(new G4VisCommandGeometrySetForceWireframe);
  RegisterMessenger(new G4VisCommandGeometrySetVisibility);

  directory = new G4UIdirectory ("/vis/scene/");
  directory -> SetGuidance ("Operations on Geant4 scenes.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandSceneCreate);
  RegisterMessenger(new G4VisCommandSceneEndOfEventAction);
  RegisterMessenger(new G4VisCommandSceneEndOfRunAction);
  RegisterMessenger(new G4VisCommandSceneList);
  RegisterMessenger(new G4VisCommandSceneNotifyHandlers);
  RegisterMessenger(new G4VisCommandSceneSelect);
  RegisterMessenger(new G4VisCommandSceneTransientsAction);

  directory = new G4UIdirectory ("/vis/scene/add/");
  directory -> SetGuidance ("Add model to current scene.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandSceneAddAxes);
  RegisterMessenger(new G4VisCommandSceneAddEventID);
  RegisterMessenger(new G4VisCommandSceneAddGhosts);
  RegisterMessenger(new G4VisCommandSceneAddHits);
  RegisterMessenger(new G4VisCommandSceneAddLogicalVolume);
  RegisterMessenger(new G4VisCommandSceneAddLogo);
  RegisterMessenger(new G4VisCommandSceneAddScale);
  RegisterMessenger(new G4VisCommandSceneAddText);
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

  directory = new G4UIdirectory ("/vis/viewer/");
  directory -> SetGuidance ("Operations on Geant4 viewers.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandViewerClear);
  RegisterMessenger(new G4VisCommandViewerCreate);
  RegisterMessenger(new G4VisCommandViewerDolly);
  RegisterMessenger(new G4VisCommandViewerFlush);
  RegisterMessenger(new G4VisCommandViewerList);
  RegisterMessenger(new G4VisCommandViewerPan);
  RegisterMessenger(new G4VisCommandViewerRebuild);
  RegisterMessenger(new G4VisCommandViewerRefresh);
  RegisterMessenger(new G4VisCommandViewerReset);
  RegisterMessenger(new G4VisCommandViewerScale);
  RegisterMessenger(new G4VisCommandViewerSelect);
  RegisterMessenger(new G4VisCommandViewerUpdate);
  RegisterMessenger(new G4VisCommandViewerZoom);

  directory = new G4UIdirectory ("/vis/viewer/set/");
  directory -> SetGuidance ("Set view parameters of current viewer.");
  fDirectoryList.push_back (directory);
  RegisterMessenger(new G4VisCommandsViewerSet);

  // List manager commands
  RegisterMessenger(new G4VisCommandListManagerList< G4VisModelManager<G4VTrajectoryModel> >
		    (fpTrajDrawModelMgr, fpTrajDrawModelMgr->Placement()));
  RegisterMessenger(new G4VisCommandListManagerSelect< G4VisModelManager<G4VTrajectoryModel> >
		    (fpTrajDrawModelMgr, fpTrajDrawModelMgr->Placement()));  

  // Filter manager commands
  RegisterMessenger(new G4VisCommandListManagerList< G4VisFilterManager<G4VTrajectory> >
                    (fpTrajFilterMgr, fpTrajFilterMgr->Placement()));
  RegisterMessenger(new G4VisCommandManagerMode< G4VisFilterManager<G4VTrajectory> >
                    (fpTrajFilterMgr, fpTrajFilterMgr->Placement()));
}

void G4VisManager::PrintAvailableGraphicsSystems () const {
  G4int nSystems = fAvailableGraphicsSystems.size ();
  G4cout << "Current available graphics systems are:";
  if (nSystems) {
    for (int i = 0; i < nSystems; i++) {
      const G4VGraphicsSystem* pSystem = fAvailableGraphicsSystems [i];
      G4cout << "\n  " << pSystem -> GetName ();
      if (pSystem -> GetNickname () != "") {
	G4cout << " (" << pSystem -> GetNickname () << ")";
      }
    }
  }
  else {
    G4cout << "\n  NONE!!!  None registered - yet!  Mmmmm!";
  }
  G4cout << G4endl;
}

void G4VisManager::PrintAvailableModels () const
{
  fpTrajDrawModelMgr->Print(G4cout);
  G4cout << G4endl;
  fpTrajFilterMgr->Print(G4cout);
}

void G4VisManager::PrintInvalidPointers () const {
  if (fVerbosity >= errors) {
    G4cout << "ERROR: G4VisManager::PrintInvalidPointers:";
    if (!fpGraphicsSystem) {
      G4cout << "\n null graphics system pointer.";
    }
    else {
      G4cout << "\n  Graphics system is " << fpGraphicsSystem -> GetName ()
	     << " but:";
      if (!fpScene)
	G4cout <<
	  "\n  Null scene pointer. Use \"/vis/drawVolume\" or"
	  " \"/vis/scene/create\".";
      if (!fpSceneHandler)
	G4cout <<
	  "\n  Null scene handler pointer. Use \"/vis/open\" or"
	  " \"/vis/sceneHandler/create\".";
      if (!fpViewer )
	G4cout <<
	  "\n  Null viewer pointer. Use \"/vis/viewer/create\".";
    }
    G4cout << G4endl;
  }
}

/************* Event copy stuff *******************************
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

G4String fLastEventRandomStatus;
G4Event* fLastEvent;
std::vector<G4String> fEventsRandomStatus;
std::vector<G4Event*> fEvents;

G4PrimaryParticle* CopyParticles(G4PrimaryParticle* pp)
{
  G4PrimaryParticle* newPP = 0;
  G4PrimaryParticle* previousPP = 0;
  while (pp) {
    G4PrimaryParticle* tmpPP = new G4PrimaryParticle(*pp);
    tmpPP->ClearNextParticlePointer();
    G4PrimaryParticle* daughters = pp->GetDaughter();
    G4PrimaryParticle* newDaughters = CopyParticles(daughters);
    tmpPP->SetDaughter(newDaughters);
    if (!newPP) newPP = tmpPP;
    else previousPP->SetNext(tmpPP);
    previousPP = tmpPP;
    pp = pp->GetNext();
  }
  // Set user info = 0???
  return newPP;
}

// Not needed?
void DeleteParticles(G4PrimaryParticle* pp)
{
  while (pp) {
    G4PrimaryParticle* daughters = pp->GetDaughter();
    DeleteParticles(daughters);
    G4PrimaryParticle* oldPP = pp;
    pp = pp->GetNext();
    delete oldPP;
  }
}

G4PrimaryVertex* CopyVertices(G4PrimaryVertex* pv)
{
  G4PrimaryVertex* newPV = 0;
  G4PrimaryVertex* previousPV = 0;
  while (pv) {
    G4PrimaryVertex* tmpPV = new G4PrimaryVertex(*pv);
    tmpPV->ClearNextVertexPointer();
    G4PrimaryParticle* primaryParticle = pv->GetPrimary();
    G4PrimaryParticle* newPrimaryParticle = CopyParticles(primaryParticle);
    tmpPV->SetPrimary(newPrimaryParticle);
    if (!newPV) newPV = tmpPV;
    else previousPV->SetNext(tmpPV);
    previousPV = tmpPV;
    pv = pv->GetNext();
  }
  // Set user info = 0???
  return newPV;
}

// Not needed?
void DeleteVertices(G4PrimaryVertex* pv)
{
  while (pv) {
    G4PrimaryParticle* primaryParticle = pv->GetPrimary();
    DeleteParticles(primaryParticle);
    G4PrimaryVertex* oldPV = pv;
    pv = pv->GetNext();
    delete oldPV;
  }
}

G4Event* CreateDeepCopy(const G4Event* event)
{
  G4Event* newEvent = new G4Event(*event);
  G4PrimaryVertex* primaryVertex = event->GetPrimaryVertex();
  G4PrimaryVertex* newPrimaryVertex = CopyVertices(primaryVertex);
  newEvent->AddPrimaryVertex(newPrimaryVertex);
  // Set user info = 0???
  return newEvent;
}

// Not needed?
void DeleteEvent(G4Event* event)
{
  if (event) {
    G4PrimaryVertex* primaryVertex = event->GetPrimaryVertex();
    DeleteVertices(primaryVertex);
    delete event;
  }
}
************ Event copy stuff ***********************************/

void G4VisManager::BeginOfRun ()
{
  //G4cout << "G4VisManager::BeginOfRun" << G4endl;
  fEventCount = 0;
  fTransientsDrawnThisRun = false;
  if (fpSceneHandler) fpSceneHandler->SetTransientsDrawnThisRun(false);
/************* Event copy stuff *******************************
  if (!fReprocessing) {
    fEventsRandomStatus.clear();
    std::vector<G4Event*>::iterator i;
    //for (i = fEvents.begin(); i != fEvents.end(); ++i) DeleteEvent(*i);
    for (i = fEvents.begin(); i != fEvents.end(); ++i) delete *i;
    fEvents.clear();
  }
************ Event copy stuff ***********************************/
}

void G4VisManager::BeginOfEvent ()
{
  //G4cout << "G4VisManager::BeginOfEvent" << G4endl;
  G4RunManager* runManager = G4RunManager::GetRunManager();
  if (runManager) {
    if (!fEventCount) {  // First event.  Do run stuff here because,
			 // curiously, currentRun is still zero in
			 // BeginOfRun.
      fBeginOfLastRunRandomStatus =
        runManager->GetRandomNumberStatusForThisRun();
      /*************************************************
      std::ostringstream oss;
      CLHEP::HepRandom::saveFullState(oss);
      fBeginOfLastRunRandomStatus = oss.str();
      ************************************************/
      const G4Run* currentRun = runManager->GetCurrentRun();
      if (currentRun) fLastRunID = currentRun->GetRunID();
      else fLastRunID++;
    }
    fBeginOfLastEventRandomStatus =
      runManager->GetRandomNumberStatusForThisEvent();
    /*************************************************
    std::ostringstream oss;
    CLHEP::HepRandom::saveFullState(oss);
    fBeginOfLastEventRandomStatus = oss.str();
    ************************************************/
    G4Event* currentEvent =
      G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
    if (currentEvent) {
      if (fReprocessingLastEvent) currentEvent->SetEventID(fLastEventID);
      else fLastEventID = currentEvent->GetEventID();
    }

    /************* Event copy stuff *******************************
    // If not triggered by transient re-computation in
    // G4VSceneHandler::ProcessScene (and re-computation requested)...
    if (!fReprocessing && fpScene && fpScene->GetRecomputeTransients()) {
      const G4Event* event =
	G4EventManager::GetEventManager()->GetConstCurrentEvent();
      std::ostringstream oss;
      CLHEP::HepRandom::saveFullState(oss);
      fLastEventRandomStatus = oss.str();
      fLastEvent = CreateDeepCopy(event);
      if (!fpScene->GetRefreshAtEndOfEvent()) {
	fEventsRandomStatus.push_back(fLastEventRandomStatus);
	fEvents.push_back(fLastEvent);
      }
    }
    ************ Event copy stuff ***********************************/

  }
  fTransientsDrawnThisEvent = false;
  if (fpSceneHandler) fpSceneHandler->SetTransientsDrawnThisEvent(false);
}

void G4VisManager::EndOfEvent ()
{
  //G4cout << "G4VisManager::EndOfEvent" << G4endl;
  ++fEventCount;

  // Don't call IsValidView unless there is a scene handler.  This
  // avoids WARNING message at end of event and run when the user has
  // not instantiated a scene handler, e.g., in batch mode.
  if (GetConcreteInstance() && fpSceneHandler && IsValidView()) {
    const std::vector<G4VModel*>& EOEModelList =
      fpScene -> GetEndOfEventModelList ();
    size_t nModels = EOEModelList.size();
    if (nModels) {
      ClearTransientStoreIfMarked();
      for (size_t i = 0; i < nModels; i++) {
	G4VModel* pModel = EOEModelList [i];
	fpSceneHandler -> SetModel (pModel);
	pModel -> DescribeYourselfTo (*fpSceneHandler);
      }
      fpSceneHandler -> SetModel (0);
    }
    if (fpScene->GetRefreshAtEndOfEvent()) {
      G4int nEvents = 0;
      G4int eventID = -2;  // (If no run manager, triggers ShowView as normal.)
      G4RunManager* runManager = G4RunManager::GetRunManager();
      if (runManager) {
	const G4Run* currentRun = runManager->GetCurrentRun();
	const G4Event* currentEvent = runManager->GetCurrentEvent();
	if (currentRun && currentEvent) {
	  nEvents = currentRun->GetNumberOfEventToBeProcessed();
	  eventID = currentEvent->GetEventID();
	}
      }
      // Unless last event (in which case wait end of run)...
      if (eventID < nEvents - 1) {
	fpViewer->ShowView();
	fpSceneHandler->SetMarkForClearingTransientStore(true);
      }
    }
  }
}

void G4VisManager::EndOfRun ()
{
  //G4cout << "G4VisManager::EndOfRun" << G4endl;
  // Don't call IsValidView unless there is a scene handler.  This
  // avoids WARNING message at end of event and run when the user has
  // not instantiated a scene handler, e.g., in batch mode.
  if (GetConcreteInstance() && fpSceneHandler && IsValidView()) {
    if (!fpSceneHandler->GetMarkForClearingTransientStore()) {
      if (fpScene->GetRefreshAtEndOfRun()) {
	fpViewer->ShowView();
	fpSceneHandler->SetMarkForClearingTransientStore(true);
      }
    }
  }
  fReprocessing = false;
  fReprocessingLastEvent = false;
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
  // transient reprocessing is not done too early.
  fTransientsDrawnThisEvent = fpSceneHandler->GetTransientsDrawnThisEvent();
  fTransientsDrawnThisRun = fpSceneHandler->GetTransientsDrawnThisRun();
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
  G4String s;
  switch (verbosity) {
  case         quiet: s = "quiet (0)"; break;
  case       startup: s = "startup (1)"; break;
  case        errors: s = "errors (2)"; break;
  case      warnings: s = "warnings (3)"; break;
  case confirmations: s = "confirmations (4)"; break;
  case    parameters: s = "parameters (5)"; break;
  case           all: s = "all (6)"; break;
  }
  return s;
}

G4VisManager::Verbosity
G4VisManager::GetVerbosityValue(const G4String& verbosityString) {
  G4String s(verbosityString); s.toLower();
  Verbosity verbosity;
  if      (s(0) == 'q') verbosity = quiet;
  else if (s(0) == 's') verbosity = startup;
  else if (s(0) == 'e') verbosity = errors;
  else if (s(0) == 'w') verbosity = warnings;
  else if (s(0) == 'c') verbosity = confirmations;
  else if (s(0) == 'p') verbosity = parameters;
  else if (s(0) == 'a') verbosity = all;
  else {
    G4int intVerbosity;
    std::istringstream is(s);
    is >> intVerbosity;
    if (!is) {
      G4cout << "ERROR: G4VisManager::GetVerbosityValue: invalid verbosity \""
	     << verbosityString << "\"";
      for (size_t i = 0; i < VerbosityGuidanceStrings.size(); ++i) {
	G4cout << '\n' << VerbosityGuidanceStrings[i];
      }
      verbosity = warnings;
      G4cout << "\n  Returning " << VerbosityString(verbosity)
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
      G4cout <<
	"ERROR: G4VisManager::IsValidView(): Current view is not valid."
	     << G4endl;
      PrintInvalidPointers ();
    }
    return false;
  }

  if (fpScene != fpSceneHandler -> GetScene ()) {
    if (fVerbosity >= errors) {
      G4cout << "ERROR: G4VisManager::IsValidView ():";
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
      G4cout <<
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
	G4cout << "ERROR: G4VisManager::IsViewValid ():";
	G4cout <<
	  "\n  Attempt at some drawing operation when scene is empty."
	  "\n  Maybe the geometry has not yet been defined."
	  "  Try /run/initialize."
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
  if (isValid) SetConcreteInstance(this);
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

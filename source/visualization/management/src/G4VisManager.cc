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
// $Id: G4VisManager.cc,v 1.48 2002-11-20 17:19:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager - John Allison 02/Jan/1996.

#include "G4VisManager.hh"

#include "G4VisCommands.hh"
#include "G4VisCommandsCompound.hh"
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

G4VisManager* G4VisManager::fpInstance = 0;

G4VisManager::G4VisManager ():
  fInitialised     (false),
  fpGraphicsSystem (0),
  fpScene          (0),
  fpSceneHandler   (0),
  fpViewer         (0),
  fVerbosity       (warnings),
  fVerbose         (1),
  fpStateDependent (0)   // All other objects use default constructors.
{
  if (fpInstance) {
    G4Exception
      ("G4VisManager: attempt to Construct more than one VisManager.");
  }
  else {

    fpInstance = this;

    fpStateDependent = new G4VisStateDependent (this);
    // No need to delete this; G4StateManager does this.

    if (fVerbosity >= startup) {
      G4cout << "Visualization Manager instantiating..." << G4endl;
    }

    // Note: The specific graphics systems must be instantiated in a
    // higher level library to avoid circular dependencies.  Also,
    // some specifically need additional external libararies that the
    // user must supply.  Therefore we ask the user to implement
    // RegisterGraphicsSystems() in a subclass.  We have to wait for
    // the subclass to instantiate so RegisterGraphicsSystems() cannot
    // be called from this constructor; it is called from
    // Initialise().  So we ask the user:
    //   (a) to write a subclass and implement
    //       RegisterGraphicsSystems().  See
    //       visualization/include/MyVisManager.hh/cc as an example.
    //   (b) instantiate the subclass.
    //   (c) invoke the Initialise() method of the subclass.
    // For example:
    //   ...
    // #ifdef G4VIS_USE
    //   // Instantiate and initialise Visualization Manager.
    //   G4VisManager* visManager = new MyVisManager;
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
    PrintAllGraphicsSystems ();
    PrintInstalledGraphicsSystems ();
    G4cout <<
      "\nYou have instantiated your own Visualization Manager, inheriting"
      "\n  G4VisManager and implementing RegisterGraphicsSystems(), in which"
      "\n  you should, normally, instantiate drivers which do not need"
      "\n  external packages or libraries, namely:"
      "\n    ASCIITree, DAWNFILE, HepRepFile, GAGTree, RayTracer, VRMLFILE"
      "\n  and, optionally, drivers under control of the following"
      "\n  environment variables:"
#ifdef G4VIS_BUILD_DAWN_DRIVER
      "\n    G4VIS_USE_DAWN"
#endif
#ifdef G4VIS_BUILD_OPACS_DRIVER
      "\n    G4VIS_USE_OPACS"
#endif
#ifdef G4VIS_BUILD_OPENGLX_DRIVER
      "\n    G4VIS_USE_OPENGLX"
#endif
#ifdef G4VIS_BUILD_OPENGLXM_DRIVER
      "\n    G4VIS_USE_OPENGLXM"
#endif
#ifdef G4VIS_BUILD_OIX_DRIVER
      "\n    G4VIS_USE_OIX"
#endif
#ifdef G4VIS_BUILD_OIWIN32_DRIVER
      "\n    G4VIS_USE_OIWIN32"
#endif
#ifdef G4VIS_BUILD_VRML_DRIVER
      "\n    G4VIS_USE_VRML"
#endif
      "\n  See visualization/include/MyVisManager.hh/cc, for example."
      "\n  In your main() you will have something like:"
      "\n  #ifdef G4VIS_USE"
      "\n    G4VisManager* visManager = new MyVisManager;"
      "\n    visManager -> SetVerboseLevel (Verbose);"
      "\n    visManager -> Initialize ();"
      "\n  #endif"
      "\n  (Don't forget to delete visManager;)"
	 << G4endl;
  }

  if (fVerbosity >= startup) {
    G4cout << "Registering graphics systems..." << G4endl;
  }

  RegisterGraphicsSystems ();

  if (fVerbosity >= confirmations) {
    G4cout <<
      "\nYou have successfully chosen to use the following graphics systems."
	 << G4endl;
    PrintAvailableGraphicsSystems ();
  }

  RegisterMessengers ();

  fInitialised = true;
}

// void G4VisManager::RegisterMessengers () - see separate file,
// G4VisManagerRegisterMessengers.cc.

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
  G4bool happy(true);
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
	     <<G4endl;
    }
    happy=false;
  }
  return happy;
}

void G4VisManager::Draw (const G4Circle& circle,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    CheckModel();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (circle);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4NURBS& nurbs,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    CheckModel();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (nurbs);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Polyhedron& polyhedron,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    CheckModel();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (polyhedron);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Polyline& line,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    CheckModel();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (line);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Polymarker& polymarker,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    CheckModel();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (polymarker);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Scale& scale,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    CheckModel();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (scale);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Square& square,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    CheckModel();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (square);
    fpSceneHandler -> EndPrimitives ();
  }
}

void G4VisManager::Draw (const G4Text& text,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    ClearTransientStoreIfMarked();
    CheckModel();
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (text);
    fpSceneHandler -> EndPrimitives ();
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
    CheckModel();
    fpSceneHandler -> PreAddThis (objectTransform, attribs);
    solid.DescribeYourselfTo (*fpSceneHandler);
    fpSceneHandler -> PostAddThis ();
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
      // Make it possible for user action code to Draw.
      SetConcreteInstance(this);

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

void G4VisManager::DeleteCurrentSceneHandler () {
  if (fVerbosity >= confirmations) {
    G4cout << "G4VisManager::DeleteCurrentSceneHandler: scene handler \""
	   << fpSceneHandler -> GetName ()
	   << "\"\n  and its viewer(s) are being deleted."
	   << G4endl;
  }
  if(fpSceneHandler) {
    fAvailableSceneHandlers.remove (fpSceneHandler);
    delete fpSceneHandler;
  }
  const G4SceneHandlerList& sceneHandlerList = fAvailableSceneHandlers;
  G4int nSH = sceneHandlerList.size ();
  G4int iSH;
  for (iSH = 0; iSH < nSH; iSH++) {
    if (sceneHandlerList [iSH] -> GetViewerList ().size ()) break;
  }
  if (iSH < nSH) {
    fpSceneHandler = sceneHandlerList [iSH];
    fpViewer = fpSceneHandler -> GetViewerList () [0];
    if (fVerbosity >= confirmations) {
      G4cout << "  Scene handler is now \""
	     << fpSceneHandler -> GetName ();
      G4cout << "\"\n  and viewer now \""
	     << fpViewer -> GetName ()
	     << "." << G4endl;
    }
    IsValidView ();  // Check.
  }
  else if (nSH) {
    fpSceneHandler = fAvailableSceneHandlers [0];
    fpViewer = 0;
    // Make it impossible for user action code to Draw.
    SetConcreteInstance(0);
    if (fVerbosity >= warnings) {
      G4cout << "WARNING: scene handler is now \""
	     << fpSceneHandler -> GetName ();
      G4cout << "\" but it has no viewers -\n  please create one."
	     << G4endl;
    }
  }
  else {
    fpSceneHandler = 0;
    fpViewer  = 0;
    // Make it impossible for user action code to Draw.
    SetConcreteInstance(0);
    if (fVerbosity >= warnings) {
      G4cout <<
	"WARNING: There are now no scene handlers left."
	"\n  /vis/sceneHandler/select to select a scene handler and"
	"\n  /vis/viewer/select to select a viewer,"
	"\n  or maybe you will have to create one."
	     << G4endl;
    }
  }
}

void G4VisManager::DeleteCurrentViewer () {
  if (fVerbosity >= confirmations) {
    G4cout << "G4VisManager::DeleteCurrentViewer: viewer \""
	   << fpViewer -> GetName () 
	   << "\" being deleted."
	   << G4endl;
  }
  if (fpViewer ) {
    fpViewer -> GetSceneHandler () -> RemoveViewerFromList (fpViewer);
    delete fpViewer;
  }
  const G4ViewerList& viewerList = fpSceneHandler -> GetViewerList ();
  if (viewerList.size () > 0) {
    fpViewer = viewerList [0];
    fpSceneHandler -> SetCurrentViewer (fpViewer);
    if (IsValidView ()) {
      if (fVerbosity >= confirmations) {
	G4cout << "  Viewer is now \""
	       << fpViewer -> GetName ()
	       << "\"." << G4endl;
      }
    }
  }
  else {
    fpViewer = 0;
    fpSceneHandler -> SetCurrentViewer (0);
    // Make it impossible for user action code to Draw.
    SetConcreteInstance(0);
    if (fVerbosity >= warnings) {
      G4cout <<
	"WARNING: There are now no viewers left for this scene handler."
	"\n  /vis/sceneHandler/select to select a different scene handler"
	"\n  and /vis/viewer/select to select a viewer, or maybe you will"
	"\n  have to create some."
	     << G4endl;
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
    G4std::vector<G4VModel*>& modelList = pScene -> SetRunDurationModelList ();

    if (modelList.size ()) {
      G4bool modelInvalid;
      do {  // Remove, if required, one at a time.
	modelInvalid = false;
	G4std::vector<G4VModel*>::iterator iterModel;
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

void G4VisManager::SetCurrentGraphicsSystemAndCreateViewer
(G4VGraphicsSystem* pSystem) {
  // This is used only by the deprecated
  // /vis~/create_view/new_graphics_system command.  Hence the ad hoc
  // creation of "scene-vis~" below as a working solution while these
  // old commands exist.
  fpGraphicsSystem = pSystem;
  if (!fpScene) {
    fpScene = new G4Scene ("scene-vis~");
    fSceneList.push_back (fpScene);
    if (fVerbosity >= confirmations) {
      G4cout << "G4VisManager::SetCurrentGraphicsSystemAndCreateViewer:"
	"\n  Empty scene \"scene-vis~\" created."
	     << G4endl;
    }
  }
  CreateSceneHandler ();
  CreateViewer ();
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
	// Make it impossible for user action code to Draw.
	SetConcreteInstance(0);
      }
    }
    else {
      fpSceneHandler = 0;
      fpViewer = 0;
      // Make it impossible for user action code to Draw.
      SetConcreteInstance(0);
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
    // Make it impossible for user action code to Draw.
    SetConcreteInstance(0);
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

  // Instantiate individual messengers/commands (often one command per
  // messenger).

  G4VVisCommand::SetVisManager (this);  // Sets shared pointer to vis manager.

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

void G4VisManager::PrintAllGraphicsSystems () const {
  G4cout <<
    "\nThe following graphics systems drivers are supported in the"
    " GEANT4 distribution:"
    "\n\n  ASCIITree (prints geometry hierarchy)"
    "\n\n  DAWN (socket connection to the Fukui Renderer DAWN) "
	 << G4VisFeaturesOfFukuiRenderer () <<
    "\n\n  DAWNFILE (file connection to the Fukui Renderer DAWN  ) "
	 << G4VisFeaturesOfDAWNFILE () <<
    "\n\n  GAGTree (prints geometry hierarchy, connectable to GAG"
    "\n  user interface)"
    "\n\n  HepRepFile (file connection to HepRep viewers, e.g., WIRED)"
    "\n\n  OPACS (the Orsay Package) "
    "\n\n  OpenGLIX (direct/immediate drawing on X Windows)\n"
       << G4VisFeaturesOfOpenGLIX () <<
    "\n\n  OpenGLSX (display List/stored drawing on X Windows)\n"
       << G4VisFeaturesOfOpenGLSX () <<
    "\n\n  OpenGLIXm (with Motif widgets)\n"
       << G4VisFeaturesOfOpenGLIXm () <<
    "\n\n  OpenGLSXm (with Motif widgets)\n"
       << G4VisFeaturesOfOpenGLSXm () <<
    "\n\n  Open Inventor"
       << G4VisFeaturesOfOpenInventor () <<
    "\n\n  RayTracer (produces JPEG file)"
    "\n\n  VRML1     (produces VRML 1 file over network)"
    "\n\n  VRML1FILE (produces VRML 1 file locally    )"
    "\n\n  VRML2     (produces VRML 2 file over network)"
    "\n\n  VRML2FILE (produces VRML 2 file locally    )"
       << G4endl;
}

void G4VisManager::PrintInstalledGraphicsSystems () const {
  G4cout << "\nThe following graphics systems drivers are installed on your"
    " system:"
       << "\n  ASCII Tree (produces ASCII file of geometry hierarchy)"
#ifdef G4VIS_BUILD_DAWN_DRIVER
       << "\n  DAWN     (socket connection to the Fukui Renderer DAWN)"
#endif
       << "\n  DAWNFILE (file connection to the Fukui Renderer DAWN)"
       << "\n  GAG Tree (produces ascii file of geometry hierarchy for GAG)"
#ifdef G4VIS_BUILD_OPACS_DRIVER
       << "\n  OPACS (the Orsay Package)"
#endif
#ifdef G4VIS_BUILD_OPENGLX_DRIVER
       << "\n  OpenGLIX (direct/immediate drawing on X Windows)"
       << "\n  OpenGLSX (display List/stored drawing on X Windows)"
#endif
#ifdef G4VIS_BUILD_OPENGLXM_DRIVER
       << "\n  OpenGLIXm (with Motif widgets)"
       << "\n  OpenGLSXm (with Motif widgets)"
#endif
#ifdef G4VIS_BUILD_OIX_DRIVER
       << "\n  Open Inventor X11"
#endif
#ifdef G4VIS_BUILD_OIWIN32_DRIVER
       << "\n  Open Inventor Win32"
#endif
       << "\n  RayTracer (produces JPEG file)"
#ifdef G4VIS_BUILD_VRML_DRIVER
       << "\n  VRML1 (produces VRML 1 file over network)"
       << "\n  VRML2 (produces VRML 2 file over network)"
#endif
       << "\n  VRML1FILE (produces VRML 1 file locally)"
       << "\n  VRML2FILE (produces VRML 2 file locally)"
       << G4endl;
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

void G4VisManager::PrintInvalidPointers () const {
  if (fVerbosity >= errors) {
    G4cout << "ERROR: G4VisManager::PrintInvalidPointers:";
    if (!fpGraphicsSystem) {
      G4cout << "\n null graphics system pointer.";
    }
    else {
      G4cout << "\n graphics system is " << fpGraphicsSystem -> GetName ()
	     << " but:";
      if (!fpScene) G4cout << "\nNull scene pointer. ";
      if (!fpSceneHandler) G4cout << "\nNull scene handler pointer. ";
      if (!fpViewer ) G4cout << "\nNull viewer pointer. ";
    }
    G4cout << G4endl;
  }
}

void G4VisManager::BeginOfRun () {
  //G4cout << "G4VisManager::BeginOfRun" << G4endl;
}

void G4VisManager::BeginOfEvent () {
  //G4cout << "G4VisManager::BeginOfEvent" << G4endl;
}

void G4VisManager::EndOfEvent () {
  //G4cout << "G4VisManager::EndOfEvent" << G4endl;
  if (GetConcreteInstance() && IsValidView ()) {
    ClearTransientStoreIfMarked();
    fVisManagerModelingParameters
      = *(fpSceneHandler -> CreateModelingParameters ());
    const G4std::vector<G4VModel*>& EOEModelList =
      fpScene -> GetEndOfEventModelList ();
    for (size_t i = 0; i < EOEModelList.size (); i++) {
      G4VModel* pModel = EOEModelList [i];
      pModel -> SetModelingParameters (&fVisManagerModelingParameters);
      fpSceneHandler -> SetModel (pModel);
      pModel -> DescribeYourselfTo (*fpSceneHandler);
    }
    if (fpScene->GetRefreshAtEndOfEvent()) {
      fpViewer->ShowView();  // ...for systems needing post processing.
      fpSceneHandler->SetMarkForClearingTransientStore(true);
    }
    fpSceneHandler -> SetModel (0);  // Flags invalid model.
  }
}

void G4VisManager::EndOfRun () {
  //G4cout << "G4VisManager::EndOfRun" << G4endl;
  if (GetConcreteInstance() && IsValidView ()) {
    fpSceneHandler->SetMarkForClearingTransientStore(true);
  }
}

void G4VisManager::ClearTransientStoreIfMarked(){
  // Assumes valid view.
  if (fpSceneHandler->GetMarkForClearingTransientStore()) {
    fpSceneHandler->SetMarkForClearingTransientStore(false);
    fpSceneHandler->ClearTransientStore();
  }
}

void G4VisManager::CheckModel () {
  G4VModel* pModel = fpSceneHandler->GetModel();
  if (!pModel) {  // provide a null model.
    pModel = &fVisManagerNullModel;
    fpSceneHandler -> SetModel (pModel);
  }
  // Ensure modeling parameters are right for this view...
  fVisManagerModelingParameters
    = *(fpSceneHandler -> CreateModelingParameters ());
  pModel->SetModelingParameters (&fVisManagerModelingParameters);
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

G4String G4VisManager::VerbosityGuidanceString
("\n  default:    warnings"
 "\nSimple graded message scheme - digit or string (1st character defines):"
 "\n  0) quiet,         // Nothing is printed."
 "\n  1) startup,       // Startup and endup messages are printed..."
 "\n  2) errors,        // ...and errors..."
 "\n  3) warnings,      // ...and warnings..."
 "\n  4) confirmations, // ...and confirming messages..."
 "\n  5) parameters,    // ...and parameters of scenes and views..."
 "\n  6) all            // ...and everything available."
 );

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
    const char* t = s;
    G4std::istrstream is((char*)t);
    is >> intVerbosity;
    if (!is) {
      G4cout << "ERROR: G4VisManager::GetVerbosityValue: invalid verbosity \""
	     << verbosityString << "\"\n"
	     << VerbosityGuidanceString
	     << "\n  Returning \"warnings\" == " << (G4int)warnings
	     << G4endl;
      verbosity = warnings;
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

  SetConcreteInstance(0);
  // Unless we survive a few preliminary tests, users must not use.

  if (!fpGraphicsSystem) return false;
  // Simply return without printing - we do not want printing if the
  // user simply does not want to use graphics, e.g., in batch mode.

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

  SetConcreteInstance(this);  // Unless we find another problem, users can use!
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
      SetConcreteInstance(0);  // Users must not use!
    }
    else {
      G4UImanager::GetUIpointer () -> ApplyCommand ("/vis/viewer/reset");
      if (fVerbosity >= warnings) {
	G4cout <<
	  "WARNING: G4VisManager: the scene was empty, \"world\" has been"
	  "\n  added and the view parameters have been reset.";
	G4cout << G4endl;
      }
    }
  }
  return isValid;
}

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManager.cc,v 1.10 1999-11-05 16:32:00 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager - John Allison 02/Jan/1996.

#include "G4VisManager.hh"

#include "G4UImanager.hh"
#include "G4StateManager.hh"
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

#include "G4VisManMessenger.hh"
// TEMPORARY /vis/ -> /vis~/ equivalence.
#include "G4VisToOldVisCommands.hh"

#include "G4ios.hh"

G4VisManager* G4VisManager::fpInstance = 0;

G4VisManager::G4VisManager ():
  fInitialised     (false),
  fpGraphicsSystem (0),
  fpScene          (0),
  fpSceneHandler   (0),
  fpViewer         (0),
  fVerbose         (0)  // All other objects use default constructors.
{
  if (fpInstance) {
    G4Exception
      ("G4VisManager: attempt to Construct more than one VisManager.");
  }
  else {
    fpInstance = this;

    G4cout << "Constructing Visualization Manager...." << endl;
    // Note: You might think that we could register graphics systems
    // and messengers here but we have to give the subclass time to
    // instantiate.  So the user has to invoke Initialise().
  }
}

G4VisManager::~G4VisManager () {
  fpInstance = 0;
  fSceneList.clearAndDestroy ();    
  fAvailableSceneHandlers.clearAndDestroy ();    
  fAvailableGraphicsSystems.clearAndDestroy ();
  G4cout << "Graphics systems deleted." << endl;
  G4cout << "VisManager deleting." << endl;
  fMessengerList.clearAndDestroy ();
  delete fpMessenger;
}

G4VisManager* G4VisManager::GetInstance () {
  if (!fpInstance) {
    G4Exception
      ("G4VisManager::GetInstance: VisManager not yet instantiated!");
  }
  return fpInstance;
}

void G4VisManager::Initialise () {

  G4cout << "Initialising Visualization Manager...." << endl;

  if (fVerbose > 0) {
    PrintAllGraphicsSystems ();
    PrintInstalledGraphicsSystems ();
    G4cout <<
      "\nYou have instantiated your own Visualization Manager, inheriting"
      "\n  G4VisManager and implementing RegisterGraphicsSystems.  An example"
      "\n  class MyVisManager is provided which is controlled by the"
      "\n  following environment variables:"
#ifdef G4VIS_BUILD_DAWN_DRIVER
      "\n    G4VIS_USE_DAWN"
#endif
#ifdef G4VIS_BUILD_DAWNFILE_DRIVER
      "\n    G4VIS_USE_DAWNFILE"
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
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
      "\n    G4VIS_USE_VRMLFILE"
#endif
      "\n  Thus, in your main() you have something like:"
      "\n    G4VisManager* visManager = new MyVisManager;"
      "\n    visManager -> Initialize ();"
      "\n  (Don't forget to delete visManager;)"
	 << endl;
  }
   
  G4cout << "Registering graphics systems...." << endl;

  RegisterGraphicsSystems ();

  RegisterMessengers ();

  // Remaining commands - not yet template-ised.
  fpMessenger = new G4VisManMessenger (this);

  // TEMPORARY /vis/ -> /vis~/ equivalence.
  fMessengerList.append (new G4VisToOldVisCommands);  

  fInitialised = true;
}

const G4GraphicsSystemList& G4VisManager::GetAvailableGraphicsSystems () {
  G4int nSystems = fAvailableGraphicsSystems.entries ();
  if (nSystems == 0) {
    G4cout << "G4VisManager::GetAvailableGraphicsSystems: WARNING: no graphics"
      " system available!"
      "\n  1) Did you have environment variables G4VIS_BUILD_xxxx_DRIVER set"
      "\n     when you compiled/built the visualization code?"
      "\n  2) Did you instantiate your own Visualization Manager and forget"
      "\n     to implement RegisterGraphicsSystems correctly?"
      "\n  3) You can register your own graphics system, e.g.,"
      "\n     G4VisManager::GetInstance () ->"
      "\n       RegisterGraphicsSystem (new MyGraphicsSystem);)"
	 << endl;
  }
  return fAvailableGraphicsSystems;
}

G4bool G4VisManager::RegisterGraphicsSystem (G4VGraphicsSystem* pSystem) {
  if (pSystem -> GetFunctionality () == G4VGraphicsSystem::noFunctionality) {
    G4cout << "G4VisManager::RegisterGraphicsSystem: WARNING: attempt to"
      "\n  register a \"no functionality\" graphics system, probably an"
      "\n  unbuilt system, i.e., a system not available locally.  Please"
      "\n  consult your computer manager.  System was "
	 << pSystem -> GetName ();
    if (pSystem -> GetNickname () != "") {
      G4cout << " (" << pSystem -> GetNickname () << ")";
    }
    G4cout << endl;
    return false;
  }
  else {
    fAvailableGraphicsSystems.append (pSystem);
    if (fVerbose > 0) {
      G4cout << "G4VisManager::RegisterGraphicsSystem: "
	   << pSystem -> GetName ();
      if (pSystem -> GetNickname () != "") {
	G4cout << " (" << pSystem -> GetNickname () << ")";
      }
      G4cout << " registered." << endl;
    }
    return true;
  }
}

void G4VisManager::CopyViewParameters () {
  if (IsValidView ()) {
    fVP = fpViewer -> GetViewParameters ();
  }
}

/***********
void G4VisManager::Clear () {

  // Clear current scene and current view, marking all its views as
  // needing refreshing.  This is a comprehensive clear which clears
  // both framebuffers of a double buffered system and clears the
  // scene's graphics databse (display lists, etc.) and clears the
  // current scene data.

  if (IsValidView ()) {

    fpViewer -> ClearView ();
    // Clears current buffer, i.e., back buffer in case of double
    // buffered system.

    fpViewer -> FinishView ();
    // Swaps buffers of double buffered systems.

    fpViewer -> ClearView ();
    // Clears the swapped buffer.

    fpViewer -> NeedKernelVisit ();
    // Informs all views of need to revisit kernel.

    fpSceneHandler -> ClearStore ();
    // Clears the graphics database (display lists) if any.

    fSD.Clear ();
    fpSceneHandler -> SetSceneData (fSD);
    // Clears current scene data - then updates scene.
  }
}
*****************/

/****************
void G4VisManager::ClearScene () {

  // Clear current scene, marking all its views as needing refreshing.
  // Clears the scene handler's graphics database (display lists, etc.)

  if (IsValidView ()) {

    fpViewer -> NeedKernelVisit ();
    // Informs all views of need to revisit kernel.

    fpSceneHandler -> ClearStore ();
    // Clears the graphics database (display lists) if any.

    fSD.Clear ();
    fpSceneHandler -> SetSceneData (fSD);
    // Clears current scene data - then updates scene.
  }
}
********************/

void G4VisManager::ClearView () {

  // Clear visible window of current view (both buffers of a double
  // buffered system).

  if (IsValidView ()) {

    fpViewer -> ClearView ();
    // Clears current buffer, i.e., back buffer in case of double
    // buffered system.

    fpViewer -> FinishView ();
    // Swaps buffers of double buffered systems.

    fpViewer -> ClearView ();
    // Clears the swapped buffer.
  }
}

void G4VisManager::Draw () {
  if (IsValidView ()) {
    fpViewer -> SetViewParameters (fVP);
    fpViewer -> DrawView ();
  }
}

void G4VisManager::Show () {
  if (IsValidView ()) {
     fpViewer -> ShowView ();
  }
}

void G4VisManager::Draw (const G4Polyline& line,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    G4ModelingParameters* pMP = fpSceneHandler -> CreateModelingParameters ();
    G4VModel* pModel = new G4NullModel (pMP);
    fpSceneHandler -> SetModel (pModel);
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (line);
    fpSceneHandler -> EndPrimitives ();
    fpSceneHandler -> SetModel (0);
    delete pModel;
    delete pMP;
  }
}

void G4VisManager::Draw (const G4Text& text,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    G4ModelingParameters* pMP = fpSceneHandler -> CreateModelingParameters ();
    G4VModel* pModel = new G4NullModel (pMP);
    fpSceneHandler -> SetModel (pModel);
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (text);
    fpSceneHandler -> EndPrimitives ();
    fpSceneHandler -> SetModel (0);
    delete pModel;
    delete pMP;
  }
}

void G4VisManager::Draw (const G4Circle& circle,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    G4ModelingParameters* pMP = fpSceneHandler -> CreateModelingParameters ();
    G4VModel* pModel = new G4NullModel (pMP);
    fpSceneHandler -> SetModel (pModel);
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (circle);
    fpSceneHandler -> EndPrimitives ();
    fpSceneHandler -> SetModel (0);
    delete pModel;
    delete pMP;
  }
}

void G4VisManager::Draw (const G4Square& Square,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    G4ModelingParameters* pMP = fpSceneHandler -> CreateModelingParameters ();
    G4VModel* pModel = new G4NullModel (pMP);
    fpSceneHandler -> SetModel (pModel);
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (Square);
    fpSceneHandler -> EndPrimitives ();
    fpSceneHandler -> SetModel (0);
    delete pModel;
    delete pMP;
  }
}

void G4VisManager::Draw (const G4Polymarker& polymarker,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    G4ModelingParameters* pMP = fpSceneHandler -> CreateModelingParameters ();
    G4VModel* pModel = new G4NullModel (pMP);
    fpSceneHandler -> SetModel (pModel);
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (polymarker);
    fpSceneHandler -> EndPrimitives ();
    fpSceneHandler -> SetModel (0);
    delete pModel;
    delete pMP;
  }
}

void G4VisManager::Draw (const G4Polyhedron& polyhedron,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    G4ModelingParameters* pMP = fpSceneHandler -> CreateModelingParameters ();
    G4VModel* pModel = new G4NullModel (pMP);
    fpSceneHandler -> SetModel (pModel);
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (polyhedron);
    fpSceneHandler -> EndPrimitives ();
    fpSceneHandler -> SetModel (0);
    delete pModel;
    delete pMP;
  }
}

void G4VisManager::Draw (const G4NURBS& nurbs,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    G4ModelingParameters* pMP = fpSceneHandler -> CreateModelingParameters ();
    G4VModel* pModel = new G4NullModel (pMP);
    fpSceneHandler -> SetModel (pModel);
    fpSceneHandler -> BeginPrimitives (objectTransform);
    fpSceneHandler -> AddPrimitive (nurbs);
    fpSceneHandler -> EndPrimitives ();
    fpSceneHandler -> SetModel (0);
    delete pModel;
    delete pMP;
  }
}

void G4VisManager::Draw (const G4VSolid& solid,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
  if (IsValidView ()) {
    G4ModelingParameters* pMP = fpSceneHandler -> CreateModelingParameters ();
    G4VModel* pModel = new G4NullModel (pMP);
    fpSceneHandler -> SetModel (pModel);
    fpSceneHandler -> PreAddThis (objectTransform, attribs);
    solid.DescribeYourselfTo (*fpSceneHandler);
    fpSceneHandler -> PostAddThis ();
    fpSceneHandler -> SetModel (0);
    delete pModel;
    delete pMP;
  }
}

void G4VisManager::Draw (const G4LogicalVolume& logicalVol,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
  // Find corresponding solid.
  G4VSolid* pSol = logicalVol.GetSolid ();
  Draw (*pSol, attribs, objectTransform);
}

void G4VisManager::Draw (const G4VPhysicalVolume& physicalVol,
			 const G4VisAttributes& attribs,
			 const G4Transform3D& objectTransform) {
  // Find corresponding logical volume and solid.
  G4LogicalVolume* pLV  = physicalVol.GetLogicalVolume ();
  G4VSolid*        pSol = pLV -> GetSolid ();
  Draw (*pSol, attribs, objectTransform);
}

G4bool G4VisManager::Notify (G4ApplicationState requestedState) {
  G4StateManager* stateManager = G4StateManager::GetStateManager ();
  if(stateManager -> GetPreviousState () == EventProc &&
     requestedState == GeomClosed) {
    if (fpConcreteInstance && IsValidView ()) {
      G4ModelingParameters* pMP =
	fpSceneHandler -> CreateModelingParameters ();
      const RWTPtrOrderedVector <G4VModel>& EOEModelList =
	fpScene -> GetEndOfEventModelList ();
      fpSceneHandler->ClearTransientStore(); //GB
      for (int i = 0; i < EOEModelList.entries (); i++) {
	G4VModel* pModel = EOEModelList [i];
	pModel -> SetModelingParameters (pMP);
	pModel -> DescribeYourselfTo (*fpSceneHandler);
	pModel -> SetModelingParameters (0);
      }
      delete pMP;
    }
  }
  return true;
}

void G4VisManager::CreateSceneHandler (G4String name) {
  if (!fInitialised) Initialise ();
  if (fpGraphicsSystem) {
    G4VSceneHandler* pSceneHandler =
      fpGraphicsSystem -> CreateSceneHandler (name);
    G4VViewer* pViewer;
    if (pSceneHandler) {
      fAvailableSceneHandlers.append (pSceneHandler);
      fpSceneHandler = pSceneHandler;                         // Make current.
    }
    else {
      G4cout << "Error in G4VisManager::CreateSceneHandler during "
	   << fpGraphicsSystem -> GetName ()
	   << " scene creation."
	   << endl;
      fpSceneHandler = 0;
      fpViewer  = 0;
      // Make it impossible for user action code to Draw.
      fpConcreteInstance = 0;
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
      fpSceneHandler -> AddViewerToList (fpViewer);
      fpSceneHandler -> SetCurrentViewer (fpViewer);
      // Make it possible for user action code to Draw.
      fpConcreteInstance = this;

      G4bool warnings = false;
      if (fVP.IsCulling () && fVP.IsCullingInvisible ()) {
	warnings = true;
	G4cout << "G4VisManager::CreateViewer: new viewer created:";
	if (fVerbose > 0) {
	  G4cout << " view parameters are:\n  " << fVP;
	}
	G4cout <<
	  "\n  NOTE: objects with visibility flag set to \"false\""
	  " will not be drawn!"
	  "\n  \"/vis~/set/culling off\" to Draw such objects.";
      }
      if (fVP.IsCullingCovered ()) {
	if (!warnings) {
	  G4cout << "G4VisManager::CreateViewer: new viewer created:";
	}
	warnings = true;
	G4cout <<
	  "\n  ALSO: covered objects in solid mode will not be part of"
	  " the scene!"
	  "\n  \"/vis~/set/cull_covered_daughters off\" to reverse this.";
      }
      if (warnings) {
	G4cout << "\n  Also see other \"/vis~/set\" commands."
	     << endl;
      }

    }
  }
  else PrintInvalidPointers ();
}

void G4VisManager::DeleteCurrentSceneHandler () {
  G4cout << "G4VisManager::DeleteCurrentSceneHandler: scene handler \""
	 << fpSceneHandler -> GetName ()
	 << "\"\n  and its viewer(s) are being deleted.";
  if(fpSceneHandler) {
    fAvailableSceneHandlers.remove (fpSceneHandler);
    delete fpSceneHandler;
  }
  const G4SceneHandlerList& sceneHandlerList = fAvailableSceneHandlers;
  G4int nSH = sceneHandlerList.entries ();
  G4int iSH;
  for (iSH = 0; iSH < nSH; iSH++) {
    if (sceneHandlerList [iSH] -> GetViewerList ().entries ()) break;
  }
  if (iSH < nSH) {
    fpSceneHandler = sceneHandlerList [iSH];
    G4cout << "\n  scene handler is now \""
	   << fpSceneHandler -> GetName ();
    fpViewer = fpSceneHandler -> GetViewerList () [0];
    G4cout << "\"\n  and viewer now \""
	   << fpViewer -> GetName ()
	   << ".";
    IsValidView ();  // Check.
  }
  else if (nSH) {
    fpSceneHandler = fAvailableSceneHandlers [0];
    G4cout << "\n  scene handler is now \""
	   << fpSceneHandler -> GetName ();
    fpViewer = 0;
    // Make it impossible for user action code to Draw.
    fpConcreteInstance = 0;
    G4cout << "\" but it has no viewers -\n  please create one.";
  }
  else {
    fpSceneHandler = 0;
    fpViewer  = 0;
    // Make it impossible for user action code to Draw.
    fpConcreteInstance = 0;
    G4cout <<
      "\n  There are now no scene handlers left.  /vis/sceneHandler/select"
      "\n  to select a scene handler and /vis/viewer/select to select"
      "\n  a viewer, or maybe you will have to create some.";
  }
  G4cout << endl;
}

void G4VisManager::DeleteCurrentViewer () {
  G4cout << "G4VisManager::DeleteCurrentViewer: viewer \""
	 << fpViewer -> GetName () 
	 << "\" being deleted.";
  if (fpViewer ) {
    fpViewer -> GetScene () -> RemoveViewerFromList (fpViewer);
    delete fpViewer;
  }
  const G4ViewerList& viewerList = fpSceneHandler -> GetViewerList ();
  if (viewerList.entries () > 0) {
    fpViewer = viewerList [0];
    fpSceneHandler -> SetCurrentViewer (fpViewer);
    if (IsValidView ()) {
      G4cout << "\n  viewer is now \""
	     << fpViewer -> GetName ()
	     << "\".";
    }
  }
  else {
    fpViewer = 0;
    fpSceneHandler -> SetCurrentViewer (0);
    // Make it impossible for user action code to Draw.
    fpConcreteInstance = 0;
    G4cout <<
      "\n  There are now no viewers left for this scene handler."
      "\n  /vis/sceneHandler/select to select a different scene handler"
      "\n  and /vis/viewer/select to select a viewer, or maybe you will"
      "\n  have to create some.";
  }
  G4cout << endl;
}

void G4VisManager::GeometryHasChanged () {
  G4cout << "G4VisManager::GeometryHasChanged() called." << endl;

  // Change the world...
  G4VPhysicalVolume* pWorld =
    G4TransportationManager::GetTransportationManager ()
    -> GetNavigatorForTracking () -> GetWorldVolume ();
  if (!pWorld) {
    G4cout << "  The world has ended!!!  (Is this serious?)" << endl;
  }

  // Check scenes.
  G4SceneList& sceneList = fSceneList;
  G4int iScene, nScenes = sceneList.entries ();
  for (iScene = 0; iScene < nScenes; iScene++) {
    G4Scene* pScene = sceneList [iScene];
    RWTPtrOrderedVector <G4VModel>& modelList =
      pScene -> SetRunDurationModelList ();

    G4int nModels, iModel;
    G4bool modelInvalid;
    do {  // Remove, if required, one at a time.
      modelInvalid = false;
      nModels = modelList.entries ();
      for (iModel = 0; iModel < nModels; iModel++) {
	if (modelInvalid = !(modelList [iModel] -> Validate ())) {
	  // Model invalid - remove and break.
	  G4cout << "  Model \""
		 << modelList [iModel] -> GetGlobalDescription ()
		 <<
	    "\" is no longer valid - being removed\n  from scene \""
		 << pScene -> GetName () << "\""
		 << endl;
	  modelList.removeAt (iModel);
	    break;
	}
      }
    } while (modelInvalid);

    nModels = modelList.entries ();
    if (nModels == 0) {
      G4cout << "  No models left in this scene \""
	     << pScene -> GetName ()
	     << "\"."
	     << endl;
    }
  }

  // Check the manager's current scene...
  if (fpScene && fpScene -> GetRunDurationModelList ().entries () == 0) {
    G4cout << "  The current scene \""
	   << fpScene -> GetName ()
	   << "\" has no models left."
	   << endl;
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
    fSceneList.append (fpScene);
    G4cout << "G4VisManager::SetCurrentGraphicsSystemAndCreateViewer:"
      "\n  Empty scene \"scene-vis~\" created."
	   << endl;
  }
  CreateSceneHandler ();
  CreateViewer ();
}

void G4VisManager::SetCurrentGraphicsSystem (G4VGraphicsSystem* pSystem) {
  fpGraphicsSystem = pSystem;
  G4cout << "G4VisManager::SetCurrentGraphicsSystem: system now "
	 << pSystem -> GetName ();
  // If current scene handler is of same graphics system, leave unchanged.
  // Else find the most recent scene handler of same graphics system.
  // Or clear pointers.
  if (!(fpSceneHandler && fpSceneHandler -> GetGraphicsSystem () == pSystem)) {
    const G4SceneHandlerList& sceneHandlerList = fAvailableSceneHandlers;
    G4int nSH = sceneHandlerList.entries ();  // No. of scene handlers.
    G4int iSH;
    for (iSH = nSH - 1; iSH >= 0; iSH--) {
      if (sceneHandlerList [iSH] -> GetGraphicsSystem () == pSystem) break;
    }
    if (iSH >= nSH) {
      fpSceneHandler = sceneHandlerList [iSH];
      G4cout << "\n  Scene Handler now "
	     << fpSceneHandler -> GetName ();
      const G4ViewerList& viewerList = fpSceneHandler -> GetViewerList ();
      if (viewerList.entries ()) {
	fpViewer = viewerList [0];
	G4cout << "\n  Viewer now " << fpViewer -> GetName ();
      }
      else {
	fpViewer = 0;
	// Make it impossible for user action code to Draw.
	fpConcreteInstance = 0;
      }
    }
    else {
      fpSceneHandler = 0;
      fpViewer = 0;
      // Make it impossible for user action code to Draw.
    fpConcreteInstance = 0;
    }
  }
  G4cout << endl;
}

void G4VisManager::SetCurrentSceneHandler (G4VSceneHandler* pSceneHandler) {
  fpSceneHandler = pSceneHandler;
  G4cout << "G4VisManager::SetCurrentSceneHandler: scene handler now \""
	 << pSceneHandler -> GetName () << "\"";
  if (fpGraphicsSystem != pSceneHandler -> GetGraphicsSystem ()) {
    fpGraphicsSystem = pSceneHandler -> GetGraphicsSystem ();
    G4cout << "\n  graphics system now \""
	   << fpGraphicsSystem -> GetName () << "\"";
  }
  const G4ViewerList& viewerList = fpSceneHandler -> GetViewerList ();
  G4int nViewers = viewerList.entries ();
  if (nViewers) {
    G4int iViewer;
    for (iViewer = 0; iViewer < nViewers; iViewer++) {
      if (fpViewer == viewerList [iViewer]) break;
    }
    if (iViewer >= nViewers) {
      fpViewer = viewerList [0];
      G4cout << "\n  Viewer now \"" << fpViewer -> GetName () << "\"";
    }
    IsValidView ();  // Checks.
  }
  else {
    fpViewer = 0;
    // Make it impossible for user action code to Draw.
    fpConcreteInstance = 0;
    G4cout << "\n  No viewers for this scene handler - please create one.";
  }
  G4cout << endl; 
}

void G4VisManager::SetCurrentViewer (G4VViewer* pViewer) {
  fpViewer  = pViewer;
  G4cout << "G4VisManager::SetCurrentViewer: viewer now "
	 << pViewer -> GetName ()
	 << endl;
  fpSceneHandler = fpViewer -> GetScene ();
  fpSceneHandler -> SetCurrentViewer (pViewer);
  fpGraphicsSystem = fpSceneHandler -> GetGraphicsSystem ();
  IsValidView ();  // Checks.
}

void G4VisManager::PrintCurrentSystems () const {
  if (fpGraphicsSystem && fpSceneHandler && fpViewer) {
    G4int nSystems = fAvailableGraphicsSystems.entries ();
    if (nSystems <= 0) {
      G4cout << "No graphics systems available yet." << endl;
    }
    else {
      PrintAvailableGraphicsSystems ();
    }
  }
  else PrintInvalidPointers ();
}

void G4VisManager::PrintCurrentSystem () const {
  if (fpGraphicsSystem && fpSceneHandler && fpViewer) {
    G4cout << "Current graphics system is: " << *fpGraphicsSystem;
    G4cout << endl;
  }
  else PrintInvalidPointers ();
}

void G4VisManager::PrintCurrentScene () const {
  if (fpGraphicsSystem && fpSceneHandler && fpViewer) {
    G4cout << "Current Scene is: " << fpSceneHandler -> GetName ();
    G4cout << '\n' << *fpSceneHandler;
    G4cout << endl;
  }
  else PrintInvalidPointers ();
}

void G4VisManager::PrintCurrentView () const {
  if (fpGraphicsSystem && fpSceneHandler && fpViewer) {
    G4cout << "Current View is: ";
    G4cout << fpGraphicsSystem -> GetName () << '-'
	   << fpSceneHandler  -> GetSceneId () << '-'
	   << fpViewer   -> GetViewId ()
	   << " selected (check: " << fpViewer -> GetName () << ").";
    G4cout << '\n' << *fpViewer;
    G4cout << "\nCurrent view parameters";
    if (fVP != fpViewer -> GetViewParameters ()) {
      G4cout << " differ in the following respect:\n";
      fVP.PrintDifferences (fpViewer -> GetViewParameters ());
    }
    else {
      G4cout << " are same."
	   << endl;
    }
  }
  else PrintInvalidPointers ();
}

void G4VisManager::PrintAllGraphicsSystems () const {
  G4cout << "\nThe following graphics systems drivers are supported in the"
    " GEANT4 distribution:"
       << "\n\n  DAWN     (socket connection to the Fukui Renderer DAWN) " << G4VisFeaturesOfFukuiRenderer ()
       << "\n\n  DAWNFILE (file connection to the Fukui Renderer DAWN  ) " << G4VisFeaturesOfDAWNFILE      ()
       << "\n\n  OPACS (the Orsay Package) "
       << "\n\n  OpenGLIX (direct/immediate drawing on X Windows)\n"
       << G4VisFeaturesOfOpenGLIX ()
       << "\n\n  OpenGLSX (display List/stored drawing on X Windows)\n"
       << G4VisFeaturesOfOpenGLSX ()
       << "\n\n  OpenGLIXm (with Motif widgets)\n"
       << G4VisFeaturesOfOpenGLIXm ()
       << "\n\n  OpenGLSXm (with Motif widgets)\n"
       << G4VisFeaturesOfOpenGLSXm ()
       << "\n\n  Open Inventor"
       << G4VisFeaturesOfOpenInventor ()
       << "\n\n  VRML1     (produces VRML 1 file over network)"
       << "\n\n  VRML1File (produces VRML 1 file locally    )"
       << "\n\n  VRML2     (produces VRML 2 file over network), in preparation"
       << "\n\n  VRML2File (produces VRML 2 file locally    ), in preparation"
       << endl;
}

void G4VisManager::PrintInstalledGraphicsSystems () const {
  G4cout << "\nThe following graphics systems drivers are installed on your"
    " system:"
#ifdef G4VIS_BUILD_DAWN_DRIVER
       << "\n  DAWN     (socket connection to the Fukui Renderer DAWN)"
#endif
#ifdef G4VIS_BUILD_DAWNFILE_DRIVER
       << "\n  DAWNFILE (file connection to the Fukui Renderer DAWN)"
#endif
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
#ifdef G4VIS_BUILD_VRML_DRIVER
       << "\n  VRML1 (produces VRML 1 file over network)"
       << "\n  VRML2 (produces VRML 2 file over network), in preparation"
#endif
#ifdef G4VIS_BUILD_VRMLFILE_DRIVER
       << "\n  VRML1File (produces VRML 1 file locally)"
       << "\n  VRML2File (produces VRML 2 file locally), in preparation"
#endif
       << endl;
}

void G4VisManager::PrintAvailableGraphicsSystems () const {
  G4int nSystems = fAvailableGraphicsSystems.entries ();
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
  G4cout << endl;
}

void G4VisManager::PrintInvalidPointers () const {
  G4cout << "Error in G4VisManager::PrintInvalidPointers:";
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
  G4cout << endl;
}

void G4VisManager::RefreshCurrentView  () {  // Soft clear, then redraw.
  if (IsValidView ()) {
    fpViewer -> ClearView ();  // Soft clear, i.e., clears buffer only.
    Draw ();
    Show ();
  }
}

G4bool G4VisManager::IsValidView () {

  if (!fInitialised) Initialise ();

  fpConcreteInstance = 0;
  // Unless we survive a few preliminary tests, users must not use.

  if (!fpGraphicsSystem) return false;
  // Simply return without printing - we do not want printing if the
  // user simply does not want to use graphics, e.g., in batch mode.

  if ((!fpScene) || (!fpSceneHandler) || (!fpViewer)) {
    G4cout << "G4VisManager::IsValidView ():"
      "\n  Current view is not valid.";
    G4cout << endl;
    PrintInvalidPointers ();
    return false;
  }

  if (fpScene != fpSceneHandler -> GetScene ()) {
    G4cout << "G4VisManager::IsValidView ():";
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
	     << endl;
    }
    else {
      G4cout << "\n  Scene handler \""
	     << fpSceneHandler -> GetName ()
	     << "\" has null scene pointer."
	"\n  Attach a scene with /vis/sceneHandler/attach [<scene-name>]"
	     << endl;
    }
    return false;
  }

  const G4ViewerList& viewerList = fpSceneHandler -> GetViewerList ();
  if (viewerList.entries () == 0) {
    G4cout << "G4VisManager::IsValidView (): the current scene handler\n  \""
	   << fpSceneHandler -> GetName ()
	   << "\" has no viewers.  Do /vis/viewer/create."
	   << endl;
    return false;
  }

  fpConcreteInstance = this;  // Unless we find another problem, users can use!
  G4bool isValid = true;
  if (fpScene -> IsEmpty ()) {
    fpScene -> AddWorldIfEmpty ();  // Add world by default if possible.
    if (fpScene -> IsEmpty ()) {    // If still empty...
      G4cout << "G4VisManager::IsViewValid ():";
      G4cout <<
	"\n  Attempt at some drawing operation when scene is empty."
	"\n  Maybe the geometry has not yet been defined."
	"  Try /run/initialize."
	     << endl;
      isValid = false;
      fpConcreteInstance = 0;  // Users must not use!
    }
    else {
      G4UImanager::GetUIpointer () -> ApplyCommand ("/vis/camera/reset");
      G4cout <<
	"G4VisManager: the scene was empty, \"world\" has been added"
	"\n  and the view parameters have been reset.";
      G4cout << endl;
    }
  }
  return isValid;
}

// "No functionality" graphics systems to trap accidental attempt to
// use unbuilt systems.

#ifndef G4VIS_BUILD_DAWN_DRIVER

class G4FukuiRenderer: public G4VGraphicsSystem {
public:
  G4FukuiRenderer ();
};
G4FukuiRenderer::G4FukuiRenderer ():
  G4VGraphicsSystem ("FukuiRenderer",
                     "DAWN",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

#ifndef G4VIS_BUILD_DAWNFILE_DRIVER

class G4DAWNFILE: public G4VGraphicsSystem {
public:
  G4DAWNFILE ();
};
G4DAWNFILE::G4DAWNFILE ():
  G4VGraphicsSystem ("FukuiRendererFile",
                     "DAWNFILE",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

#ifndef G4VIS_BUILD_OPACS_DRIVER

class G4Wo: public G4VGraphicsSystem {
public:
  G4Wo ();
};
G4Wo::G4Wo ():
  G4VGraphicsSystem ("Wo",
		     G4VGraphicsSystem::noFunctionality) {}

class G4Xo: public G4VGraphicsSystem {
public:
  G4Xo ();
};
G4Xo::G4Xo ():
  G4VGraphicsSystem ("Xo",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

#ifndef G4VIS_BUILD_OPENGLX_DRIVER

class G4OpenGLImmediateX: public G4VGraphicsSystem {
public:
  G4OpenGLImmediateX ();
};
G4OpenGLImmediateX::G4OpenGLImmediateX ():
  G4VGraphicsSystem ("OpenGLImmediateX",
                     "OGLIX",
		     G4VGraphicsSystem::noFunctionality) {}

class G4OpenGLStoredX: public G4VGraphicsSystem {
public:
  G4OpenGLStoredX ();
};
G4OpenGLStoredX::G4OpenGLStoredX ():
  G4VGraphicsSystem ("OpenGLStoredX",
                     "OGLSX",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

#ifndef G4VIS_BUILD_OPENGLXM_DRIVER

class G4OpenGLImmediateXm: public G4VGraphicsSystem {
public:
  G4OpenGLImmediateXm ();
};
G4OpenGLImmediateXm::G4OpenGLImmediateXm ():
  G4VGraphicsSystem ("OpenGLImmediateXm",
                     "OGLIXm",
		     G4VGraphicsSystem::noFunctionality) {}

class G4OpenGLStoredXm: public G4VGraphicsSystem {
public:
  G4OpenGLStoredXm ();
};
G4OpenGLStoredXm::G4OpenGLStoredXm ():
  G4VGraphicsSystem ("OpenGLStoredXm",
                     "OGLSXm",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

#ifndef G4VIS_BUILD_OIX_DRIVER

class G4OpenInventorX: public G4VGraphicsSystem {
public:
  G4OpenInventorX ();
};
G4OpenInventorX::G4OpenInventorX ():
  G4VGraphicsSystem ("OpenInventorX",
		     "OIX",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

#ifndef G4VIS_BUILD_OIWIN32_DRIVER

class G4OpenInventorWin32: public G4VGraphicsSystem {
public:
  G4OpenInventorWin32 ();
};
G4OpenInventorWin32::G4OpenInventorWin32 ():
  G4VGraphicsSystem ("OpenInventorWin32",
		     "OIWIN32",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

#ifndef G4VIS_BUILD_VRML_DRIVER

class G4VRML1: public G4VGraphicsSystem {
public:
  G4VRML1 ();
};
G4VRML1::G4VRML1 ():
  G4VGraphicsSystem ("VRML1.0",
		     "VRML1",
		     G4VGraphicsSystem::noFunctionality) {}

class G4VRML2: public G4VGraphicsSystem {
public:
  G4VRML2 ();
};
G4VRML2::G4VRML2 ():
  G4VGraphicsSystem ("VRML2.0",
		     "VRML2",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

#ifndef G4VIS_BUILD_VRMLFILE_DRIVER

class G4VRML1File: public G4VGraphicsSystem {
public:
  G4VRML1File ();
};
G4VRML1File::G4VRML1File ():
  G4VGraphicsSystem ("VRML1.0File",
		     "VRML1File",
		     G4VGraphicsSystem::noFunctionality) {}

class G4VRML2File: public G4VGraphicsSystem {
public:
  G4VRML2File ();
};
G4VRML2File::G4VRML2File ():
  G4VGraphicsSystem ("VRML2.0File",
		     "VRML2File",
		     G4VGraphicsSystem::noFunctionality) {}

#endif

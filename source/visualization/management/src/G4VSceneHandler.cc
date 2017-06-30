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
//
// $Id: G4VSceneHandler.cc 102802 2017-02-22 15:19:00Z gcosmo $
//
// 
// John Allison  19th July 1996
// Abstract interface class for graphics scenes.

#include "G4VSceneHandler.hh"

#include "G4ios.hh"
#include <sstream>

#include "G4VisManager.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VViewer.hh"
#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4Polyline.hh"
#include "G4Scale.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polymarker.hh"
#include "G4Polyhedron.hh"
#include "G4Visible.hh"
#include "G4VisAttributes.hh"
#include "G4VModel.hh"
#include "G4TrajectoriesModel.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4DisplacedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ModelingParameters.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4HitsModel.hh"
#include "G4VHit.hh"
#include "G4VDigi.hh"
#include "G4ScoringManager.hh"
#include "G4DefaultLinearColorMap.hh"
#include "Randomize.hh"
#include "G4StateManager.hh"
#include "G4RunManager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif
#include "G4Run.hh"
#include "G4Transform3D.hh"
#include "G4AttHolder.hh"
#include "G4AttDef.hh"
#include "G4VVisCommand.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4VSceneHandler::G4VSceneHandler (G4VGraphicsSystem& system, G4int id, const G4String& name):
  fSystem                (system),
  fSceneHandlerId        (id),
  fViewCount             (0),
  fpViewer               (0),
  fpScene                (0),
  fMarkForClearingTransientStore (true),  // Ready for first
					  // ClearTransientStoreIfMarked(),
					  // e.g., at end of run (see
					  // G4VisManager.cc).
  fReadyForTransients    (true),  // Only false while processing scene.
  fProcessingSolid       (false),
  fProcessing2D          (false),
  fpModel                (0),
  fNestingDepth          (0),
  fpVisAttribs           (0)
{
  G4VisManager* pVMan = G4VisManager::GetInstance ();
  fpScene = pVMan -> GetCurrentScene ();
  if (name == "") {
    std::ostringstream ost;
    ost << fSystem.GetName () << '-' << fSceneHandlerId;
    fName = ost.str();
  }
  else {
    fName = name;
  }
  fTransientsDrawnThisEvent = pVMan->GetTransientsDrawnThisEvent();
  fTransientsDrawnThisRun = pVMan->GetTransientsDrawnThisRun();
}

G4VSceneHandler::~G4VSceneHandler () {
  G4VViewer* last;
  while( ! fViewerList.empty() ) {
    last = fViewerList.back();
    fViewerList.pop_back();
    delete last;
  }
}

const G4VisExtent& G4VSceneHandler::GetExtent() const
{
  if (fpScene) {
    return fpScene->GetExtent();
  } else {
    static const G4VisExtent defaultExtent = G4VisExtent();
    return defaultExtent;
  }
}

void G4VSceneHandler::PreAddSolid (const G4Transform3D& objectTransformation,
				   const G4VisAttributes& visAttribs) {
  fObjectTransformation = objectTransformation;
  fpVisAttribs = &visAttribs;
  fProcessingSolid = true;
}

void G4VSceneHandler::PostAddSolid () {
  fpVisAttribs = 0;
  fProcessingSolid = false;
  if (fReadyForTransients) {
    fTransientsDrawnThisEvent = true;
    fTransientsDrawnThisRun = true;
  }
}

void G4VSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation) {
  //static G4int count = 0;
  //G4cout << "G4VSceneHandler::BeginPrimitives: " << count++ << G4endl;
  fNestingDepth++;
  if (fNestingDepth > 1)
    G4Exception
      ("G4VSceneHandler::BeginPrimitives",
       "visman0101", FatalException,
       "Nesting detected. It is illegal to nest Begin/EndPrimitives.");
  fObjectTransformation = objectTransformation;
}

void G4VSceneHandler::EndPrimitives () {
  if (fNestingDepth <= 0)
    G4Exception("G4VSceneHandler::EndPrimitives",
		"visman0102", FatalException, "Nesting error.");
  fNestingDepth--;
  if (fReadyForTransients) {
    fTransientsDrawnThisEvent = true;
    fTransientsDrawnThisRun = true;
  }
}

void G4VSceneHandler::BeginPrimitives2D
(const G4Transform3D& objectTransformation) {
  fNestingDepth++;
  if (fNestingDepth > 1)
    G4Exception
      ("G4VSceneHandler::BeginPrimitives2D",
       "visman0103", FatalException,
       "Nesting detected. It is illegal to nest Begin/EndPrimitives.");
  fObjectTransformation = objectTransformation;
  fProcessing2D = true;
}

void G4VSceneHandler::EndPrimitives2D () {
  if (fNestingDepth <= 0)
    G4Exception("G4VSceneHandler::EndPrimitives2D",
		"visman0104", FatalException, "Nesting error.");
  fNestingDepth--;
  if (fReadyForTransients) {
    fTransientsDrawnThisEvent = true;
    fTransientsDrawnThisRun = true;
  }
  fProcessing2D = false;
}

void G4VSceneHandler::BeginModeling () {
}

void G4VSceneHandler::EndModeling ()
{
  fpModel = 0;
}

void G4VSceneHandler::ClearStore () {}

void G4VSceneHandler::ClearTransientStore () {}

template <class T> void G4VSceneHandler::AddSolidT
(const T& solid)
{
  // Get and check applicable vis attributes.
  fpVisAttribs = fpViewer->GetApplicableVisAttributes(fpVisAttribs);
  RequestPrimitives (solid);
}

template <class T> void G4VSceneHandler::AddSolidWithAuxiliaryEdges
(const T& solid)
{
  // Get and check applicable vis attributes.
  fpVisAttribs = fpViewer->GetApplicableVisAttributes(fpVisAttribs);
  // Draw with auxiliary edges unless otherwise specified.
  if (!fpVisAttribs->IsForceAuxEdgeVisible()) {
    // Create a vis atts object for the modified vis atts.
    // It is static so that we may return a reliable pointer to it.
    static G4VisAttributes visAttsWithAuxEdges;
    // Initialise it with the current vis atts and reset the pointer.
    visAttsWithAuxEdges = *fpVisAttribs;
    // Force auxiliary edges visible.
    visAttsWithAuxEdges.SetForceAuxEdgeVisible();
    fpVisAttribs = &visAttsWithAuxEdges;
  }
  RequestPrimitives (solid);
}

void G4VSceneHandler::AddSolid (const G4Box& box) {
  AddSolidT (box);
  // If your graphics system is sophisticated enough to handle a
  //  particular solid shape as a primitive, in your derived class write a
  //  function to override this.
  // Your function might look like this...
  // void G4MySceneHandler::AddSolid (const G4Box& box) {
  // Get and check applicable vis attributes.
  //   fpVisAttribs = fpViewer->GetApplicableVisAttributes(fpVisAttribs);
  // Do not draw if not visible.
  //   if (fpVisAttribs->IsVisible()) {
  //   Get parameters of appropriate object, e.g.:
  //     G4double dx = box.GetXHalfLength ();
  //     G4double dy = box.GetYHalfLength ();
  //     G4double dz = box.GetZHalfLength ();
  //     ...
  //     and Draw or Store in your display List.
}

void G4VSceneHandler::AddSolid (const G4Cons& cons) {
  AddSolidT (cons);
}

void G4VSceneHandler::AddSolid (const G4Orb& orb) {
  AddSolidWithAuxiliaryEdges (orb);
}

void G4VSceneHandler::AddSolid (const G4Para& para) {
  AddSolidT (para);
}

void G4VSceneHandler::AddSolid (const G4Sphere& sphere) {
  AddSolidWithAuxiliaryEdges (sphere);
}

void G4VSceneHandler::AddSolid (const G4Torus& torus) {
  AddSolidWithAuxiliaryEdges (torus);
}

void G4VSceneHandler::AddSolid (const G4Trap& trap) {
  AddSolidT (trap);
}

void G4VSceneHandler::AddSolid (const G4Trd& trd) {
  AddSolidT (trd);
}

void G4VSceneHandler::AddSolid (const G4Tubs& tubs) {
  AddSolidT (tubs);
}

void G4VSceneHandler::AddSolid (const G4Ellipsoid& ellipsoid) {
  AddSolidWithAuxiliaryEdges (ellipsoid);
}

void G4VSceneHandler::AddSolid (const G4Polycone& polycone) {
  AddSolidT (polycone);
}

void G4VSceneHandler::AddSolid (const G4Polyhedra& polyhedra) {
  AddSolidT (polyhedra);
}

void G4VSceneHandler::AddSolid (const G4VSolid& solid) {
  AddSolidT (solid);
}

void G4VSceneHandler::AddCompound (const G4VTrajectory& traj) {
  G4TrajectoriesModel* trajectoriesModel =
    dynamic_cast<G4TrajectoriesModel*>(fpModel);
  if (trajectoriesModel)
    traj.DrawTrajectory();
  else {
    G4Exception
    ("G4VSceneHandler::AddCompound(const G4VTrajectory&)",
     "visman0105", FatalException, "Not a G4TrajectoriesModel.");
  }
}

void G4VSceneHandler::AddCompound (const G4VHit& hit) {
  // Cast away const because Draw is non-const!!!!
  const_cast<G4VHit&>(hit).Draw();
}

void G4VSceneHandler::AddCompound (const G4VDigi& digi) {
  // Cast away const because Draw is non-const!!!!
  const_cast<G4VDigi&>(digi).Draw();
}

void G4VSceneHandler::AddCompound (const G4THitsMap<G4double>& hits) {
  //G4cout << "AddCompound: hits: " << &hits << G4endl;
  G4bool scoreMapHits = false;
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManagerIfExist();
  if (scoringManager) {
    size_t nMeshes = scoringManager->GetNumberOfMesh();
    for (size_t iMesh = 0; iMesh < nMeshes; ++iMesh) {
      G4VScoringMesh* mesh = scoringManager->GetMesh(iMesh);
      if (mesh && mesh->IsActive()) {
	MeshScoreMap scoreMap = mesh->GetScoreMap();
        const G4String& mapNam = const_cast<G4THitsMap<G4double>&>(hits).GetName();
	for(MeshScoreMap::const_iterator i = scoreMap.begin();
	    i != scoreMap.end(); ++i) {
	  const G4String& scoreMapName = i->first;
	  if (scoreMapName == mapNam) {
	    G4DefaultLinearColorMap colorMap("G4VSceneHandlerColorMap");
	    scoreMapHits = true;
	    mesh->DrawMesh(scoreMapName, &colorMap);
	  }
	}
      }
    }
  }
  if (scoreMapHits) {
    static G4bool first = true;
    if (first) {
      first = false;
      G4cout <<
	"Scoring map drawn with default parameters."
	"\n  To get gMocren file for gMocren browser:"
	"\n    /vis/open gMocrenFile"
	"\n    /vis/viewer/flush"
	"\n  Many other options available with /score/draw... commands."
	"\n  You might want to \"/vis/viewer/set/autoRefresh false\"."
	     << G4endl;
    }
  } else {  // Not score map hits.  Just call DrawAllHits.
    // Cast away const because DrawAllHits is non-const!!!!
    const_cast<G4THitsMap<G4double>&>(hits).DrawAllHits();
  }
}

void G4VSceneHandler::AddCompound (const G4THitsMap<G4StatDouble>& hits) {
  //G4cout << "AddCompound: hits: " << &hits << G4endl;
  G4bool scoreMapHits = false;
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManagerIfExist();
  if (scoringManager) {
    size_t nMeshes = scoringManager->GetNumberOfMesh();
    for (size_t iMesh = 0; iMesh < nMeshes; ++iMesh) {
      G4VScoringMesh* mesh = scoringManager->GetMesh(iMesh);
      if (mesh && mesh->IsActive()) {
	MeshScoreMap scoreMap = mesh->GetScoreMap();
	for(MeshScoreMap::const_iterator i = scoreMap.begin();
	    i != scoreMap.end(); ++i) {
	  const G4String& scoreMapName = i->first;
	  const G4THitsMap<G4StatDouble>* foundHits = i->second;
	  if (foundHits == &hits) {
	    G4DefaultLinearColorMap colorMap("G4VSceneHandlerColorMap");
	    scoreMapHits = true;
	    mesh->DrawMesh(scoreMapName, &colorMap);
	  }
	}
      }
    }
  }
  if (scoreMapHits) {
    static G4bool first = true;
    if (first) {
      first = false;
      G4cout <<
	"Scoring map drawn with default parameters."
	"\n  To get gMocren file for gMocren browser:"
	"\n    /vis/open gMocrenFile"
	"\n    /vis/viewer/flush"
	"\n  Many other options available with /score/draw... commands."
	"\n  You might want to \"/vis/viewer/set/autoRefresh false\"."
	     << G4endl;
    }
  } else {  // Not score map hits.  Just call DrawAllHits.
    // Cast away const because DrawAllHits is non-const!!!!
    const_cast<G4THitsMap<G4StatDouble>&>(hits).DrawAllHits();
  }
}

void G4VSceneHandler::AddViewerToList (G4VViewer* pViewer) {
  fViewerList.push_back (pViewer);
}

void G4VSceneHandler::AddPrimitive (const G4Scale& scale) {

  const G4double margin(0.01);
  // Fractional margin - ensures scale is comfortably inside viewing
  // volume.
  const G4double oneMinusMargin (1. - margin);

  const G4VisExtent& sceneExtent = fpScene->GetExtent();

  // Useful constants...
  const G4double length(scale.GetLength());
  const G4double halfLength(length / 2.);
  const G4double tickLength(length / 20.);
  const G4double piBy2(halfpi);

  // Get size of scene...
  const G4double xmin = sceneExtent.GetXmin();
  const G4double xmax = sceneExtent.GetXmax();
  const G4double ymin = sceneExtent.GetYmin();
  const G4double ymax = sceneExtent.GetYmax();
  const G4double zmin = sceneExtent.GetZmin();
  const G4double zmax = sceneExtent.GetZmax();

  // Create (empty) polylines having the same vis attributes...
  G4Polyline scaleLine, tick11, tick12, tick21, tick22;
  G4VisAttributes visAtts(*scale.GetVisAttributes());  // Long enough life.
  scaleLine.SetVisAttributes(&visAtts);
  tick11.SetVisAttributes(&visAtts);
  tick12.SetVisAttributes(&visAtts);
  tick21.SetVisAttributes(&visAtts);
  tick22.SetVisAttributes(&visAtts);

  // Add points to the polylines to represent an scale parallel to the
  // x-axis centred on the origin...
  G4Point3D r1(G4Point3D(-halfLength, 0., 0.));
  G4Point3D r2(G4Point3D( halfLength, 0., 0.));
  scaleLine.push_back(r1);
  scaleLine.push_back(r2);
  G4Point3D ticky(0., tickLength, 0.);
  G4Point3D tickz(0., 0., tickLength);
  tick11.push_back(r1 + ticky);
  tick11.push_back(r1 - ticky);
  tick12.push_back(r1 + tickz);
  tick12.push_back(r1 - tickz);
  tick21.push_back(r2 + ticky);
  tick21.push_back(r2 - ticky);
  tick22.push_back(r2 + tickz);
  tick22.push_back(r2 - tickz);
  G4Point3D textPosition(0., tickLength, 0.);

  // Transform appropriately...

  G4Transform3D transformation;
  if (scale.GetAutoPlacing()) {
    G4Transform3D rotation;
    switch (scale.GetDirection()) {
    case G4Scale::x:
      break;
    case G4Scale::y:
      rotation = G4RotateZ3D(piBy2);
      break;
    case G4Scale::z:
      rotation = G4RotateY3D(piBy2);
      break;
    }
    G4double sxmid;
    G4double symid;
    G4double szmid;
    sxmid = xmin + oneMinusMargin * (xmax - xmin);
    symid = ymin + margin * (ymax - ymin);
    szmid = zmin + oneMinusMargin * (zmax - zmin);
    switch (scale.GetDirection()) {
    case G4Scale::x:
      sxmid -= halfLength;
      break;
    case G4Scale::y:
      symid += halfLength;
      break;
    case G4Scale::z:
      szmid -= halfLength;
      break;
    }
    G4Translate3D translation(sxmid, symid, szmid);
    transformation = translation * rotation;
  } else {
    if (fpModel) transformation = fpModel->GetTransformation();
  }

  // Draw...
  // We would like to call BeginPrimitives(transformation) here but
  // calling BeginPrimitives from within an AddPrimitive is not
  // allowed!  So we have to do our own transformation...
  AddPrimitive(scaleLine.transform(transformation));
  AddPrimitive(tick11.transform(transformation));
  AddPrimitive(tick12.transform(transformation));
  AddPrimitive(tick21.transform(transformation));
  AddPrimitive(tick22.transform(transformation));
  G4Text text(scale.GetAnnotation(),textPosition.transform(transformation));
  G4VisAttributes va(G4VVisCommand::GetCurrentTextColour());
  text.SetVisAttributes(va);
  text.SetScreenSize(scale.GetAnnotationSize());
  AddPrimitive(text);
}

void G4VSceneHandler::AddPrimitive (const G4Polymarker& polymarker) {
  switch (polymarker.GetMarkerType()) {
  default:
  case G4Polymarker::dots:
    {
      for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
	G4Circle dot (polymarker);
        dot.SetPosition (polymarker[iPoint]);
	dot.SetWorldSize  (0.);
	dot.SetScreenSize (0.1);  // Very small circle.
	AddPrimitive (dot);
      }
    }
    break;
  case G4Polymarker::circles:
    {
      for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
	G4Circle circle (polymarker);
	circle.SetPosition (polymarker[iPoint]);
	AddPrimitive (circle);
      }
    }
    break;
  case G4Polymarker::squares:
    {
      for (size_t iPoint = 0; iPoint < polymarker.size (); iPoint++) {
	G4Square square (polymarker);
	square.SetPosition (polymarker[iPoint]);
	AddPrimitive (square);
      }
    }
    break;
  }
}

void G4VSceneHandler::RemoveViewerFromList (G4VViewer* pViewer) {
  fViewerList.remove(pViewer);
}

void G4VSceneHandler::SetScene (G4Scene* pScene) {
  fpScene = pScene;
  // Notify all viewers that a kernel visit is required.
  G4ViewerListIterator i;
  for (i = fViewerList.begin(); i != fViewerList.end(); i++) {
    (*i) -> SetNeedKernelVisit (true);
  }
}

void G4VSceneHandler::RequestPrimitives (const G4VSolid& solid) {
  G4Polyhedron::SetNumberOfRotationSteps (GetNoOfSides (fpVisAttribs));
  G4Polyhedron* pPolyhedron = solid.GetPolyhedron ();
  G4Polyhedron::ResetNumberOfRotationSteps ();
  if (pPolyhedron) {
    pPolyhedron -> SetVisAttributes (fpVisAttribs);
    BeginPrimitives (fObjectTransformation);
    AddPrimitive (*pPolyhedron);
    EndPrimitives ();
  }
  else {
    G4VisManager::Verbosity verbosity = G4VisManager::GetVerbosity();
    if (verbosity >= G4VisManager::errors) {
      G4cerr <<
      "ERROR: G4VSceneHandler::RequestPrimitives"
      "\n  Polyhedron not available for " << solid.GetName () <<
      "\n  Touchable path: ";
      G4PhysicalVolumeModel* pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
      if (pPVModel) {
        G4cerr << pPVModel->GetFullPVPath();
      }
      G4cerr <<
      "\n  This means it cannot be visualized on most systems (try RayTracer)."
      "\n  1) The solid may not have implemented the CreatePolyhedron method."
      "\n  2) For Boolean solids, the BooleanProcessor, which attempts to create"
      "\n     the resultant polyhedron, may have failed."
      << G4endl;
    }
  }
}

void G4VSceneHandler::ProcessScene () {

  // Assumes graphics database store has already been cleared if
  // relevant for the particular scene handler.

  if (!fpScene) return;

  G4VisManager* visManager = G4VisManager::GetInstance();

  if (!visManager->GetConcreteInstance()) return;

  G4VisManager::Verbosity verbosity = visManager->GetVerbosity();

  fReadyForTransients = false;

  // Reset fMarkForClearingTransientStore. (Leaving
  // fMarkForClearingTransientStore true causes problems with
  // recomputing transients below.)  Restore it again at end...
  G4bool tmpMarkForClearingTransientStore = fMarkForClearingTransientStore;
  fMarkForClearingTransientStore = false;

  // Traverse geometry tree and send drawing primitives to window(s).

  const std::vector<G4Scene::Model>& runDurationModelList =
    fpScene -> GetRunDurationModelList ();

  if (runDurationModelList.size ()) {
    if (verbosity >= G4VisManager::confirmations) {
      G4cout << "Traversing scene data..." << G4endl;
    }

    BeginModeling ();

    // Create modeling parameters from view parameters...
    G4ModelingParameters* pMP = CreateModelingParameters ();

    for (size_t i = 0; i < runDurationModelList.size (); i++) {
      if (runDurationModelList[i].fActive) {
	G4VModel* pModel = runDurationModelList[i].fpModel;
	// Note: this is not the place to take action on
	// pModel->GetTransformation().  The model must take care of
	// this in pModel->DescribeYourselfTo(*this).  See, for example,
	// G4PhysicalVolumeModel and /vis/scene/add/logo.
	pModel -> SetModelingParameters (pMP);
	SetModel (pModel);  // Store for use by derived class.
	pModel -> DescribeYourselfTo (*this);
	pModel -> SetModelingParameters (0);
      }
    }

    delete pMP;
    EndModeling ();
  }

  fReadyForTransients = true;

  // Refresh event from end-of-event model list.
  // Allow only in Idle or GeomClosed state...
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState state = stateManager->GetCurrentState();
  if (state == G4State_Idle || state == G4State_GeomClosed) {

    visManager->SetEventRefreshing(true);

    if (visManager->GetRequestedEvent()) {
      DrawEvent(visManager->GetRequestedEvent());

    } else {

      G4RunManager* runManager = G4RunManager::GetRunManager();
#ifdef G4MULTITHREADED
      if(G4Threading::IsMultithreadedApplication())
      { runManager = G4MTRunManager::GetMasterRunManager(); }
#endif
      if (runManager) {
	const G4Run* run = runManager->GetCurrentRun();
        const std::vector<const G4Event*>* events =
	  run? run->GetEventVector(): 0;
	size_t nKeptEvents = 0;
	if (events) nKeptEvents = events->size();
	if (nKeptEvents) {

	  if (fpScene->GetRefreshAtEndOfEvent()) {

	    if (verbosity >= G4VisManager::confirmations) {
	      G4cout << "Refreshing event..." << G4endl;
	    }
	    const G4Event* event = 0;
	    if (events && events->size()) event = events->back();
	    if (event) DrawEvent(event);

	  } else {  // Accumulating events.

	    if (verbosity >= G4VisManager::confirmations) {
	      G4cout << "Refreshing events in run..." << G4endl;
	    }
            for (const auto& event: *events) {
              if (event) DrawEvent(event);
            }

	    if (!fpScene->GetRefreshAtEndOfRun()) {
	      if (verbosity >= G4VisManager::warnings) {
		G4cout <<
		  "WARNING: Cannot refresh events accumulated over more"
		  "\n  than one runs.  Refreshed just the last run."
		       << G4endl;
	      }
	    }
	  }
	}
      }
    }
    visManager->SetEventRefreshing(false);
  }

  // Refresh end-of-run model list.
  // Allow only in Idle or GeomClosed state...
  if (state == G4State_Idle || state == G4State_GeomClosed) {
    DrawEndOfRunModels();
  }

  fMarkForClearingTransientStore = tmpMarkForClearingTransientStore;
}

void G4VSceneHandler::DrawEvent(const G4Event* event)
{
  const std::vector<G4Scene::Model>& EOEModelList =
    fpScene -> GetEndOfEventModelList ();
  size_t nModels = EOEModelList.size();
  if (nModels) {
    G4ModelingParameters* pMP = CreateModelingParameters();
    pMP->SetEvent(event);
    for (size_t i = 0; i < nModels; i++) {
      if (EOEModelList[i].fActive) {
	G4VModel* pModel = EOEModelList[i].fpModel;
	pModel -> SetModelingParameters(pMP);
	SetModel (pModel);
	pModel -> DescribeYourselfTo (*this);
	pModel -> SetModelingParameters(0);
      }
    }
    delete pMP;
    SetModel (0);
  }
}

void G4VSceneHandler::DrawEndOfRunModels()
{
  const std::vector<G4Scene::Model>& EORModelList =
    fpScene -> GetEndOfRunModelList ();
  size_t nModels = EORModelList.size();
  if (nModels) {
    G4ModelingParameters* pMP = CreateModelingParameters();
    pMP->SetEvent(0);
    for (size_t i = 0; i < nModels; i++) {
      if (EORModelList[i].fActive) {
	G4VModel* pModel = EORModelList[i].fpModel;
	pModel -> SetModelingParameters(pMP);
	SetModel (pModel);
	pModel -> DescribeYourselfTo (*this);
	pModel -> SetModelingParameters(0);
      }
    }
    delete pMP;
    SetModel (0);
  }
}

G4ModelingParameters* G4VSceneHandler::CreateModelingParameters ()
{
  // Create modeling parameters from View Parameters...
  if (!fpViewer) return NULL;

  const G4ViewParameters& vp = fpViewer -> GetViewParameters ();

  // Convert drawing styles...
  G4ModelingParameters::DrawingStyle modelDrawingStyle =
    G4ModelingParameters::wf;
  switch (vp.GetDrawingStyle ()) {
  default:
  case G4ViewParameters::wireframe:
    modelDrawingStyle = G4ModelingParameters::wf;
    break;
  case G4ViewParameters::hlr:
    modelDrawingStyle = G4ModelingParameters::hlr;
    break;
  case G4ViewParameters::hsr:
    modelDrawingStyle = G4ModelingParameters::hsr;
    break;
  case G4ViewParameters::hlhsr:
    modelDrawingStyle = G4ModelingParameters::hlhsr;
    break;
  }

  // Decide if covered daughters are really to be culled...
  G4bool reallyCullCovered =
    vp.IsCullingCovered()   // Culling daughters depends also on...
    && !vp.IsSection ()     // Sections (DCUT) not requested.
    && !vp.IsCutaway ()     // Cutaways not requested.
    ;

  G4ModelingParameters* pModelingParams = new G4ModelingParameters
    (vp.GetDefaultVisAttributes (),
     modelDrawingStyle,
     vp.IsCulling (),
     vp.IsCullingInvisible (),
     vp.IsDensityCulling (),
     vp.GetVisibleDensity (),
     reallyCullCovered,
     vp.GetNoOfSides ()
     );

  pModelingParams->SetWarning
    (G4VisManager::GetVerbosity() >= G4VisManager::warnings);

  pModelingParams->SetExplodeFactor(vp.GetExplodeFactor());
  pModelingParams->SetExplodeCentre(vp.GetExplodeCentre());

  pModelingParams->SetSectionSolid(CreateSectionSolid());
  pModelingParams->SetCutawaySolid(CreateCutawaySolid());
  // The polyhedron objects are deleted in the modeling parameters destructor.
  
  pModelingParams->SetVisAttributesModifiers(vp.GetVisAttributesModifiers());

  return pModelingParams;
}

G4VSolid* G4VSceneHandler::CreateSectionSolid()
{
  G4VSolid* sectioner = 0;
  const G4ViewParameters& vp = fpViewer->GetViewParameters();
  if (vp.IsSection () ) {
    G4double radius = fpScene->GetExtent().GetExtentRadius();
    G4double safe = radius + fpScene->GetExtent().GetExtentCentre().mag();
    G4VSolid* sectionBox =
      new G4Box("_sectioner", safe, safe, 1.e-5 * radius);  // Thin in z-plane.
    const G4Plane3D& sp = vp.GetSectionPlane ();
    G4double a = sp.a();
    G4double b = sp.b();
    G4double c = sp.c();
    G4double d = sp.d();
    G4Transform3D transform = G4TranslateZ3D(-d);
    const G4Normal3D normal(a,b,c);
    if (normal != G4Normal3D(0,0,1)) {
      const G4double angle = std::acos(normal.dot(G4Normal3D(0,0,1)));
      const G4Vector3D axis = G4Normal3D(0,0,1).cross(normal);
      transform = G4Rotate3D(angle, axis) * transform;
    }
    sectioner = new G4DisplacedSolid
      ("_displaced_sectioning_box", sectionBox, transform);
  }
  return sectioner;
}

G4VSolid* G4VSceneHandler::CreateCutawaySolid()
{
  // To be reviewed.
  return 0;
  /*** An alternative way of getting a cutaway is to use
  Command /vis/scene/add/volume
  Guidance :
  Adds a physical volume to current scene, with optional clipping volume.
  If physical-volume-name is "world" (the default), the top of the
  main geometry tree (material world) is added.  If "worlds", the
  top of all worlds - material world and parallel worlds, if any - are
    added.  Otherwise a search of all worlds is made, taking the first
    matching occurence only.  To see a representation of the geometry
    hierarchy of the worlds, try "/vis/drawTree [worlds]" or one of the
    driver/browser combinations that have the required functionality, e.g., HepRep.
    If clip-volume-type is specified, the subsequent parameters are used to
    to define a clipping volume.  For example,
    "/vis/scene/add/volume ! ! ! -box km 0 1 0 1 0 1" will draw the world
    with the positive octant cut away.  (If the Boolean Processor issues
    warnings try replacing 0 by 0.000000001 or something.)
    If clip-volume-type is prepended with '-', the clip-volume is subtracted
    (cutaway). (This is the default if there is no prepended character.)
    If '*' is prepended, the intersection of the physical-volume and the
    clip-volume is made. (You can make a section/DCUT with a thin box, for
    example).
    For "box", the parameters are xmin,xmax,ymin,ymax,zmin,zmax.
    Only "box" is programmed at present.
   ***/
}

void G4VSceneHandler::LoadAtts(const G4Visible& visible, G4AttHolder* holder)
{
  // Load G4Atts from G4VisAttributes, if any...
  const G4VisAttributes* va = visible.GetVisAttributes();
  if (va) {
    const std::map<G4String,G4AttDef>* vaDefs =
      va->GetAttDefs();
    if (vaDefs) {
      holder->AddAtts(visible.GetVisAttributes()->CreateAttValues(), vaDefs);
    }
  }

  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (pPVModel) {
    // Load G4Atts from G4PhysicalVolumeModel...
    const std::map<G4String,G4AttDef>* pvDefs = pPVModel->GetAttDefs();
    if (pvDefs) {
      holder->AddAtts(pPVModel->CreateCurrentAttValues(), pvDefs);
    }
  }

  G4TrajectoriesModel* trajModel = dynamic_cast<G4TrajectoriesModel*>(fpModel);
  if (trajModel) {
    // Load G4Atts from trajectory model...
    const std::map<G4String,G4AttDef>* trajModelDefs = trajModel->GetAttDefs();
    if (trajModelDefs) {
      holder->AddAtts(trajModel->CreateCurrentAttValues(), trajModelDefs);
    }
    // Load G4Atts from trajectory...
    const G4VTrajectory* traj = trajModel->GetCurrentTrajectory();
    if (traj) {
      const std::map<G4String,G4AttDef>* trajDefs = traj->GetAttDefs();
      if (trajDefs) {
        holder->AddAtts(traj->CreateAttValues(), trajDefs);
      }
      G4int nPoints = traj->GetPointEntries();
      for (G4int i = 0; i < nPoints; ++i) {
        G4VTrajectoryPoint* trajPoint = traj->GetPoint(i);
        if (trajPoint) {
          const std::map<G4String,G4AttDef>* pointDefs = trajPoint->GetAttDefs();
          if (pointDefs) {
            holder->AddAtts(trajPoint->CreateAttValues(), pointDefs);
          }
        }
      }
    }
  }

  G4HitsModel* hitsModel = dynamic_cast<G4HitsModel*>(fpModel);
  if (hitsModel) {
    // Load G4Atts from hit...
    const G4VHit* hit = hitsModel->GetCurrentHit();
    const std::map<G4String,G4AttDef>* hitsDefs = hit->GetAttDefs();
    if (hitsDefs) {
      holder->AddAtts(hit->CreateAttValues(), hitsDefs);
    }
  }
}

const G4Colour& G4VSceneHandler::GetTextColour (const G4Text& text) {
  const G4VisAttributes* pVA = text.GetVisAttributes ();
  if (!pVA) {
    return G4VVisCommand::GetCurrentTextColour();
  }
  const G4Colour& colour = pVA -> GetColour ();
  return colour;
}

G4double G4VSceneHandler::GetLineWidth(const G4VisAttributes* pVisAttribs)
{
  G4double lineWidth = pVisAttribs->GetLineWidth();
  if (lineWidth < 1.) lineWidth = 1.;
  lineWidth *= fpViewer -> GetViewParameters().GetGlobalLineWidthScale();
  if (lineWidth < 1.) lineWidth = 1.;
  return lineWidth;
}

G4ViewParameters::DrawingStyle G4VSceneHandler::GetDrawingStyle
(const G4VisAttributes* pVisAttribs) {
  // Drawing style is normally determined by the view parameters, but
  // it can be overriddden by the ForceDrawingStyle flag in the vis
  // attributes.
  G4ViewParameters::DrawingStyle style = 
    fpViewer->GetViewParameters().GetDrawingStyle();
  if (pVisAttribs -> IsForceDrawingStyle ()) {
    G4VisAttributes::ForcedDrawingStyle forcedStyle =
      pVisAttribs -> GetForcedDrawingStyle ();
    // This is complicated because if hidden line and surface removal
    // has been requested we wish to preserve this sometimes.
    switch (forcedStyle) {
    case (G4VisAttributes::solid):
      switch (style) {
      case (G4ViewParameters::hlr):
	style = G4ViewParameters::hlhsr;
	break;
      case (G4ViewParameters::wireframe):
	style = G4ViewParameters::hsr;
	break;
      case (G4ViewParameters::hlhsr):
      case (G4ViewParameters::hsr):
      default:
	break;
      }	
      break;
    case (G4VisAttributes::wireframe):
    default:
      // But if forced style is wireframe, do it, because one of its
      // main uses is in displaying the consituent solids of a Boolean
      // solid and their surfaces overlap with the resulting Booean
      // solid, making a mess if hlr is specified.
      style = G4ViewParameters::wireframe;
      break;
    }
  }
  return style;
}

G4bool G4VSceneHandler::GetAuxEdgeVisible (const G4VisAttributes* pVisAttribs) {
  G4bool isAuxEdgeVisible = fpViewer->GetViewParameters().IsAuxEdgeVisible ();
  if (pVisAttribs -> IsForceAuxEdgeVisible()) {
    isAuxEdgeVisible = pVisAttribs->IsForcedAuxEdgeVisible();
  }
  return isAuxEdgeVisible;
}

G4double G4VSceneHandler::GetMarkerSize
(const G4VMarker& marker, 
 G4VSceneHandler::MarkerSizeType& markerSizeType)
{
  G4bool userSpecified = marker.GetWorldSize() || marker.GetScreenSize();
  const G4VMarker& defaultMarker =
    fpViewer -> GetViewParameters().GetDefaultMarker();
  G4double size = userSpecified ?
    marker.GetWorldSize() : defaultMarker.GetWorldSize();
  if (size) {
    // Draw in world coordinates.
    markerSizeType = world;
  }
  else {
    size = userSpecified ?
      marker.GetScreenSize() : defaultMarker.GetScreenSize();
    // Draw in screen coordinates.
    markerSizeType = screen;
  }
  size *= fpViewer -> GetViewParameters().GetGlobalMarkerScale();
  if (markerSizeType == screen && size < 1.) size = 1.;
  return size;
}

G4int G4VSceneHandler::GetNoOfSides(const G4VisAttributes* pVisAttribs)
{
  // No. of sides (lines segments per circle) is normally determined
  // by the view parameters, but it can be overriddden by the
  // ForceLineSegmentsPerCircle in the vis attributes.
  G4int lineSegmentsPerCircle = fpViewer->GetViewParameters().GetNoOfSides();
  if (pVisAttribs) {
    if (pVisAttribs->IsForceLineSegmentsPerCircle())
      lineSegmentsPerCircle = pVisAttribs->GetForcedLineSegmentsPerCircle();
    if (lineSegmentsPerCircle < pVisAttribs->GetMinLineSegmentsPerCircle()) {
      lineSegmentsPerCircle = pVisAttribs->GetMinLineSegmentsPerCircle();
      G4cout <<
	"G4VSceneHandler::GetNoOfSides: attempt to set the"
	"\nnumber of line segements per circle < " << lineSegmentsPerCircle
	     << "; forced to " << pVisAttribs->GetMinLineSegmentsPerCircle() << G4endl;
    }
  }
  return lineSegmentsPerCircle;
}

std::ostream& operator << (std::ostream& os, const G4VSceneHandler& sh) {

  os << "Scene handler " << sh.fName << " has "
     << sh.fViewerList.size () << " viewer(s):";
  for (size_t i = 0; i < sh.fViewerList.size (); i++) {
    os << "\n  " << *(sh.fViewerList [i]);
  }

  if (sh.fpScene) {
    os << "\n  " << *sh.fpScene;
  }
  else {
    os << "\n  This scene handler currently has no scene.";
  }

  return os;
}

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
#include "G4Tet.hh"
#include "G4DisplacedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ModelingParameters.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4HitsModel.hh"
#include "G4VHit.hh"
#include "G4VDigi.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"
#include "G4Mesh.hh"
#include "G4DefaultLinearColorMap.hh"
#include "G4QuickRand.hh"
#include "G4StateManager.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Run.hh"
#include "G4Transform3D.hh"
#include "G4AttHolder.hh"
#include "G4AttDef.hh"
#include "G4VVisCommand.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <set>

#define G4warn G4cout

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

void G4VSceneHandler::AddSolid (const G4TessellatedSolid& tess) {
  AddSolidT (tess);
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
  using MeshScoreMap = G4VScoringMesh::MeshScoreMap;
  //G4cout << "AddCompound: hits: " << &hits << G4endl;
  G4bool scoreMapHits = false;
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManagerIfExist();
  if (scoringManager) {
    std::size_t nMeshes = scoringManager->GetNumberOfMesh();
    for (std::size_t iMesh = 0; iMesh < nMeshes; ++iMesh) {
      G4VScoringMesh* mesh = scoringManager->GetMesh((G4int)iMesh);
      if (mesh && mesh->IsActive()) {
	MeshScoreMap scoreMap = mesh->GetScoreMap();
        const G4String& mapNam = const_cast<G4THitsMap<G4double>&>(hits).GetName();
	for(MeshScoreMap::const_iterator i = scoreMap.cbegin();
	    i != scoreMap.cend(); ++i) {
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
  using MeshScoreMap = G4VScoringMesh::MeshScoreMap;
  //G4cout << "AddCompound: hits: " << &hits << G4endl;
  G4bool scoreMapHits = false;
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManagerIfExist();
  if (scoringManager) {
    std::size_t nMeshes = scoringManager->GetNumberOfMesh();
    for (std::size_t iMesh = 0; iMesh < nMeshes; ++iMesh) {
      G4VScoringMesh* mesh = scoringManager->GetMesh((G4int)iMesh);
      if (mesh && mesh->IsActive()) {
	MeshScoreMap scoreMap = mesh->GetScoreMap();
	for(MeshScoreMap::const_iterator i = scoreMap.cbegin();
	    i != scoreMap.cend(); ++i) {
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

void G4VSceneHandler::AddCompound(const G4Mesh& mesh)
{
  G4warn <<
  "There has been an attempt to draw a mesh with option \""
  << fpViewer->GetViewParameters().GetSpecialMeshRenderingOption()
  << "\":\n" << mesh
  << "but it is not of a recognised type or is not implemented"
  "\nby the current graphics driver. Instead we draw its"
  "\ncontainer \"" << mesh.GetContainerVolume()->GetName() << "\"."
  << G4endl;
  const auto& pv = mesh.GetContainerVolume();
  const auto& lv = pv->GetLogicalVolume();
  const auto& solid = lv->GetSolid();
  const auto& transform = mesh.GetTransform();
  // Make sure container is visible
  G4VisAttributes tmpVisAtts;  // Visible, white, not forced.
  const auto& saveVisAtts = lv->GetVisAttributes();
  if (saveVisAtts) {
    tmpVisAtts = *saveVisAtts;
    tmpVisAtts.SetVisibility(true);
    auto colour = saveVisAtts->GetColour();
    colour.SetAlpha(1.);
    tmpVisAtts.SetColour(colour);
  }
  // Draw container
  PreAddSolid(transform,tmpVisAtts);
  solid->DescribeYourselfTo(*this);
  PostAddSolid();
  // Restore vis attributes
  lv->SetVisAttributes(saveVisAtts);
}

void G4VSceneHandler::AddViewerToList (G4VViewer* pViewer) {
  fViewerList.push_back (pViewer);
}

void G4VSceneHandler::AddPrimitive (const G4Polymarker& polymarker) {
  switch (polymarker.GetMarkerType()) {
    default:
    case G4Polymarker::dots:
    {
      G4Circle dot (polymarker);
      dot.SetWorldSize  (0.);
      dot.SetScreenSize (0.1);  // Very small circle.
      for (std::size_t iPoint = 0; iPoint < polymarker.size (); ++iPoint) {
        dot.SetPosition (polymarker[iPoint]);
        AddPrimitive (dot);
      }
    }
      break;
    case G4Polymarker::circles:
    {
      G4Circle circle (polymarker);  // Default circle
      for (std::size_t iPoint = 0; iPoint < polymarker.size (); ++iPoint) {
        circle.SetPosition (polymarker[iPoint]);
        AddPrimitive (circle);
      }
    }
      break;
    case G4Polymarker::squares:
    {
      G4Square square (polymarker);  // Default square
      for (std::size_t iPoint = 0; iPoint < polymarker.size (); ++iPoint) {
        square.SetPosition (polymarker[iPoint]);
        AddPrimitive (square);
      }
    }
      break;
  }
}

void G4VSceneHandler::RemoveViewerFromList (G4VViewer* pViewer) {
  fViewerList.remove(pViewer);  // Does nothing if already removed
  // And reset current viewer
  auto visManager = G4VisManager::GetInstance();
  visManager->SetCurrentViewer(nullptr);
}


void G4VSceneHandler::AddPrimitive (const G4Plotter&) {
  G4warn << "WARNING: Plotter not implemented for " << fSystem.GetName() << G4endl;
  G4warn << "  Open a plotter-aware graphics system or remove plotter with" << G4endl;
  G4warn << "  /vis/scene/removeModel Plotter" << G4endl;
}

void G4VSceneHandler::SetScene (G4Scene* pScene) {
  fpScene = pScene;
  // Notify all viewers that a kernel visit is required.
  G4ViewerListIterator i;
  for (i = fViewerList.begin(); i != fViewerList.end(); i++) {
    (*i) -> SetNeedKernelVisit (true);
  }
}

void G4VSceneHandler::RequestPrimitives (const G4VSolid& solid)
{
  // Sometimes solids that have no substance get requested. They may
  // be part of the geometry tree but have been "spirited away", for
  // example by a Boolean subtraction in wich the original volume
  // is entirely inside the subtractor.
  // The problem is that the Boolean Processor still returns a
  // polyhedron in these cases (IMHO it should not), so the
  // workaround is to return before the damage is done.
  auto pSolid = &solid;
  auto pBooleanSolid = dynamic_cast<const G4BooleanSolid*>(pSolid);
  if (pBooleanSolid) {
    G4ThreeVector bmin, bmax;
    pBooleanSolid->BoundingLimits(bmin, bmax);
    G4bool isGood = false;
    for (G4int i=0; i<100000; ++i) {
      G4double x = bmin.x() + (bmax.x() - bmin.x())*G4QuickRand();
      G4double y = bmin.y() + (bmax.y() - bmin.y())*G4QuickRand();
      G4double z = bmin.z() + (bmax.z() - bmin.z())*G4QuickRand();
      if (pBooleanSolid->Inside(G4ThreeVector(x,y,z)) == kInside) {
        isGood = true;
        break;
      }
    }
    if (!isGood) return;
  }
  
  const G4ViewParameters::DrawingStyle style = GetDrawingStyle(fpVisAttribs);
  const G4ViewParameters& vp = fpViewer->GetViewParameters();

  switch (style) {
    default:
    case G4ViewParameters::wireframe:
    case G4ViewParameters::hlr:
    case G4ViewParameters::hsr:
    case G4ViewParameters::hlhsr:
    {
      // Use polyhedral representation
      G4Polyhedron::SetNumberOfRotationSteps (GetNoOfSides (fpVisAttribs));
      G4Polyhedron* pPolyhedron = solid.GetPolyhedron ();
      G4Polyhedron::ResetNumberOfRotationSteps ();
      if (pPolyhedron) {
        pPolyhedron -> SetVisAttributes (fpVisAttribs);
        BeginPrimitives (fObjectTransformation);
        AddPrimitive (*pPolyhedron);
        EndPrimitives ();
        break;
      } else {  // Print warnings and drop through to cloud
        G4VisManager::Verbosity verbosity = G4VisManager::GetVerbosity();
        static std::set<const G4VSolid*> problematicSolids;
        if (verbosity >= G4VisManager::errors &&
            problematicSolids.find(&solid) == problematicSolids.end()) {
          problematicSolids.insert(&solid);
          G4warn <<
          "ERROR: G4VSceneHandler::RequestPrimitives"
          "\n  Polyhedron not available for " << solid.GetName ();
          G4PhysicalVolumeModel* pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
          if (pPVModel) {
            G4warn << "\n  Touchable path: " << pPVModel->GetFullPVPath();
          }
          static G4bool explanation = false;
          if (!explanation) {
            explanation = true;
            G4warn <<
            "\n  This means it cannot be visualized in the usual way on most systems."
            "\n  1) The solid may not have implemented the CreatePolyhedron method."
            "\n  2) For Boolean solids, the BooleanProcessor, which attempts to create"
            "\n     the resultant polyhedron, may have failed."
            "\n  Try RayTracer. It uses Geant4's tracking algorithms instead.";
          }
          G4warn << "\n  Drawing solid with cloud of points.";
          G4warn << G4endl;
        }
      }
    }
      [[fallthrough]];

    case G4ViewParameters::cloud:
    {
      // Form solid out of cloud of dots on surface of solid
      G4Polymarker dots;
      // Note: OpenGL has a fast implementation of polymarker so it's better
      // to build a polymarker rather than add a succession of circles.
      // And anyway, in Qt, in the latter case each circle would be a scene-tree
      // entry, something we would want to avoid.
      dots.SetVisAttributes(fpVisAttribs);
      dots.SetMarkerType(G4Polymarker::dots);
      dots.SetSize(G4VMarker::screen,1.);
      G4int numberOfCloudPoints = GetNumberOfCloudPoints(fpVisAttribs);
      if (numberOfCloudPoints <= 0) numberOfCloudPoints = vp.GetNumberOfCloudPoints();
      for (G4int i = 0; i < numberOfCloudPoints; ++i) {
	G4ThreeVector p = solid.GetPointOnSurface();
	dots.push_back(p);
      }
      BeginPrimitives (fObjectTransformation);
      AddPrimitive(dots);
      EndPrimitives ();
      break;
    }
  }
}

//namespace {
//  void DrawExtent(const G4VModel* pModel)
//  {
//    // Show extent boxes - debug only, OGLSX only (OGLSQt problem?)
//    if (pModel->GetExtent() != G4VisExtent::GetNullExtent()) {
//      const auto& extent = pModel->GetExtent();
//      const auto& centre = extent.GetExtentCenter();
//      const auto& position = G4Translate3D(centre);
//      const auto& dx = (extent.GetXmax()-extent.GetXmin())/2.;
//      const auto& dy = (extent.GetYmax()-extent.GetYmin())/2.;
//      const auto& dz = (extent.GetZmax()-extent.GetZmin())/2.;
//      auto visAtts = G4VisAttributes();
//      visAtts.SetForceWireframe();
//      G4Box extentBox("Extent",dx,dy,dz);
//      G4VisManager::GetInstance()->Draw(extentBox,visAtts,position);
//    }
//  }
//}

void G4VSceneHandler::ProcessScene()
{
  // Assumes graphics database store has already been cleared if
  // relevant for the particular scene handler.

  if(!fpScene)
    return;

  if(fpScene->GetExtent() == G4VisExtent::GetNullExtent())
  {
    G4Exception("G4VSceneHandler::ProcessScene", "visman0106", JustWarning,
                "The scene has no extent.");
  }

  G4VisManager* visManager = G4VisManager::GetInstance();

  if(!visManager->GetConcreteInstance())
    return;

  G4VisManager::Verbosity verbosity = visManager->GetVerbosity();

  fReadyForTransients = false;

  // Reset fMarkForClearingTransientStore. (Leaving
  // fMarkForClearingTransientStore true causes problems with
  // recomputing transients below.)  Restore it again at end...
  G4bool tmpMarkForClearingTransientStore = fMarkForClearingTransientStore;
  fMarkForClearingTransientStore          = false;

  // Traverse geometry tree and send drawing primitives to window(s).

  const std::vector<G4Scene::Model>& runDurationModelList =
    fpScene->GetRunDurationModelList();

  if(runDurationModelList.size())
  {
    if(verbosity >= G4VisManager::confirmations)
    {
      G4cout << "Traversing scene data..." << G4endl;
    }

    BeginModeling();

    // Create modeling parameters from view parameters...
    G4ModelingParameters* pMP = CreateModelingParameters();

    for(std::size_t i = 0; i < runDurationModelList.size(); ++i)
    {
      if(runDurationModelList[i].fActive)
      {
	fpModel = runDurationModelList[i].fpModel;
	fpModel->SetModelingParameters(pMP);
	fpModel->DescribeYourselfTo(*this);
	// To see the extents of each model represented as wireframe boxes,
	// uncomment the next line and DrawExtent in namespace above
	// DrawExtent(fpModel);
	fpModel->SetModelingParameters(0);
      }
    }

    fpModel = 0;
    delete pMP;

    EndModeling();
  }

  fReadyForTransients = true;

  // Refresh event from end-of-event model list.
  // Allow only in Idle or GeomClosed state...
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState state     = stateManager->GetCurrentState();
  if(state == G4State_Idle || state == G4State_GeomClosed)
  {
    visManager->SetEventRefreshing(true);

    if(visManager->GetRequestedEvent())
    {
      DrawEvent(visManager->GetRequestedEvent());
    }
    else
    {
      G4RunManager* runManager = G4RunManagerFactory::GetMasterRunManager();
      if(runManager)
      {
        const G4Run* run = runManager->GetCurrentRun();
        const std::vector<const G4Event*>* events =
          run ? run->GetEventVector() : 0;
        std::size_t nKeptEvents = 0;
        if(events)
          nKeptEvents = events->size();
        if(nKeptEvents)
        {
          if(fpScene->GetRefreshAtEndOfEvent())
          {
            if(verbosity >= G4VisManager::confirmations)
            {
              G4cout << "Refreshing event..." << G4endl;
            }
            const G4Event* event = 0;
            if(events && events->size())
              event = events->back();
            if(event)
              DrawEvent(event);
          }
          else
          {  // Accumulating events.

            if(verbosity >= G4VisManager::confirmations)
            {
              G4cout << "Refreshing events in run..." << G4endl;
            }
            for(const auto& event : *events)
            {
              if(event)
                DrawEvent(event);
            }

            if(!fpScene->GetRefreshAtEndOfRun())
            {
              if(verbosity >= G4VisManager::warnings)
              {
                G4warn << "WARNING: Cannot refresh events accumulated over more"
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
  if(state == G4State_Idle || state == G4State_GeomClosed)
  {
    DrawEndOfRunModels();
  }

  fMarkForClearingTransientStore = tmpMarkForClearingTransientStore;
}

void G4VSceneHandler::DrawEvent(const G4Event* event)
{
  const std::vector<G4Scene::Model>& EOEModelList =
    fpScene -> GetEndOfEventModelList ();
  std::size_t nModels = EOEModelList.size();
  if (nModels) {
    G4ModelingParameters* pMP = CreateModelingParameters();
    pMP->SetEvent(event);
    for (std::size_t i = 0; i < nModels; ++i) {
      if (EOEModelList[i].fActive) {
	fpModel = EOEModelList[i].fpModel;
	fpModel -> SetModelingParameters(pMP);
	fpModel -> DescribeYourselfTo (*this);
	fpModel -> SetModelingParameters(0);
      }
    }
    fpModel = 0;
    delete pMP;
  }
}

void G4VSceneHandler::DrawEndOfRunModels()
{
  const std::vector<G4Scene::Model>& EORModelList =
    fpScene -> GetEndOfRunModelList ();
  std::size_t nModels = EORModelList.size();
  if (nModels) {
    G4ModelingParameters* pMP = CreateModelingParameters();
    pMP->SetEvent(0);
    for (std::size_t i = 0; i < nModels; ++i) {
      if (EORModelList[i].fActive) {
        fpModel = EORModelList[i].fpModel;
	fpModel -> SetModelingParameters(pMP);
	fpModel -> DescribeYourselfTo (*this);
	fpModel -> SetModelingParameters(0);
      }
    }
    fpModel = 0;
    delete pMP;
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
    case G4ViewParameters::cloud:
      modelDrawingStyle = G4ModelingParameters::cloud;
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

  pModelingParams->SetNumberOfCloudPoints(vp.GetNumberOfCloudPoints());
  pModelingParams->SetWarning
    (G4VisManager::GetVerbosity() >= G4VisManager::warnings);

  pModelingParams->SetCBDAlgorithmNumber(vp.GetCBDAlgorithmNumber());
  pModelingParams->SetCBDParameters(vp.GetCBDParameters());

  pModelingParams->SetExplodeFactor(vp.GetExplodeFactor());
  pModelingParams->SetExplodeCentre(vp.GetExplodeCentre());

  pModelingParams->SetSectionSolid(CreateSectionSolid());
  pModelingParams->SetCutawaySolid(CreateCutawaySolid());
  // The polyhedron objects are deleted in the modeling parameters destructor.
  
  pModelingParams->SetVisAttributesModifiers(vp.GetVisAttributesModifiers());

  pModelingParams->SetSpecialMeshRendering(vp.IsSpecialMeshRendering());
  pModelingParams->SetSpecialMeshVolumes(vp.GetSpecialMeshVolumes());

  return pModelingParams;
}

G4DisplacedSolid* G4VSceneHandler::CreateSectionSolid()
{
  G4DisplacedSolid* sectioner = 0;

  const G4ViewParameters& vp = fpViewer->GetViewParameters();
  if (vp.IsSection () ) {

    G4double radius = fpScene->GetExtent().GetExtentRadius();
    G4double safe = radius + fpScene->GetExtent().GetExtentCentre().mag();
    G4VSolid* sectionBox =
      new G4Box("_sectioner", safe, safe, 1.e-5 * radius);  // Thin in z-plane...
    const G4Normal3D originalNormal(0,0,1);  // ...so this is original normal.

    const G4Plane3D& sp = vp.GetSectionPlane ();
    const G4double& a = sp.a();
    const G4double& b = sp.b();
    const G4double& c = sp.c();
    const G4double& d = sp.d();
    const G4Normal3D newNormal(a,b,c);

    G4Transform3D requiredTransform;
    // Rotate
    if (newNormal != originalNormal) {
      const G4double& angle = std::acos(newNormal.dot(originalNormal));
      const G4Vector3D& axis = originalNormal.cross(newNormal);
      requiredTransform = G4Rotate3D(angle, axis);
    }
    // Translate
    requiredTransform = requiredTransform * G4TranslateZ3D(-d);

    sectioner = new G4DisplacedSolid
      ("_displaced_sectioning_box", sectionBox, requiredTransform);
  }
  
  return sectioner;
}

G4DisplacedSolid* G4VSceneHandler::CreateCutawaySolid()
{
  const G4ViewParameters& vp = fpViewer->GetViewParameters();
  if (vp.IsCutaway()) {

    std::vector<G4DisplacedSolid*> cutaway_solids;

    G4double radius = fpScene->GetExtent().GetExtentRadius();
    G4double safe = radius + fpScene->GetExtent().GetExtentCentre().mag();
    G4VSolid* cutawayBox =
    new G4Box("_cutaway_box", safe, safe, safe);  // world box...

    for (int plane_no = 0; plane_no < int(vp.GetCutawayPlanes().size()); plane_no++){

      const G4Normal3D originalNormal(0,0,1);  // ...so this is original normal.

      const G4Plane3D& sp = vp.GetCutawayPlanes()[plane_no]; //];
      const G4double& a = sp.a();
      const G4double& b = sp.b();
      const G4double& c = sp.c();
      const G4double& d = sp.d();
      const G4Normal3D newNormal(-a,-b,-c);  // Convention: keep a*x+b*y+c*z+d>=0
      // Not easy to see why the above gives the right convention, but it has been
      // arrived at by trial and error to agree with the OpenGL implementation
      // of clipping planes.

      G4Transform3D requiredTransform;  // Null transform
      // Calculate the rotation
      // If newNormal is (0,0,1), no need to do anything
      // Treat (0,0,-1) as a special case, since cannot define axis in this case
      if (newNormal == G4Normal3D(0,0,-1)) {
        requiredTransform = G4Rotate3D(pi,G4Vector3D(1,0,0));
      } else if (newNormal != originalNormal) {
        const G4double& angle = std::acos(newNormal.dot(originalNormal));
        const G4Vector3D& axis = originalNormal.cross(newNormal);
        requiredTransform = G4Rotate3D(angle, axis);
      }
      // Translation
      requiredTransform = requiredTransform * G4TranslateZ3D(d + safe);
      cutaway_solids.push_back
      (new G4DisplacedSolid("_displaced_cutaway_box", cutawayBox, requiredTransform));
    }

    if (cutaway_solids.size() == 1){
      return (G4DisplacedSolid*) cutaway_solids[0];
    } else if (vp.GetCutawayMode() == G4ViewParameters::cutawayUnion) {
      G4UnionSolid* union2 =
      new G4UnionSolid("_union_2", cutaway_solids[0], cutaway_solids[1]);
      if (cutaway_solids.size() == 2)
        return (G4DisplacedSolid*)union2;
      else
        return (G4DisplacedSolid*)
        new G4UnionSolid("_union_3", union2, cutaway_solids[2]);
    } else if (vp.GetCutawayMode() == G4ViewParameters::cutawayIntersection){
      G4IntersectionSolid* intersection2 =
      new G4IntersectionSolid("_intersection_2", cutaway_solids[0], cutaway_solids[1]);
      if (cutaway_solids.size() == 2)
        return (G4DisplacedSolid*)intersection2;
      else
        return (G4DisplacedSolid*)
        new G4IntersectionSolid("_intersection_3", intersection2, cutaway_solids[2]);
    }
  }

  return 0;
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

const G4Colour& G4VSceneHandler::GetColour () {
  fpVisAttribs = fpViewer->GetApplicableVisAttributes(fpVisAttribs);
  const G4Colour& colour = fpVisAttribs -> GetColour ();
  return colour;
}

const G4Colour& G4VSceneHandler::GetColour (const G4Visible& visible) {
  auto pVA = visible.GetVisAttributes();
  if (!pVA) pVA = fpViewer->GetViewParameters().GetDefaultVisAttributes();
  return pVA->GetColour();
}

const G4Colour& G4VSceneHandler::GetTextColour (const G4Text& text) {
  auto pVA = text.GetVisAttributes();
  if (!pVA) pVA = fpViewer->GetViewParameters().GetDefaultTextVisAttributes();
  return pVA->GetColour();
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
  const G4ViewParameters& vp = fpViewer->GetViewParameters();
  const G4ViewParameters::DrawingStyle viewerStyle = vp.GetDrawingStyle();
  G4ViewParameters::DrawingStyle resultantStyle = viewerStyle;
  if (pVisAttribs -> IsForceDrawingStyle ()) {
    G4VisAttributes::ForcedDrawingStyle forcedStyle =
    pVisAttribs -> GetForcedDrawingStyle ();
    // This is complicated because if hidden line and surface removal
    // has been requested we wish to preserve this sometimes.
    switch (forcedStyle) {
      case (G4VisAttributes::solid):
        switch (viewerStyle) {
          case (G4ViewParameters::hlr):
            resultantStyle = G4ViewParameters::hlhsr;
            break;
          case (G4ViewParameters::wireframe):
            resultantStyle = G4ViewParameters::hsr;
            break;
          case (G4ViewParameters::cloud):
            resultantStyle = G4ViewParameters::hsr;
            break;
          case (G4ViewParameters::hlhsr):
          case (G4ViewParameters::hsr):
            break;
        }
        break;
      case (G4VisAttributes::cloud):
        resultantStyle = G4ViewParameters::cloud;
        break;
      case (G4VisAttributes::wireframe):
      default:
        // But if forced style is wireframe, do it, because one of its
        // main uses is in displaying the consituent solids of a Boolean
        // solid and their surfaces overlap with the resulting Booean
        // solid, making a mess if hlr is specified.
        resultantStyle = G4ViewParameters::wireframe;
        break;
    }
  }
  return resultantStyle;
}

G4int G4VSceneHandler::GetNumberOfCloudPoints
(const G4VisAttributes* pVisAttribs) const {
  // Returns no of cloud points from current view parameters, unless the user
  // has forced through the vis attributes, thereby over-riding the
  // current view parameter.
  G4int numberOfCloudPoints = fpViewer->GetViewParameters().GetNumberOfCloudPoints();
  if (pVisAttribs -> IsForceDrawingStyle() &&
      pVisAttribs -> GetForcedDrawingStyle() == G4VisAttributes::cloud &&
      pVisAttribs -> GetForcedNumberOfCloudPoints() > 0) {
    numberOfCloudPoints = pVisAttribs -> GetForcedNumberOfCloudPoints();
  }
  return numberOfCloudPoints;
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
      G4warn <<
	"G4VSceneHandler::GetNoOfSides: attempt to set the"
	"\nnumber of line segments per circle < " << lineSegmentsPerCircle
	     << "; forced to " << pVisAttribs->GetMinLineSegmentsPerCircle() << G4endl;
    }
  }
  return lineSegmentsPerCircle;
}

std::ostream& operator << (std::ostream& os, const G4VSceneHandler& sh) {

  os << "Scene handler " << sh.fName << " has "
     << sh.fViewerList.size () << " viewer(s):";
  for (std::size_t i = 0; i < sh.fViewerList.size (); ++i) {
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

void G4VSceneHandler::PseudoSceneFor3DRectMeshPositions::AddSolid(const G4Box&) {
  if (fpPVModel->GetCurrentDepth() == fDepth) {  // Leaf-level cells only
    const auto& material = fpPVModel->GetCurrentLV()->GetMaterial();
    const auto& name = material->GetName();
    const auto* pVisAtts = fpPVModel->GetCurrentLV()->GetVisAttributes();
    // Get position in world coordinates
    // As a parameterisation the box is transformed by the current transformation
    // and its centre, originally by definition at (0,0,0), is now translated.
    const G4ThreeVector& position = fpCurrentObjectTransformation->getTranslation();
    fPositionByMaterial.insert(std::make_pair(material,position));
    if (fNameAndVisAttsByMaterial.find(material) == fNameAndVisAttsByMaterial.end())
      // Store name and vis attributes of first encounter with this material
      fNameAndVisAttsByMaterial[material] = NameAndVisAtts(name,*pVisAtts);
  }
}

void G4VSceneHandler::PseudoSceneForTetVertices::AddSolid(const G4VSolid& solid) {
  if (fpPVModel->GetCurrentDepth() == fDepth) {  // Leaf-level cells only
    // Need to know it's a tet !!!! or implement G4VSceneHandler::AddSolid (const G4Tet&) !!!!
    try {
      const G4Tet& tet = dynamic_cast<const G4Tet&>(solid);
      const auto& material = fpPVModel->GetCurrentLV()->GetMaterial();
      const auto& name = material->GetName();
      const auto* pVisAtts = fpPVModel->GetCurrentLV()->GetVisAttributes();
      // Transform into world coordinates if necessary
      if (fpCurrentObjectTransformation->xx() == 1. &&
          fpCurrentObjectTransformation->yy() == 1. &&
          fpCurrentObjectTransformation->zz() == 1.) { // No transformation necessary
        const auto& vertices = tet.GetVertices();
        fVerticesByMaterial.insert(std::make_pair(material,vertices));
      } else {
        auto vertices = tet.GetVertices();
        for (auto&& vertex: vertices) {
          vertex = G4Point3D(vertex).transform(*fpCurrentObjectTransformation);
        }
        fVerticesByMaterial.insert(std::make_pair(material,vertices));
      }
      if (fNameAndVisAttsByMaterial.find(material) == fNameAndVisAttsByMaterial.end())
        // Store name and vis attributes of first encounter with this material
        fNameAndVisAttsByMaterial[material] = NameAndVisAtts(name,*pVisAtts);
    }
    catch (const std::bad_cast&) {
      G4ExceptionDescription ed;
      ed << "Called for a mesh that is not a tetrahedron mesh: " << solid.GetName();
      G4Exception("PseudoSceneForTetVertices","visman0108",JustWarning,ed);
    }
  }
}

void G4VSceneHandler::StandardSpecialMeshRendering(const G4Mesh& mesh)
// Standard way of special mesh rendering.
// MySceneHandler::AddCompound(const G4Mesh& mesh) may use this if
// appropriate or implement its own special mesh rendereing.
{
  G4bool implemented = false;
  switch (mesh.GetMeshType()) {
    case G4Mesh::rectangle: [[fallthrough]];
    case G4Mesh::nested3DRectangular:
      switch (fpViewer->GetViewParameters().GetSpecialMeshRenderingOption()) {
        case G4ViewParameters::meshAsDots:
          Draw3DRectMeshAsDots(mesh);  // Rectangular 3-deep mesh as dots
          implemented = true;
          break;
        case G4ViewParameters::meshAsSurfaces:
          Draw3DRectMeshAsSurfaces(mesh);  // Rectangular 3-deep mesh as surfaces
          implemented = true;
          break;
      }
      break;
    case G4Mesh::tetrahedron:
      switch (fpViewer->GetViewParameters().GetSpecialMeshRenderingOption()) {
        case G4ViewParameters::meshAsDots:
          DrawTetMeshAsDots(mesh);  // Tetrahedron mesh as dots
          implemented = true;
          break;
        case G4ViewParameters::meshAsSurfaces:
          DrawTetMeshAsSurfaces(mesh);  // Tetrahedron mesh as surfaces
          implemented = true;
          break;
      }
      break;
    case G4Mesh::cylinder: [[fallthrough]];
    case G4Mesh::sphere: [[fallthrough]];
    case G4Mesh::invalid: break;
  }
  if (implemented) {
    // Draw container if not marked invisible...
    auto container = mesh.GetContainerVolume();
    auto containerLogical = container->GetLogicalVolume();
    auto containerVisAtts = containerLogical->GetVisAttributes();
    if (containerVisAtts == nullptr || containerVisAtts->IsVisible()) {
      auto solid = containerLogical->GetSolid();
      auto polyhedron = solid->GetPolyhedron();
      // Always draw as wireframe
      G4VisAttributes tmpVisAtts;
      if (containerVisAtts != nullptr) tmpVisAtts = *containerVisAtts;
      tmpVisAtts.SetForceWireframe();
      polyhedron->SetVisAttributes(tmpVisAtts);
      BeginPrimitives(mesh.GetTransform());
      AddPrimitive(*polyhedron);
      EndPrimitives();
    }
  } else {
    // Invoke base class function
    G4VSceneHandler::AddCompound(mesh);
  }
  return;
}

void G4VSceneHandler::Draw3DRectMeshAsDots(const G4Mesh& mesh)
// For a rectangular 3-D mesh, draw as coloured dots by colour and material,
// one dot randomly placed in each visible mesh cell.
{
  // Check
  if (mesh.GetMeshType() != G4Mesh::rectangle &&
      mesh.GetMeshType() != G4Mesh::nested3DRectangular) {
    G4ExceptionDescription ed;
    ed << "Called with a mesh that is not rectangular:" << mesh;
    G4Exception("G4VSceneHandler::Draw3DRectMeshAsDots","visman0108",JustWarning,ed);
    return;
  }

  static G4bool firstPrint = true;
  const auto& verbosity = G4VisManager::GetVerbosity();
  G4bool print = firstPrint && verbosity >= G4VisManager::errors;
  if (print) {
    G4cout
    << "Special case drawing of 3D rectangular G4VNestedParameterisation as dots:"
    << '\n' << mesh
    << G4endl;
  }

  const auto& container = mesh.GetContainerVolume();

  // This map is static so that once filled it stays filled.
  static std::map<G4String,std::map<const G4Material*,G4Polymarker>> dotsByMaterialAndMesh;
  auto& dotsByMaterial = dotsByMaterialAndMesh[mesh.GetContainerVolume()->GetName()];

  // Fill map if not already filled
  if (dotsByMaterial.empty()) {

    // Get positions and material one cell at a time (using PseudoSceneFor3DRectMeshPositions).
    // The pseudo scene allows a "private" descent into the parameterisation.
    // Instantiate a temporary G4PhysicalVolumeModel
    G4ModelingParameters tmpMP;
    tmpMP.SetCulling(true);  // This avoids drawing transparent...
    tmpMP.SetCullingInvisible(true);  // ... or invisble volumes.
    const G4bool useFullExtent = true;  // To avoid calculating the extent
    G4PhysicalVolumeModel tmpPVModel
    (container,
     G4PhysicalVolumeModel::UNLIMITED,
     G4Transform3D(),  // so that positions are in local coordinates
     &tmpMP,
     useFullExtent);
    // Accumulate information in temporary maps by material
    std::multimap<const G4Material*,const G4ThreeVector> positionByMaterial;
    std::map<const G4Material*,G4VSceneHandler::NameAndVisAtts> nameAndVisAttsByMaterial;
    // Instantiate the pseudo scene
    PseudoSceneFor3DRectMeshPositions pseudoScene
    (&tmpPVModel,mesh.GetMeshDepth(),positionByMaterial,nameAndVisAttsByMaterial);
    // Make private descent into the parameterisation
    tmpPVModel.DescribeYourselfTo(pseudoScene);
    // Now we have a map of positions by material.
    // Also a map of name and colour by material.

    const auto& prms = mesh.GetThreeDRectParameters();
    const auto& halfX = prms.fHalfX;
    const auto& halfY = prms.fHalfY;
    const auto& halfZ = prms.fHalfZ;

    // Fill the permanent (static) map of dots by material
    G4int nDotsTotal = 0;
    for (const auto& entry: nameAndVisAttsByMaterial) {
      G4int nDots = 0;
      const auto& material = entry.first;
      const auto& nameAndVisAtts = nameAndVisAttsByMaterial[material];
      const auto& name = nameAndVisAtts.fName;
      const auto& visAtts = nameAndVisAtts.fVisAtts;
      G4Polymarker dots;
      dots.SetInfo(name);
      dots.SetVisAttributes(visAtts);
      dots.SetMarkerType(G4Polymarker::dots);
      dots.SetSize(G4VMarker::screen,1.);
      // Enter empty polymarker into the map
      dotsByMaterial[material] = dots;
      // Now fill it in situ
      auto& dotsInMap = dotsByMaterial[material];
      const auto& range = positionByMaterial.equal_range(material);
      for (auto posByMat = range.first; posByMat != range.second; ++posByMat) {
        dotsInMap.push_back(GetPointInBox(posByMat->second, halfX, halfY, halfZ));
        ++nDots;
      }

      if (print) {
        G4cout
        << std::setw(30) << std::left << name.substr(0,30) << std::right
        << ": " << std::setw(7) << nDots << " dots"
        << ": colour " << std::fixed << std::setprecision(2)
        << visAtts.GetColour() << std::defaultfloat
        << G4endl;
      }

      nDotsTotal += nDots;
    }

    if (print) {
      G4cout << "Total number of dots: " << nDotsTotal << G4endl;
    }
  }

  // Some subsequent expressions apply only to G4PhysicalVolumeModel
  auto pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);

  G4String parameterisationName;
  if (pPVModel) {
    parameterisationName = pPVModel->GetFullPVPath().back().GetPhysicalVolume()->GetName();
  }

  // Draw the dots by material
  // Ensure they are "hidden", i.e., use the z-buffer as non-marker primitives do
  auto keepVP = fpViewer->GetViewParameters();
  auto vp = fpViewer->GetViewParameters();
  vp.SetMarkerHidden();
  fpViewer->SetViewParameters(vp);
  // Now we transform to world coordinates
  BeginPrimitives (mesh.GetTransform());
  for (const auto& entry: dotsByMaterial) {
    const auto& dots = entry.second;
    // The current "leaf" node in the PVPath is the parameterisation. Here it has
    // been converted into polymarkers by material. So...temporarily...change
    // its name to that of the material (whose name has been stored in Info)
    // so that its appearance in the scene tree of, e.g., G4OpenGLQtViewer, has
    // an appropriate name and its visibility and colour may be changed.
    if (pPVModel) {
      const auto& fullPVPath = pPVModel->GetFullPVPath();
      auto leafPV = fullPVPath.back().GetPhysicalVolume();
      leafPV->SetName(dots.GetInfo());
    }
    // Add dots to the scene
    AddPrimitive(dots);
  }
  EndPrimitives ();
  // Restore view parameters
  fpViewer->SetViewParameters(keepVP);
  // Restore parameterisation name
  if (pPVModel) {
    pPVModel->GetFullPVPath().back().GetPhysicalVolume()->SetName(parameterisationName);
  }

  firstPrint = false;
  return;
}

void G4VSceneHandler::Draw3DRectMeshAsSurfaces(const G4Mesh& mesh)
// For a rectangular 3-D mesh, draw as surfaces by colour and material
// with inner shared faces removed.
{
  // Check
  if (mesh.GetMeshType() != G4Mesh::rectangle &&
      mesh.GetMeshType() != G4Mesh::nested3DRectangular) {
    G4ExceptionDescription ed;
    ed << "Called with a mesh that is not rectangular:" << mesh;
    G4Exception("G4VSceneHandler::Draw3DRectMeshAsSurfaces","visman0108",JustWarning,ed);
    return;
  }

  static G4bool firstPrint = true;
  const auto& verbosity = G4VisManager::GetVerbosity();
  G4bool print = firstPrint && verbosity >= G4VisManager::errors;
  if (print) {
    G4cout
    << "Special case drawing of 3D rectangular G4VNestedParameterisation as surfaces:"
    << '\n' << mesh
    << G4endl;
  }

  const auto& container = mesh.GetContainerVolume();

  // This map is static so that once filled it stays filled.
  static std::map<G4String,std::map<const G4Material*,G4Polyhedron>> boxesByMaterialAndMesh;
  auto& boxesByMaterial = boxesByMaterialAndMesh[mesh.GetContainerVolume()->GetName()];

  // Fill map if not already filled
  if (boxesByMaterial.empty()) {

    // Get positions and material one cell at a time (using PseudoSceneFor3DRectMeshPositions).
    // The pseudo scene allows a "private" descent into the parameterisation.
    // Instantiate a temporary G4PhysicalVolumeModel
    G4ModelingParameters tmpMP;
    tmpMP.SetCulling(true);  // This avoids drawing transparent...
    tmpMP.SetCullingInvisible(true);  // ... or invisble volumes.
    const G4bool useFullExtent = true;  // To avoid calculating the extent
    G4PhysicalVolumeModel tmpPVModel
    (container,
     G4PhysicalVolumeModel::UNLIMITED,
     G4Transform3D(),  // so that positions are in local coordinates
     &tmpMP,
     useFullExtent);
    // Accumulate information in temporary maps by material
    std::multimap<const G4Material*,const G4ThreeVector> positionByMaterial;
    std::map<const G4Material*,G4VSceneHandler::NameAndVisAtts> nameAndVisAttsByMaterial;
    // Instantiate the pseudo scene
    PseudoSceneFor3DRectMeshPositions pseudoScene
    (&tmpPVModel,mesh.GetMeshDepth(),positionByMaterial,nameAndVisAttsByMaterial);
    // Make private descent into the parameterisation
    tmpPVModel.DescribeYourselfTo(pseudoScene);
    // Now we have a map of positions by material.
    // Also a map of name and colour by material.

    const auto& prms = mesh.GetThreeDRectParameters();
    const auto& sizeX = 2.*prms.fHalfX;
    const auto& sizeY = 2.*prms.fHalfY;
    const auto& sizeZ = 2.*prms.fHalfZ;

    // Fill the permanent (static) map of boxes by material
    G4int nBoxesTotal = 0, nFacetsTotal = 0;
    for (const auto& entry: nameAndVisAttsByMaterial) {
      G4int nBoxes = 0;
      const auto& material = entry.first;
      const auto& nameAndVisAtts = nameAndVisAttsByMaterial[material];
      const auto& name = nameAndVisAtts.fName;
      const auto& visAtts = nameAndVisAtts.fVisAtts;
      // Transfer positions into a vector ready for creating polyhedral surface
      std::vector<G4ThreeVector> positionsForPolyhedron;
      const auto& range = positionByMaterial.equal_range(material);
      for (auto posByMat = range.first; posByMat != range.second; ++posByMat) {
        const auto& position = posByMat->second;
        positionsForPolyhedron.push_back(position);
        ++nBoxes;
      }
      // The polyhedron will be in local coordinates
      // Add an empty place-holder to the map and get a reference to it
      auto& polyhedron = boxesByMaterial[material];
      // Replace with the desired polyhedron (uses efficient "move assignment")
      polyhedron = G4PolyhedronBoxMesh(sizeX,sizeY,sizeZ,positionsForPolyhedron);
      polyhedron.SetVisAttributes(visAtts);
      polyhedron.SetInfo(name);

      if (print) {
        G4cout
        << std::setw(30) << std::left << name.substr(0,30) << std::right
        << ": " << std::setw(7) << nBoxes << " boxes"
        << " (" << std::setw(7) << 6*nBoxes << " faces)"
        << ": reduced to " << std::setw(7) << polyhedron.GetNoFacets() << " facets ("
        << std::setw(2) << std::fixed << std::setprecision(2) << 100*polyhedron.GetNoFacets()/(6*nBoxes)
        << "%): colour " << std::fixed << std::setprecision(2)
        << visAtts.GetColour() << std::defaultfloat
        << G4endl;
      }

      nBoxesTotal += nBoxes;
      nFacetsTotal += polyhedron.GetNoFacets();
    }

    if (print) {
      G4cout << "Total number of boxes: " << nBoxesTotal << " (" << 6*nBoxesTotal << " faces)"
      << ": reduced to " << nFacetsTotal << " facets ("
      << std::setw(2) << std::fixed << std::setprecision(2) << 100*nFacetsTotal/(6*nBoxesTotal) << "%)"
      << G4endl;
    }
  }

  // Some subsequent expressions apply only to G4PhysicalVolumeModel
  auto pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);

  G4String parameterisationName;
  if (pPVModel) {
    parameterisationName = pPVModel->GetFullPVPath().back().GetPhysicalVolume()->GetName();
  }

  // Draw the boxes by material
  // Now we transform to world coordinates
  BeginPrimitives (mesh.GetTransform());
  for (const auto& entry: boxesByMaterial) {
    const auto& poly = entry.second;
    // The current "leaf" node in the PVPath is the parameterisation. Here it has
    // been converted into polyhedra by material. So...temporarily...change
    // its name to that of the material (whose name has been stored in Info)
    // so that its appearance in the scene tree of, e.g., G4OpenGLQtViewer, has
    // an appropriate name and its visibility and colour may be changed.
    if (pPVModel) {
      const auto& fullPVPath = pPVModel->GetFullPVPath();
      auto leafPV = fullPVPath.back().GetPhysicalVolume();
      leafPV->SetName(poly.GetInfo());
    }
    AddPrimitive(poly);
  }
  EndPrimitives ();
  // Restore parameterisation name
  if (pPVModel) {
    pPVModel->GetFullPVPath().back().GetPhysicalVolume()->SetName(parameterisationName);
  }

  firstPrint = false;
  return;
}

void G4VSceneHandler::DrawTetMeshAsDots(const G4Mesh& mesh)
// For a tetrahedron mesh, draw as coloured dots by colour and material,
// one dot randomly placed in each visible mesh cell.
{
  // Check
  if (mesh.GetMeshType() != G4Mesh::tetrahedron) {
    G4ExceptionDescription ed;
    ed << "Called with mesh that is not a tetrahedron mesh:" << mesh;
    G4Exception("G4VSceneHandler::DrawTetMeshAsDots","visman0108",JustWarning,ed);
    return;
  }

  static G4bool firstPrint = true;
  const auto& verbosity = G4VisManager::GetVerbosity();
  G4bool print = firstPrint && verbosity >= G4VisManager::errors;

  if (print) {
    G4cout
    << "Special case drawing of tetrahedron mesh as dots"
    << '\n' << mesh
    << G4endl;
  }

  const auto& container = mesh.GetContainerVolume();

  // This map is static so that once filled it stays filled.
  static std::map<G4String,std::map<const G4Material*,G4Polymarker>> dotsByMaterialAndMesh;
  auto& dotsByMaterial = dotsByMaterialAndMesh[mesh.GetContainerVolume()->GetName()];

  // Fill map if not already filled
  if (dotsByMaterial.empty()) {

    // Get vertices and colour one cell at a time (using PseudoSceneForTetVertices).
    // The pseudo scene allows a "private" descent into the parameterisation.
    // Instantiate a temporary G4PhysicalVolumeModel
    G4ModelingParameters tmpMP;
    tmpMP.SetCulling(true);  // This avoids drawing transparent...
    tmpMP.SetCullingInvisible(true);  // ... or invisble volumes.
    const G4bool useFullExtent = true;  // To avoid calculating the extent
    G4PhysicalVolumeModel tmpPVModel
    (container,
     G4PhysicalVolumeModel::UNLIMITED,
     G4Transform3D(),  // so that positions are in local coordinates
     &tmpMP,
     useFullExtent);
    // Accumulate information in temporary maps by material
    std::multimap<const G4Material*,std::vector<G4ThreeVector>> verticesByMaterial;
    std::map<const G4Material*,G4VSceneHandler::NameAndVisAtts> nameAndVisAttsByMaterial;
    // Instantiate a pseudo scene
    PseudoSceneForTetVertices pseudoScene
    (&tmpPVModel,mesh.GetMeshDepth(),verticesByMaterial,nameAndVisAttsByMaterial);
    // Make private descent into the parameterisation
    tmpPVModel.DescribeYourselfTo(pseudoScene);
    // Now we have a map of vertices by material.
    // Also a map of name and colour by material.

    // Fill the permanent (static) map of dots by material
    G4int nDotsTotal = 0;
    for (const auto& entry: nameAndVisAttsByMaterial) {
      G4int nDots = 0;
      const auto& material = entry.first;
      const auto& nameAndVisAtts = nameAndVisAttsByMaterial[material];
      const auto& name = nameAndVisAtts.fName;
      const auto& visAtts = nameAndVisAtts.fVisAtts;
      G4Polymarker dots;
      dots.SetVisAttributes(visAtts);
      dots.SetMarkerType(G4Polymarker::dots);
      dots.SetSize(G4VMarker::screen,1.);
      dots.SetInfo(name);
      // Enter empty polymarker into the map
      dotsByMaterial[material] = dots;
      // Now fill it in situ
      auto& dotsInMap = dotsByMaterial[material];
      const auto& range = verticesByMaterial.equal_range(material);
      for (auto vByMat = range.first; vByMat != range.second; ++vByMat) {
        dotsInMap.push_back(GetPointInTet(vByMat->second));
        ++nDots;
      }

      if (print) {
        G4cout
        << std::setw(30) << std::left << name.substr(0,30) << std::right
        << ": " << std::setw(7) << nDots << " dots"
        << ": colour " << std::fixed << std::setprecision(2)
        << visAtts.GetColour() << std::defaultfloat
        << G4endl;
      }

      nDotsTotal += nDots;
    }

    if (print) {
      G4cout << "Total number of dots: " << nDotsTotal << G4endl;
    }
  }

  // Some subsequent expressions apply only to G4PhysicalVolumeModel
  auto pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);

  G4String parameterisationName;
  if (pPVModel) {
    parameterisationName = pPVModel->GetFullPVPath().back().GetPhysicalVolume()->GetName();
  }

  // Draw the dots by material
  // Ensure they are "hidden", i.e., use the z-buffer as non-marker primitives do
  auto keepVP = fpViewer->GetViewParameters();
  auto vp = fpViewer->GetViewParameters();
  vp.SetMarkerHidden();
  fpViewer->SetViewParameters(vp);

  // Now we transform to world coordinates
  BeginPrimitives (mesh.GetTransform());
  for (const auto& entry: dotsByMaterial) {
    const auto& dots = entry.second;
    // The current "leaf" node in the PVPath is the parameterisation. Here it has
    // been converted into polymarkers by material. So...temporarily...change
    // its name to that of the material (whose name has been stored in Info)
    // so that its appearance in the scene tree of, e.g., G4OpenGLQtViewer, has
    // an appropriate name and its visibility and colour may be changed.
    if (pPVModel) {
      const auto& fullPVPath = pPVModel->GetFullPVPath();
      auto leafPV = fullPVPath.back().GetPhysicalVolume();
      leafPV->SetName(dots.GetInfo());
    }
    AddPrimitive(dots);
  }
  EndPrimitives ();

  // Restore view parameters
  fpViewer->SetViewParameters(keepVP);
  // Restore parameterisation name
  if (pPVModel) {
    pPVModel->GetFullPVPath().back().GetPhysicalVolume()->SetName(parameterisationName);
  }

  firstPrint = false;
  return;
}

void G4VSceneHandler::DrawTetMeshAsSurfaces(const G4Mesh& mesh)
// For a tetrahedron mesh, draw as surfaces by colour and material
// with inner shared faces removed.
{
  // Check
  if (mesh.GetMeshType() != G4Mesh::tetrahedron) {
    G4ExceptionDescription ed;
    ed << "Called with mesh that is not a tetrahedron mesh:" << mesh;
    G4Exception("G4VSceneHandler::DrawTetMeshAsSurfaces","visman0108",JustWarning,ed);
    return;
  }

  static G4bool firstPrint = true;
  const auto& verbosity = G4VisManager::GetVerbosity();
  G4bool print = firstPrint && verbosity >= G4VisManager::errors;

  if (print) {
    G4cout
    << "Special case drawing of tetrahedron mesh as surfaces"
    << '\n' << mesh
    << G4endl;
  }

  // This map is static so that once filled it stays filled.
  static std::map<G4String,std::map<const G4Material*,G4Polyhedron>> surfacesByMaterialAndMesh;
  auto& surfacesByMaterial = surfacesByMaterialAndMesh[mesh.GetContainerVolume()->GetName()];

  // Fill map if not already filled
  if (surfacesByMaterial.empty()) {

    // Get vertices and colour one cell at a time (using PseudoSceneForTetVertices).
    // The pseudo scene allows a "private" descent into the parameterisation.
    // Instantiate a temporary G4PhysicalVolumeModel
    G4ModelingParameters tmpMP;
    tmpMP.SetCulling(true);  // This avoids drawing transparent...
    tmpMP.SetCullingInvisible(true);  // ... or invisble volumes.
    const G4bool useFullExtent = true;  // To avoid calculating the extent
    G4PhysicalVolumeModel tmpPVModel
    (mesh.GetContainerVolume(),
     G4PhysicalVolumeModel::UNLIMITED,
     G4Transform3D(),  // so that positions are in local coordinates
     &tmpMP,
     useFullExtent);
    // Accumulate information in temporary maps by material
    std::multimap<const G4Material*,std::vector<G4ThreeVector>> verticesByMaterial;
    std::map<const G4Material*,G4VSceneHandler::NameAndVisAtts> nameAndVisAttsByMaterial;
    // Instantiate a pseudo scene
    PseudoSceneForTetVertices pseudoScene
    (&tmpPVModel,mesh.GetMeshDepth(),verticesByMaterial,nameAndVisAttsByMaterial);
    // Make private descent into the parameterisation
    tmpPVModel.DescribeYourselfTo(pseudoScene);
    // Now we have a map of vertices by material.
    // Also a map of name and colour by material.

    // Fill the permanent (static) map of surfaces by material
    G4int nTetsTotal = 0, nFacetsTotal = 0;
    for (const auto& entry: nameAndVisAttsByMaterial) {
      G4int nTets = 0;
      const auto& material = entry.first;
      const auto& nameAndVisAtts = nameAndVisAttsByMaterial[material];
      const auto& name = nameAndVisAtts.fName;
      const auto& visAtts = nameAndVisAtts.fVisAtts;
      // Transfer vertices into a vector ready for creating polyhedral surface
      std::vector<G4ThreeVector> verticesForPolyhedron;
      const auto& range = verticesByMaterial.equal_range(material);
      for (auto vByMat = range.first; vByMat != range.second; ++vByMat) {
        const std::vector<G4ThreeVector>& vertices = vByMat->second;
        for (const auto& vertex: vertices)
          verticesForPolyhedron.push_back(vertex);
        ++nTets;
      }
      // The polyhedron will be in local coordinates
      // Add an empty place-holder to the map and get a reference to it
      auto& polyhedron = surfacesByMaterial[material];
      // Replace with the desired polyhedron (uses efficient "move assignment")
      polyhedron = G4PolyhedronTetMesh(verticesForPolyhedron);
      polyhedron.SetVisAttributes(visAtts);
      polyhedron.SetInfo(name);

      if (print) {
        G4cout
        << std::setw(30) << std::left << name.substr(0,30) << std::right
        << ": " << std::setw(7) << nTets << " tetrahedra"
        << " (" << std::setw(7) << 4*nTets << " faces)"
        << ": reduced to " << std::setw(7) << polyhedron.GetNoFacets() << " facets ("
        << std::setw(2) << std::fixed << std::setprecision(2) << 100*polyhedron.GetNoFacets()/(4*nTets)
        << "%): colour " << std::fixed << std::setprecision(2)
        << visAtts.GetColour() << std::defaultfloat
        << G4endl;
     }

      nTetsTotal += nTets;
      nFacetsTotal += polyhedron.GetNoFacets();
    }

    if (print) {
      G4cout << "Total number of tetrahedra: " << nTetsTotal << " (" << 4*nTetsTotal << " faces)"
      << ": reduced to " << nFacetsTotal << " facets ("
      << std::setw(2) << std::fixed << std::setprecision(2) << 100*nFacetsTotal/(4*nTetsTotal) << "%)"
      << G4endl;
    }
  }

  // Some subsequent expressions apply only to G4PhysicalVolumeModel
  auto pPVModel = dynamic_cast<G4PhysicalVolumeModel*>(fpModel);

  G4String parameterisationName;
  if (pPVModel) {
    parameterisationName = pPVModel->GetFullPVPath().back().GetPhysicalVolume()->GetName();
  }

  // Draw the surfaces by material
  // Now we transform to world coordinates
  BeginPrimitives (mesh.GetTransform());
  for (const auto& entry: surfacesByMaterial) {
    const auto& poly = entry.second;
    // The current "leaf" node in the PVPath is the parameterisation. Here it has
    // been converted into polyhedra by material. So...temporarily...change
    // its name to that of the material (whose name has been stored in Info)
    // so that its appearance in the scene tree of, e.g., G4OpenGLQtViewer, has
    // an appropriate name and its visibility and colour may be changed.
    if (pPVModel) {
      const auto& fullPVPath = pPVModel->GetFullPVPath();
      auto leafPV = fullPVPath.back().GetPhysicalVolume();
      leafPV->SetName(poly.GetInfo());
    }
    AddPrimitive(poly);
  }
  EndPrimitives ();

  // Restore parameterisation name
  if (pPVModel) {
    pPVModel->GetFullPVPath().back().GetPhysicalVolume()->SetName(parameterisationName);
  }

  firstPrint = false;
  return;
}

G4ThreeVector
G4VSceneHandler::GetPointInBox(const G4ThreeVector& pos,
                               G4double halfX,
                               G4double halfY,
                               G4double halfZ) const
{
  G4double x = pos.getX() + (2.*G4QuickRand() - 1.)*halfX;
  G4double y = pos.getY() + (2.*G4QuickRand() - 1.)*halfY;
  G4double z = pos.getZ() + (2.*G4QuickRand() - 1.)*halfZ;
  return G4ThreeVector(x, y, z);
}

G4ThreeVector
G4VSceneHandler::GetPointInTet(const std::vector<G4ThreeVector>& vertices) const
{
  G4double p = G4QuickRand();
  G4double q = G4QuickRand();
  G4double r = G4QuickRand();
  if (p + q > 1.)
  {
    p = 1. - p;
    q = 1. - q;
  }
  if (q + r > 1.)
  {
    G4double tmp = r;
    r = 1. - p - q;
    q = 1. - tmp;
  }
  else if (p + q + r > 1.)
  {
    G4double tmp = r;
    r = p + q + r - 1.;
    p = 1. - q - tmp;
  }
  G4double a = 1. - p - q - r;
  return vertices[0]*a + vertices[1]*p + vertices[2]*q + vertices[3]*r;
}

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSceneHandler.cc,v 1.7 1999-11-11 15:38:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996
// Abstract interface class for graphics scenes.

#include "G4VSceneHandler.hh"

#include "G4ios.hh"
#ifdef WIN32
#include <strstrea.h>
#else
#include <strstream.h>
#endif

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
#include "G4NURBS.hh"
#include "G4Visible.hh"
#include "G4VisAttributes.hh"
#include "G4VModel.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4ModelingParameters.hh"

G4VSceneHandler::G4VSceneHandler (G4VGraphicsSystem& system, G4int id, const G4String& name):
  fSystem                (system),
  fSceneId               (id),
  fViewCount             (0),
  fpViewer               (0),
  fReadyForTransients    (false),
  fpModel                (0),
  fpObjectTransformation (&G4Transform3D::Identity),
  fpVisAttribs           (0),
  fCurrentDepth          (0),
  fpCurrentPV            (0),
  fpCurrentLV            (0)
{
  G4VisManager* pVMan = G4VisManager::GetInstance ();
  fpScene = pVMan -> GetCurrentScene ();
  if (name == "") {
    char charname [50];
    ostrstream ost (charname, 50);
    ost << fSystem.GetName () << '-' << fSceneId << ends;
    fName = charname;
  }
  else {
    fName = name;
  }
}

G4VSceneHandler::~G4VSceneHandler () {
  fViewerList.clearAndDestroy ();
}

void G4VSceneHandler::EndModeling () {}

void G4VSceneHandler::PreAddThis (const G4Transform3D& objectTransformation,
                                  const G4VisAttributes& visAttribs) {
  fpObjectTransformation = &objectTransformation;
  fpVisAttribs = &visAttribs;
}

void G4VSceneHandler::PostAddThis () {
  fpObjectTransformation = &G4Transform3D::Identity;
  fpVisAttribs = 0;
}

void G4VSceneHandler::ClearStore () {
  fpViewer -> NeedKernelVisit ();
}

void G4VSceneHandler::ClearTransientStore () {}

void G4VSceneHandler::AddThis (const G4Box& box) {
  RequestPrimitives (box);
// If your graphics system is sophisticated enough to handle a
//  particular solid shape as a primitive, in your derived class write a
//  function to override this.  (Note: some compilers warn that your
//  function "hides" this one.  That's OK.)
// Your function might look like this...
// void G4MyScene::AddThis (const G4Box& box) {
// Get parameters of appropriate object, e.g.:
//   G4double dx = box.GetXHalfLength ();
//   G4double dy = box.GetYHalfLength ();
//   G4double dz = box.GetZHalfLength ();
// and Draw or Store in your display List.
}

void G4VSceneHandler::AddThis (const G4Tubs& tubs) {
  RequestPrimitives (tubs);
}

void G4VSceneHandler::AddThis (const G4Cons& cons) {
  RequestPrimitives (cons);
}

void G4VSceneHandler::AddThis (const G4Trd& trd) {
  RequestPrimitives (trd);
}

void G4VSceneHandler::AddThis (const G4Trap& trap) {
  RequestPrimitives (trap);
}

void G4VSceneHandler::AddThis (const G4Sphere& sphere) {
  RequestPrimitives (sphere );
}

void G4VSceneHandler::AddThis (const G4Para& para) {
  RequestPrimitives (para);
}

void G4VSceneHandler::AddThis (const G4Torus& torus) {
  RequestPrimitives (torus);
}

void G4VSceneHandler::AddThis (const G4Polycone& polycone) {
  RequestPrimitives (polycone);
}

void G4VSceneHandler::AddThis (const G4Polyhedra& polyhedra) {
  RequestPrimitives (polyhedra);
}

void G4VSceneHandler::AddThis (const G4VSolid& solid) {
  RequestPrimitives (solid);
}

void G4VSceneHandler::AddViewerToList (G4VViewer* pViewer) {
  fViewerList.append (pViewer);
}

void G4VSceneHandler::EstablishSpecials (G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace (&fCurrentDepth,
					&fpCurrentPV,
					&fpCurrentLV);
}

void G4VSceneHandler::BeginModeling () {
  if (!GetModel ()) G4Exception ("G4VSceneHandler::BeginModeling: NO MODEL!!!");
}

void G4VSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation) {
  if (!GetModel ()) G4Exception ("G4VSceneHandler::BeginPrimitives: NO MODEL!!!");
  fpObjectTransformation = &objectTransformation;
}

void G4VSceneHandler::EndPrimitives () {}

void G4VSceneHandler::AddPrimitive (const G4Polymarker& polymarker) {
  switch (polymarker.GetMarkerType()) {
  default:
  case G4Polymarker::line:
    {
      G4Polyline polyline (polymarker);
      for (int iPoint = 0; iPoint < polymarker.entries (); iPoint++) {
	polyline.append (polymarker[iPoint]);
      }
      AddPrimitive (polyline);
    }
    break;
  case G4Polymarker::dots:
    {
      for (int iPoint = 0; iPoint < polymarker.entries (); iPoint++) {
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
      for (int iPoint = 0; iPoint < polymarker.entries (); iPoint++) {
	G4Circle circle (polymarker);
	circle.SetPosition (polymarker[iPoint]);
	AddPrimitive (circle);
      }
    }
    break;
  case G4Polymarker::squares:
    {
      for (int iPoint = 0; iPoint < polymarker.entries (); iPoint++) {
	G4Square Square (polymarker);
	Square.SetPosition (polymarker[iPoint]);
	AddPrimitive (Square);
      }
    }
    break;
  }
}

void G4VSceneHandler::RemoveViewerFromList (G4VViewer* pViewer) {
  fViewerList.remove (pViewer);
}

void G4VSceneHandler::SetScene (G4Scene* pScene) {
  fpScene = pScene;
  // Notify all viewers that a kernel visit is required.
  for (int i = 0; i < fViewerList.entries (); i++) {
    fViewerList [i] -> SetNeedKernelVisit ();
  }
}

void G4VSceneHandler::RequestPrimitives (const G4VSolid& solid) {
  G4Polyhedron* pPolyhedron;
  G4NURBS*      pNURBS;
  BeginPrimitives (*fpObjectTransformation);
  switch (GetModel () -> GetModelingParameters () -> GetRepStyle ()) {
  case G4ModelingParameters::nurbs:
    pNURBS = solid.CreateNURBS ();
    if (pNURBS) {
      pNURBS -> SetVisAttributes
	(fpViewer -> GetApplicableVisAttributes (fpVisAttribs));
      AddPrimitive (*pNURBS);
      delete pNURBS;
      break;
    }
    else {
      G4cerr << "NURBS not available for " << solid.GetName () << endl;
      G4cerr << "Trying polyhedron." << endl;
    }
    // Dropping through to polyhedron...
  case G4ModelingParameters::polyhedron:
  default:
    G4Polyhedron::SetNumberOfRotationSteps
	(GetModel () -> GetModelingParameters () -> GetNoOfSides ());
    pPolyhedron = solid.CreatePolyhedron ();
    G4Polyhedron::ResetNumberOfRotationSteps ();
    if (pPolyhedron) {
      pPolyhedron -> SetVisAttributes
	(fpViewer -> GetApplicableVisAttributes (fpVisAttribs));
      AddPrimitive (*pPolyhedron);
      delete pPolyhedron;
    }
    else {
      G4cerr << "Polyhedron not available for " << solid.GetName ()
	   <<".\nThis means it cannot be visualized on most systems."
	"\nContact the Visualization Coordinator." << endl;
    }
    break;
  case G4ModelingParameters::hierarchy:
    G4cout << "Geometry hierarchy requested: tag " 
	 << GetModel () -> GetCurrentTag ()
	 << ", depth "
	 << fCurrentDepth
	 << endl;
    break;
  }
  EndPrimitives ();
}

void G4VSceneHandler::ProcessScene (G4VViewer& view) {

  fReadyForTransients = false;

  // Clear stored scene, if any, i.e., display lists, scene graphs.
  ClearStore ();

  // Traverse geometry tree and send drawing primitives to window(s).

  const G4RWTPtrOrderedVector <G4VModel>& runDurationModelList =
    fpScene -> GetRunDurationModelList ();

  if (runDurationModelList.entries ()) {
    G4cout << "Traversing scene data..." << endl;
    G4ModelingParameters* pMP = CreateModelingParameters ();
    for (int i = 0; i < runDurationModelList.entries (); i++) {
      G4VModel* pModel = runDurationModelList[i];
      const G4ModelingParameters* tempMP =
	pModel -> GetModelingParameters ();
      // NOTE THAT tempMP COULD BE ZERO.
      // (Not sure the above is necessary; but in future we might
      // want to take notice of the modeling parameters with which
      // the model was created.  For the time being we are ignoring
      // them and simply using the view parameters.  When the time
      // comes to do this, then perhaps there should be a default
      // set of modeling parameters in the view parameters for the
      // case of a zero modeling parameters pointer.)
      // (I think for the G4 Vis System we'll rely on view parameters
      // and convert using pMP = CreateModelingParameters () as above.)
      pModel -> SetModelingParameters (pMP);
      SetModel (pModel);  // Store for use by derived class.
      BeginModeling ();
      pModel -> DescribeYourselfTo (*this);
      EndModeling ();
      pModel -> SetModelingParameters (tempMP);
    }
    delete pMP;
    SetModel (0);  // Flags invalid model.
  }
  else {
    G4cerr << "No run-duration models in scene data." << endl;
  }

  fReadyForTransients = true;
}

G4ModelingParameters* G4VSceneHandler::CreateModelingParameters () {
  // Create modeling parameters from View Parameters...
  const G4ViewParameters& vp = fpViewer -> GetViewParameters ();
  // Convert rep styles...
  G4ModelingParameters::RepStyle modelRepStyle =
    G4ModelingParameters::wireframe;
  if (vp.GetDrawingStyle () != G4ViewParameters::wireframe) {
    switch (vp.GetRepStyle ()) {
    default:
    case G4ViewParameters::polyhedron:
      modelRepStyle = G4ModelingParameters::polyhedron;
      break;
    case G4ViewParameters::nurbs:
      modelRepStyle = G4ModelingParameters::nurbs;
      break;
    }
  }

  // Decide if covered daughters are really to be culled...
  G4bool reallyCullCovered =
    vp.IsCullingCovered()   // Culling daughters depends also on...
    && !vp.IsSection ()     // Sections (DCUT) not requested.
    && !vp.IsCutaway ()     // Cutaways not requested.
    && (                    // Surface drawing in operation.
	vp.GetDrawingStyle () == G4ViewParameters::hsr ||
	vp.GetDrawingStyle () == G4ViewParameters::hlhsr
	);

  /////////////// TESTING TESTING
  //modelRepStyle = G4ModelingParameters::hierarchy;
  ///////////////////////////////

  G4ModelingParameters* pModelingParams = new G4ModelingParameters
    (vp.GetDefaultVisAttributes (),
     modelRepStyle,
     vp.IsCulling (),
     vp.IsCullingInvisible (),
     vp.IsDensityCulling (),
     vp.GetVisibleDensity (),
     reallyCullCovered,
     vp.GetNoOfSides (),
     vp.IsViewGeom (),
     vp.IsViewHits (),
     vp.IsViewDigis ()
     );

  return pModelingParams;
}

const G4Colour& G4VSceneHandler::GetColour (const G4Visible& visible) {
  // Colour is determined by the applicable (real) vis attributes.
  const G4VisAttributes* pVA = visible.GetVisAttributes ();
  pVA = fpViewer -> GetApplicableVisAttributes (pVA);
  return pVA -> GetColour ();
}

const G4Colour& G4VSceneHandler::GetTextColour (const G4Text& text) {
  const G4VisAttributes* pVA = text.GetVisAttributes ();
  if (!pVA) {
    pVA = fpViewer -> GetViewParameters (). GetDefaultTextVisAttributes ();
  }
  return pVA -> GetColour ();
}

G4ViewParameters::DrawingStyle G4VSceneHandler::GetDrawingStyle
(const G4Visible& visible) {
  // Drawing style is determined by the applicable (real) vis
  // attributes, except when overridden - see GetDrawingStyle (const
  // G4VisAttributes* pVisAttribs).
  const G4VisAttributes* pVA = visible.GetVisAttributes ();
  pVA = fpViewer -> GetApplicableVisAttributes (pVA);
  return GetDrawingStyle (pVA);
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
    // This is complicated because is hidden line removal has been
    // requested we wish to preserve this.
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
      switch (style) {
      case (G4ViewParameters::hlhsr):
	style = G4ViewParameters::hlr;
	break;
      case (G4ViewParameters::hsr):
	style = G4ViewParameters::wireframe;
	break;
      case (G4ViewParameters::hlr):
      case (G4ViewParameters::wireframe):
      default:
	break;
      }	
      break;
    }
  }
  return style;
}

G4double G4VSceneHandler::GetMarkerSize (const G4VMarker& marker, 
				  G4VSceneHandler::MarkerSizeType& markerSizeType) {
  G4bool userSpecified = marker.GetWorldSize() || marker.GetScreenSize();
  const G4VMarker& defaultMarker =
    fpViewer -> GetViewParameters().GetDefaultMarker();
  G4double size;
  if (size = // Assignment intentional.
      userSpecified ? marker.GetWorldSize() : defaultMarker.GetWorldSize()) {
    // Draw in world coordinates.
    markerSizeType = world;
  }
  else {
    size = // Assignment intentional.
      userSpecified ? marker.GetScreenSize() : defaultMarker.GetScreenSize();
    // Draw in screen coordinates.
    markerSizeType = screen;
  }
  if (size <= 0.) size = 1.;
  size *= fpViewer -> GetViewParameters().GetGlobalMarkerScale();
  return size;
}

ostream& operator << (ostream& os, const G4VSceneHandler& s) {

  G4VisManager* pVMan = G4VisManager::GetInstance ();

  os << "Scene " << s.fName << " has "
     << s.fViewerList.entries () << " viewers:";
  for (int i = 0; i < s.fViewerList.entries (); i++) {
    os << "\n  " << s.fViewerList [i];
  }

  os << "\n  " << *s.fpScene;

  return os;
}

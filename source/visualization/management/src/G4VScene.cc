// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VScene.cc,v 1.1 1999-01-07 16:15:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996
// Abstract interface class for graphics scenes.

#include "G4VScene.hh"

#include "G4ios.hh"
#ifdef WIN32
#include <strstrea.h>
#else
#include <strstream.h>
#endif

#include "G4VisManager.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VView.hh"
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

G4VScene::G4VScene (G4VGraphicsSystem& system, G4int id, const G4String& name):
  fSystem                (system),
  fSceneId               (id),
  fViewCount             (0),
  fpView                 (0),
  fReadyForTransients    (false),
  fpModel                (0),
  fpObjectTransformation (0),
  fpVisAttribs           (0),
  fCurrentDepth          (0),
  fpCurrentPV            (0),
  fpCurrentLV            (0)
{
  G4VisManager* pVMan = G4VisManager::GetInstance ();
  fSD = pVMan -> GetCurrentSceneData ();
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

G4VScene::~G4VScene () {
  fViewList.clearAndDestroy ();
}

void G4VScene::AddViewToList (G4VView* pView) {
  fViewList.append (pView);
}

void G4VScene::EstablishSpecials (G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace (&fCurrentDepth,
					&fpCurrentPV,
					&fpCurrentLV);
}

void G4VScene::BeginModeling () {
  if (!GetModel ()) G4Exception ("G4VScene::BeginModeling: NO MODEL!!!");
}

void G4VScene::BeginPrimitives
(const G4Transform3D& objectTransformation) {
  if (!GetModel ()) G4Exception ("G4VScene::BeginPrimitives: NO MODEL!!!");
  fpObjectTransformation = &objectTransformation;
}

void G4VScene::EndPrimitives () {}

void G4VScene::AddPrimitive (const G4Polymarker& polymarker) {
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

void G4VScene::RemoveViewFromList (G4VView* pView) {
  fViewList.remove (pView);
}

void G4VScene::SetSceneData (const G4SceneData& sd) {
  fSD = sd;
  // Notify all views that a kernel visit is required.
  for (int i = 0; i < fViewList.entries (); i++) {
    fViewList [i] -> SetNeedKernelVisit ();
  }
}

void G4VScene::RequestPrimitives (const G4VSolid& solid) {
  G4Polyhedron* pPolyhedron;
  G4NURBS*      pNURBS;
  BeginPrimitives (*fpObjectTransformation);
  switch (GetModel () -> GetModelingParameters () -> GetRepStyle ()) {
  case G4ModelingParameters::nurbs:
    pNURBS = solid.CreateNURBS ();
    if (pNURBS) {
      pNURBS -> SetVisAttributes
	(fpView -> GetApplicableVisAttributes (fpVisAttribs));
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
	(fpView -> GetApplicableVisAttributes (fpVisAttribs));
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

void G4VScene::ProcessScene (G4VView& view) {

  fReadyForTransients = false;

  // Clear stored scene, if any, i.e., display lists, scene graphs.
  ClearStore ();

  if (fpView -> GetViewParameters ().IsViewGeom ()) {

    // Traverse geometry tree and send drawing primitives to window(s).

    const RWTPtrOrderedVector <G4VModel>& runDurationModelList =
      fSD.GetRunDurationModelList ();

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
  }

  fReadyForTransients = true;
}

G4ModelingParameters* G4VScene::CreateModelingParameters () {
  // Create modeling parameters from View Parameters...
  const G4ViewParameters& vp = fpView -> GetViewParameters ();
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
     vp.GetNoOfSides ());

  return pModelingParams;
}

const G4Colour& G4VScene::GetColour (const G4Visible& visible) {
  // Colour is determined by the applicable (real) vis attributes.
  const G4VisAttributes* pVA = visible.GetVisAttributes ();
  pVA = fpView -> GetApplicableVisAttributes (pVA);
  return pVA -> GetColour ();
}

const G4Colour& G4VScene::GetTextColour (const G4Text& text) {
  const G4VisAttributes* pVA = text.GetVisAttributes ();
  if (!pVA) {
    pVA = fpView -> GetViewParameters (). GetDefaultTextVisAttributes ();
  }
  return pVA -> GetColour ();
}

G4ViewParameters::DrawingStyle G4VScene::GetDrawingStyle
(const G4Visible& visible) {
  // Drawing style is determined by the applicable (real) vis
  // attributes, except when overridden - see GetDrawingStyle (const
  // G4VisAttributes* pVisAttribs).
  const G4VisAttributes* pVA = visible.GetVisAttributes ();
  pVA = fpView -> GetApplicableVisAttributes (pVA);
  return GetDrawingStyle (pVA);
}

G4ViewParameters::DrawingStyle G4VScene::GetDrawingStyle
(const G4VisAttributes* pVisAttribs) {
  // Drawing style is normally determined by the view parameters, but
  // it can be overriddden by the ForceDrawingStyle flag in the vis
  // attributes.
  G4ViewParameters::DrawingStyle style =
    fpView->GetViewParameters().GetDrawingStyle();
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

G4double G4VScene::GetMarkerSize (const G4VMarker& marker, 
				  G4VScene::MarkerSizeType& markerSizeType) {
  G4bool userSpecified = marker.GetWorldSize() || marker.GetScreenSize();
  const G4VMarker& defaultMarker =
    fpView -> GetViewParameters().GetDefaultMarker();
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
  size *= fpView -> GetViewParameters().GetGlobalMarkerScale();
  return size;
}

ostream& operator << (ostream& os, const G4VScene& s) {

  G4VisManager* pVMan = G4VisManager::GetInstance ();

  os << "Scene " << s.fName << " has "
     << s.fViewList.entries () << " views:";
  for (int i = 0; i < s.fViewList.entries (); i++) {
    os << "\n  " << s.fViewList [i];
  }

  os << "\n  " << s.fSD;

  return os;
}

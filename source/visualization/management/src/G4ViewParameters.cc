// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ViewParameters.cc,v 1.2 1999-01-09 16:31:17 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996
// View parameters and options.

#include "G4ViewParameters.hh"

#include "G4VisManager.hh"
#include "G4ios.hh"

G4ViewParameters::G4ViewParameters ():
  fDrawingStyle (wireframe),
  fRepStyle (polyhedron),
  fCulling (true),
  fCullInvisible (true),
  fDensityCulling (false),
  fVisibleDensity (0.01 * g / cm3),
  fCullCovered (false),
  fSection (false),
  fCutaway (false),
  fExplode (false),
  fExplodeFactor (1.),
  fNoOfSides (24),
  fViewpointDirection (G4Vector3D (0., 0., 1.)),  // On z-axis.
  fUpVector (G4Vector3D (0., 1., 0.)),            // y-axis up.
  fFieldHalfAngle (0.),                           // Orthogonal projection.
  fZoomFactor (1.),
  fDolly (0.),
  fRelativeLightpointDirection (G4Vector3D (1., 1., 1.)),
  fActualLightpointDirection (G4Vector3D (1., 1., 1.)),
  fLightsMoveWithCamera (false),
  fViewGeom (true),
  fViewHits (false),
  fViewDigis (false),
  fDefaultTextVisAttributes (G4Colour (0., 0., 1.)),
  fGlobalMarkerScale (1.),
  fMarkerNotHidden (true),
  fWindowSizeHintX (600),
  fWindowSizeHintY (600)
// Note: all other parameters use default constructors.
{
  fDefaultMarker.SetScreenSize (5.);
  // Markers are 5 pixels radius, 10 pixels diameter.
}

G4ViewParameters::~G4ViewParameters () {
// Clear cutaway planes?  Rogue Wave probably destroys OK.
}

G4Vector3D& G4ViewParameters::GetActualLightpointDirection () {
  SetViewAndLights (fViewpointDirection);
  return fActualLightpointDirection;
}

// Useful quantities - begin snippet.
// Here Follow functions to evaluate the above algorithms as a
// function of the radius of the Bounding Sphere of the object being
// viewed.  Call them in the order given - for efficiency, later
// functions depend on the results of earlier ones (Store the
// results of earlier functions in your own temporary variables -
// see, for example, G4OpenGLView::SetView ().)

G4double G4ViewParameters::GetCameraDistance (G4double radius) const {
  G4double cameraDistance;
  if (fFieldHalfAngle == 0.) {
    cameraDistance = radius;
  }
  else {
    cameraDistance = radius / sin (fFieldHalfAngle) - fDolly;
  }
  return cameraDistance;
}

G4double G4ViewParameters::GetNearDistance (G4double cameraDistance,
					    G4double radius) const {
  const G4double small = 1.e-6 * radius;
  G4double nearDistance = cameraDistance - radius;
  if (nearDistance < small) nearDistance = small;
  return nearDistance;
}

G4double G4ViewParameters::GetFarDistance (G4double cameraDistance,
					   G4double nearDistance,
					   G4double radius) const {
  G4double farDistance = cameraDistance + radius;
  if (farDistance < nearDistance) farDistance = nearDistance;
  return farDistance;
}

G4double G4ViewParameters::GetFrontHalfHeight (G4double nearDistance,
					       G4double radius) const {
  G4double frontHalfHeight;
  if (fFieldHalfAngle == 0.) {
    frontHalfHeight = radius / fZoomFactor;
  }
  else {
    frontHalfHeight = nearDistance * tan (fFieldHalfAngle) / fZoomFactor;
  }
  return frontHalfHeight;
}
// Useful quantities - end snippet.

void G4ViewParameters::AddCutawayPlane (const G4Plane3D& cutawayPlane) {
  fCutaway = true;
  if (fCutawayPlanes.entries () < 3 ) {
    fCutawayPlanes.insert (cutawayPlane);
  }
  else {
    G4cerr << "A maximum of 3 cutaway planes supported." << endl;
  }
}

void G4ViewParameters::ClearCutawayPlanes () {
  fCutaway = false;
  fCutawayPlanes.clear ();
}

void G4ViewParameters::SetVisibleDensity (G4double visibleDensity) {
  const G4double reasonableMaximum = 10.0 * g / cm3;
  if (visibleDensity < 0) {
    G4cout << "G4ViewParameters::SetVisibleDensity: attempt to set negative "
      "density - ignored." << endl;
  }
  else {
    if (fVisibleDensity > reasonableMaximum) {
      G4cout << "G4ViewParameters::SetVisibleDensity: density > "
	   << reasonableMaximum
	   << " g / cm3 - did you mean this?"
	   << endl;
    }
    fVisibleDensity = visibleDensity;
  }
}

void G4ViewParameters::SetNoOfSides (G4int nSides) {
  const G4int  nSidesMin = 3;
  if (nSides < nSidesMin) {
    nSides = nSidesMin;
    G4cout << "G4ViewParameters::SetNoOfSides: attempt to set the"
      "\nnumber of sides per circle < " << nSidesMin
	 << "; forced to" << nSides << endl;
  }
  fNoOfSides = nSides;
}

void G4ViewParameters::SetViewAndLights
(const G4Vector3D& viewpointDirection) {
  fViewpointDirection = viewpointDirection;
  if (fLightsMoveWithCamera) {
  G4Vector3D zprime = fViewpointDirection.unit ();
  G4Vector3D xprime = (fUpVector.cross (zprime)).unit ();
  G4Vector3D yprime = zprime.cross (xprime);
   fActualLightpointDirection =
     fRelativeLightpointDirection.x () * xprime +
     fRelativeLightpointDirection.y () * yprime +
     fRelativeLightpointDirection.x () * zprime;     
  }
  else {
    fActualLightpointDirection = fRelativeLightpointDirection;
  }
}

void G4ViewParameters::SetLightpointDirection
(const G4Vector3D& lightpointDirection) {
  fRelativeLightpointDirection = lightpointDirection;
  SetViewAndLights (fViewpointDirection);
}

void G4ViewParameters::Pan (G4double right, G4double up) {
  G4Vector3D unitRight = (fUpVector.cross (fViewpointDirection)).unit();
  G4Vector3D unitUp    = (fViewpointDirection.cross (unitRight)).unit();
  fCurrentTargetPoint += right * unitRight + up * unitUp;
}

void G4ViewParameters::PrintDifferences (const G4ViewParameters& v) const {

  // Put performance-sensitive parameters first.
  if (
      // This first to optimise spin, etc.
      (fViewpointDirection   != v.fViewpointDirection)   ||

      // No particular order from here on.
      (fDrawingStyle         != v.fDrawingStyle)         ||
      (fRepStyle             != v.fRepStyle)             ||
      (fCulling              != v.fCulling)              ||
      (fCullInvisible        != v.fCullInvisible)        ||
      (fDensityCulling       != v.fDensityCulling)       ||
      (fVisibleDensity       != v.fVisibleDensity)       ||
      (fCullCovered          != v.fCullCovered)          ||
      (fSection              != v.fSection)              ||
      (fCutaway              != v.fCutaway)              || 
      (fExplode              != v.fExplode)              ||
      (fNoOfSides            != v.fNoOfSides)            ||
      (fUpVector             != v.fUpVector)             ||
      (fFieldHalfAngle       != v.fFieldHalfAngle)       ||
      (fZoomFactor           != v.fZoomFactor)           ||
      (fCurrentTargetPoint   != v.fCurrentTargetPoint)   ||
      (fDolly                != v.fDolly)                ||
      (fRelativeLightpointDirection != v.fRelativeLightpointDirection)  ||
      (fLightsMoveWithCamera != v.fLightsMoveWithCamera) ||
      (fViewGeom             != v.fViewGeom)             ||
      (fViewHits             != v.fViewHits)             ||
      (fViewDigis            != v.fViewDigis)            ||
      (fDefaultVisAttributes != v.fDefaultVisAttributes) ||
      (fDefaultTextVisAttributes != v.fDefaultTextVisAttributes) ||
      (fDefaultMarker        != v.fDefaultMarker)        ||
      (fGlobalMarkerScale    != v.fGlobalMarkerScale)    ||
      (fMarkerNotHidden      != v.fMarkerNotHidden)      ||
      (fWindowSizeHintY      != v.fWindowSizeHintY))
    G4cout << "Difference in 1st batch." << endl;

  if (fSection) {
    if (!(fSectionPlane == v.fSectionPlane))
      G4cout << "Difference in section planes batch." << endl;
  }

  if (fCutaway) {
    if (fCutawayPlanes.entries () != v.fCutawayPlanes.entries ()) {
      G4cout << "Difference in no of cutaway planes." << endl;
    }
    else {
      for (int i = 0; i < fCutawayPlanes.entries (); i++) {
	if (!(fCutawayPlanes (i) == v.fCutawayPlanes (i)))
	  G4cout << "Difference in cutaway plane no. " << i << endl;
      }
    }
  }

  if (fExplode) {
    if (fExplodeFactor != v.fExplodeFactor)
      G4cout << "Difference in explode factor." << endl;
  }
}

ostream& operator << (ostream& os, const G4ViewParameters& v) {
  os << "View parameters and options:";

  os << "\n  Drawing style: ";
  switch (v.fDrawingStyle) {
  case G4ViewParameters::wireframe:
    os << "wireframe"; break;
  case G4ViewParameters::hlr:
    os << "hlr - hidden lines removed"; break;
  case G4ViewParameters::hsr:
    os << "hsr - hidden surfaces removed"; break;
  case G4ViewParameters::hlhsr:
    os << "hlhsr - hidden line, hidden surface removed"; break;
  default: os << "unrecognised"; break;
  }

  os << "\n  Representation style: ";
  switch (v.fRepStyle) {
  case G4ViewParameters::polyhedron:
    os << "polyhedron"; break;
  case G4ViewParameters::nurbs:
    os << "nurbs"; break;
  default: os << "unrecognised"; break;
  }

  os << "\n  Culling: ";
  if (v.fCulling) os << "on";
  else            os << "off";

  os << "\n  Culling invisible objects: ";
  if (v.fCullInvisible) os << "on";
  else                  os << "off";

  os << "\n  Density culling: ";
  if (v.fDensityCulling) {
    os << "on - invisible if density less than "
       << v.fVisibleDensity / (1. * g / cm3) << " g cm^-3";
  }
  else os << "off";

  os << "\n  Culling daughters covered by opaque mothers: ";
  if (v.fCullCovered) os << "on";
  else                os << "off";

  os << "\n  Section flag: ";
  if (v.fSection) os << "true, section/cut plane: " << v.fSectionPlane;
  else            os << "false";

  os << "\n  Cutaway flag: ";
  if (v.fCutaway) {
    os << "true, cutaway planes: ";
    for (int i = 0; i < v.fCutawayPlanes.entries (); i++) {
      os << ' ' << v.fCutawayPlanes[i];
    }
  }
  else {
    os << "false";
  }

  os << "\n  Explode flag: ";
  if (v.fExplode) os << "true, explode factor: " << v.fExplodeFactor;
  else            os << "false";

  os << "\n  No. of sides used in circle polygon approximation: "
     << v.fNoOfSides;

  os << "\n  Viewpoint direction:  " << v.fViewpointDirection;

  os << "\n  Up vector:            " << v.fUpVector;

  os << "\n  Field half angle:     " << v.fFieldHalfAngle;

  os << "\n  Zoom factor:          " << v.fZoomFactor;

  os << "\n  Current target point: " << v.fCurrentTargetPoint;

  os << "\n  Dolly distance:       " << v.fDolly;

  os << "\n  Light ";
  if (v.fLightsMoveWithCamera) os << "moves";
  else                         os << "does not move";
  os << " with camera";

  os << "\n  Relative lightpoint direction: "
     << v.fRelativeLightpointDirection;

  os << "\n  Actual lightpoint direction: "
     << v.fActualLightpointDirection;

  os << "\n  Derived parameters for standard view of object of unit radius:";
  G4ViewParameters tempVP = v;
  tempVP.fDolly = 0.;
  tempVP.fZoomFactor = 1.;
  const G4double radius = 1.;
  const G4double cameraDistance = tempVP.GetCameraDistance (radius);
  const G4double nearDistance =
    tempVP.GetNearDistance (cameraDistance, radius);
  const G4double farDistance =
    tempVP.GetFarDistance  (cameraDistance, nearDistance, radius);
  const G4double right  = tempVP.GetFrontHalfHeight (nearDistance, radius);
  os << "\n    Camera distance:   " << cameraDistance;
  os << "\n    Near distance:     " << nearDistance;
  os << "\n    Far distance:      " << farDistance;
  os << "\n    Front half height: " << right;

  os << "\n  View geometry: ";
  if (v.fViewGeom) os << "true";
  else os << "false";

  os << "\n  View hits    : ";
  if (v.fViewHits) os << "true";
  else os << "false";

  os << "\n  View digits  : ";
  if (v.fViewDigis) os << "true";
  else os << "false";

  os << "\n  Default VisAttributes:\n  " << v.fDefaultVisAttributes;

  os << "\n  Default TextVisAttributes:\n  " << v.fDefaultTextVisAttributes;

  os << "\n  Default marker: " << v.fDefaultMarker;

  os << "\n  Global marker scale: " << v.fGlobalMarkerScale;

  os << "\n  Marker ";
  if (v.fMarkerNotHidden) os << "not";
  os << " hidden by surfaces.";

  os << "\n  Window size hint: "
     << v.fWindowSizeHintX << 'x'<< v.fWindowSizeHintX;

  return os;
}

G4bool operator != (const G4ViewParameters& v1,
		    const G4ViewParameters& v2) {

  // Put performance-sensitive parameters first.
  if (
      // This first to optimise spin, etc.
      (v1.fViewpointDirection   != v2.fViewpointDirection)   ||

      // No particular order from here on.
      (v1.fDrawingStyle         != v2.fDrawingStyle)         ||
      (v1.fRepStyle             != v2.fRepStyle)             ||
      (v1.fCulling              != v2.fCulling)              ||
      (v1.fCullInvisible        != v2.fCullInvisible)        ||
      (v1.fDensityCulling       != v2.fDensityCulling)       ||
      (v1.fCullCovered          != v2.fCullCovered)          ||
      (v1.fSection              != v2.fSection)              ||
      (v1.fCutaway              != v2.fCutaway)              || 
      (v1.fExplode              != v2.fExplode)              ||
      (v1.fNoOfSides            != v2.fNoOfSides)            ||
      (v1.fUpVector             != v2.fUpVector)             ||
      (v1.fFieldHalfAngle       != v2.fFieldHalfAngle)       ||
      (v1.fZoomFactor           != v2.fZoomFactor)           ||
      (v1.fCurrentTargetPoint   != v2.fCurrentTargetPoint)   ||
      (v1.fDolly                != v2.fDolly)                ||
      (v1.fRelativeLightpointDirection != v2.fRelativeLightpointDirection)  ||
      (v1.fLightsMoveWithCamera != v2.fLightsMoveWithCamera) ||
      (v1.fViewGeom             != v2.fViewGeom)             ||
      (v1.fViewHits             != v2.fViewHits)             ||
      (v1.fViewDigis            != v2.fViewDigis)            ||
      (v1.fDefaultVisAttributes != v2.fDefaultVisAttributes) ||
      (v1.fDefaultTextVisAttributes != v2.fDefaultTextVisAttributes) ||
      (v1.fDefaultMarker        != v2.fDefaultMarker)        ||
      (v1.fGlobalMarkerScale    != v2.fGlobalMarkerScale)    ||
      (v1.fMarkerNotHidden      != v2.fMarkerNotHidden)      ||
      (v1.fWindowSizeHintX      != v2.fWindowSizeHintX)      ||
      (v1.fWindowSizeHintY      != v2.fWindowSizeHintY))
    return true;

  if (v1.fDensityCulling &&
      (v1.fVisibleDensity != v2.fVisibleDensity)) return true;

  if (v1.fSection &&
      (!(v1.fSectionPlane == v2.fSectionPlane))) return true;

  if (v1.fCutaway) {
    if (v1.fCutawayPlanes.entries () != v2.fCutawayPlanes.entries ())
      return true;
    else {
      for (int i = 0; i < v1.fCutawayPlanes.entries (); i++) {
	if (!(v1.fCutawayPlanes (i) == v2.fCutawayPlanes (i))) return true;
      }
    }
  }

  if (v1.fExplode &&
      (v1.fExplodeFactor != v2.fExplodeFactor)) return true;

  return false;
}

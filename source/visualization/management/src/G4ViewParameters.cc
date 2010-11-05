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
// $Id: G4ViewParameters.cc,v 1.38 2010-11-05 16:00:11 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996
// View parameters and options.

#include <sstream>

#include "G4ViewParameters.hh"

#include "G4VisManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

G4ViewParameters::G4ViewParameters ():
  fNoValue(0x0000),
  fXValue(0x0001),
  fYValue(0x0002),
  fWidthValue(0x0004),
  fHeightValue(0x0008),
  fAllValues(0x000F),
  fXNegative(0x0010),
  fYNegative(0x0020),
  fGeometryMask(0),
  fDrawingStyle (wireframe),
  fAuxEdgeVisible (false),
  fRepStyle (polyhedron),
  fCulling (true),
  fCullInvisible (true),
  fDensityCulling (false),
  fVisibleDensity (0.01 * g / cm3),
  fCullCovered (false),
  fSection (false),
  fSectionPlane (),
  fCutawayMode (cutawayUnion),
  fCutawayPlanes (),
  fExplodeFactor (1.),
  fNoOfSides (24),
  fViewpointDirection (G4Vector3D (0., 0., 1.)),  // On z-axis.
  fUpVector (G4Vector3D (0., 1., 0.)),            // y-axis up.
  fFieldHalfAngle (0.),                           // Orthogonal projection.
  fZoomFactor (1.),
  fScaleFactor (G4Vector3D (1., 1., 1.)),
  fCurrentTargetPoint (),
  fDolly (0.),
  fLightsMoveWithCamera (false),
  fRelativeLightpointDirection (G4Vector3D (1., 1., 1.)),
  fActualLightpointDirection (G4Vector3D (1., 1., 1.)),
  fDefaultVisAttributes (),
  fDefaultTextVisAttributes (G4Colour (0., 0., 1.)),
  fDefaultMarker (),
  fGlobalMarkerScale (1.),
  fGlobalLineWidthScale (1.),
  fMarkerNotHidden (true),
  fWindowSizeHintX (600),
  fWindowSizeHintY (600),
  fWindowLocationHintX(0),
  fWindowLocationHintY(0),
  fWindowLocationHintXNegative(true),
  fWindowLocationHintYNegative(false),
  fAutoRefresh (false),
  fBackgroundColour (G4Colour(0.,0.,0.)),         // Black
  fPicking (false)
{
  fDefaultMarker.SetScreenSize (5.);
  // Markers are 5 pixels "overall" size, i.e., diameter.
}

G4ViewParameters::~G4ViewParameters () {}

void G4ViewParameters::MultiplyScaleFactor
(const G4Vector3D& scaleFactorMultiplier) {
  fScaleFactor.setX(fScaleFactor.x() * scaleFactorMultiplier.x());
  fScaleFactor.setY(fScaleFactor.y() * scaleFactorMultiplier.y());
  fScaleFactor.setZ(fScaleFactor.z() * scaleFactorMultiplier.z());
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
    cameraDistance = radius / std::sin (fFieldHalfAngle) - fDolly;
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
    frontHalfHeight = nearDistance * std::tan (fFieldHalfAngle) / fZoomFactor;
  }
  return frontHalfHeight;
}
// Useful quantities - end snippet.

void G4ViewParameters::AddCutawayPlane (const G4Plane3D& cutawayPlane) {
  if (fCutawayPlanes.size () < 3 ) {
    fCutawayPlanes.push_back (cutawayPlane);
  }
  else {
    G4cout <<
      "ERROR: G4ViewParameters::AddCutawayPlane:"
      "\n  A maximum of 3 cutaway planes supported." << G4endl;
  }
}

void G4ViewParameters::ChangeCutawayPlane
(size_t index, const G4Plane3D& cutawayPlane) {
  if (index >= fCutawayPlanes.size()) {
    G4cout <<
      "ERROR: G4ViewParameters::ChangeCutawayPlane:"
      "\n  Plane " << index << " does not exist." << G4endl;
  } else {
    fCutawayPlanes[index] = cutawayPlane;
  }
}

void G4ViewParameters::SetVisibleDensity (G4double visibleDensity) {
  const G4double reasonableMaximum = 10.0 * g / cm3;
  if (visibleDensity < 0) {
    G4cout << "G4ViewParameters::SetVisibleDensity: attempt to set negative "
      "density - ignored." << G4endl;
  }
  else {
    if (visibleDensity > reasonableMaximum) {
      G4cout << "G4ViewParameters::SetVisibleDensity: density > "
	     << G4BestUnit (reasonableMaximum, "Volumic Mass")
	     << " - did you mean this?"
	     << G4endl;
    }
    fVisibleDensity = visibleDensity;
  }
}

G4int G4ViewParameters::SetNoOfSides (G4int nSides) {
  const G4int  nSidesMin = 12;
  if (nSides < nSidesMin) {
    nSides = nSidesMin;
    G4cout << "G4ViewParameters::SetNoOfSides: attempt to set the"
      "\nnumber of sides per circle < " << nSidesMin
	 << "; forced to " << nSides << G4endl;
  }
  fNoOfSides = nSides;
  return fNoOfSides;
}

void G4ViewParameters::SetViewAndLights
(const G4Vector3D& viewpointDirection) {

  fViewpointDirection = viewpointDirection;

  // If the requested viewpoint direction is parallel to the up
  // vector, the orientation of the view is undefined...
  if (fViewpointDirection.unit() * fUpVector.unit() > .9999) {
    G4cout <<
      "WARNING: Viewpoint direction is very close to the up vector direction."
      "\n  Consider setting the up vector to obtain definable behaviour."
	   << G4endl;
  }

  // Move the lights too if requested...
  if (fLightsMoveWithCamera) {
    G4Vector3D zprime = fViewpointDirection.unit ();
    G4Vector3D xprime = (fUpVector.cross (zprime)).unit ();
    G4Vector3D yprime = zprime.cross (xprime);
    fActualLightpointDirection =
      fRelativeLightpointDirection.x () * xprime +
      fRelativeLightpointDirection.y () * yprime +
      fRelativeLightpointDirection.x () * zprime;     
  } else {
    fActualLightpointDirection = fRelativeLightpointDirection;
  }
}

void G4ViewParameters::SetLightpointDirection
(const G4Vector3D& lightpointDirection) {
  fRelativeLightpointDirection = lightpointDirection;
  SetViewAndLights (fViewpointDirection);
}

void G4ViewParameters::SetPan (G4double right, G4double up) {
  G4Vector3D unitRight = (fUpVector.cross (fViewpointDirection)).unit();
  G4Vector3D unitUp    = (fViewpointDirection.cross (unitRight)).unit();
  fCurrentTargetPoint  = right * unitRight + up * unitUp;
}

void G4ViewParameters::IncrementPan (G4double right, G4double up) {
  IncrementPan (right,up, 0);
}

void G4ViewParameters::IncrementPan (G4double right, G4double up, G4double distance) {
  G4Vector3D unitRight = (fUpVector.cross (fViewpointDirection)).unit();
  G4Vector3D unitUp    = (fViewpointDirection.cross (unitRight)).unit();
  fCurrentTargetPoint += right * unitRight + up * unitUp + distance * fViewpointDirection;
}

void G4ViewParameters::PrintDifferences (const G4ViewParameters& v) const {

  // Put performance-sensitive parameters first.
  if (
      // This first to optimise spin, etc.
      (fViewpointDirection   != v.fViewpointDirection)   ||

      // No particular order from here on.
      (fDrawingStyle         != v.fDrawingStyle)         ||
      (fAuxEdgeVisible       != v.fAuxEdgeVisible)       ||
      (fRepStyle             != v.fRepStyle)             ||
      (fCulling              != v.fCulling)              ||
      (fCullInvisible        != v.fCullInvisible)        ||
      (fDensityCulling       != v.fDensityCulling)       ||
      (fVisibleDensity       != v.fVisibleDensity)       ||
      (fCullCovered          != v.fCullCovered)          ||
      (fSection              != v.fSection)              ||
      (fNoOfSides            != v.fNoOfSides)            ||
      (fUpVector             != v.fUpVector)             ||
      (fFieldHalfAngle       != v.fFieldHalfAngle)       ||
      (fZoomFactor           != v.fZoomFactor)           ||
      (fScaleFactor          != v.fScaleFactor)          ||
      (fCurrentTargetPoint   != v.fCurrentTargetPoint)   ||
      (fDolly                != v.fDolly)                ||
      (fRelativeLightpointDirection != v.fRelativeLightpointDirection)  ||
      (fLightsMoveWithCamera != v.fLightsMoveWithCamera) ||
      (fDefaultVisAttributes != v.fDefaultVisAttributes) ||
      (fDefaultTextVisAttributes != v.fDefaultTextVisAttributes) ||
      (fDefaultMarker        != v.fDefaultMarker)        ||
      (fGlobalMarkerScale    != v.fGlobalMarkerScale)    ||
      (fGlobalLineWidthScale != v.fGlobalLineWidthScale) ||
      (fMarkerNotHidden      != v.fMarkerNotHidden)      ||
      (fWindowSizeHintX      != v.fWindowSizeHintX)      ||
      (fWindowSizeHintY      != v.fWindowSizeHintY)      ||
      (fXGeometryString      != v.fXGeometryString)      ||
      (fAutoRefresh          != v.fAutoRefresh)          ||
      (fBackgroundColour     != v.fBackgroundColour)     || 
      (fPicking              != v.fPicking)
      )
    G4cout << "Difference in 1st batch." << G4endl;

  if (fSection) {
    if (!(fSectionPlane == v.fSectionPlane))
      G4cout << "Difference in section planes batch." << G4endl;
  }

  if (IsCutaway()) {
    if (fCutawayPlanes.size () != v.fCutawayPlanes.size ()) {
      G4cout << "Difference in no of cutaway planes." << G4endl;
    }
    else {
      for (size_t i = 0; i < fCutawayPlanes.size (); i++) {
	if (!(fCutawayPlanes[i] == v.fCutawayPlanes[i]))
	  G4cout << "Difference in cutaway plane no. " << i << G4endl;
      }
    }
  }

  if (IsExplode()) {
    if (fExplodeFactor != v.fExplodeFactor)
      G4cout << "Difference in explode factor." << G4endl;
    if (fExplodeCentre != v.fExplodeCentre)
      G4cout << "Difference in explode centre." << G4endl;
  }
}

std::ostream& operator << (std::ostream& os,
			     const G4ViewParameters::DrawingStyle& style) {
  switch (style) {
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
  return os;
}

std::ostream& operator << (std::ostream& os, const G4ViewParameters& v) {
  os << "View parameters and options:";

  os << "\n  Drawing style: " << v.fDrawingStyle;

  os << "\n  Auxiliary edges: ";
  if (!v.fAuxEdgeVisible) os << "in";
  os << "visible";

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

  if (v.IsCutaway()) {
    os << "\n  Cutaway planes: ";
    for (size_t i = 0; i < v.fCutawayPlanes.size (); i++) {
      os << ' ' << v.fCutawayPlanes[i];
    }
  }
  else {
    os << "\n  No cutaway planes";
  }

  os << "\n  Explode factor: " << v.fExplodeFactor
     << " about centre: " << v.fExplodeCentre;

  os << "\n  No. of sides used in circle polygon approximation: "
     << v.fNoOfSides;

  os << "\n  Viewpoint direction:  " << v.fViewpointDirection;

  os << "\n  Up vector:            " << v.fUpVector;

  os << "\n  Field half angle:     " << v.fFieldHalfAngle;

  os << "\n  Zoom factor:          " << v.fZoomFactor;

  os << "\n  Scale factor:         " << v.fScaleFactor;

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

  os << "\n  Default VisAttributes:\n  " << v.fDefaultVisAttributes;

  os << "\n  Default TextVisAttributes:\n  " << v.fDefaultTextVisAttributes;

  os << "\n  Default marker: " << v.fDefaultMarker;

  os << "\n  Global marker scale: " << v.fGlobalMarkerScale;

  os << "\n  Global lineWidth scale: " << v.fGlobalLineWidthScale;

  os << "\n  Marker ";
  if (v.fMarkerNotHidden) os << "not ";
  os << "hidden by surfaces.";

  os << "\n  Window size hint: "
     << v.fWindowSizeHintX << 'x'<< v.fWindowSizeHintX;

  os << "\n  X geometry string: " << v.fXGeometryString;

  os << "\n  Auto refresh: ";
  if (v.fAutoRefresh) os << "true";
  else os << "false";

  os << "\n  Background colour: " << v.fBackgroundColour;

  os << "\n  Picking requested: ";
  if (v.fPicking) os << "true";
  else os << "false";

  return os;
}

G4bool G4ViewParameters::operator != (const G4ViewParameters& v) const {

  // Put performance-sensitive parameters first.
  if (
      // This first to optimise spin, etc.
      (fViewpointDirection   != v.fViewpointDirection)   ||

      // No particular order from here on.
      (fDrawingStyle         != v.fDrawingStyle)         ||
      (fAuxEdgeVisible       != v.fAuxEdgeVisible)       ||
      (fRepStyle             != v.fRepStyle)             ||
      (fCulling              != v.fCulling)              ||
      (fCullInvisible        != v.fCullInvisible)        ||
      (fDensityCulling       != v.fDensityCulling)       ||
      (fCullCovered          != v.fCullCovered)          ||
      (fSection              != v.fSection)              ||
      (IsCutaway()           != v.IsCutaway())           ||
      (IsExplode()           != v.IsExplode())           ||
      (fNoOfSides            != v.fNoOfSides)            ||
      (fUpVector             != v.fUpVector)             ||
      (fFieldHalfAngle       != v.fFieldHalfAngle)       ||
      (fZoomFactor           != v.fZoomFactor)           ||
      (fScaleFactor          != v.fScaleFactor)          ||
      (fCurrentTargetPoint   != v.fCurrentTargetPoint)   ||
      (fDolly                != v.fDolly)                ||
      (fRelativeLightpointDirection != v.fRelativeLightpointDirection)  ||
      (fLightsMoveWithCamera != v.fLightsMoveWithCamera) ||
      (fDefaultVisAttributes != v.fDefaultVisAttributes) ||
      (fDefaultTextVisAttributes != v.fDefaultTextVisAttributes) ||
      (fDefaultMarker        != v.fDefaultMarker)        ||
      (fGlobalMarkerScale    != v.fGlobalMarkerScale)    ||
      (fGlobalLineWidthScale != v.fGlobalLineWidthScale) ||
      (fMarkerNotHidden      != v.fMarkerNotHidden)      ||
      (fWindowSizeHintX      != v.fWindowSizeHintX)      ||
      (fWindowSizeHintY      != v.fWindowSizeHintY)      ||
      (fXGeometryString      != v.fXGeometryString)      ||
      (fAutoRefresh          != v.fAutoRefresh)          ||
      (fBackgroundColour     != v.fBackgroundColour)     ||
      (fPicking              != v.fPicking)
      )
    return true;

  if (fDensityCulling &&
      (fVisibleDensity != v.fVisibleDensity)) return true;

  if (fSection &&
      (!(fSectionPlane == v.fSectionPlane))) return true;

  if (IsCutaway()) {
    if (fCutawayPlanes.size () != v.fCutawayPlanes.size ())
      return true;
    else {
      for (size_t i = 0; i < fCutawayPlanes.size (); i++) {
	if (!(fCutawayPlanes[i] == v.fCutawayPlanes[i])) return true;
      }
    }
  }

  if (IsExplode() &&
      ((fExplodeFactor != v.fExplodeFactor) ||
       (fExplodeCentre != v.fExplodeCentre))) return true;

  return false;
}


void G4ViewParameters::SetXGeometryString (const G4String& geomStringArg) {


  G4int x,y = 0;
  unsigned int w,h = 0;
  G4String geomString = geomStringArg;
  // Parse windowSizeHintString for backwards compatibility...
  const G4String delimiters("xX+-");
  G4String::size_type i = geomString.find_first_of(delimiters);
  if (i == G4String::npos) {  // Does not contain "xX+-".  Assume single number
    std::istringstream iss(geomString);
    G4int size;
    iss >> size;
    if (!iss) {
      size = 600;
      G4cout << "Unrecognised windowSizeHint string: \""
	     << geomString
	     << "\".  Asuuming " << size << G4endl;
    }
    std::ostringstream oss;
    oss << size << 'x' << size;
    geomString = oss.str();
  }
 
  fGeometryMask = ParseGeometry( geomString, &x, &y, &w, &h );

  // Handle special case :
  if ((fGeometryMask & fYValue) == 0)
    {  // Using default
      y =  fWindowLocationHintY;
    }
  if ((fGeometryMask & fXValue) == 0)
    {  // Using default
      x =  fWindowLocationHintX;
    }

  // Check errors
  // if there is no Width and Height
  if ( ((fGeometryMask & fHeightValue) == 0 ) &&
       ((fGeometryMask & fWidthValue)  == 0 )) {
    h = fWindowSizeHintY;
    w = fWindowSizeHintX;
  } else  if ((fGeometryMask & fHeightValue) == 0 ) {

    // if there is only Width. Special case to be backward compatible
    // We set Width and Height the same to obtain a square windows.
    
    G4cout << "Unrecognised geometry string \""
           << geomString
           << "\".  No Height found. Using Width value instead"
           << G4endl;
    h = w;
  }
  if ( ((fGeometryMask & fXValue) == 0 ) ||
       ((fGeometryMask & fYValue)  == 0 )) {
    //Using defaults
    x = fWindowLocationHintX;
    y = fWindowLocationHintY;
  }
  // Set the string
  fXGeometryString = geomString;
  
  // Set values
  fWindowSizeHintX = w;
  fWindowSizeHintY = h;
  fWindowLocationHintX = x;
  fWindowLocationHintY = y;

  if ( ((fGeometryMask & fXValue)) &&
       ((fGeometryMask & fYValue))) {

    if ( (fGeometryMask & fXNegative) ) {
      fWindowLocationHintXNegative = true;
    } else {
      fWindowLocationHintXNegative = false;
    }
    if ( (fGeometryMask & fYNegative) ) {
    fWindowLocationHintYNegative = true;
    } else {
      fWindowLocationHintYNegative = false;
    }
  }
}

G4int G4ViewParameters::GetWindowAbsoluteLocationHintX (G4int sizeX ) const {
  if ( fWindowLocationHintXNegative ) {
    return  sizeX  + fWindowLocationHintX - fWindowSizeHintX;
  }
  return fWindowLocationHintX;
}

G4int G4ViewParameters::GetWindowAbsoluteLocationHintY (G4int sizeY ) const {
  if (  fWindowLocationHintYNegative ) {
    return  sizeY  + fWindowLocationHintY - fWindowSizeHintY;
  }
  return fWindowLocationHintY;
}

/* Keep from :
 * ftp://ftp.trolltech.com/qt/source/qt-embedded-free-3.0.6.tar.gz/qt-embedded-free-3.0.6/src/kernel/qapplication_qws.cpp
 *
 *    ParseGeometry parses strings of the form
 *   "=<width>x<height>{+-}<xoffset>{+-}<yoffset>", where
 *   width, height, xoffset, and yoffset are unsigned integers.
 *   Example:  "=80x24+300-49"
 *   The equal sign is optional.
 *   It returns a bitmask that indicates which of the four values
 *   were actually found in the string. For each value found,
 *   the corresponding argument is updated;  for each value
 *   not found, the corresponding argument is left unchanged.
 */

int G4ViewParameters::ParseGeometry (
 const char *string,
 G4int *x,
 G4int *y,
 unsigned int *width,
 unsigned int *height)
{

  G4int mask = fNoValue;
  register char *strind;
  unsigned int tempWidth  = 0;
  unsigned int tempHeight = 0;
  G4int tempX = 0;
  G4int tempY = 0;
  char *nextCharacter;
  if ( (string == NULL) || (*string == '\0')) {
    return(mask);
  }
  if (*string == '=')
    string++;  /* ignore possible '=' at beg of geometry spec */
  strind = (char *)string;
  if (*strind != '+' && *strind != '-' && *strind != 'x') {
    tempWidth = ReadInteger(strind, &nextCharacter);
    if (strind == nextCharacter)
      return (0);
    strind = nextCharacter;
    mask |= fWidthValue;
  }
  if (*strind == 'x' || *strind == 'X') {
    strind++;
    tempHeight = ReadInteger(strind, &nextCharacter);
    if (strind == nextCharacter)
      return (0);
    strind = nextCharacter;
    mask |= fHeightValue;
  }

  if ((*strind == '+') || (*strind == '-')) {
    if (*strind == '-') {
      strind++;
      tempX = -ReadInteger(strind, &nextCharacter);
      if (strind == nextCharacter)
        return (0);
      strind = nextCharacter;
      mask |= fXNegative;

    }
    else
      {	strind++;
        tempX = ReadInteger(strind, &nextCharacter);
        if (strind == nextCharacter)
          return(0);
        strind = nextCharacter;
      }
    mask |= fXValue;
    if ((*strind == '+') || (*strind == '-')) {
      if (*strind == '-') {
        strind++;
        tempY = -ReadInteger(strind, &nextCharacter);
        if (strind == nextCharacter)
          return(0);
        strind = nextCharacter;
        mask |= fYNegative;
      }
      else
        {
          strind++;
          tempY = ReadInteger(strind, &nextCharacter);
          if (strind == nextCharacter)
            return(0);
          strind = nextCharacter;
        }
      mask |= fYValue;
    }
  }
  /* If strind isn't at the end of the string the it's an invalid
     geometry specification. */
  if (*strind != '\0') return (0);
  if (mask & fXValue)
    *x = tempX;
  if (mask & fYValue)
    *y = tempY;
  if (mask & fWidthValue)
    *width = tempWidth;
  if (mask & fHeightValue)
    *height = tempHeight;
  return (mask);
}

/* Keep from :
 * ftp://ftp.trolltech.com/qt/source/qt-embedded-free-3.0.6.tar.gz/qt-embedded-free-3.0.6/src/kernel/qapplication_qws.cpp
 *
 */
G4int G4ViewParameters::ReadInteger(char *string, char **NextString)
{
    register G4int Result = 0;
    G4int Sign = 1;

    if (*string == '+')
	string++;
    else if (*string == '-')
    {
	string++;
	Sign = -1;
    }
    for (; (*string >= '0') && (*string <= '9'); string++)
    {
	Result = (Result * 10) + (*string - '0');
    }
    *NextString = string;
    if (Sign >= 0)
	return (Result);
    else
 	return (-Result);
}
   

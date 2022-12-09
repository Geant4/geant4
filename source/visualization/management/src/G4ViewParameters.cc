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
// View parameters and options.

#include "G4ViewParameters.hh"

#include "G4VisManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedron.hh"

#include <sstream>
#include <cmath>

#define G4warn G4cout

G4ViewParameters::G4ViewParameters ():
  fDrawingStyle (wireframe),
  fNumberOfCloudPoints(10000),
  fAuxEdgeVisible (false),
  fCulling (true),
  fCullInvisible (true),
  fDensityCulling (false),
  fVisibleDensity (0.01 * g / cm3),
  fCullCovered (false),
  fCBDAlgorithmNumber (0),
  fSection (false),
  fSectionPlane (),
  fCutawayMode (cutawayUnion),
  fCutawayPlanes (),
  fExplodeFactor (1.),
  fNoOfSides (),
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
  fGeometryMask(0),
  fAutoRefresh (false),
  fBackgroundColour (G4Colour(0.,0.,0.)),         // Black
  fPicking (false),
  fRotationStyle (constrainUpDirection),
  fStartTime(-G4VisAttributes::fVeryLongTime),
  fEndTime(G4VisAttributes::fVeryLongTime),
  fFadeFactor(0.),
  fDisplayHeadTime(false),
  fDisplayHeadTimeX(-0.9),
  fDisplayHeadTimeY(-0.9),
  fDisplayHeadTimeSize(24.),
  fDisplayHeadTimeRed(0.),
  fDisplayHeadTimeGreen(1.),
  fDisplayHeadTimeBlue(1.),
  fDisplayLightFront(false),
  fDisplayLightFrontX(0.),
  fDisplayLightFrontY(0.),
  fDisplayLightFrontZ(0.),
  fDisplayLightFrontT(0.),
  fDisplayLightFrontRed(0.),
  fDisplayLightFrontGreen(1.),
  fDisplayLightFrontBlue(0.),
  fSpecialMeshRendering(false),
  fSpecialMeshRenderingOption(meshAsDots)
{
  // Pick up default no of sides from G4Polyhedron.
  // Note that this parameter is variously called:
  //   No of sides
  //   NumberOfRotationSteps
  //   Line segments per circle
  // It refers to the approximation of a circle by a polygon of
  // stated number of sides.
  fNoOfSides = G4Polyhedron::GetNumberOfRotationSteps();
  
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
    G4warn <<
      "ERROR: G4ViewParameters::AddCutawayPlane:"
      "\n  A maximum of 3 cutaway planes supported." << G4endl;
  }
}

void G4ViewParameters::ChangeCutawayPlane
(size_t index, const G4Plane3D& cutawayPlane) {
  if (index >= fCutawayPlanes.size()) {
    G4warn <<
      "ERROR: G4ViewParameters::ChangeCutawayPlane:"
      "\n  Plane " << index << " does not exist." << G4endl;
  } else {
    fCutawayPlanes[index] = cutawayPlane;
  }
}

void G4ViewParameters::SetVisibleDensity (G4double visibleDensity) {
  const G4double reasonableMaximum = 10.0 * g / cm3;
  if (visibleDensity < 0) {
    G4warn << "G4ViewParameters::SetVisibleDensity: attempt to set negative "
      "density - ignored." << G4endl;
  }
  else {
    if (visibleDensity > reasonableMaximum) {
      G4warn << "G4ViewParameters::SetVisibleDensity: density > "
	     << G4BestUnit (reasonableMaximum, "Volumic Mass")
	     << " - did you mean this?"
	     << G4endl;
    }
    fVisibleDensity = visibleDensity;
  }
}

G4int G4ViewParameters::SetNoOfSides (G4int nSides) {
  const G4int nSidesMin = fDefaultVisAttributes.GetMinLineSegmentsPerCircle();
  if (nSides < nSidesMin) {
    nSides = nSidesMin;
    G4warn << "G4ViewParameters::SetNoOfSides: attempt to set the"
    "\nnumber of sides per circle < " << nSidesMin
    << "; forced to " << nSides << G4endl;
  }
  fNoOfSides = nSides;
  return fNoOfSides;
}

G4int G4ViewParameters::SetNumberOfCloudPoints(G4int nPoints) {
  const G4int nPointsMin = 100;
  if (nPoints < nPointsMin) {
    nPoints = nPointsMin;
    G4warn << "G4ViewParameters::SetNumberOfCloudPoints:"
    "\nnumber of points per cloud set to minimum " << nPoints
    << G4endl;
  }
  fNumberOfCloudPoints = nPoints;
  return fNumberOfCloudPoints;
}

void G4ViewParameters::SetViewAndLights
(const G4Vector3D& viewpointDirection) {

  fViewpointDirection = viewpointDirection;

  // If the requested viewpoint direction is parallel to the up
  // vector, the orientation of the view is undefined...
  if (fViewpointDirection.unit() * fUpVector.unit() > .9999) {
    static G4bool firstTime = true;
    if (firstTime) {
      firstTime = false;
      G4warn <<
      "WARNING: Viewpoint direction is very close to the up vector direction."
      "\n  Change the up vector or \"/vis/viewer/set/rotationStyle freeRotation\"."
      << G4endl;
    }
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

void G4ViewParameters::AddVisAttributesModifier
(const G4ModelingParameters::VisAttributesModifier& vam) {
  // If target exists with same signifier just change vis attributes.
  G4bool duplicateTarget = false;
  auto i = fVisAttributesModifiers.begin();
  for (; i < fVisAttributesModifiers.end(); ++i) {
    if (vam.GetPVNameCopyNoPath() == (*i).GetPVNameCopyNoPath() &&
        vam.GetVisAttributesSignifier() == (*i).GetVisAttributesSignifier()) {
      duplicateTarget = true;
      break;
    }
  }
  if (duplicateTarget) (*i).SetVisAttributes(vam.GetVisAttributes());
  else fVisAttributesModifiers.push_back(vam);
}

G4String G4ViewParameters::CameraAndLightingCommands
(const G4Point3D standardTargetPoint) const
{
  std::ostringstream oss;

  oss << "#\n# Camera and lights commands";
  
  oss << "\n/vis/viewer/set/viewpointVector "
  << fViewpointDirection.x()
  << ' ' << fViewpointDirection.y()
  << ' ' << fViewpointDirection.z();
  
  oss << "\n/vis/viewer/set/upVector "
  << fUpVector.x()
  << ' ' << fUpVector.y()
  << ' ' << fUpVector.z();
  
  oss << "\n/vis/viewer/set/projection ";
    if (fFieldHalfAngle == 0.) {
    oss
    << "orthogonal";
  } else {
    oss
    << "perspective "
    << fFieldHalfAngle/deg
    << " deg";
  }
  
  oss << "\n/vis/viewer/zoomTo "
  << fZoomFactor;
  
  oss << "\n/vis/viewer/scaleTo "
  << fScaleFactor.x()
  << ' ' << fScaleFactor.y()
  << ' ' << fScaleFactor.z();
  
  oss << "\n/vis/viewer/set/targetPoint "
  << G4BestUnit(standardTargetPoint+fCurrentTargetPoint,"Length")
  << "\n# Note that if you have not set a target point, the vis system sets"
  << "\n# a target point based on the scene - plus any panning and dollying -"
  << "\n# so don't be alarmed by strange coordinates here.";
  
  oss << "\n/vis/viewer/dollyTo "
  << G4BestUnit(fDolly,"Length");
  
  oss << "\n/vis/viewer/set/lightsMove ";
  if (fLightsMoveWithCamera) {
    oss << "camera";
  } else {
    oss << "object";
  }
  
  oss << "\n/vis/viewer/set/lightsVector "
  << fRelativeLightpointDirection.x()
  << ' ' << fRelativeLightpointDirection.y()
  << ' ' << fRelativeLightpointDirection.z();
  
  oss << "\n/vis/viewer/set/rotationStyle ";
  if (fRotationStyle == constrainUpDirection) {
    oss << "constrainUpDirection";
  } else {
    oss << "freeRotation";
  }
  
  G4Colour c = fBackgroundColour;
  oss << "\n/vis/viewer/set/background "
  << c.GetRed()
  << ' ' << c.GetGreen()
  << ' ' << c.GetBlue()
  << ' ' << c.GetAlpha();
  
  c = fDefaultVisAttributes.GetColour();
  oss << "\n/vis/viewer/set/defaultColour "
  << c.GetRed()
  << ' ' << c.GetGreen()
  << ' ' << c.GetBlue()
  << ' ' << c.GetAlpha();
  
  c = fDefaultTextVisAttributes.GetColour();
  oss << "\n/vis/viewer/set/defaultTextColour "
  << c.GetRed()
  << ' ' << c.GetGreen()
  << ' ' << c.GetBlue()
  << ' ' << c.GetAlpha();
  
  oss << std::endl;
  
  return oss.str();
}

G4String G4ViewParameters::DrawingStyleCommands() const
{
  std::ostringstream oss;
  
  oss << "#\n# Drawing style commands";
  
  oss << "\n/vis/viewer/set/style ";
  switch (fDrawingStyle) {
    case wireframe:
    case hlr:
      oss << "wireframe";
      break;
    case hsr:
    case hlhsr:
      oss << "surface";
      break;
    case cloud:
      oss << "cloud";
      break;
  }

  oss << "\n/vis/viewer/set/hiddenEdge ";
  if (fDrawingStyle == hlr || fDrawingStyle == hlhsr) {
    oss << "true";
  } else {
    oss << "false";
  }
  
  oss << "\n/vis/viewer/set/auxiliaryEdge ";
  if (fAuxEdgeVisible) {
    oss << "true";
  } else {
    oss << "false";
  }
  
  oss << "\n/vis/viewer/set/hiddenMarker ";
  if (fMarkerNotHidden) {
    oss << "false";
  } else {
    oss << "true";
  }
  
  oss << "\n/vis/viewer/set/globalLineWidthScale "
  << fGlobalLineWidthScale;
  
  oss << "\n/vis/viewer/set/globalMarkerScale "
  << fGlobalMarkerScale;

  oss << "\n/vis/viewer/set/numberOfCloudPoints "
  << fNumberOfCloudPoints;

  oss << "\n/vis/viewer/set/specialMeshRendering ";
  if (fSpecialMeshRendering) {
    oss << "true";
  } else {
    oss << "false";
  }

  oss << "\n/vis/viewer/set/specialMeshRenderingOption "
  << fSpecialMeshRenderingOption;

  oss << "\n/vis/viewer/set/specialMeshVolumes";
  for (const auto& volume : fSpecialMeshVolumes) {
    oss << ' ' << volume.GetName() << ' ' << volume.GetCopyNo();
  }

  oss << std::endl;
  
  return oss.str();
}

G4String G4ViewParameters::SceneModifyingCommands() const
{
  std::ostringstream oss;
  
  oss << "#\n# Scene-modifying commands";
  
  oss << "\n/vis/viewer/set/culling global ";
  if (fCulling) {
    oss << "true";
  } else {
    oss << "false";
  }

  oss << "\n/vis/viewer/set/culling invisible ";
  if (fCullInvisible) {
    oss << "true";
  } else {
    oss << "false";
  }
  
  oss << "\n/vis/viewer/set/culling density ";
  if (fDensityCulling) {
    oss << "true " << fVisibleDensity/(g/cm3) << " g/cm3";
  } else {
    oss << "false";
  }
  
  oss << "\n/vis/viewer/set/culling coveredDaughters ";
  if (fCullCovered) {
    oss << "true";
  } else {
    oss << "false";
  }

  oss << "\n/vis/viewer/colourByDensity "
  << fCBDAlgorithmNumber << " g/cm3";
  for (auto p: fCBDParameters) {
    oss << ' ' << p/(g/cm3);
  }

  oss << "\n/vis/viewer/set/sectionPlane ";
  if (fSection) {
    oss << "on "
    << G4BestUnit(fSectionPlane.point(),"Length")
    << fSectionPlane.normal().x()
    << ' ' << fSectionPlane.normal().y()
    << ' ' << fSectionPlane.normal().z();
  } else {
    oss << "off";
  }
  
  oss << "\n/vis/viewer/set/cutawayMode ";
  if (fCutawayMode == cutawayUnion) {
    oss << "union";
  } else {
    oss << "intersection";
  }
  
  oss << "\n/vis/viewer/clearCutawayPlanes";
  if (fCutawayPlanes.size()) {
    for (size_t i = 0; i < fCutawayPlanes.size(); i++) {
      oss << "\n/vis/viewer/addCutawayPlane "
      << G4BestUnit(fCutawayPlanes[i].point(),"Length")
      << fCutawayPlanes[i].normal().x()
      << ' ' << fCutawayPlanes[i].normal().y()
      << ' ' << fCutawayPlanes[i].normal().z();
    }
  } else {
    oss << "\n# No cutaway planes defined.";
  }
  
  oss << "\n/vis/viewer/set/explodeFactor "
  << fExplodeFactor
  << ' ' << G4BestUnit(fExplodeCentre,"Length");
  
  oss << "\n/vis/viewer/set/lineSegmentsPerCircle "
  << fNoOfSides;
  
  oss << std::endl;
  
  return oss.str();
}

G4String G4ViewParameters::TouchableCommands() const
{
  std::ostringstream oss;
  
  oss << "#\n# Touchable commands";

  const std::vector<G4ModelingParameters::VisAttributesModifier>& vams =
    fVisAttributesModifiers;

  if (vams.empty()) {
    oss
    << "\n# None"
    << "\n/vis/viewer/clearVisAttributesModifiers";
    oss << std::endl;
    return oss.str();
  }

  oss
  << "\n/vis/viewer/clearVisAttributesModifiers";

  G4ModelingParameters::PVNameCopyNoPath lastPath;
  std::vector<G4ModelingParameters::VisAttributesModifier>::const_iterator
    iModifier;
  for (iModifier = vams.begin();
       iModifier != vams.end();
       ++iModifier) {
    const G4ModelingParameters::PVNameCopyNoPath& vamPath =
      iModifier->GetPVNameCopyNoPath();
    if (vamPath != lastPath) {
      lastPath = vamPath;
      oss << "\n/vis/set/touchable";
      G4ModelingParameters::PVNameCopyNoPathConstIterator iVAM;
      for (iVAM = vamPath.begin();
           iVAM != vamPath.end();
           ++iVAM) {
        oss << ' ' << iVAM->GetName() << ' ' << iVAM->GetCopyNo();
      }
    }
    const G4VisAttributes& vamVisAtts = iModifier->GetVisAttributes();
    const G4Colour& c = vamVisAtts.GetColour();
    switch (iModifier->GetVisAttributesSignifier()) {
      case G4ModelingParameters::VASVisibility:
        oss << "\n/vis/touchable/set/visibility ";
        if (vamVisAtts.IsVisible()) {
          oss << "true";
        } else {
          oss << "false";
        }
        break;
      case G4ModelingParameters::VASDaughtersInvisible:
        oss << "\n/vis/touchable/set/daughtersInvisible ";
        if (vamVisAtts.IsDaughtersInvisible()) {
          oss << "true";
        } else {
          oss << "false";
        }
        break;
      case G4ModelingParameters::VASColour:
        oss << "\n/vis/touchable/set/colour "
        << c.GetRed()
        << ' ' << c.GetGreen()
        << ' ' << c.GetBlue()
        << ' ' << c.GetAlpha();
        break;
      case G4ModelingParameters::VASLineStyle:
        oss << "\n/vis/touchable/set/lineStyle ";
        switch (vamVisAtts.GetLineStyle()) {
          case G4VisAttributes::unbroken:
            oss << "unbroken";
            break;
          case G4VisAttributes::dashed:
            oss << "dashed";
            break;
          case G4VisAttributes::dotted:
          oss << "dotted";
        }
        break;
      case G4ModelingParameters::VASLineWidth:
        oss << "\n/vis/touchable/set/lineWidth "
        << vamVisAtts.GetLineWidth();
        break;
      case G4ModelingParameters::VASForceWireframe:
        if (vamVisAtts.IsForceDrawingStyle()) {
          if (vamVisAtts.GetForcedDrawingStyle() == G4VisAttributes::wireframe) {
            oss << "\n/vis/touchable/set/forceWireframe ";
            if (vamVisAtts.IsForceDrawingStyle()) {
              oss << "true";
            } else {
              oss << "false";
            }
          }
        }
        break;
      case G4ModelingParameters::VASForceSolid:
        if (vamVisAtts.IsForceDrawingStyle()) {
          if (vamVisAtts.GetForcedDrawingStyle() == G4VisAttributes::solid) {
            oss << "\n/vis/touchable/set/forceSolid ";
            if (vamVisAtts.IsForceDrawingStyle()) {
              oss << "true";
            } else {
              oss << "false";
            }
          }
        }
        break;
      case G4ModelingParameters::VASForceCloud:
        if (vamVisAtts.IsForceDrawingStyle()) {
          if (vamVisAtts.GetForcedDrawingStyle() == G4VisAttributes::cloud) {
            oss << "\n/vis/touchable/set/forceCloud ";
            if (vamVisAtts.IsForceDrawingStyle()) {
              oss << "true";
            } else {
              oss << "false";
            }
          }
        }
        break;
      case G4ModelingParameters::VASForceAuxEdgeVisible:
        if (vamVisAtts.IsForceAuxEdgeVisible()) {
          oss << "\n/vis/touchable/set/forceAuxEdgeVisible ";
          if (vamVisAtts.IsForcedAuxEdgeVisible()) {
            oss << "true";
          } else {
            oss << "false";
          }
        }
        break;
      case G4ModelingParameters::VASForceLineSegmentsPerCircle:
        oss << "\n/vis/touchable/set/lineSegmentsPerCircle "
        << vamVisAtts.GetForcedLineSegmentsPerCircle();
        break;
      case G4ModelingParameters::VASForceNumberOfCloudPoints:
        oss << "\n/vis/touchable/set/numberOfCloudPoints "
        << vamVisAtts.GetForcedNumberOfCloudPoints();
        break;
    }
  }
  
  oss << std::endl;
  
  return oss.str();
}

G4String G4ViewParameters::TimeWindowCommands() const
{
  std::ostringstream oss;

  oss <<  "#\n# Time window commands";

  oss
  << "\n/vis/viewer/set/timeWindow/startTime "
  << fStartTime/ns << " ns ";

  oss
  << "\n/vis/viewer/set/timeWindow/endTime "
  << fEndTime/ns << " ns ";

  oss << "\n/vis/viewer/set/timeWindow/fadeFactor "
  << fFadeFactor;

  oss
  << "\n/vis/viewer/set/timeWindow/displayHeadTime ";
  if (!fDisplayHeadTime) {
    oss << "false";
  } else {
    oss
    << "true"
    << ' ' << fDisplayHeadTimeX
    << ' ' << fDisplayHeadTimeY
    << ' ' << fDisplayHeadTimeSize
    << ' ' << fDisplayHeadTimeRed
    << ' ' << fDisplayHeadTimeGreen
    << ' ' << fDisplayHeadTimeBlue;
  }

  oss
  << "\n/vis/viewer/set/timeWindow/displayLightFront ";
  if (!fDisplayLightFront) {
    oss << "false";
  } else {
    oss
    << "true"
    << ' ' << fDisplayLightFrontX/mm
    << ' ' << fDisplayLightFrontY/mm
    << ' ' << fDisplayLightFrontZ/mm
    << " mm"
    << ' ' << fDisplayLightFrontT/ns
    << " ns"
    << ' ' << fDisplayLightFrontRed
    << ' ' << fDisplayLightFrontGreen
    << ' ' << fDisplayLightFrontBlue;
  }

  oss << std::endl;

  return oss.str();
}

void G4ViewParameters::PrintDifferences (const G4ViewParameters& v) const {

  // Put performance-sensitive parameters first.
  if (
      // This first to optimise spin, etc.
      (fViewpointDirection   != v.fViewpointDirection)   ||

      // No particular order from here on.
      (fDrawingStyle         != v.fDrawingStyle)         ||
      (fNumberOfCloudPoints  != v.fNumberOfCloudPoints)  ||
      (fAuxEdgeVisible       != v.fAuxEdgeVisible)       ||
      (fCulling              != v.fCulling)              ||
      (fCullInvisible        != v.fCullInvisible)        ||
      (fDensityCulling       != v.fDensityCulling)       ||
      (fVisibleDensity       != v.fVisibleDensity)       ||
      (fCullCovered          != v.fCullCovered)          ||
      (fCBDAlgorithmNumber   != v.fCBDAlgorithmNumber)   ||
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
      (fGeometryMask         != v.fGeometryMask)         ||
      (fAutoRefresh          != v.fAutoRefresh)          ||
      (fBackgroundColour     != v.fBackgroundColour)     || 
      (fPicking              != v.fPicking)              ||
      (fRotationStyle        != v.fRotationStyle)
      )
    G4cout << "Difference in 1st batch." << G4endl;

  if (fCBDAlgorithmNumber > 0) {
    if (fCBDParameters.size() != v.fCBDParameters.size()) {
      G4cout << "Difference in number of colour by density parameters." << G4endl;
    } else if (fCBDParameters != v.fCBDParameters) {
      G4cout << "Difference in values of colour by density parameters." << G4endl;
    }
  }

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

  if (fVisAttributesModifiers != v.fVisAttributesModifiers) {
    G4cout << "Difference in vis attributes modifiers." << G4endl;
  }

  if (fStartTime != v.fStartTime ||
      fEndTime   != v.fEndTime)  {
    G4cout << "Difference in time window." << G4endl;
  }

  if (fFadeFactor != v.fFadeFactor) {
    G4cout << "Difference in time window fade factor." << G4endl;
  }

  if (fDisplayHeadTime != v.fDisplayHeadTime) {
    G4cout << "Difference in display head time flag." << G4endl;
  } else {
    if (fDisplayHeadTimeX     != v.fDisplayHeadTimeX     ||
        fDisplayHeadTimeY     != v.fDisplayHeadTimeY     ||
        fDisplayHeadTimeSize  != v.fDisplayHeadTimeSize  ||
        fDisplayHeadTimeRed   != v.fDisplayHeadTimeRed   ||
        fDisplayHeadTimeGreen != v.fDisplayHeadTimeGreen ||
        fDisplayHeadTimeBlue  != v.fDisplayHeadTimeBlue) {
      G4cout << "Difference in display head time parameters." << G4endl;
    }
  }

  if (fDisplayLightFront != v.fDisplayLightFront) {
    G4cout << "Difference in display light front flag." << G4endl;
  } else {
    if (fDisplayLightFrontX     != v.fDisplayLightFrontX     ||
        fDisplayLightFrontY     != v.fDisplayLightFrontY     ||
        fDisplayLightFrontZ     != v.fDisplayLightFrontZ     ||
        fDisplayLightFrontT     != v.fDisplayLightFrontT     ||
        fDisplayLightFrontRed   != v.fDisplayLightFrontRed   ||
        fDisplayLightFrontGreen != v.fDisplayLightFrontGreen ||
        fDisplayLightFrontBlue  != v.fDisplayLightFrontBlue) {
      G4cout << "Difference in display light front parameters." << G4endl;
    }
  }
}

std::ostream& operator <<
 (std::ostream& os, G4ViewParameters::DrawingStyle style)
{
  switch (style) {
    case G4ViewParameters::wireframe:
      os << "wireframe"; break;
    case G4ViewParameters::hlr:
      os << "hlr - hidden lines removed"; break;
    case G4ViewParameters::hsr:
      os << "hsr - hidden surfaces removed"; break;
    case G4ViewParameters::hlhsr:
      os << "hlhsr - hidden line, hidden surface removed"; break;
    case G4ViewParameters::cloud:
      os << "cloud - draw volume as a cloud of dots"; break;
    default: os << "unrecognised"; break;
  }
  return os;
}

std::ostream& operator <<
(std::ostream& os, G4ViewParameters::SMROption option)
{
  switch (option) {
    case G4ViewParameters::meshAsDots:
      os << "dots"; break;
    case G4ViewParameters::meshAsSurfaces:
      os << "surfaces"; break;
  }
  return os;
}

std::ostream& operator << (std::ostream& os, const G4ViewParameters& v) {
  os << "View parameters and options:";

  os << "\n  Drawing style: " << v.fDrawingStyle;

  os << "\n  Number of cloud points: " << v.fNumberOfCloudPoints;

  os << "\n  Auxiliary edges: ";
  if (!v.fAuxEdgeVisible) os << "in";
  os << "visible";

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

  os << "\n  Colour by density: ";
  if (v.fCBDAlgorithmNumber <= 0) {
    os << "inactive";
  } else {
    os << "Algorithm " << v.fCBDAlgorithmNumber << ", Parameters:";
    for (auto p: v.fCBDParameters) {
      os << ' ' << G4BestUnit(p,"Volumic Mass");
    }
  }

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
  os << "\n  X geometry mask: "
     << std::showbase << std::hex << v.fGeometryMask
     << std::noshowbase << std::dec;

  os << "\n  Auto refresh: ";
  if (v.fAutoRefresh) os << "true";
  else os << "false";

  os << "\n  Background colour: " << v.fBackgroundColour;

  os << "\n  Picking requested: ";
  if (v.fPicking) os << "true";
  else os << "false";

  os << "\n  Rotation style: ";
  switch (v.fRotationStyle) {
  case G4ViewParameters::constrainUpDirection:
    os << "constrainUpDirection (conventional HEP view)"; break;
  case G4ViewParameters::freeRotation:
    os << "freeRotation (Google-like rotation, using mouse-grab)"; break;
  default: os << "unrecognised"; break;
  }

  os << "\n  Vis attributes modifiers: ";
  const std::vector<G4ModelingParameters::VisAttributesModifier>& vams =
    v.fVisAttributesModifiers;
  if (vams.empty()) {
    os << "None";
  } else {
    os << vams;
  }

  os << "\n  Time window parameters:"
  << "\n  Start time:  " << v.fStartTime/ns << " ns"
  << "\n  End time:    " << v.fEndTime/ns << " ns"
  << "\n  Fade factor: " << v.fFadeFactor;
  if (!v.fDisplayHeadTime) {
    os << "\n  Head time display not requested.";
  } else {
    os
    << "\n  Head time position: "
    << v.fDisplayHeadTimeX << ' ' << v.fDisplayHeadTimeY
    << "\n  Head time size:     " << v.fDisplayHeadTimeSize
    << "\n  Head time colour:   " << v.fDisplayHeadTimeRed
    << ' ' << v.fDisplayHeadTimeGreen << ' ' << v.fDisplayHeadTimeBlue;
  }
  if (!v.fDisplayLightFront) {
    os << "\n  Light front display not requested.";
  } else {
    os
    << "\n  Light front position: "
    << v.fDisplayLightFrontX/mm << ' ' << v.fDisplayLightFrontY/mm
    << ' ' << v.fDisplayLightFrontZ/mm << " mm"
    << "\n  Light front time:     " << v.fDisplayLightFrontT/ns << " ns"
    << "\n  Light front colour:   " << v.fDisplayLightFrontRed
    << ' ' << v.fDisplayLightFrontGreen << ' ' << v.fDisplayLightFrontBlue;
  }

  os << "\n  Special Mesh Rendering";
  if (v.fSpecialMeshRendering) {
    os << " requested with option \"" << v.fSpecialMeshRenderingOption;
    os << "\" for ";
    if (v.fSpecialMeshVolumes.empty()) {
      os << "any mesh";
    } else {
      os << "selected meshes";
      for (const auto& vol: v.fSpecialMeshVolumes) {
	os << "\n    " << vol.GetName() << ':' << vol.GetCopyNo();
      }
    }
  } else os << ": off";
  return os;
}

G4bool G4ViewParameters::operator != (const G4ViewParameters& v) const {

  // Put performance-sensitive parameters first.
  if (
      // This first to optimise spin, etc.
      (fViewpointDirection   != v.fViewpointDirection)   ||

      // No particular order from here on.
      (fDrawingStyle         != v.fDrawingStyle)         ||
      (fNumberOfCloudPoints  != v.fNumberOfCloudPoints)  ||
      (fAuxEdgeVisible       != v.fAuxEdgeVisible)       ||
      (fCulling              != v.fCulling)              ||
      (fCullInvisible        != v.fCullInvisible)        ||
      (fDensityCulling       != v.fDensityCulling)       ||
      (fCullCovered          != v.fCullCovered)          ||
      (fCBDAlgorithmNumber   != v.fCBDAlgorithmNumber)   ||
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
      (fGeometryMask         != v.fGeometryMask)         ||
      (fAutoRefresh          != v.fAutoRefresh)          ||
      (fBackgroundColour     != v.fBackgroundColour)     ||
      (fPicking              != v.fPicking)              ||
      (fRotationStyle        != v.fRotationStyle)        ||
      (fSpecialMeshRendering != v.fSpecialMeshRendering) ||
      (fSpecialMeshRenderingOption != v.fSpecialMeshRenderingOption)
      )
    return true;

  if (fDensityCulling &&
      (fVisibleDensity != v.fVisibleDensity)) return true;

  if (fCBDAlgorithmNumber > 0) {
    if (fCBDParameters.size() != v.fCBDParameters.size()) return true;
    else if (fCBDParameters != v.fCBDParameters) return true;
  }

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

  if (fVisAttributesModifiers != v.fVisAttributesModifiers) return true;

  if (fStartTime  != v.fStartTime ||
      fEndTime    != v.fEndTime   ||
      fFadeFactor != v.fFadeFactor) return true;

  if (fDisplayHeadTime != v.fDisplayHeadTime) return true;
  if (fDisplayHeadTime) {
    if (fDisplayHeadTimeX     != v.fDisplayHeadTimeX     ||
        fDisplayHeadTimeY     != v.fDisplayHeadTimeY     ||
        fDisplayHeadTimeSize  != v.fDisplayHeadTimeSize  ||
        fDisplayHeadTimeRed   != v.fDisplayHeadTimeRed   ||
        fDisplayHeadTimeGreen != v.fDisplayHeadTimeGreen ||
        fDisplayHeadTimeBlue  != v.fDisplayHeadTimeBlue) {
      return true;
    }
  }

  if (fDisplayLightFront != v.fDisplayLightFront) return true;
  if (fDisplayLightFront) {
    if (fDisplayLightFrontX     != v.fDisplayLightFrontX     ||
        fDisplayLightFrontY     != v.fDisplayLightFrontY     ||
        fDisplayLightFrontZ     != v.fDisplayLightFrontZ     ||
        fDisplayLightFrontT     != v.fDisplayLightFrontT     ||
        fDisplayLightFrontRed   != v.fDisplayLightFrontRed   ||
        fDisplayLightFrontGreen != v.fDisplayLightFrontGreen ||
        fDisplayLightFrontBlue  != v.fDisplayLightFrontBlue) {
      return true;
    }
  }

  if (fSpecialMeshRendering) {
    if (fSpecialMeshVolumes != v.fSpecialMeshVolumes)
      return true;;
  }

  return false;
}

void G4ViewParameters::SetXGeometryString (const G4String& geomString)
{
  const G4String delimiters("xX+-");
  G4String::size_type i = geomString.find_first_of(delimiters);
  if (i == G4String::npos) {
    // Does not contain "xX+-".
    // Is it a single number?
    std::istringstream iss(geomString);
    G4int size;
    iss >> size;
    if (iss) {
      // It is a number
      fWindowSizeHintX = size;
      fWindowSizeHintY = size;
    }
    // Accept other or all defaults (in G4ViewParameters constructor)
    // Reconstruct a geometry string coherent with the above
    char signX, signY;
    if (fWindowLocationHintXNegative) signX = '-'; else signX ='+';
    if (fWindowLocationHintYNegative) signY = '-'; else signY ='+';
    std::ostringstream oss;
    oss << fWindowSizeHintX << 'x' << fWindowSizeHintY
    << signX << fWindowLocationHintX << signY << fWindowLocationHintY;
    fXGeometryString = oss.str();
    return;
  }

  // Assume it's a parseable X geometry string
  G4int x = 0, y = 0;
  unsigned int w = 0, h = 0;
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
    
    G4warn << "Unrecognised geometry string \""
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
  char *strind;
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
    G4int Result = 0;
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

G4ViewParameters* G4ViewParameters::CatmullRomCubicSplineInterpolation
(const std::vector<G4ViewParameters>& views,
 G4int nInterpolationPoints)  // No of interpolations points per interval
{
  // Returns a null pointer when no more to be done.  For example:
  // do {
  //   G4ViewParameters* vp =
  //   G4ViewParameters::CatmullRomCubicSplineInterpolation(viewVector,nInterpolationPoints);
  //   if (!vp) break;
  //     ...
  // } while (true);

  // See https://en.wikipedia.org/wiki/Cubic_Hermite_spline

  // Assumes equal intervals

  if (views.size() < 2) {
    G4Exception
    ("G4ViewParameters::CatmullRomCubicSplineInterpolation",
     "visman0301", JustWarning,
     "There must be at least two views.");
    return 0;
  }

  if (nInterpolationPoints < 1) {
    G4Exception
    ("G4ViewParameters::CatmullRomCubicSplineInterpolation",
     "visman0302", JustWarning,
     "Number of interpolation points cannot be zero or negative.");
    return 0;
  }

  const size_t nIntervals = views.size() - 1;
  const G4double dt = 1./nInterpolationPoints;

  static G4ViewParameters holdingValues;
  static G4double t = 0.;  // 0. <= t <= 1.
  static G4int iInterpolationPoint = 0;
  static size_t iInterval = 0;

//  G4cout << "Interval " << iInterval << ", t = " << t << G4endl;

  // Hermite polynomials.
  const G4double h00 = 2.*t*t*t - 3.*t*t +1;
  const G4double h10 = t*t*t -2.*t*t + t;
  const G4double h01 = -2.*t*t*t + 3.*t*t;
  const G4double h11 = t*t*t - t*t;

  // Aliases (to simplify code)
  const size_t& n = nIntervals;
  size_t& i = iInterval;
  const std::vector<G4ViewParameters>& v = views;

  // The Catmull-Rom cubic spline prescription is as follows:
  // Slope at first way point is v[1] - v[0].
  // Slope at last way point is v[n] - v[n-1].
  // Otherwise slope at way point i is 0.5*(v[i+1] - v[i-1]).
  // Result = h00*v[i] + h10*m[i] + h01*v[i+1] + h11*m[i+1],
  // where m[i] amd m[i+1] are the slopes at the start and end
  // of the interval for the particular value.
  // If (n == 1), linear interpolation results.
  // If (n == 2), quadratic interpolation results.

  // Working variables
  G4double mi, mi1, real, x, y, z;
  
  // First, a crude interpolation of all parameters.  Then, below, a
  // smooth interpolation of those for which it makes sense.
  holdingValues = t < 0.5? v[i]: v[i+1];

  // Catmull-Rom cubic spline interpolation
#define INTERPOLATE(param) \
/* This works out the interpolated param in i'th interval */ \
/* Assumes n >= 1 */ \
if (i == 0) { \
/* First interval */ \
mi = v[1].param - v[0].param; \
/* If there is only one interval, make start and end slopes equal */ \
/* (This results in a linear interpolation) */ \
if (n == 1) mi1 = mi; \
/* else the end slope of the interval takes account of the next waypoint along */ \
else mi1 = 0.5 * (v[2].param - v[0].param); \
} else if (i >= n - 1) { \
/* Similarly for last interval */ \
mi1 = v[i+1].param - v[i].param; \
/* If there is only one interval, make start and end slopes equal */ \
if (n == 1) mi = mi1; \
/* else the start slope of the interval takes account of the previous waypoint */ \
else mi = 0.5 * (v[i+1].param - v[i-1].param); \
} else { \
/* Full Catmull-Rom slopes use previous AND next waypoints */ \
mi  = 0.5 * (v[i+1].param - v[i-1].param); \
mi1 = 0.5 * (v[i+2].param - v[i  ].param); \
} \
real = h00 * v[i].param + h10 * mi + h01 * v[i+1].param + h11 * mi1;

#define INTERPOLATELOG(param) \
if (i == 0) { \
mi = std::log(v[1].param) - std::log(v[0].param); \
if (n == 1) mi1 = mi; \
else mi1 = 0.5 * (std::log(v[2].param) - std::log(v[0].param)); \
} else if (i >= n - 1) { \
mi1 = std::log(v[i+1].param) - std::log(v[i].param); \
if (n == 1) mi = mi1; \
else mi = 0.5 * (std::log(v[i+1].param) - std::log(v[i-1].param)); \
} else { \
mi  = 0.5 * (std::log(v[i+1].param) - std::log(v[i-1].param)); \
mi1 = 0.5 * (std::log(v[i+2].param) - std::log(v[i  ].param)); \
} \
real = std::exp(h00 * std::log(v[i].param) + h10 * mi + h01 * std::log(v[i+1].param) + h11 * mi1);

  // Real parameters
  INTERPOLATE(fVisibleDensity);
  if (real < 0.) real = 0.;
  holdingValues.fVisibleDensity = real;
  INTERPOLATELOG(fExplodeFactor);
  holdingValues.fExplodeFactor = real;
  INTERPOLATE(fFieldHalfAngle);
  if (real < 0.) real = 0.;
  holdingValues.fFieldHalfAngle = real;
  INTERPOLATELOG(fZoomFactor);
  holdingValues.fZoomFactor = real;
  INTERPOLATE(fDolly);
  holdingValues.fDolly = real;
  INTERPOLATE(fGlobalMarkerScale);
  if (real < 0.) real = 0.;
  holdingValues.fGlobalMarkerScale = real;
  INTERPOLATE(fGlobalLineWidthScale);
  if (real < 0.) real = 0.;
  holdingValues.fGlobalLineWidthScale = real;

  // Unit vectors
#define INTERPOLATEUNITVECTOR(vector) \
INTERPOLATE(vector.x()); x = real; \
INTERPOLATE(vector.y()); y = real; \
INTERPOLATE(vector.z()); z = real;
  INTERPOLATEUNITVECTOR(fViewpointDirection);
  holdingValues.fViewpointDirection          = G4Vector3D(x,y,z).unit();
  INTERPOLATEUNITVECTOR(fUpVector);
  holdingValues.fUpVector                    = G4Vector3D(x,y,z).unit();
  INTERPOLATEUNITVECTOR(fRelativeLightpointDirection);
  holdingValues.fRelativeLightpointDirection = G4Vector3D(x,y,z).unit();
  INTERPOLATEUNITVECTOR(fActualLightpointDirection);
  holdingValues.fActualLightpointDirection   = G4Vector3D(x,y,z).unit();

  // Un-normalised vectors
#define INTERPOLATEVECTOR(vector) \
INTERPOLATE(vector.x()); x = real; \
INTERPOLATE(vector.y()); y = real; \
INTERPOLATE(vector.z()); z = real;
  INTERPOLATEVECTOR(fScaleFactor);
  holdingValues.fScaleFactor = G4Vector3D(x,y,z);

  // Points
#define INTERPOLATEPOINT(point) \
INTERPOLATE(point.x()); x = real; \
INTERPOLATE(point.y()); y = real; \
INTERPOLATE(point.z()); z = real;
  INTERPOLATEPOINT(fExplodeCentre);
  holdingValues.fExplodeCentre      = G4Point3D(x,y,z);
  INTERPOLATEPOINT(fCurrentTargetPoint);
  holdingValues.fCurrentTargetPoint = G4Point3D(x,y,z);

  // Colour
  G4double red, green, blue, alpha;
#define INTERPOLATECOLOUR(colour) \
INTERPOLATE(colour.GetRed());   red   = real; \
INTERPOLATE(colour.GetGreen()); green = real; \
INTERPOLATE(colour.GetBlue());  blue  = real; \
INTERPOLATE(colour.GetAlpha()); alpha = real;
  INTERPOLATECOLOUR(fBackgroundColour);
  // Components are clamped to 0. <= component <= 1.
  holdingValues.fBackgroundColour = G4Colour(red,green,blue,alpha);

  // For some parameters we need to check some continuity
  G4bool continuous;
#define CONTINUITY(quantity) \
  continuous = false; \
  /* This follows the logic of the INTERPOLATE macro above; see comments therein */ \
  if (i == 0) { \
    if (v[1].quantity == v[0].quantity) { \
       if (n == 1) continuous = true; \
       else if (v[2].quantity == v[0].quantity) \
       continuous = true; \
    } \
  } else if (i >= n - 1) { \
    if (v[i+1].quantity == v[i].quantity) { \
      if (n == 1) continuous = true; \
      else if (v[i+1].quantity == v[i-1].quantity) \
      continuous = true; \
    } \
  } else { \
    if (v[i-1].quantity == v[i].quantity && \
        v[i+1].quantity == v[i].quantity && \
        v[i+2].quantity == v[i].quantity) \
    continuous = true; \
  }

  G4double a, b, c, d;
#define INTERPOLATEPLANE(plane) \
INTERPOLATE(plane.a()); a = real; \
INTERPOLATE(plane.b()); b = real; \
INTERPOLATE(plane.c()); c = real; \
INTERPOLATE(plane.d()); d = real;

  // Section plane
  CONTINUITY(fSection);
  if (continuous) {
    INTERPOLATEPLANE(fSectionPlane);
    holdingValues.fSectionPlane = G4Plane3D(a,b,c,d);
  }

  // Cutaway planes
  if (v[i].fCutawayPlanes.size()) {
    CONTINUITY(fCutawayPlanes.size());
    if (continuous) {
      for (size_t j = 0; j < v[i].fCutawayPlanes.size(); ++j) {
        INTERPOLATEPLANE(fCutawayPlanes[j]);
        holdingValues.fCutawayPlanes[j] = G4Plane3D(a,b,c,d);
      }
    }
  }

  // Vis attributes modifiers
  // Really, we are only interested in colour - other attributes can follow
  // the "crude" interpolation that is guaranteed above.
  static G4VisAttributes workingVA;
  if  (v[i].fVisAttributesModifiers.size()) {
    CONTINUITY(fVisAttributesModifiers.size());
    if (continuous) {
      for (size_t j = 0; j < v[i].fVisAttributesModifiers.size(); ++j) {
        CONTINUITY(fVisAttributesModifiers[j].GetPVNameCopyNoPath());
        if (continuous) {
          CONTINUITY(fVisAttributesModifiers[j].GetVisAttributesSignifier());
          if (continuous) {
            if (v[i].fVisAttributesModifiers[j].GetVisAttributesSignifier() ==
                G4ModelingParameters::VASColour) {
              INTERPOLATECOLOUR(fVisAttributesModifiers[j].GetVisAttributes().GetColour());
              workingVA = v[i].fVisAttributesModifiers[j].GetVisAttributes();
              workingVA.SetColour(G4Colour(red,green,blue,alpha));
              holdingValues.fVisAttributesModifiers[j].SetVisAttributes(workingVA);
            }
          }
        }
      }
    }
  }

  // Time window parameters (for showing particles in flight)
  // Only two parameters are interpolated. The others are usually chosen
  // once and for all by the user for a given series of views - or at least,
  // if not, they will be interpolated by the default "crude" method above.
  INTERPOLATE(fStartTime)
  holdingValues.fStartTime = real;
  INTERPOLATE(fEndTime)
  holdingValues.fEndTime = real;

  // Increment counters
  iInterpolationPoint++;
  t += dt;
  if (iInterpolationPoint > nInterpolationPoints) {
    iInterpolationPoint = 1;  // Ready for next interval.
    t = dt;
    iInterval++;
  }
  if (iInterval >= nIntervals) {
    iInterpolationPoint = 0;  // Ready for a complete restart.
    t = 0.;
    iInterval = 0;
    return 0;
  }

  return &holdingValues;
}

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
// Michael Kelsey  31 January 2019
//
// Class Description:
//
// Abstract base class to implement drawing vector field geometries
// (e.g., electric, magnetic or gravity).  Implementation extracted
// from G4MagneticFieldModel, with field-value access left pure
// virtual for implementation by base classes.

#include "G4VFieldModel.hh"

#include "G4ArrowModel.hh"
#include "G4Colour.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Point3D.hh"
#include "G4Polyline.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include <sstream>
#include <limits>
#include <vector>

#define G4warn G4cout

// Constructor and destructor

G4VFieldModel::~G4VFieldModel() {;}

G4VFieldModel::G4VFieldModel
(const G4String& typeOfField, const G4String& symbol,
 const G4VisExtent& extentForField,
 const std::vector<G4PhysicalVolumesSearchScene::Findings>& pvFindings,
 G4int nDataPointsPerMaxHalfExtent,
 Representation representation,
 G4int arrow3DLineSegmentsPerCircle)
: fExtentForField(extentForField)
, fPVFindings(pvFindings)
, fNDataPointsPerMaxHalfExtent(nDataPointsPerMaxHalfExtent)
, fRepresentation(representation)
, fArrow3DLineSegmentsPerCircle(arrow3DLineSegmentsPerCircle)
, fTypeOfField(typeOfField)
, fArrowPrefix(symbol)
{
  fType = "G4"+typeOfField+"FieldModel";
  fGlobalTag = fType;

  std::ostringstream oss;
  oss << ':' << fNDataPointsPerMaxHalfExtent
  << ':' << fArrow3DLineSegmentsPerCircle;
  if (fExtentForField == G4VisExtent::GetNullExtent()) {
    oss << " whole scene";
  } else {
    oss
    << ':' << fExtentForField.GetXmin()
    << ':' << fExtentForField.GetXmax()
    << ':' << fExtentForField.GetYmin()
    << ':' << fExtentForField.GetYmax()
    << ':' << fExtentForField.GetZmin()
    << ':' << fExtentForField.GetZmax();
  }
  for (const auto& findings: fPVFindings) {
    oss
    << ',' << findings.fpFoundPV->GetName()
    << ':' << findings.fFoundPVCopyNo;
  }
  if (fRepresentation == Representation::fullArrow) {
    oss << " full arrow";
  } else if (fRepresentation == Representation::lightArrow) {
    oss << " light arrow";
  }

  fGlobalDescription = fType + oss.str();
}


// The main task of a model is to describe itself to the graphics scene.

void G4VFieldModel::DescribeYourselfTo(G4VGraphicsScene& sceneHandler) {
//  G4cout << "G4VFieldModel::DescribeYourselfTo" << G4endl;

  G4TransportationManager* tMgr =
  G4TransportationManager::GetTransportationManager();
  assert(tMgr);
  G4Navigator* navigator = tMgr->GetNavigatorForTracking();
  assert(navigator);

  G4FieldManager* globalFieldMgr = tMgr->GetFieldManager();
  const G4Field* globalField = 0;
  const G4String intro = "G4VFieldModel::DescribeYourselfTo: ";
  if (globalFieldMgr) {
    if (globalFieldMgr->DoesFieldExist()) {
      globalField = globalFieldMgr->GetDetectorField();
      if (!globalField) {
        static G4bool warned = false;
        if (!warned) {
          G4warn << intro << "Null global field pointer." << G4endl;
          warned = true;
        }
      }
    }
  } else {
    static G4bool warned = false;
    if (!warned) {
      G4warn << intro << "No global field manager." << G4endl;
      warned = true;
    }
  }

  G4VisExtent extent = sceneHandler.GetExtent();
  if (fExtentForField == G4VisExtent::GetNullExtent()) {
    extent = sceneHandler.GetExtent();
  } else {
    extent = fExtentForField;
  }
  const G4double& xMin = extent.GetXmin();
  const G4double& yMin = extent.GetYmin();
  const G4double& zMin = extent.GetZmin();
  const G4double& xMax = extent.GetXmax();
  const G4double& yMax = extent.GetYmax();
  const G4double& zMax = extent.GetZmax();
  const G4double xHalfScene = 0.5 * (xMax - xMin);
  const G4double yHalfScene = 0.5 * (yMax - yMin);
  const G4double zHalfScene = 0.5 * (zMax - zMin);
  const G4double xSceneCentre = 0.5 * (xMax + xMin);
  const G4double ySceneCentre = 0.5 * (yMax + yMin);
  const G4double zSceneCentre = 0.5 * (zMax + zMin);
  const G4double maxHalfScene =
  std::max(xHalfScene,std::max(yHalfScene,zHalfScene));
  if (maxHalfScene <= 0.) {
    G4warn << "Scene extent non-positive." << G4endl;
    return;
  }

  // Constants
  const G4double interval = maxHalfScene / fNDataPointsPerMaxHalfExtent;
  const G4int nDataPointsPerXHalfScene = G4int(xHalfScene / interval);
  const G4int nDataPointsPerYHalfScene = G4int(yHalfScene / interval);
  const G4int nDataPointsPerZHalfScene = G4int(zHalfScene / interval);
  const G4int nXSamples = 2 * nDataPointsPerXHalfScene + 1;
  const G4int nYSamples = 2 * nDataPointsPerYHalfScene + 1;
  const G4int nZSamples = 2 * nDataPointsPerZHalfScene + 1;
  const G4int nSamples = nXSamples * nYSamples * nZSamples;
  const G4double arrowLengthMax = 0.8 * interval;

  // Working vectors for field values, etc.
  std::vector<G4Point3D> Field(nSamples);          // Initialises to (0,0,0)
  std::vector<G4Point3D> xyz(nSamples);            // Initialises to (0,0,0)
  G4double FieldMagnitudeMax = -std::numeric_limits<G4double>::max();

  // Get field values and ascertain maximum field.
  for (G4int i = 0; i < nXSamples; i++) {
    G4double x = xSceneCentre + (i - nDataPointsPerXHalfScene) * interval;

    for (G4int j = 0; j < nYSamples; j++) {
      G4double y = ySceneCentre + (j - nDataPointsPerYHalfScene) * interval;

      for (G4int k = 0; k < nZSamples; k++) {
        G4double z = zSceneCentre + (k - nDataPointsPerZHalfScene) * interval;

        // Calculate indices into working vectors
        const G4int ijk = i * nYSamples * nZSamples + j * nZSamples + k;
	xyz[ijk].set(x,y,z);

        G4ThreeVector pos(x,y,z);

        // Check if point is in findings
        if (!fPVFindings.empty()) {
          G4bool isInPV = false;
          for (const auto& findings: fPVFindings) {
            G4VPhysicalVolume* pv = findings.fpFoundPV;
            G4int copyNo = findings.fFoundPVCopyNo;
            G4VSolid* solid = pv->GetLogicalVolume()->GetSolid();
            G4PVParameterised* pvParam = dynamic_cast<G4PVParameterised*>(pv);
            if (pvParam) {
              auto* param = pvParam->GetParameterisation();
              solid = param->ComputeSolid(copyNo,pvParam);
              solid->ComputeDimensions(param,copyNo,pvParam);
            }
            // Transform point to local coordinate system
            const auto& transform = findings.fFoundObjectTransformation;
            auto rotation = transform.getRotation();
            auto translation = transform.getTranslation();
            G4ThreeVector lPos = pos; lPos -= translation; lPos.transform(rotation.invert());
            if (solid->Inside(lPos)==kInside) {
              isInPV = true;
              break;
            }
          }
          if (!isInPV) continue;
        }
        // Point is in findings - or there were no findings

        // Find volume and field at this location.
        const G4VPhysicalVolume* pPV =
        navigator->LocateGlobalPointAndSetup(pos,0,false,true);
        const G4Field* field = globalField;
        if (pPV) {
          // Get logical volume.
          const G4LogicalVolume* pLV = pPV->GetLogicalVolume();
          if (pLV) {
            // Value for Region, if any, overrides
            G4Region* pRegion = pLV->GetRegion();
            if (pRegion) {
              G4FieldManager* pRegionFieldMgr = pRegion->GetFieldManager();
              if (pRegionFieldMgr) {
                field = pRegionFieldMgr->GetDetectorField();
                // G4cout << "Region with field" << G4endl;
              }
            }
            // 'Local' value from logical volume, if any, overrides
            G4FieldManager* pLVFieldMgr = pLV->GetFieldManager();
            if (pLVFieldMgr) {
              field = pLVFieldMgr->GetDetectorField();
              // G4cout << "Logical volume with field" << G4endl;
            }
          }
        }

	G4double time = 0.;	// FIXME:  Can we get event time in some way?

	// Subclasses will have implemented this for their own field
	GetFieldAtLocation(field, xyz[ijk], time, Field[ijk]);

	G4double mag = Field[ijk].mag();
	if (mag > FieldMagnitudeMax) FieldMagnitudeMax = mag;
      }	// for (k, z
    }	// for (j, y
  }	// for (i, x

  if (FieldMagnitudeMax <= 0.) {
    G4warn << "No " << fTypeOfField << " field in this extent." << G4endl;
    return;
  }

  for (G4int i = 0; i < nSamples; i++) {
    const G4double Fmag = Field[i].mag();
    const G4double f = Fmag / FieldMagnitudeMax;
    if (f <= 0.) continue;  // Skip zero field locations

    G4double red = 0., green = 0., blue = 0., alpha = 1.;
    if (f < 0.5) {  // Linear colour scale: 0->0.5->1 is red->green->blue.
      green = 2. * f;
      red = 2. * (0.5 - f);
    } else {
      blue = 2. * (f - 0.5);
      green = 2. * (1.0 - f);
    }
    const G4Colour arrowColour(red,green,blue,alpha);

    // Very small arrows are difficult to see. Better to draw a line.
    G4bool drawAsLine = false;
    switch (fRepresentation) {
      case Representation::fullArrow:
        if (f < 0.1) {
          drawAsLine = true;
        }
        break;
      case Representation::lightArrow:
        drawAsLine = true;
        break;
      default:
        break;
    }

    // Head of arrow depends on field direction and strength...
    G4double arrowLength = arrowLengthMax * f;
    // ...but limit the length so it's visible.
    if (f < 0.01) arrowLength = arrowLengthMax * 0.01;
    const G4Point3D head = xyz[i] + arrowLength*Field[i]/Fmag;

    if (drawAsLine) {
      G4Polyline FArrowLite;
      G4VisAttributes va(arrowColour);
      va.SetLineWidth(2.);
      FArrowLite.SetVisAttributes(va);
      FArrowLite.push_back(xyz[i]);
      FArrowLite.push_back(head);
      sceneHandler.BeginPrimitives();
      sceneHandler.AddPrimitive(FArrowLite);
      sceneHandler.EndPrimitives();
    } else {
      G4ArrowModel FArrow(xyz[i].x(), xyz[i].y(), xyz[i].z(),
                          head.x(), head.y(), head.z(),
                          arrowLength/5, arrowColour,
                          fArrowPrefix+"Field",
                          fArrow3DLineSegmentsPerCircle);
      FArrow.DescribeYourselfTo(sceneHandler);
    }
  }	// for (i, nSamples
}

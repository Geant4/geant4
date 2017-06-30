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
// $Id: G4MagneticFieldModel.cc 101958 2016-12-12 08:04:35Z gcosmo $
//
// 
// John Allison  17th August 2013
// Model that knows how to draw the magnetic field.

#include "G4MagneticFieldModel.hh"

#include "G4VGraphicsScene.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4Colour.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ArrowModel.hh"
#include "G4Polyline.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"

#include <sstream>
#include <limits>
#include <vector>

G4MagneticFieldModel::~G4MagneticFieldModel () {}

G4MagneticFieldModel::G4MagneticFieldModel
(G4int nDataPointsPerMaxHalfScene,
 Representation representation,
 G4int arrow3DLineSegmentsPerCircle)
: fNDataPointsPerMaxHalfScene(nDataPointsPerMaxHalfScene)
, fRepresentation(representation)
, fArrow3DLineSegmentsPerCircle(arrow3DLineSegmentsPerCircle)
{
  fType = "G4MagneticFieldModel";
  fGlobalTag = fType;
  std::ostringstream oss;
  oss << ':' << fNDataPointsPerMaxHalfScene
  << ':' << fArrow3DLineSegmentsPerCircle;
  fGlobalDescription = fType + oss.str();
}

void G4MagneticFieldModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
//  G4cout << "G4MagneticFieldModel::DescribeYourselfTo" << G4endl;

  const G4VisExtent& extent = sceneHandler.GetExtent();
  const G4double xMin = extent.GetXmin();
  const G4double yMin = extent.GetYmin();
  const G4double zMin = extent.GetZmin();
  const G4double xMax = extent.GetXmax();
  const G4double yMax = extent.GetYmax();
  const G4double zMax = extent.GetZmax();
  const G4double xHalfScene = 0.5 * (xMax - xMin);
  const G4double yHalfScene = 0.5 * (yMax - yMin);
  const G4double zHalfScene = 0.5 * (zMax - zMin);
  const G4double xSceneCentre = 0.5 * (xMax + xMin);
  const G4double ySceneCentre = 0.5 * (yMax + yMin);
  const G4double zSceneCentre = 0.5 * (zMax + zMin);
  const G4double maxHalfScene =
  std::max(xHalfScene,std::max(yHalfScene,zHalfScene));
  if (maxHalfScene <= 0.) {
    G4cout
    << "Extent non-positive."
    << G4endl;
    return;
  }

  G4TransportationManager* tMgr =
  G4TransportationManager::GetTransportationManager();
  assert(tMgr);
  G4Navigator* navigator = tMgr->GetNavigatorForTracking();
  assert(navigator);

  G4FieldManager* globalFieldMgr = tMgr->GetFieldManager();
  const G4Field* globalField = 0;
  const G4String intro = "G4MagneticFieldModel::DescribeYourselfTo: ";
  if (globalFieldMgr) {
    if (globalFieldMgr->DoesFieldExist()) {
      globalField = globalFieldMgr->GetDetectorField();
      if (!globalField) {
        static G4bool warned = false;
        if (!warned) {
          G4cout  << intro
          << "Null global field pointer."
          << G4endl;
          warned = true;
        }
      }
    }
  } else {
    static G4bool warned = false;
    if (!warned) {
      G4cout  << intro
      << "No global field manager."
      << G4endl;
      warned = true;
    }
  }

  // Constants
  const G4double interval = maxHalfScene / fNDataPointsPerMaxHalfScene;
  const G4int nDataPointsPerXHalfScene = G4int(xHalfScene / interval);
  const G4int nDataPointsPerYHalfScene = G4int(yHalfScene / interval);
  const G4int nDataPointsPerZHalfScene = G4int(zHalfScene / interval);
  const G4int nXSamples = 2 * nDataPointsPerXHalfScene + 1;
  const G4int nYSamples = 2 * nDataPointsPerYHalfScene + 1;
  const G4int nZSamples = 2 * nDataPointsPerZHalfScene + 1;
  const G4int nSamples = nXSamples * nYSamples * nZSamples;
  const G4int nSamples3 = nSamples * 3;
  const G4double arrowLengthMax = 0.8 * interval;
  const G4int nResults = 6;  // 3 B-field + 3 E-field.

  // Working space for GetFieldValue.
  double position_time[4] = {0,0,0,0};
  double result[nResults];

  // Working vectors for field values, etc.
  std::vector<G4double> BField(nSamples3);          // Initialises to zero.
  std::vector<G4double> BFieldMagnitude(nSamples);  // Initialises to zero.
  std::vector<G4double> xyz(nSamples3);             // Initialises to zero.

  // Get field values and ascertain maximum field.
  G4double BFieldMagnitudeMax = -std::numeric_limits<G4double>::max();
  for (G4int i = 0; i < nXSamples; i++) {
    G4double x = xSceneCentre + (i - nDataPointsPerXHalfScene) * interval;
    position_time[0] = x;
    for (G4int j = 0; j < nYSamples; j++) {
      G4double y = ySceneCentre + (j - nDataPointsPerYHalfScene) * interval;
      position_time[1] = y;
      for (G4int k = 0; k < nZSamples; k++) {
        G4double z = zSceneCentre + (k - nDataPointsPerZHalfScene) * interval;
        position_time[2] = z;
        // Calculate indices into working vectors
        const G4int ijk = i * nYSamples * nZSamples + j * nZSamples + k;
        const G4int ijk3 = ijk * 3;
        // Find volume at this location.
        G4ThreeVector pos(x,y,z);
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
        // If field found, get values and store in working vectors.
        if (field) {
          // Get field values in result array.
          field->GetFieldValue(position_time,result);
          //                G4cout
          //                << "BField/T:"
          //                << " " << result[0]/tesla
          //                << " " << result[1]/tesla
          //                << " " << result[2]/tesla
          //                << G4endl;
          // Store B-field components.
          for (G4int l = 0; l < 3; l++) {
            BField[ijk3 + l] = result[l];
          }
          // Calculate magnitude and store.
          G4double mag = sqrt
          (result[0]*result[0]+result[1]*result[1]+result[2]*result[2]);
          BFieldMagnitude[ijk] = mag;
          // Store position.
          xyz[ijk3] = x;
          xyz[ijk3 + 1] = y;
          xyz[ijk3 + 2] = z;
          // Find maximum field magnitude.
          if (mag > BFieldMagnitudeMax) {
            BFieldMagnitudeMax = mag;
          }
        }
      }
    }
  }

  if (BFieldMagnitudeMax <= 0) {
    G4cout
    << "No field in this scene."
    << G4endl;
    return;
  }

  if (fRepresentation == Representation::lightArrow) sceneHandler.BeginPrimitives();
  for (G4int i = 0; i < nSamples; i++) {
    if (BFieldMagnitude[i] > 0) {
      const G4int i3 = i * 3;
      // Field (Bx,By,Bz) at (x,y,z).
      const G4double Bx = BField[i3];
      const G4double By = BField[i3 + 1];
      const G4double Bz = BField[i3 + 2];
      const G4double x = xyz[i3];
      const G4double y = xyz[i3 + 1];
      const G4double z = xyz[i3 + 2];
      const G4double B = BFieldMagnitude[i];
//      G4cout
//      << "Position/mm, BField/T unpacked:"
//      << ' ' << x/mm
//      << ' ' << y/mm
//      << ' ' << z/mm
//      << " " << Bx/tesla
//      << " " << By/tesla
//      << " " << Bz/tesla
//      << G4endl;
      if (B > 0.) {
        const G4double f = B / BFieldMagnitudeMax;
        G4double red = 0., green = 0., blue = 0., alpha = 1.;
        if (f < 0.5) {  // Linear colour scale: 0->0.5->1 is blue->green->red.
          green = 2. * f;
          blue = 2. * (0.5 - f);
        } else {
          red = 2. * (f - 0.5);
          green = 2. * (1.0 - f);
        }
        const G4Colour arrowColour(red,green,blue,alpha);
        const G4double arrowLength = arrowLengthMax * f;
        // Base of arrow is at (x,y,z).
        const G4double& x1 = x;
        const G4double& y1 = y;
        const G4double& z1 = z;
        // Head of arrow depends on field direction and strength.
        const G4double x2 = x1 + arrowLength * Bx / B;
        const G4double y2 = y1 + arrowLength * By / B;
        const G4double z2 = z1 + arrowLength * Bz / B;
        if (fRepresentation == Representation::fullArrow) {
          G4ArrowModel BArrow(x1,y1,z1,x2,y2,z2,arrowLength/5,arrowColour,
                              "BField",
                              fArrow3DLineSegmentsPerCircle);
          BArrow.DescribeYourselfTo(sceneHandler);
        } else if (fRepresentation == Representation::lightArrow) {
          G4Polyline BArrowLite;
          G4VisAttributes va(arrowColour);
          BArrowLite.SetVisAttributes(va);
          BArrowLite.push_back(G4Point3D(x1,y1,z1));
          BArrowLite.push_back(G4Point3D(x2,y2,z2));
          sceneHandler.AddPrimitive(BArrowLite);
        }
      }
    }
  }
  if (fRepresentation == Representation::lightArrow) sceneHandler.EndPrimitives();
}

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
// John Allison  26th July 1999.
// Model for logical volumes.

#include "G4LogicalVolumeModel.hh"

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4DrawVoxels.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VReadOutGeometry.hh"
#include "G4Circle.hh"

#include <vector>
#include <utility>

G4LogicalVolumeModel::G4LogicalVolumeModel
(G4LogicalVolume*            pLV,
 G4int                       soughtDepth,
 G4bool                      booleans,
 G4bool                      voxels,
 G4bool                      readout,
 G4bool                      checkOverlaps,
 const G4Transform3D&        modelTransformation,
 const G4ModelingParameters* pMP):
  // Instantiate a G4PhysicalVolumeModel with a G4PVPlacement to
  // represent this logical volume.  It has no rotation and a null
  // translation so that the logical volume will be seen in its own
  // reference system.  It will be added to the physical volume store
  // but it will not be part of the normal geometry heirarchy so it
  // has no mother.
  G4PhysicalVolumeModel
(new G4PVPlacement (0,                   // No rotation.
		    G4ThreeVector(),     // Null traslation.
		    "PhysVol representation of LogVol " + pLV -> GetName (),
		    pLV,
		    0,                   // No mother.
		    false,               // Not "MANY".
		    0),                  // Copy number.
 soughtDepth,
 modelTransformation,
 pMP,
 true),                                  // Use full extent.
  fpLV (pLV),
  fBooleans (booleans),
  fVoxels (voxels),
  fReadout (readout),
  fCheckOverlaps(checkOverlaps),
  fOverlapsPrinted(false)
{
  fType = "G4LogicalVolumeModel";
  fGlobalTag = fpLV -> GetName ();
  fGlobalDescription = "G4LogicalVolumeModel " + fGlobalTag;
}

G4LogicalVolumeModel::~G4LogicalVolumeModel () {}

namespace {
  // Vis attributes
  const G4Colour highlightSolidColour(1.0,0.8,0.8);
  const G4double highlightSolidLineWidth(10./*pixels*/);
  const G4Colour highlightPointColour(0.5,0.5,1.0);
  const G4double highlightPointDiameter(20./*pixels*/);
  // Keep a vector of solid-copy number pairs to avoid duplication.
  typedef std::pair<G4VSolid*,G4int> solidCopyNoPair;
  std::vector<solidCopyNoPair> solidCopyNoVector;
  void DrawSolid
  (G4VGraphicsScene& sceneHandler,
   G4VSolid* sol, G4int copyNo, const G4Transform3D& t) {
    // Avoid duplication.
    std::pair<G4VSolid*,G4int> pair(sol,copyNo);
    auto iter = solidCopyNoVector.begin();
    for ( ; iter != solidCopyNoVector.end(); ++iter) {
      if (*iter == pair) break;
    }
    if (iter == solidCopyNoVector.end()) {
      solidCopyNoVector.push_back(pair);
      G4VisAttributes highlightSolidVisAtts(highlightSolidColour);
      highlightSolidVisAtts.SetLineWidth(highlightSolidLineWidth);
      sceneHandler.PreAddSolid(t,highlightSolidVisAtts);
      sceneHandler.AddSolid(*sol);
      sceneHandler.PostAddSolid();
    }
  }
  void DrawPoint
  (G4VGraphicsScene& sceneHandler,
   const G4ThreeVector& point) {
    G4VisAttributes highlightPointVisAtts(highlightPointColour);
    G4Circle overlapPoint;
    overlapPoint.SetVisAttributes(highlightPointVisAtts);
    overlapPoint.SetPosition(point);
    overlapPoint.SetDiameter(G4VMarker::SizeType::screen,highlightPointDiameter);
    overlapPoint.SetFillStyle(G4VMarker::FillStyle::filled);
    sceneHandler.BeginPrimitives();
    sceneHandler.AddPrimitive(overlapPoint);
    sceneHandler.EndPrimitives();
  }
}

void G4LogicalVolumeModel::DescribeYourselfTo
(G4VGraphicsScene& sceneHandler) {

  // Store current modeling parameters and ensure nothing is culled.
  const G4ModelingParameters* tmpMP = fpMP;
  G4ModelingParameters nonCulledMP;
  if (fpMP) nonCulledMP = *fpMP;
  nonCulledMP.SetCulling (false);
  fpMP = &nonCulledMP;    
  G4PhysicalVolumeModel::DescribeYourselfTo (sceneHandler);
  fpMP = tmpMP;

  if (fVoxels) {
    if (fpTopPV->GetLogicalVolume()->GetVoxelHeader()) {
      // Add Voxels.
      G4DrawVoxels dv;
      G4PlacedPolyhedronList* pPPL =
	dv.CreatePlacedPolyhedra (fpTopPV -> GetLogicalVolume ());
      for (size_t i = 0; i < pPPL -> size (); i++) {
	const G4Transform3D& transform = (*pPPL)[i].GetTransform ();
	const G4Polyhedron& polyhedron = (*pPPL)[i].GetPolyhedron ();
	sceneHandler.BeginPrimitives (transform);
	sceneHandler.AddPrimitive (polyhedron);
	sceneHandler.EndPrimitives ();
      }
      delete pPPL;
    }
  }

  if (fReadout) {
    // Draw readout geometry...
    G4VSensitiveDetector* sd = fpLV->GetSensitiveDetector();
    if (sd) {
      G4VReadOutGeometry* roGeom = sd->GetROgeometry();
      if (roGeom) {
	G4VPhysicalVolume* roWorld = roGeom->GetROWorld();
//	G4cout << "Readout geometry \"" << roGeom->GetName()
//	       << "\" with top physical volume \""
//	       << roWorld->GetName()
//	       << "\"" << G4endl;
	G4PhysicalVolumeModel pvModel(roWorld);
	pvModel.SetModelingParameters(fpMP);
	pvModel.DescribeYourselfTo(sceneHandler);
      }
    }
  }

  if (fCheckOverlaps) {
    G4LogicalVolume* motherLog = fpTopPV->GetLogicalVolume();
    G4VSolid* motherSolid = motherLog->GetSolid();
    G4int nDaughters = (G4int)motherLog->GetNoDaughters();

    // Models are called repeatedly by the scene handler so be careful...
    // Print overlaps - but only the first time for a given instantiation of G4LogicalVolume
    if (!fOverlapsPrinted) {
      for (G4int iDaughter = 0; iDaughter < nDaughters; ++iDaughter) {
        G4VPhysicalVolume* daughterPhys = motherLog->GetDaughter(iDaughter);
        daughterPhys->CheckOverlaps();
      }
      fOverlapsPrinted = true;
    }

    // Draw overlaps
    solidCopyNoVector.clear();
    for (G4int iDaughter = 0; iDaughter < nDaughters; ++iDaughter) {
      G4VPhysicalVolume* daughterPhys = motherLog->GetDaughter(iDaughter);
      G4PVPlacement* daughterPVPlace = dynamic_cast<G4PVPlacement*>(daughterPhys);
      G4PVParameterised* daughterPVParam = dynamic_cast<G4PVParameterised*>(daughterPhys);
      const G4int nPoints = 1000;

      if (daughterPVPlace) {

        // This algorithm is based on G4PVPlacement::CheckOverlaps.
        G4AffineTransform tDaughter(daughterPhys->GetRotation(),daughterPhys->GetTranslation());
        G4VSolid* daughterSolid = daughterPhys->GetLogicalVolume()->GetSolid();
        for (G4int i = 0; i < nPoints; ++i) {
          G4ThreeVector point = daughterSolid->GetPointOnSurface();
          // Transform to mother's coordinate system
          G4ThreeVector motherPoint = tDaughter.TransformPoint(point);
          // Check overlaps with the mother volume
          if (motherSolid->Inside(motherPoint)==kOutside) {
            // Draw mother and daughter and point
            DrawSolid(sceneHandler,motherSolid,0,G4Transform3D());
            DrawSolid(sceneHandler,daughterSolid,daughterPhys->GetCopyNo(),tDaughter);
            DrawPoint(sceneHandler,motherPoint);
          }
          // Check other daughters
          for (G4int iSister = 0; iSister < nDaughters; ++iSister) {
            if (iSister == iDaughter) continue;
            G4VPhysicalVolume* sisterPhys = motherLog->GetDaughter(iSister);
            G4AffineTransform tSister(sisterPhys->GetRotation(),sisterPhys->GetTranslation());
            // Transform to sister's coordinate system
            G4ThreeVector sisterPoint = tSister.InverseTransformPoint(motherPoint);
            G4LogicalVolume* sisterLog = sisterPhys->GetLogicalVolume();
            G4VSolid* sisterSolid = sisterLog->GetSolid();
            if (sisterSolid->Inside(sisterPoint)==kInside) {
              // Draw daughter and sister and point
              DrawSolid(sceneHandler,daughterSolid,daughterPhys->GetCopyNo(),tDaughter);
              DrawSolid(sceneHandler,sisterSolid,sisterPhys->GetCopyNo(),tSister);
              DrawPoint(sceneHandler,motherPoint);
            }
          }
        }

      } else if (daughterPVParam) {

        // This algorithm is based on G4PVParameterised::CheckOverlaps
        const G4int multiplicity = daughterPVParam->GetMultiplicity();
        auto* param = daughterPVParam->GetParameterisation();
        // Cache points for later checking against other parameterisations
        std::vector<G4ThreeVector> motherPoints;
        for (G4int iP = 0; iP < multiplicity; iP++) {
          G4VSolid* daughterSolid = param->ComputeSolid(iP, daughterPhys);
          daughterSolid->ComputeDimensions(param, iP, daughterPhys);
          param->ComputeTransformation(iP, daughterPhys);
          G4AffineTransform tDaughter(daughterPVParam->GetRotation(),daughterPVParam->GetTranslation());
          for (G4int i = 0; i < nPoints; ++i) {
            G4ThreeVector point = daughterSolid->GetPointOnSurface();
            // Transform to mother's coordinate system
            G4ThreeVector motherPoint = tDaughter.TransformPoint(point);
            // Check overlaps with the mother volume
            if (motherSolid->Inside(motherPoint)==kOutside) {
              // Draw mother and daughter and point
              DrawSolid(sceneHandler,motherSolid,0,G4Transform3D());
              DrawSolid(sceneHandler,daughterSolid,iP,tDaughter);
              DrawPoint(sceneHandler,motherPoint);
            }
            motherPoints.push_back(motherPoint);
          }
          // Check sister parameterisations
          for (G4int iPP = iP + 1; iPP < multiplicity; iPP++) {
            G4VSolid* sisterSolid = param->ComputeSolid(iPP, daughterPhys);
            sisterSolid->ComputeDimensions(param, iPP, daughterPhys);
            param->ComputeTransformation(iPP, daughterPhys);
            G4AffineTransform tSister
            (daughterPVParam->GetRotation(),daughterPVParam->GetTranslation());
            for (const auto& motherPoint: motherPoints) {
              // Transform each point into daughter's frame
              G4ThreeVector sisterPoint = tSister.InverseTransformPoint(motherPoint);
              if (sisterSolid->Inside(sisterPoint)==kInside) {
                // Draw sister
                DrawSolid(sceneHandler,sisterSolid,iPP,tSister);
                // Recompute daughter parameterisation before drawing
                daughterSolid->ComputeDimensions(param, iP, daughterPhys);
                param->ComputeTransformation(iP, daughterPhys);
                tDaughter = G4AffineTransform
                (daughterPVParam->GetRotation(),daughterPVParam->GetTranslation());
                DrawSolid(sceneHandler,daughterSolid,iP,tDaughter);
                DrawPoint(sceneHandler,motherPoint);
              }
            }
          }
        }
      }
    }
  }
}

// This called from G4PhysicalVolumeModel::DescribeAndDescend by the
// virtual function mechanism.
void G4LogicalVolumeModel::DescribeSolid
(const G4Transform3D& theAT,
 G4VSolid* pSol,
 const G4VisAttributes* pVisAttribs,
 G4VGraphicsScene& sceneHandler) {

  if (fBooleans) {
    // Look for "constituents".  Could be a Boolean solid.
    G4VSolid* pSol0 = pSol -> GetConstituentSolid (0);
    if (pSol0) {  // Composite solid...
      G4VSolid* pSol1 = pSol -> GetConstituentSolid (1);
      if (!pSol1) {
	G4Exception
	  ("G4PhysicalVolumeModel::DescribeSolid",
	   "modeling0001", FatalException,
	   "2nd component solid in Boolean is missing.");
      }
      // Draw these constituents white and "forced wireframe"...
      G4VisAttributes constituentAttributes;
      constituentAttributes.SetForceWireframe(true);
      DescribeSolid (theAT, pSol0, &constituentAttributes, sceneHandler);
      DescribeSolid (theAT, pSol1, &constituentAttributes, sceneHandler);
    }
  }

  // In any case draw the original/resultant solid...
  sceneHandler.PreAddSolid (theAT, *pVisAttribs);
  pSol -> DescribeYourselfTo (sceneHandler);
  sceneHandler.PostAddSolid ();
}

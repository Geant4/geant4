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
// John Allison  31st December 1997.
// Model for physical volumes.

#include "G4PhysicalVolumeModel.hh"

#include "G4VGraphicsScene.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPVParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4BoundingExtentScene.hh"
#include "G4TransportationManager.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4Vector3D.hh"
#include "G4Mesh.hh"

#include <sstream>
#include <iomanip>

#define G4warn G4cout

namespace {
  G4int volumeCount = 0;
}

G4PhysicalVolumeModel::G4PhysicalVolumeModel
(G4VPhysicalVolume*            pVPV
 , G4int                       requestedDepth
 , const G4Transform3D&        modelTransform
 , const G4ModelingParameters* pMP
 , G4bool                      useFullExtent
 , const std::vector<G4PhysicalVolumeNodeID>& baseFullPVPath)
: G4VModel           (pMP)
, fpTopPV            (pVPV)
, fTopPVCopyNo       (pVPV? pVPV->GetCopyNo(): 0)
, fRequestedDepth    (requestedDepth)
, fUseFullExtent     (useFullExtent)
, fTransform         (modelTransform)
, fCurrentDepth      (0)
, fpCurrentPV        (fpTopPV)
, fCurrentPVCopyNo   (fpTopPV? fpTopPV->GetCopyNo(): 0)
, fpCurrentLV        (fpTopPV? fpTopPV->GetLogicalVolume(): 0)
, fpCurrentMaterial  (fpCurrentLV? fpCurrentLV->GetMaterial(): 0)
, fCurrentTransform  (modelTransform)
, fBaseFullPVPath    (baseFullPVPath)
, fAbort             (false)
, fCurtailDescent    (false)
, fpClippingSolid    (0)
, fClippingMode      (subtraction)
{
  fType = "G4PhysicalVolumeModel";

  if (!fpTopPV) {

    // In some circumstances creating an "empty" G4PhysicalVolumeModel is
    // allowed, so I have supressed the G4Exception below.  If it proves to
    // be a problem we might have to re-instate it, but it is unlikley to
    // be used except by visualisation experts.  See, for example, /vis/list,
    // where it is used simply to get a list of G4AttDefs.
    //    G4Exception
    //    ("G4PhysicalVolumeModel::G4PhysicalVolumeModel",
    //     "modeling0010", FatalException, "Null G4PhysicalVolumeModel pointer.");

    fTopPVName = "NULL";
    fGlobalTag = "Empty";
    fGlobalDescription = "G4PhysicalVolumeModel " + fGlobalTag;

  } else {

    fTopPVName = fpTopPV -> GetName ();
    std::ostringstream oss;
    oss << fpTopPV->GetName() << ':' << fpTopPV->GetCopyNo()
    << " BasePath:" << fBaseFullPVPath;
    fGlobalTag = oss.str();
    fGlobalDescription = "G4PhysicalVolumeModel " + fGlobalTag;
    CalculateExtent ();
  }
}

G4PhysicalVolumeModel::~G4PhysicalVolumeModel ()
{
  delete fpClippingSolid;
}

G4ModelingParameters::PVNameCopyNoPath G4PhysicalVolumeModel::GetPVNameCopyNoPath
(const std::vector<G4PhysicalVolumeNodeID>& path)
{
  G4ModelingParameters::PVNameCopyNoPath PVNameCopyNoPath;
  for (const auto& node: path) {
    PVNameCopyNoPath.push_back
    (G4ModelingParameters::PVNameCopyNo
     (node.GetPhysicalVolume()->GetName(),node.GetCopyNo()));
  }
  return PVNameCopyNoPath;
}

void G4PhysicalVolumeModel::CalculateExtent ()
{
  // To handle paramaterisations, set copy number and compute dimensions
  // to get extent right
  G4VPVParameterisation* pP = fpTopPV -> GetParameterisation ();
  if (pP) {
    fpTopPV -> SetCopyNo (fTopPVCopyNo);
    G4VSolid* solid = pP -> ComputeSolid (fTopPVCopyNo, fpTopPV);
    solid -> ComputeDimensions (pP, fTopPVCopyNo, fpTopPV);
  }
  if (fUseFullExtent) {
    fExtent = fpTopPV -> GetLogicalVolume () -> GetSolid () -> GetExtent ();
  } else {
    // Calculate extent of *drawn* volumes, i.e., ignoring culled, e.g.,
    // invisible volumes, by traversing the whole geometry hierarchy below
    // this physical volume.
    G4BoundingExtentScene beScene(this);
    const G4int tempRequestedDepth = fRequestedDepth;
    const G4Transform3D tempTransform = fTransform;
    const G4ModelingParameters* tempMP = fpMP;
    fRequestedDepth = -1;  // Always search to all depths to define extent.
    fTransform = G4Transform3D();  // Extent is in local cooridinates
    G4ModelingParameters mParams
      (0,      // No default vis attributes needed.
       G4ModelingParameters::wf,  // wireframe (not relevant for this).
       true,   // Global culling.
       true,   // Cull invisible volumes.
       false,  // Density culling.
       0.,     // Density (not relevant if density culling false).
       true,   // Cull daughters of opaque mothers.
       24);    // No of sides (not relevant for this operation).
    mParams.SetSpecialMeshRendering(true);  // Avoids traversing parameterisations
    fpMP = &mParams;
    DescribeYourselfTo (beScene);
    fExtent = beScene.GetBoundingExtent();
    fpMP = tempMP;
    fTransform = tempTransform;
    fRequestedDepth = tempRequestedDepth;
  }
  G4double radius = fExtent.GetExtentRadius();
  if (radius < 0.) {  // Nothing in the scene - revert to top extent
    fExtent = fpTopPV -> GetLogicalVolume () -> GetSolid () -> GetExtent ();
  }
  fExtent.Transform(fTransform);
}

void G4PhysicalVolumeModel::DescribeYourselfTo
(G4VGraphicsScene& sceneHandler)
{
  if (!fpTopPV) G4Exception
    ("G4PhysicalVolumeModel::DescribeYourselfTo",
     "modeling0012", FatalException, "No model.");

  if (!fpMP) G4Exception
    ("G4PhysicalVolumeModel::DescribeYourselfTo",
     "modeling0003", FatalException, "No modeling parameters.");

  G4Transform3D startingTransformation = fTransform;

  volumeCount = 0;

  VisitGeometryAndGetVisReps
    (fpTopPV,
     fRequestedDepth,
     startingTransformation,
     sceneHandler);

//  G4cout
//  << "G4PhysicalVolumeModel::DescribeYourselfTo: volume count: "
//  << volumeCount
//  << G4endl;

  // Reset or clear data...
  fCurrentDepth     = 0;
  fpCurrentPV       = fpTopPV;
  fCurrentPVCopyNo  = fpTopPV->GetCopyNo();
  fpCurrentLV       = fpTopPV->GetLogicalVolume();
  fpCurrentMaterial = fpCurrentLV? fpCurrentLV->GetMaterial(): 0;
  fFullPVPath       = fBaseFullPVPath;
  fDrawnPVPath.clear();
  fAbort            = false;
  fCurtailDescent   = false;
}

G4String G4PhysicalVolumeModel::GetCurrentTag () const
{
  if (fpCurrentPV) {
    std::ostringstream o;
    o << fpCurrentPV -> GetCopyNo ();
    return fpCurrentPV -> GetName () + ":" + o.str();
  }
  else {
    return "WARNING: NO CURRENT VOLUME - global tag is " + fGlobalTag;
  }
}

G4String G4PhysicalVolumeModel::GetCurrentDescription () const
{
  return "G4PhysicalVolumeModel " + GetCurrentTag ();
}

void G4PhysicalVolumeModel::VisitGeometryAndGetVisReps
(G4VPhysicalVolume* pVPV,
 G4int requestedDepth,
 const G4Transform3D& theAT,
 G4VGraphicsScene& sceneHandler)
{
  // Visits geometry structure to a given depth (requestedDepth), starting
  //   at given physical volume with given starting transformation and
  //   describes volumes to the scene handler.
  // requestedDepth < 0 (default) implies full visit.
  // theAT is the Accumulated Transformation.

  // Find corresponding logical volume and (later) solid, storing in
  // local variables to preserve re-entrancy.
  G4LogicalVolume* pLV  = pVPV -> GetLogicalVolume ();
  G4VSolid* pSol = nullptr;
  G4Material* pMaterial = nullptr;

  if (!(pVPV -> IsReplicated ())) {
    // Non-replicated physical volume.
    pSol = pLV -> GetSolid ();
    pMaterial = pLV -> GetMaterial ();
    DescribeAndDescend (pVPV, requestedDepth, pLV, pSol, pMaterial,
			theAT, sceneHandler);
  }
  else {
    // Replicated or parametrised physical volume.
    EAxis axis;
    G4int nReplicas;
    G4double width;
    G4double offset;
    G4bool consuming;
    pVPV -> GetReplicationData (axis, nReplicas, width,  offset, consuming);
    G4int nBegin = 0;
    G4int nEnd = nReplicas;
    if (fCurrentDepth == 0) { // i.e., top volume
      nBegin = fTopPVCopyNo;  // Describe only one volume, namely the one
      nEnd = nBegin + 1;      // specified by the given copy number.
    }
    G4VPVParameterisation* pP = pVPV -> GetParameterisation ();
    if (pP) {  // Parametrised volume.
      for (int n = nBegin; n < nEnd; n++) {
	pSol = pP -> ComputeSolid (n, pVPV);
	pP -> ComputeTransformation (n, pVPV);
	pSol -> ComputeDimensions (pP, n, pVPV);
	pVPV -> SetCopyNo (n);
        fCurrentPVCopyNo = n;
	// Create a touchable of current parent for ComputeMaterial.
	// fFullPVPath has not been updated yet so at this point it
	// corresponds to the parent.
	G4PhysicalVolumeModelTouchable parentTouchable(fFullPVPath);
	pMaterial = pP -> ComputeMaterial (n, pVPV, &parentTouchable);
	DescribeAndDescend (pVPV, requestedDepth, pLV, pSol, pMaterial,
			    theAT, sceneHandler);
      }
    }
    else {  // Plain replicated volume.  From geometry_guide.txt...
      // The replica's positions are claculated by means of a linear formula.
      // Replication may occur along:
      // 
      // o Cartesian axes (kXAxis,kYAxis,kZAxis)
      // 
      //   The replications, of specified width have coordinates of
      //   form (-width*(nReplicas-1)*0.5+n*width,0,0) where n=0.. nReplicas-1
      //   for the case of kXAxis, and are unrotated.
      // 
      // o Radial axis (cylindrical polar) (kRho)
      // 
      //   The replications are cons/tubs sections, centred on the origin
      //   and are unrotated.
      //   They have radii of width*n+offset to width*(n+1)+offset
      //                      where n=0..nReplicas-1
      // 
      // o Phi axis (cylindrical polar) (kPhi)
      //   The replications are `phi sections' or wedges, and of cons/tubs form
      //   They have phi of offset+n*width to offset+(n+1)*width where
      //   n=0..nReplicas-1
      // 
      pSol = pLV -> GetSolid ();
      pMaterial = pLV -> GetMaterial ();
      G4ThreeVector originalTranslation = pVPV -> GetTranslation ();
      G4RotationMatrix* pOriginalRotation = pVPV -> GetRotation ();
      G4double originalRMin = 0., originalRMax = 0.;
      if (axis == kRho && pSol->GetEntityType() == "G4Tubs") {
	originalRMin = ((G4Tubs*)pSol)->GetInnerRadius();
	originalRMax = ((G4Tubs*)pSol)->GetOuterRadius();
      }
      G4bool visualisable = true;
      for (int n = nBegin; n < nEnd; n++) {
	G4ThreeVector translation;  // Identity.
	G4RotationMatrix rotation;  // Identity - life enough for visualizing.
	G4RotationMatrix* pRotation = 0;
	switch (axis) {
	default:
	case kXAxis:
	  translation = G4ThreeVector (-width*(nReplicas-1)*0.5+n*width,0,0);
	  break;
	case kYAxis:
	  translation = G4ThreeVector (0,-width*(nReplicas-1)*0.5+n*width,0);
	  break;
	case kZAxis:
	  translation = G4ThreeVector (0,0,-width*(nReplicas-1)*0.5+n*width);
	  break;
	case kRho:
	  if (pSol->GetEntityType() == "G4Tubs") {
	    ((G4Tubs*)pSol)->SetInnerRadius(width*n+offset);
	    ((G4Tubs*)pSol)->SetOuterRadius(width*(n+1)+offset);
	  } else {
	    if (fpMP->IsWarning())
	      G4warn <<
		"G4PhysicalVolumeModel::VisitGeometryAndGetVisReps: WARNING:"
		"\n  built-in replicated volumes replicated in radius for "
		     << pSol->GetEntityType() <<
		"-type\n  solids (your solid \""
		     << pSol->GetName() <<
		"\") are not visualisable."
		     << G4endl;
	    visualisable = false;
	  }
	  break;
	case kPhi:
	  rotation.rotateZ (-(offset+(n+0.5)*width));
	  // Minus Sign because for the physical volume we need the
	  // coordinate system rotation.
	  pRotation = &rotation;
	  break;
	} 
	pVPV -> SetTranslation (translation);
	pVPV -> SetRotation    (pRotation);
	pVPV -> SetCopyNo (n);
        fCurrentPVCopyNo = n;
	if (visualisable) {
	  DescribeAndDescend (pVPV, requestedDepth, pLV, pSol, pMaterial,
			    theAT, sceneHandler);
	}
      }
      // Restore originals...
      pVPV -> SetTranslation (originalTranslation);
      pVPV -> SetRotation    (pOriginalRotation);
      if (axis == kRho && pSol->GetEntityType() == "G4Tubs") {
	((G4Tubs*)pSol)->SetInnerRadius(originalRMin);
	((G4Tubs*)pSol)->SetOuterRadius(originalRMax);
      }
    }
  }
}

void G4PhysicalVolumeModel::DescribeAndDescend
(G4VPhysicalVolume* pVPV,
 G4int requestedDepth,
 G4LogicalVolume* pLV,
 G4VSolid* pSol,
 G4Material* pMaterial,
 const G4Transform3D& theAT,
 G4VGraphicsScene& sceneHandler)
{
  // Maintain useful data members...
  fpCurrentPV = pVPV;
  fCurrentPVCopyNo = pVPV->GetCopyNo();
  fpCurrentLV = pLV;
  fpCurrentMaterial = pMaterial;

  // Create a nodeID for use below - note the "drawn" flag is true
  G4int copyNo = fpCurrentPV->GetCopyNo();
  auto nodeID = G4PhysicalVolumeNodeID
  (fpCurrentPV,copyNo,fCurrentDepth,fCurrentTransform);

  // Update full path of physical volumes...
  fFullPVPath.push_back(nodeID);

  const G4RotationMatrix objectRotation = pVPV -> GetObjectRotationValue ();
  const G4ThreeVector&  translation     = pVPV -> GetTranslation ();
  G4Transform3D theLT (G4Transform3D (objectRotation, translation));

  // Compute the accumulated transformation...
  // Note that top volume's transformation relative to the world
  // coordinate system is specified in theAT == startingTransformation
  // = fTransform (see DescribeYourselfTo), so first time through the
  // volume's own transformation, which is only relative to its
  // mother, i.e., not relative to the world coordinate system, should
  // not be accumulated.
  G4Transform3D theNewAT (theAT);
  if (fCurrentDepth != 0) theNewAT = theAT * theLT;
  fCurrentTransform = theNewAT;

  const G4VisAttributes* pVisAttribs = pLV->GetVisAttributes();
  //  If the volume does not have any vis attributes, create it.
  G4VisAttributes* tempVisAtts = nullptr;
  if (!pVisAttribs) {
    if (fpMP->GetDefaultVisAttributes()) {
      tempVisAtts = new G4VisAttributes(*fpMP->GetDefaultVisAttributes());
    } else {
      tempVisAtts = new G4VisAttributes;
    }
    // The user may request /vis/viewer/set/colourByDensity.
    if (fpMP->GetCBDAlgorithmNumber() == 1) {
      // Algorithm 1: 3 parameters: Simple rainbow mapping.
      if (fpMP->GetCBDParameters().size() != 3) {
        G4Exception("G4PhysicalVolumeModelTouchable::DescribeAndDescend",
                    "modeling0014",
                    FatalErrorInArgument,
                    "Algorithm-parameter mismatch for Colour By Density");
      } else {
        const G4double d = pMaterial? pMaterial->GetDensity(): 0.;
        const G4double d0 = fpMP->GetCBDParameters()[0]; // Invisible d < d0.
        const G4double d1 = fpMP->GetCBDParameters()[1]; // Rainbow d0->d1->d2.
        const G4double d2 = fpMP->GetCBDParameters()[2]; // Blue d > d2.
        if (d < d0) { // Density < d0 is invisible.
          tempVisAtts->SetVisibility(false);
        } else { // Intermediate densities are on a spectrum.
          G4double red, green, blue;
          if (d < d1) {
            red = (d1-d)/(d1-d0); green = (d-d0)/(d1-d0); blue = 0.;
          } else if (d < d2) {
            red = 0.; green = (d2-d)/(d2-d1); blue = (d-d1)/(d2-d1);
          } else {  // Density >= d2 is blue.
            red = 0.; green = 0.; blue = 1.;
          }
          tempVisAtts->SetColour(G4Colour(red,green,blue));
        }
      }
    } else if (fpMP->GetCBDAlgorithmNumber() == 2) {
      // Algorithm 2
      // ...etc.
    }
    pVisAttribs = tempVisAtts;
  }
  // From here, can assume pVisAttribs is a valid pointer.  This is necessary
  // because PreAddSolid needs a vis attributes object.

  // Check if vis attributes are to be modified by a /vis/touchable/set/ command.
  const auto& vams = fpMP->GetVisAttributesModifiers();
  if (vams.size()) {
    // OK, we have some VAMs (Vis Attributes Modifiers).
    for (const auto& vam: vams) {
      const auto& vamPath = vam.GetPVNameCopyNoPath();
      if (vamPath.size() == fFullPVPath.size()) {
        // OK, we have a size match.
        // Check the volume name/copy number path.
        auto iVAMNameCopyNo = vamPath.begin();
        auto iPVNodeId = fFullPVPath.begin();
        for (; iVAMNameCopyNo != vamPath.end(); ++iVAMNameCopyNo, ++iPVNodeId) {
          if (!(
                iVAMNameCopyNo->GetName() ==
                iPVNodeId->GetPhysicalVolume()->GetName() &&
                iVAMNameCopyNo->GetCopyNo() ==
                iPVNodeId->GetPhysicalVolume()->GetCopyNo()
                )) {
            // This path element does NOT match.
            break;
          }
        }
        if (iVAMNameCopyNo == vamPath.end()) {
          // OK, the paths match (the above loop terminated normally).
          // Create a vis atts object for the modified vis atts.
          // It is static so that we may return a reliable pointer to it.
          static G4VisAttributes modifiedVisAtts;
          // Initialise it with the current vis atts and reset the pointer.
          modifiedVisAtts = *pVisAttribs;
          pVisAttribs = &modifiedVisAtts;
          const G4VisAttributes& transVisAtts = vam.GetVisAttributes();
          switch (vam.GetVisAttributesSignifier()) {
            case G4ModelingParameters::VASVisibility:
              modifiedVisAtts.SetVisibility(transVisAtts.IsVisible());
              break;
            case G4ModelingParameters::VASDaughtersInvisible:
              modifiedVisAtts.SetDaughtersInvisible
              (transVisAtts.IsDaughtersInvisible());
              break;
            case G4ModelingParameters::VASColour:
              modifiedVisAtts.SetColour(transVisAtts.GetColour());
              break;
            case G4ModelingParameters::VASLineStyle:
              modifiedVisAtts.SetLineStyle(transVisAtts.GetLineStyle());
              break;
            case G4ModelingParameters::VASLineWidth:
              modifiedVisAtts.SetLineWidth(transVisAtts.GetLineWidth());
              break;
            case G4ModelingParameters::VASForceWireframe:
              if (transVisAtts.IsForceDrawingStyle()) {
                if (transVisAtts.GetForcedDrawingStyle() ==
                    G4VisAttributes::wireframe) {
                  modifiedVisAtts.SetForceWireframe(true);
                }
              }
              break;
            case G4ModelingParameters::VASForceSolid:
              if (transVisAtts.IsForceDrawingStyle()) {
                if (transVisAtts.GetForcedDrawingStyle() ==
                    G4VisAttributes::solid) {
                  modifiedVisAtts.SetForceSolid(true);
                }
              }
              break;
            case G4ModelingParameters::VASForceCloud:
              if (transVisAtts.IsForceDrawingStyle()) {
                if (transVisAtts.GetForcedDrawingStyle() ==
                    G4VisAttributes::cloud) {
                  modifiedVisAtts.SetForceCloud(true);
                }
              }
              break;
            case G4ModelingParameters::VASForceNumberOfCloudPoints:
              modifiedVisAtts.SetForceNumberOfCloudPoints
              (transVisAtts.GetForcedNumberOfCloudPoints());
              break;
            case G4ModelingParameters::VASForceAuxEdgeVisible:
              if (transVisAtts.IsForceAuxEdgeVisible()) {
                modifiedVisAtts.SetForceAuxEdgeVisible
                (transVisAtts.IsForcedAuxEdgeVisible());
              }
              break;
            case G4ModelingParameters::VASForceLineSegmentsPerCircle:
              modifiedVisAtts.SetForceLineSegmentsPerCircle
              (transVisAtts.GetForcedLineSegmentsPerCircle());
              break;
          }
        }
      }
    }
  }

  // Check for special mesh rendering
  if (fpMP->IsSpecialMeshRendering()) {
    G4bool potentialG4Mesh = false;
    if (fpMP->GetSpecialMeshVolumes().empty()) {
      // No volumes specified - all are potentially possible
      potentialG4Mesh = true;
    } else {
      // Name and (optionally) copy number of container volume is specified
      for (const auto& pvNameCopyNo: fpMP->GetSpecialMeshVolumes()) {
        if (pVPV->GetName() == pvNameCopyNo.GetName()) {
          // We have a name match
          if (pvNameCopyNo.GetCopyNo() < 0) {
            // Any copy number is OK
            potentialG4Mesh = true;
          } else {
            if (pVPV->GetCopyNo() == pvNameCopyNo.GetCopyNo()) {
              // We have a name and copy number match
              potentialG4Mesh = true;
            }
          }
        }
      }
    }
    if (potentialG4Mesh) {
      // Create - or at least attempt to create - a mesh. If it cannot be created
      // out of this pVPV the type will be "invalid".
      G4Mesh mesh(pVPV,theNewAT);
      if (mesh.GetMeshType() != G4Mesh::invalid) {
        // Create "artificial" nodeID to represent the replaced volumes
        G4int artCopyNo = 0;
        auto artPV = mesh.GetParameterisedVolume();
        auto artDepth = fCurrentDepth + 1;
        auto artNodeID = G4PhysicalVolumeNodeID(artPV,artCopyNo,artDepth);
        fFullPVPath.push_back(artNodeID);
        fDrawnPVPath.push_back(artNodeID);
        sceneHandler.AddCompound(mesh);
        fFullPVPath.pop_back();
        fDrawnPVPath.pop_back();
        delete tempVisAtts;  // Needs cleaning up (Coverity warning!!)
        return;  // Mesh found and processed - nothing more to do.
      }  // else continue processing
    }
  }

  // Make decision to draw...
  G4bool thisToBeDrawn = true;

  // There are various reasons why this volume
  // might not be drawn...
  G4bool culling = fpMP->IsCulling();
  G4bool cullingInvisible = fpMP->IsCullingInvisible();
  G4bool markedVisible
  = pVisAttribs->IsVisible() && pVisAttribs->GetColour().GetAlpha() > 0;
  G4bool cullingLowDensity = fpMP->IsDensityCulling();
  G4double density = pMaterial? pMaterial->GetDensity(): 0;
  G4double densityCut = fpMP -> GetVisibleDensity ();

  // 1) Global culling is on....
  if (culling) {
    // 2) Culling of invisible volumes is on...
    if (cullingInvisible) {
      // 3) ...and the volume is marked not visible...
      if (!markedVisible) thisToBeDrawn = false;
    }
    // 4) Or culling of low density volumes is on...
    if (cullingLowDensity) {
      // 5) ...and density is less than cut value...
      if (density < densityCut) thisToBeDrawn = false;
    }
  }
  // 6) The user has asked for all further traversing to be aborted...
  if (fAbort) thisToBeDrawn = false;

  // Set "drawn" flag (it was true by default) - thisToBeDrawn may be false
  nodeID.SetDrawn(thisToBeDrawn);

  if (thisToBeDrawn) {

    // Update path of drawn physical volumes...
    fDrawnPVPath.push_back(nodeID);

    if (fpMP->IsExplode() && fDrawnPVPath.size() == 1) {
      // For top-level drawn volumes, explode along radius...
      G4Transform3D centering = G4Translate3D(fpMP->GetExplodeCentre());
      G4Transform3D centred = centering.inverse() * theNewAT;
      G4Scale3D oldScale;
      G4Rotate3D oldRotation;
      G4Translate3D oldTranslation;
      centred.getDecomposition(oldScale, oldRotation, oldTranslation);
      G4double explodeFactor = fpMP->GetExplodeFactor();
      G4Translate3D newTranslation =
	G4Translate3D(explodeFactor * oldTranslation.dx(),
		      explodeFactor * oldTranslation.dy(),
		      explodeFactor * oldTranslation.dz());
      theNewAT = centering * newTranslation * oldRotation * oldScale;
    }

    volumeCount++;
    DescribeSolid (theNewAT, pSol, pVisAttribs, sceneHandler);

  }

  // Make decision to draw daughters, if any.  There are various
  // reasons why daughters might not be drawn...

  // First, reasons that do not depend on culling policy...
  G4int nDaughters = (G4int)pLV->GetNoDaughters();
  G4bool daughtersToBeDrawn = true;
  // 1) There are no daughters...
  if (!nDaughters) daughtersToBeDrawn = false;
  // 2) We are at the limit if requested depth...
  else if (requestedDepth == 0) daughtersToBeDrawn = false;
  // 3) The user has asked for all further traversing to be aborted...
  else if (fAbort) daughtersToBeDrawn = false;
  // 4) The user has asked that the descent be curtailed...
  else if (fCurtailDescent) daughtersToBeDrawn = false;

  // Now, reasons that depend on culling policy...
  else {
    G4bool daughtersInvisible = pVisAttribs->IsDaughtersInvisible();
    // Culling of covered daughters request.  This is computed in
    // G4VSceneHandler::CreateModelingParameters() depending on view
    // parameters...
    G4bool cullingCovered = fpMP->IsCullingCovered();
    G4bool surfaceDrawing =
      fpMP->GetDrawingStyle() == G4ModelingParameters::hsr ||
      fpMP->GetDrawingStyle() == G4ModelingParameters::hlhsr;    
    if (pVisAttribs->IsForceDrawingStyle()) {
      switch (pVisAttribs->GetForcedDrawingStyle()) {
      default:
      case G4VisAttributes::wireframe: surfaceDrawing = false; break;
      case G4VisAttributes::solid: surfaceDrawing = true; break;
      }
    }
    G4bool opaque = pVisAttribs->GetColour().GetAlpha() >= 1.;
    // 5) Global culling is on....
    if (culling) {
      // 6) ..and culling of invisible volumes is on...
      if (cullingInvisible) {
	// 7) ...and the mother requests daughters invisible
	if (daughtersInvisible) daughtersToBeDrawn = false;
      }
      // 8) Or culling of covered daughters is requested...
      if (cullingCovered) {
	// 9) ...and surface drawing is operating...
	if (surfaceDrawing) {
	  // 10) ...but only if mother is visible...
	  if (thisToBeDrawn) {
	    // 11) ...and opaque...
	      if (opaque) daughtersToBeDrawn = false;
	  }
	}
      }
    }
  }

  if (daughtersToBeDrawn) {
    for (G4int iDaughter = 0; iDaughter < nDaughters; iDaughter++) {
      // Store daughter pVPV in local variable ready for recursion...
      G4VPhysicalVolume* pDaughterVPV = pLV -> GetDaughter (iDaughter);
      // Descend the geometry structure recursively...
      fCurrentDepth++;
      VisitGeometryAndGetVisReps
	(pDaughterVPV, requestedDepth - 1, theNewAT, sceneHandler);
      fCurrentDepth--;
    }
  }

  // Clean up
  delete tempVisAtts;

  // Reset for normal descending of next volume at this level...
  fCurtailDescent = false;

  // Pop item from paths physical volumes...
  fFullPVPath.pop_back();
  if (thisToBeDrawn) {
    fDrawnPVPath.pop_back();
  }
}

void G4PhysicalVolumeModel::DescribeSolid
(const G4Transform3D& theAT,
 G4VSolid* pSol,
 const G4VisAttributes* pVisAttribs,
 G4VGraphicsScene& sceneHandler)
{
  G4DisplacedSolid* pSectionSolid = fpMP->GetSectionSolid();
  G4DisplacedSolid* pCutawaySolid = fpMP->GetCutawaySolid();

  if (!fpClippingSolid && !pSectionSolid && !pCutawaySolid) {

    sceneHandler.PreAddSolid (theAT, *pVisAttribs);
    pSol -> DescribeYourselfTo (sceneHandler);  // Standard treatment.
    sceneHandler.PostAddSolid ();

  } else {

    G4VSolid* pResultantSolid = nullptr;

    if (fpClippingSolid) {
      switch (fClippingMode) {
        default:
        case subtraction:
          pResultantSolid = new G4SubtractionSolid
          ("subtracted_clipped_solid", pSol, fpClippingSolid, theAT.inverse());
          break;
        case intersection:
          pResultantSolid = new G4IntersectionSolid
          ("intersected_clipped_solid", pSol, fpClippingSolid, theAT.inverse());
          break;
      }
    }

    if (pSectionSolid) {
      pResultantSolid = new G4IntersectionSolid
      ("sectioned_solid", pSol, pSectionSolid, theAT.inverse());
    }

    if (pCutawaySolid) {
      pResultantSolid = new G4SubtractionSolid
      ("cutaway_solid", pSol, pCutawaySolid, theAT.inverse());
    }

    sceneHandler.PreAddSolid (theAT, *pVisAttribs);
    pResultantSolid -> DescribeYourselfTo (sceneHandler);
    sceneHandler.PostAddSolid ();

    delete pResultantSolid;
  }
}

G4bool G4PhysicalVolumeModel::Validate (G4bool warn)
{
// Not easy to see how to validate this sort of model. Previously there was
// a check that a volume of the same name (fTopPVName) existed somewhere in
// the geometry tree but under some circumstances this consumed lots of CPU
// time. Instead, let us simply check that the volume (fpTopPV) exists in the
// physical volume store.
  const auto& pvStore = G4PhysicalVolumeStore::GetInstance();
  auto iterator = find(pvStore->begin(),pvStore->end(),fpTopPV);
  if (iterator == pvStore->end()) {
    if (warn) {
      G4ExceptionDescription ed;
      ed << "Attempt to validate a volume that is no longer in the physical volume store.";
      G4Exception("G4PhysicalVolumeModel::Validate", "modeling0015", JustWarning, ed);
    }
    return false;
  } else {
    return true;
  }
}

const std::map<G4String,G4AttDef>* G4PhysicalVolumeModel::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store
      = G4AttDefStore::GetInstance("G4PhysicalVolumeModel", isNew);
    if (isNew) {
      (*store)["PVPath"] =
      G4AttDef("PVPath","Physical Volume Path","Physics","","G4String");
      (*store)["BasePVPath"] =
      G4AttDef("BasePVPath","Base Physical Volume Path","Physics","","G4String");
      (*store)["LVol"] =
      G4AttDef("LVol","Logical Volume","Physics","","G4String");
      (*store)["Solid"] =
      G4AttDef("Solid","Solid Name","Physics","","G4String");
      (*store)["EType"] =
      G4AttDef("EType","Entity Type","Physics","","G4String");
      (*store)["DmpSol"] =
      G4AttDef("DmpSol","Dump of Solid properties","Physics","","G4String");
      (*store)["LocalTrans"] =
      G4AttDef("LocalTrans","Local transformation of volume","Physics","","G4String");
      (*store)["LocalExtent"] =
      G4AttDef("LocalExtent","Local extent of volume","Physics","","G4String");
      (*store)["GlobalTrans"] =
      G4AttDef("GlobalTrans","Global transformation of volume","Physics","","G4String");
      (*store)["GlobalExtent"] =
      G4AttDef("GlobalExtent","Global extent of volume","Physics","","G4String");
      (*store)["Material"] =
      G4AttDef("Material","Material Name","Physics","","G4String");
      (*store)["Density"] =
      G4AttDef("Density","Material Density","Physics","G4BestUnit","G4double");
      (*store)["State"] =
      G4AttDef("State","Material State (enum undefined,solid,liquid,gas)","Physics","","G4String");
      (*store)["Radlen"] =
      G4AttDef("Radlen","Material Radiation Length","Physics","G4BestUnit","G4double");
      (*store)["Region"] =
      G4AttDef("Region","Cuts Region","Physics","","G4String");
      (*store)["RootRegion"] =
      G4AttDef("RootRegion","Root Region (0/1 = false/true)","Physics","","G4bool");
    }
  return store;
}

static std::ostream& operator<< (std::ostream& o, const G4Transform3D t)
{
  using namespace std;

  G4Scale3D sc;
  G4Rotate3D r;
  G4Translate3D tl;
  t.getDecomposition(sc, r, tl);

  const int w = 10;

  // Transformation itself
  o << setw(w) << t.xx() << setw(w) << t.xy() << setw(w) << t.xz() << setw(w) << t.dx() << endl;
  o << setw(w) << t.yx() << setw(w) << t.yy() << setw(w) << t.yz() << setw(w) << t.dy() << endl;
  o << setw(w) << t.zx() << setw(w) << t.zy() << setw(w) << t.zz() << setw(w) << t.dz() << endl;

  // Translation
  o << "= translation:" << endl;
  o << setw(w) << tl.dx() << setw(w) << tl.dy() << setw(w) << tl.dz() << endl;

  // Rotation
  o << "* rotation:" << endl;
  o << setw(w) << r.xx() << setw(w) << r.xy() << setw(w) << r.xz() << endl;
  o << setw(w) << r.yx() << setw(w) << r.yy() << setw(w) << r.yz() << endl;
  o << setw(w) << r.zx() << setw(w) << r.zy() << setw(w) << r.zz() << endl;

  // Scale
  o << "* scale:" << endl;
  o << setw(w) << sc.xx() << setw(w) << sc.yy() << setw(w) << sc.zz() << endl;

  // Transformed axes
  o << "Transformed axes:" << endl;
  o << "x': " << r * G4Vector3D(1., 0., 0.) << endl;
  o << "y': " << r * G4Vector3D(0., 1., 0.) << endl;
  o << "z': " << r * G4Vector3D(0., 0., 1.) << endl;

  return o;
}

std::vector<G4AttValue>* G4PhysicalVolumeModel::CreateCurrentAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  if (!fpCurrentLV) {
     G4Exception
        ("G4PhysicalVolumeModel::CreateCurrentAttValues",
         "modeling0004",
         JustWarning,
         "Current logical volume not defined.");
     return values;
  }

  std::ostringstream oss; oss << fFullPVPath;
  values->push_back(G4AttValue("PVPath", oss.str(),""));

  oss.str(""); oss << fBaseFullPVPath;
  values->push_back(G4AttValue("BasePVPath", oss.str(),""));

  values->push_back(G4AttValue("LVol", fpCurrentLV->GetName(),""));
  G4VSolid* pSol = fpCurrentLV->GetSolid();

  values->push_back(G4AttValue("Solid", pSol->GetName(),""));

  values->push_back(G4AttValue("EType", pSol->GetEntityType(),""));

  oss.str(""); oss << '\n' << *pSol;
  values->push_back(G4AttValue("DmpSol", oss.str(),""));

  const G4RotationMatrix localRotation = fpCurrentPV->GetObjectRotationValue();
  const G4ThreeVector& localTranslation = fpCurrentPV->GetTranslation();
  oss.str(""); oss << '\n' << G4Transform3D(localRotation,localTranslation);
  values->push_back(G4AttValue("LocalTrans", oss.str(),""));

  oss.str(""); oss << '\n' << pSol->GetExtent() << std::endl;
  values->push_back(G4AttValue("LocalExtent", oss.str(),""));

  oss.str(""); oss << '\n' << fCurrentTransform;
  values->push_back(G4AttValue("GlobalTrans", oss.str(),""));

  oss.str(""); oss << '\n' << (pSol->GetExtent()).Transform(fCurrentTransform) << std::endl;
  values->push_back(G4AttValue("GlobalExtent", oss.str(),""));

  G4String matName = fpCurrentMaterial? fpCurrentMaterial->GetName(): G4String("No material");
  values->push_back(G4AttValue("Material", matName,""));

  G4double matDensity = fpCurrentMaterial? fpCurrentMaterial->GetDensity(): 0.;
  values->push_back(G4AttValue("Density", G4BestUnit(matDensity,"Volumic Mass"),""));

  G4State matState = fpCurrentMaterial? fpCurrentMaterial->GetState(): kStateUndefined;
  oss.str(""); oss << matState;
  values->push_back(G4AttValue("State", oss.str(),""));

  G4double matRadlen = fpCurrentMaterial? fpCurrentMaterial->GetRadlen(): 0.;
  values->push_back(G4AttValue("Radlen", G4BestUnit(matRadlen,"Length"),""));

  G4Region* region = fpCurrentLV->GetRegion();
  G4String regionName = region? region->GetName(): G4String("No region");
  values->push_back(G4AttValue("Region", regionName,""));

  oss.str(""); oss << fpCurrentLV->IsRootRegion();
  values->push_back(G4AttValue("RootRegion", oss.str(),""));

  return values;
}

G4bool G4PhysicalVolumeModel::G4PhysicalVolumeNodeID::operator<
  (const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID& right) const
{
  if (fpPV < right.fpPV) return true;
  if (fpPV == right.fpPV) {
    if (fCopyNo < right.fCopyNo) return true;
    if (fCopyNo == right.fCopyNo)
      return fNonCulledDepth < right.fNonCulledDepth;
  }
  return false;
}

G4bool G4PhysicalVolumeModel::G4PhysicalVolumeNodeID::operator!=
  (const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID& right) const
{
  if (fpPV            != right.fpPV ||
      fCopyNo         != right.fCopyNo ||
      fNonCulledDepth != right.fNonCulledDepth ||
      fTransform      != right.fTransform ||
      fDrawn          != right.fDrawn) return true;
  return false;
}

std::ostream& operator<<
  (std::ostream& os, const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID& node)
{
  G4VPhysicalVolume* pPV = node.GetPhysicalVolume();
  if (pPV) {
    os << pPV->GetName()
       << ' ' << node.GetCopyNo()
//       << '[' << node.GetNonCulledDepth() << ']'
//       << ':' << node.GetTransform()
    ;
//    os << " (";
//    if (!node.GetDrawn()) os << "not ";
//    os << "drawn)";
  } else {
    os << " (Null PV node)";
  }
  return os;
}

std::ostream& operator<<
(std::ostream& os, const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& path)
{
  if (path.empty()) {
    os << " TOP";
  } else {
    for (const auto& nodeID: path) {
      os << ' ' << nodeID;
    }
  }
  return os;
}

G4PhysicalVolumeModel::G4PhysicalVolumeModelTouchable::G4PhysicalVolumeModelTouchable
(const std::vector<G4PhysicalVolumeNodeID>& fullPVPath):
  fFullPVPath(fullPVPath) {}

const G4ThreeVector& G4PhysicalVolumeModel::G4PhysicalVolumeModelTouchable::GetTranslation(G4int depth) const
{
  size_t i = fFullPVPath.size() - depth - 1;
  if (i >= fFullPVPath.size()) {
    G4Exception("G4PhysicalVolumeModelTouchable::GetTranslation",
		"modeling0005",
		FatalErrorInArgument,
		"Index out of range. Asking for non-existent depth");
  }
  static G4ThreeVector tempTranslation;
  tempTranslation = fFullPVPath[i].GetTransform().getTranslation();
  return tempTranslation;
}

const G4RotationMatrix* G4PhysicalVolumeModel::G4PhysicalVolumeModelTouchable::GetRotation(G4int depth) const
{
  size_t i = fFullPVPath.size() - depth - 1;
  if (i >= fFullPVPath.size()) {
    G4Exception("G4PhysicalVolumeModelTouchable::GetRotation",
		"modeling0006",
		FatalErrorInArgument,
		"Index out of range. Asking for non-existent depth");
  }
  static G4RotationMatrix tempRotation;
  tempRotation = fFullPVPath[i].GetTransform().getRotation();
  return &tempRotation;
}

G4VPhysicalVolume* G4PhysicalVolumeModel::G4PhysicalVolumeModelTouchable::GetVolume(G4int depth) const
{
  size_t i = fFullPVPath.size() - depth - 1;
  if (i >= fFullPVPath.size()) {
    G4Exception("G4PhysicalVolumeModelTouchable::GetVolume",
		"modeling0007",
		FatalErrorInArgument,
		"Index out of range. Asking for non-existent depth");
  }
  return fFullPVPath[i].GetPhysicalVolume();
}

G4VSolid* G4PhysicalVolumeModel::G4PhysicalVolumeModelTouchable::GetSolid(G4int depth) const
{
  size_t i = fFullPVPath.size() - depth - 1;
  if (i >= fFullPVPath.size()) {
    G4Exception("G4PhysicalVolumeModelTouchable::GetSolid",
		"modeling0008",
		FatalErrorInArgument,
		"Index out of range. Asking for non-existent depth");
  }
  return fFullPVPath[i].GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
}

G4int G4PhysicalVolumeModel::G4PhysicalVolumeModelTouchable::GetReplicaNumber(G4int depth) const
{
  size_t i = fFullPVPath.size() - depth - 1;
  if (i >= fFullPVPath.size()) {
    G4Exception("G4PhysicalVolumeModelTouchable::GetReplicaNumber",
		"modeling0009",
		FatalErrorInArgument,
		"Index out of range. Asking for non-existent depth");
  }
  return fFullPVPath[i].GetCopyNo();
}

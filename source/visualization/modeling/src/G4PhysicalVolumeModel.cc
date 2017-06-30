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
// $Id: G4PhysicalVolumeModel.cc 103926 2017-05-03 13:43:27Z gcosmo $
//
// 
// John Allison  31st December 1997.
// Model for physical volumes.

#include "G4PhysicalVolumeModel.hh"

#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4BoundingSphereScene.hh"
#include "G4PhysicalVolumeSearchScene.hh"
#include "G4TransportationManager.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4Vector3D.hh"

#include <sstream>

G4PhysicalVolumeModel::G4PhysicalVolumeModel
(G4VPhysicalVolume*          pVPV
 , G4int                       requestedDepth
 , const G4Transform3D& modelTransformation
 , const G4ModelingParameters* pMP
 , G4bool useFullExtent)
: G4VModel           (modelTransformation, pMP)
, fpTopPV            (pVPV)
, fTopPVCopyNo       (0)
, fRequestedDepth    (requestedDepth)
, fUseFullExtent     (useFullExtent)
, fCurrentDepth      (0)
, fpCurrentPV        (0)
, fpCurrentLV        (0)
, fpCurrentMaterial  (0)
, fpCurrentTransform (0)
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
    fTopPVCopyNo = fpTopPV -> GetCopyNo ();
    std::ostringstream o;
    o << fpTopPV -> GetCopyNo ();
    fGlobalTag = fpTopPV -> GetName () + "." + o.str();
    fGlobalDescription = "G4PhysicalVolumeModel " + fGlobalTag;

    fpCurrentPV = fpTopPV;
    if (fpCurrentPV) fpCurrentLV = fpCurrentPV->GetLogicalVolume();
    if (fpCurrentLV) fpCurrentMaterial = fpCurrentLV->GetMaterial();
    fpCurrentTransform = const_cast<G4Transform3D*>(&modelTransformation);

    CalculateExtent ();
  }
}

G4PhysicalVolumeModel::~G4PhysicalVolumeModel ()
{
  delete fpClippingSolid;
}

void G4PhysicalVolumeModel::CalculateExtent ()
{
  if (fUseFullExtent) {
    fExtent = fpTopPV -> GetLogicalVolume () -> GetSolid () -> GetExtent ();
  }
  else {
    G4BoundingSphereScene bsScene(this);
    const G4int tempRequestedDepth = fRequestedDepth;
    fRequestedDepth = -1;  // Always search to all depths to define extent.
    const G4ModelingParameters* tempMP = fpMP;
    G4ModelingParameters mParams
      (0,      // No default vis attributes needed.
       G4ModelingParameters::wf,  // wireframe (not relevant for this).
       true,   // Global culling.
       true,   // Cull invisible volumes.
       false,  // Density culling.
       0.,     // Density (not relevant if density culling false).
       true,   // Cull daughters of opaque mothers.
       72);    // No of sides (not relevant for this operation).
    fpMP = &mParams;
    DescribeYourselfTo (bsScene);
    G4double radius = bsScene.GetRadius();
    if (radius < 0.) {  // Nothing in the scene.
      fExtent = fpTopPV -> GetLogicalVolume () -> GetSolid () -> GetExtent ();
    } else {
      // Transform back to coordinates relative to the top
      // transformation, which is in G4VModel::fTransform.  This makes
      // it conform to all models, which are defined by a
      // transformation and an extent relative to that
      // transformation...
      G4Point3D centre = bsScene.GetCentre();
      centre.transform(fTransform.inverse());
      fExtent = G4VisExtent(centre, radius);
    }
    fpMP = tempMP;
    fRequestedDepth = tempRequestedDepth;
  }
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

  // For safety...
  fCurrentDepth = 0;

  G4Transform3D startingTransformation = fTransform;

  VisitGeometryAndGetVisReps
    (fpTopPV,
     fRequestedDepth,
     startingTransformation,
     sceneHandler);

  // Clear data...
  fCurrentDepth     = 0;
  fpCurrentPV       = 0;
  fpCurrentLV       = 0;
  fpCurrentMaterial = 0;
  if (fFullPVPath.size() != fBaseFullPVPath.size()) {
    // They should be equal if pushing and popping is happening properly.
    G4Exception
    ("G4PhysicalVolumeModel::DescribeYourselfTo",
     "modeling0013",
     FatalException,
     "Path at start of modeling not equal to base path.  Something badly"
     "\nwrong.  Please contact visualisation coordinator.");
  }
  fDrawnPVPath.clear();
  fAbort = false;
  fCurtailDescent = false;
}

G4String G4PhysicalVolumeModel::GetCurrentTag () const
{
  if (fpCurrentPV) {
    std::ostringstream o;
    o << fpCurrentPV -> GetCopyNo ();
    return fpCurrentPV -> GetName () + "." + o.str();
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

  G4VSolid* pSol;
  G4Material* pMaterial;

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
    if (fCurrentDepth == 0) nReplicas = 1;  // Just draw first
    G4VPVParameterisation* pP = pVPV -> GetParameterisation ();
    if (pP) {  // Parametrised volume.
      for (int n = 0; n < nReplicas; n++) {
	pSol = pP -> ComputeSolid (n, pVPV);
	pP -> ComputeTransformation (n, pVPV);
	pSol -> ComputeDimensions (pP, n, pVPV);
	pVPV -> SetCopyNo (n);
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
      for (int n = 0; n < nReplicas; n++) {
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
	      G4cout <<
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
  fpCurrentLV = pLV;
  fpCurrentMaterial = pMaterial;

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
  fpCurrentTransform = &theNewAT;

  const G4VisAttributes* pVisAttribs = pLV->GetVisAttributes();
  if (!pVisAttribs) pVisAttribs = fpMP->GetDefaultVisAttributes();
  // Beware - pVisAttribs might still be zero - probably will, since that's
  // the default for G4ModelingParameters.  So create one if necessary...
  if (!pVisAttribs) {
    static G4VisAttributes defaultVisAttribs;
    pVisAttribs = &defaultVisAttribs;
  }
  // From here, can assume pVisAttribs is a valid pointer.  This is necessary
  // because PreAddSolid needs a vis attributes object.

  // Make decision to draw...
  G4bool thisToBeDrawn = true;
  
  // Update full path of physical volumes...
  G4int copyNo = fpCurrentPV->GetCopyNo();
  fFullPVPath.push_back
  (G4PhysicalVolumeNodeID
   (fpCurrentPV,copyNo,fCurrentDepth,*fpCurrentTransform));
  
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

  // There are various reasons why this volume
  // might not be drawn...
  G4bool culling = fpMP->IsCulling();
  G4bool cullingInvisible = fpMP->IsCullingInvisible();
  G4bool markedVisible = pVisAttribs->IsVisible();
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

  // Record thisToBeDrawn in path...
  fFullPVPath.back().SetDrawn(thisToBeDrawn);

  if (thisToBeDrawn) {

    // Update path of drawn physical volumes...
    fDrawnPVPath.push_back
      (G4PhysicalVolumeNodeID
       (fpCurrentPV,copyNo,fCurrentDepth,*fpCurrentTransform,thisToBeDrawn));

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

    DescribeSolid (theNewAT, pSol, pVisAttribs, sceneHandler);

  }

  // Make decision to draw daughters, if any.  There are various
  // reasons why daughters might not be drawn...

  // First, reasons that do not depend on culling policy...
  G4int nDaughters = pLV->GetNoDaughters();
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
  sceneHandler.PreAddSolid (theAT, *pVisAttribs);

  G4VSolid* pSectionSolid = fpMP->GetSectionSolid();
  G4VSolid* pCutawaySolid = fpMP->GetCutawaySolid();

  if (!fpClippingSolid && !pSectionSolid && !pCutawaySolid) {

    pSol -> DescribeYourselfTo (sceneHandler);  // Standard treatment.

  } else {

    // Clipping, etc., performed by Boolean operations.

    // First, get polyhedron for current solid...
    if (pVisAttribs->IsForceLineSegmentsPerCircle())
      G4Polyhedron::SetNumberOfRotationSteps
	(pVisAttribs->GetForcedLineSegmentsPerCircle());
    else
      G4Polyhedron::SetNumberOfRotationSteps(fpMP->GetNoOfSides());
    const G4Polyhedron* pOriginal = pSol->GetPolyhedron();
    G4Polyhedron::ResetNumberOfRotationSteps();

    if (!pOriginal) {

      if (fpMP->IsWarning())
	G4cout <<
	  "WARNING: G4PhysicalVolumeModel::DescribeSolid: solid\n  \""
	       << pSol->GetName() <<
	  "\" has no polyhedron.  Cannot by clipped."
	       << G4endl;
      pSol -> DescribeYourselfTo (sceneHandler);  // Standard treatment.

    } else {

      G4Polyhedron resultant(*pOriginal);
      G4VisAttributes resultantVisAttribs(*pVisAttribs);
      G4VSolid* resultantSolid = 0;

      if (fpClippingSolid) {
	switch (fClippingMode) {
	default:
	case subtraction:
	  resultantSolid = new G4SubtractionSolid
	    ("resultant_solid", pSol, fpClippingSolid, theAT.inverse());
	  break;
	case intersection:
	  resultantSolid = new G4IntersectionSolid
	    ("resultant_solid", pSol, fpClippingSolid, theAT.inverse());
	  break;
	}
      }

      if (pSectionSolid) {
	resultantSolid = new G4IntersectionSolid
	  ("sectioned_solid", pSol, pSectionSolid, theAT.inverse());
      }

      if (pCutawaySolid) {
	resultantSolid = new G4SubtractionSolid
	  ("cutaway_solid", pSol, pCutawaySolid, theAT.inverse());
      }

      G4Polyhedron* tmpResultant = resultantSolid->GetPolyhedron();
      if (tmpResultant) resultant = *tmpResultant;
      else {
	if (fpMP->IsWarning())
	  G4cout <<
	    "WARNING: G4PhysicalVolumeModel::DescribeSolid: resultant polyhedron for"
	    "\n  solid \"" << pSol->GetName() <<
	    "\" not defined due to error during Boolean processing."
	    "\n  Original will be drawn in red."
		 << G4endl;
	resultantVisAttribs.SetColour(G4Colour::Red());
      }

      delete resultantSolid;

      // Finally, force polyhedron drawing...
      resultant.SetVisAttributes(resultantVisAttribs);
      sceneHandler.BeginPrimitives(theAT);
      sceneHandler.AddPrimitive(resultant);
      sceneHandler.EndPrimitives();
    }
  }
  sceneHandler.PostAddSolid ();
}

G4bool G4PhysicalVolumeModel::Validate (G4bool warn)
{
  G4TransportationManager* transportationManager =
    G4TransportationManager::GetTransportationManager ();

  size_t nWorlds = transportationManager->GetNoWorlds();

  G4bool found = false;

  std::vector<G4VPhysicalVolume*>::iterator iterWorld =
    transportationManager->GetWorldsIterator();
  for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4VPhysicalVolume* world = (*iterWorld);
    if (!world) break; // This can happen if geometry has been cleared/destroyed.
    // The idea now is to seek a PV with the same name and copy no
    // in the hope it's the same one!!
    G4PhysicalVolumeModel searchModel (world);
    G4int verbosity = 0;  // Suppress messages from G4PhysicalVolumeSearchScene.
    G4PhysicalVolumeSearchScene searchScene
      (&searchModel, fTopPVName, fTopPVCopyNo, verbosity);
    G4ModelingParameters mp;  // Default modeling parameters for this search.
    mp.SetDefaultVisAttributes(fpMP? fpMP->GetDefaultVisAttributes(): 0);
    searchModel.SetModelingParameters (&mp);
    searchModel.DescribeYourselfTo (searchScene);
    G4VPhysicalVolume* foundVolume = searchScene.GetFoundVolume ();
    if (foundVolume) {
      if (foundVolume != fpTopPV && warn) {
	G4cout <<
	  "G4PhysicalVolumeModel::Validate(): A volume of the same name and"
	  "\n  copy number (\""
	       << fTopPVName << "\", copy " << fTopPVCopyNo
	       << ") still exists and is being used."
	  "\n  But it is not the same volume you originally specified"
	  "\n  in /vis/scene/add/."
	       << G4endl;
      }
      fpTopPV = foundVolume;
      CalculateExtent ();
      found = true;
    }
  }
  if (found) return true;
  else {
    if (warn) {
      G4cout <<
	"G4PhysicalVolumeModel::Validate(): No volume of name and"
	"\n  copy number (\""
	     << fTopPVName << "\", copy " << fTopPVCopyNo
	     << ") exists."
	     << G4endl;
    }
    return false;
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
      (*store)["GlobalTrans"] =
    G4AttDef("GlobalTrans","Global transformation of volume","Physics","","G4String");
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

#include <iomanip>

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
  std::ostringstream oss;
  for (size_t i = 0; i < fFullPVPath.size(); ++i) {
    oss << fFullPVPath[i].GetPhysicalVolume()->GetName()
	<< ':' << fFullPVPath[i].GetCopyNo();
    if (i != fFullPVPath.size() - 1) oss << '/';
  }

  if (!fpCurrentLV) {
     G4Exception
        ("G4PhysicalVolumeModel::CreateCurrentAttValues",
         "modeling0004",
         JustWarning,
         "Current logical volume not defined.");
     return values;
  }

  values->push_back(G4AttValue("PVPath", oss.str(),""));
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
  oss.str(""); oss << '\n' << *fpCurrentTransform;
  values->push_back(G4AttValue("GlobalTrans", oss.str(),""));
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
    os << " (Null node)";
  }
  return os;
}

std::ostream& operator<<
(std::ostream& os, const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& path)
{
  for (const auto& nodeID: path) {
    os << nodeID << ' ';
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

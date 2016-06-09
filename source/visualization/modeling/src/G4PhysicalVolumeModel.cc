//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PhysicalVolumeModel.cc,v 1.36 2005/11/24 11:15:21 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4BoundingSphereScene.hh"
#include "G4PhysicalVolumeSearchScene.hh"
#include "G4TransportationManager.hh"
#include "G4VVisManager.hh"
#include "G4Polyhedron.hh"

#include <sstream>

G4PhysicalVolumeModel::G4PhysicalVolumeModel
(G4VPhysicalVolume*          pVPV,
 G4int                       requestedDepth,
 const G4Transform3D& modelTransformation,
 const G4ModelingParameters* pMP,
 G4bool useFullExtent):
  G4VModel        (modelTransformation, pMP),
  fpTopPV         (pVPV),
  fTopPVName      (pVPV -> GetName ()),
  fTopPVCopyNo    (pVPV -> GetCopyNo ()),
  fRequestedDepth (requestedDepth),
  fUseFullExtent  (useFullExtent),
  fCurrentDepth   (0),
  fpCurrentPV     (0),
  fpCurrentLV     (0),
  fCurtailDescent (false),
  fpClippingPolyhedron (0),
  fpCurrentDepth  (0),
  fppCurrentPV    (0),
  fppCurrentLV    (0),
  fppCurrentMaterial (0)
{
  std::ostringstream o;
  o << fpTopPV -> GetCopyNo ();
  fGlobalTag = fpTopPV -> GetName () + "." + o.str();
  fGlobalDescription = "G4PhysicalVolumeModel " + fGlobalTag;

  CalculateExtent ();
}

G4PhysicalVolumeModel::~G4PhysicalVolumeModel () {}

void G4PhysicalVolumeModel::CalculateExtent () {
  if (fUseFullExtent) {
    fExtent = fpTopPV -> GetLogicalVolume () -> GetSolid () -> GetExtent ();
  }
  else {
    G4BoundingSphereScene bsScene(this);
    const G4int tempRequestedDepth = fRequestedDepth;
    fRequestedDepth = -1;  // Always search to all depths to define extent.
    const G4ModelingParameters* tempMP = fpMP;
    G4ModelingParameters mParams
      (0,      // No default vis attributes.
       G4ModelingParameters::wireframe,
       true,   // Global culling.
       true,   // Cull invisible volumes.
       false,  // Density culling.
       0.,     // Density (not relevant if density culling false).
       true,   // Cull daughters of opaque mothers.
       24,     // No of sides (not relevant for this operation).
       true,   // View geometry.
       false,  // View hits - not relevant for physical volume model.
       false); // View digis - not relevant for physical volume model.
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
(G4VGraphicsScene& sceneHandler) {

  if (fpMP && fpMP -> IsViewGeom ()) {

    sceneHandler.EstablishSpecials (*this);
    // See .hh file for explanation of this mechanism.

    fCurrentDepth = 0;
    // Store in working space (via pointer to working space).
    if (fpCurrentDepth) *fpCurrentDepth = fCurrentDepth;

    G4Transform3D startingTransformation = fTransform;

    VisitGeometryAndGetVisReps (fpTopPV,
				fRequestedDepth,
				startingTransformation,
				sceneHandler);

    // Clear current data and working space (via pointers to working space).
    fCurrentDepth = 0;
    fpCurrentPV   = 0;
    fpCurrentLV   = 0;
    if (fpCurrentDepth) *fpCurrentDepth = fCurrentDepth;
    if (fppCurrentPV)   *fppCurrentPV   = fpCurrentPV;
    if (fppCurrentLV)   *fppCurrentLV   = fpCurrentLV;
    if (fppCurrentMaterial) *fppCurrentMaterial = fpCurrentMaterial;
    sceneHandler.DecommissionSpecials (*this);

    // Clear pointers to working space.
    fpCurrentDepth = 0;
    fppCurrentPV   = 0;
    fppCurrentLV   = 0;
    fppCurrentMaterial = 0;
  }
}

G4String G4PhysicalVolumeModel::GetCurrentTag () const {
  if (fpCurrentPV) {
    std::ostringstream o;
    o << fpCurrentPV -> GetCopyNo ();
    return fpCurrentPV -> GetName () + "." + o.str();
  }
  else {
    return "WARNING: NO CURRENT VOLUME - global tag is " + fGlobalTag;
  }
}

G4String G4PhysicalVolumeModel::GetCurrentDescription () const {
  return "G4PhysicalVolumeModel " + GetCurrentTag ();
}

void G4PhysicalVolumeModel::DefinePointersToWorkingSpace
(G4int*              pCurrentDepth,
 G4VPhysicalVolume** ppCurrentPV,
 G4LogicalVolume**   ppCurrentLV,
 G4Material**        ppCurrentMaterial) {
  fpCurrentDepth = pCurrentDepth;
  fppCurrentPV   = ppCurrentPV;
  fppCurrentLV   = ppCurrentLV;
  fppCurrentMaterial = ppCurrentMaterial;
}

void G4PhysicalVolumeModel::VisitGeometryAndGetVisReps
(G4VPhysicalVolume* pVPV,
 G4int requestedDepth,
 const G4Transform3D& theAT,
 G4VGraphicsScene& sceneHandler) {

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
    G4VPVParameterisation* pP = pVPV -> GetParameterisation ();
    if (pP) {  // Parametrised volume.
      for (int n = 0; n < nReplicas; n++) {
	pSol = pP -> ComputeSolid (n, pVPV);
	pMaterial = pP -> ComputeMaterial (n, pVPV);
	pP -> ComputeTransformation (n, pVPV);
	pSol -> ComputeDimensions (pP, n, pVPV);
	pVPV -> SetCopyNo (n);
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
	G4ThreeVector translation;  // Null.
	G4RotationMatrix rotation;  // Null - life long enough for visualizing.
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

  return;
}

void G4PhysicalVolumeModel::DescribeAndDescend
(G4VPhysicalVolume* pVPV,
 G4int requestedDepth,
 G4LogicalVolume* pLV,
 G4VSolid* pSol,
 G4Material* pMaterial,
 const G4Transform3D& theAT,
 G4VGraphicsScene& sceneHandler) {

  // Maintain data members and store in working space (via pointers to
  // working space).
  if (fpCurrentDepth) *fpCurrentDepth = fCurrentDepth;
  fpCurrentPV = pVPV;
  fpCurrentLV = pLV;
  fpCurrentMaterial = pMaterial;
  if (fppCurrentPV) *fppCurrentPV = fpCurrentPV;
  if (fppCurrentLV) *fppCurrentLV = fpCurrentLV;
  if (fppCurrentMaterial) *fppCurrentMaterial = fpCurrentMaterial;

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

  /********************************************************
  G4cout << "G4PhysicalVolumeModel::DescribeAndDescend: "
	 << pVPV -> GetName () << "." << pVPV -> GetCopyNo ();
  G4cout << "\n  theAT: ";
  G4cout << "\n    Rotation: ";
  G4RotationMatrix rotation = theAT.getRotation ();
  G4cout << rotation.thetaX() << ", "
	 << rotation.phiX() << ", "
	 << rotation.thetaY() << ", "
	 << rotation.phiY() << ", "
	 << rotation.thetaZ() << ", "
	 << rotation.phiZ();
  G4cout << "\n    Translation: " << theAT.getTranslation();
  G4cout << "\n  theNewAT: ";
  G4cout << "\n    Rotation: ";
  rotation = theNewAT.getRotation ();
  G4cout << rotation.thetaX() << ", "
	 << rotation.phiX() << ", "
	 << rotation.thetaY() << ", "
	 << rotation.phiY() << ", "
	 << rotation.thetaZ() << ", "
	 << rotation.phiZ();
  G4cout << "\n    Translation: " << theNewAT.getTranslation();
  G4cout << G4endl;
  **********************************************************/

  // Make decision to Draw.
  G4bool thisToBeDrawn = !IsThisCulled (pLV, pMaterial);
  if (thisToBeDrawn) {
    const G4VisAttributes* pVisAttribs = pLV -> GetVisAttributes ();
    if (!pVisAttribs) pVisAttribs = fpMP -> GetDefaultVisAttributes ();
    DescribeSolid (theNewAT, pSol, pVisAttribs, sceneHandler);
  }

  if (fCurtailDescent) {
    fCurtailDescent = false;
    // Reset for normal descending of next volume at this level.
    return;
  }

  // First check if mother covers...

  // This is only effective in surface drawing style, and then only if
  // the volumes are visible and opaque, and then only if no sections
  // or cutways are in operation.
  G4bool cullDaughter = thisToBeDrawn && IsDaughterCulled (pLV);
  if (!cullDaughter) {
    // OK, now let's check for daughters...
    if (requestedDepth != 0) {
      int nDaughters = pLV -> GetNoDaughters ();
      if (nDaughters) {
	for (int iDaughter = 0; iDaughter < nDaughters; iDaughter++) {
	  G4VPhysicalVolume* pVPV = pLV -> GetDaughter (iDaughter);
	  // Descend the geometry structure recursively...
	  fCurrentDepth++;
	  VisitGeometryAndGetVisReps
	    (pVPV, requestedDepth - 1, theNewAT, sceneHandler);
	  fCurrentDepth--;
	}
      }
    }
  }
}

void G4PhysicalVolumeModel::DescribeSolid
(const G4Transform3D& theAT,
 G4VSolid* pSol,
 const G4VisAttributes* pVisAttribs,
 G4VGraphicsScene& sceneHandler) {

  sceneHandler.PreAddSolid (theAT, *pVisAttribs);
  if (fpClippingPolyhedron)  // Clip and force polyhedral representation...
    {
      G4Polyhedron clipper(*fpClippingPolyhedron);  // Local copy.
      clipper.Transform(theAT.inverse());

      G4Polyhedron::SetNumberOfRotationSteps (fpMP->GetNoOfSides());
      G4Polyhedron* pClippee = pSol->GetPolyhedron();
      G4Polyhedron::ResetNumberOfRotationSteps ();
      if (pClippee)  // Solid can provide.
	{
	  G4Polyhedron clipped(pClippee->subtract(clipper));
	  if(clipped.IsErrorBooleanProcess())
	    {
	      G4cout <<
 "WARNING: G4PhysicalVolumeModel::DescribeSolid: polyhedron for solid\n  \""
		     << pSol->GetName() <<
 "\" skipped due to error during Boolean processing."
		     << G4endl;
	    }
	  else
	    {
	      clipped.SetVisAttributes(pVisAttribs);
	      G4VVisManager::GetConcreteInstance()->Draw(clipped,theAT);
	    }
	}
    }
  else  // Standard treatment...
    {
      pSol -> DescribeYourselfTo (sceneHandler);
      
    }
  sceneHandler.PostAddSolid ();
}

G4bool G4PhysicalVolumeModel::IsThisCulled (const G4LogicalVolume* pLV,
					    const G4Material* pMaterial) {
  // If true, cull, i.e., do not Draw.
  G4double density = 0.;
  if (pMaterial) density = pMaterial -> GetDensity ();
  const G4VisAttributes* pVisAttribs = pLV -> GetVisAttributes ();
  if (!pVisAttribs) pVisAttribs = fpMP -> GetDefaultVisAttributes ();
  if (fpMP) {
    return
      fpMP -> IsCulling () &&       // Global culling flag.
      (   
       // Invisible volumes...
       (fpMP -> IsCullingInvisible () &&
	!(pVisAttribs ? pVisAttribs -> IsVisible () : true)) ||

       // Low density volumes...
       (fpMP -> IsDensityCulling () &&
	(density < fpMP -> GetVisibleDensity ()))
       )
      ;
  }
  else {
    return false;
  }
}

G4bool G4PhysicalVolumeModel::IsDaughterCulled
(const G4LogicalVolume* pMotherLV) {
  // If true, cull, i.e., do not Draw.
  const G4VisAttributes* pVisAttribs = pMotherLV -> GetVisAttributes ();
  if (!pVisAttribs) pVisAttribs = fpMP -> GetDefaultVisAttributes ();
  if (fpMP) {
    return
      fpMP -> IsCulling ()           // Global culling flag.
      &&
      (
       // Does mother request daughters not to be drawn?
       (pVisAttribs ? pVisAttribs -> IsDaughtersInvisible () : false)
       ||
       (
	// Global covered daughter flag.  This is affected by drawing
	// style, etc.  The enforcing of this is done in
	// G4VScene::CreateModelingParameters ()
	fpMP -> IsCullingCovered ()
	&&
	// Cull only if mother is visible...
	(pVisAttribs ? pVisAttribs -> IsVisible () : true)
	&&
	// ...and opaque...
	(pVisAttribs ? (pVisAttribs -> GetColour ().GetAlpha() >= 1.) : true)
	)
       )
      ;
  }
  else {
    return false;
  }
}

G4bool G4PhysicalVolumeModel::Validate (G4bool warn) {
  G4VPhysicalVolume* world =
    G4TransportationManager::GetTransportationManager ()
    -> GetNavigatorForTracking () -> GetWorldVolume ();
  // The idea now is to seek a PV with the same name and copy no
  // in the hope it's the same one!!
  if (warn) {
    G4cout << "G4PhysicalVolumeModel::Validate() called." << G4endl;
  }
  G4PhysicalVolumeSearchScene searchScene (fTopPVName, fTopPVCopyNo);
  G4PhysicalVolumeModel searchModel (world);
  G4ModelingParameters mp;  // Default modeling parameters for this search.
  searchModel.SetModelingParameters (&mp);
  searchModel.DescribeYourselfTo (searchScene);
  G4VPhysicalVolume* foundVolume = searchScene.GetFoundVolume ();
  if (foundVolume) {
    if (warn) {
      G4cout << "  Volume of the same name and copy number (\""
	     << fTopPVName << "\", copy " << fTopPVCopyNo
	     << ") still exists and is being used."
	"\n  Be warned that this does not necessarily guarantee it's the same"
	"\n  volume you originally specified in /vis/scene/add/."
	     << G4endl;
    }
    fpTopPV = foundVolume;
    CalculateExtent ();
    return true;
  }
  else {
    if (warn) {
      G4cout << "  A volume of the same name and copy number (\""
	     << fTopPVName << "\", copy " << fTopPVCopyNo
	     << ") no longer exists."
	     << G4endl;
    }
    return false;
  }
}

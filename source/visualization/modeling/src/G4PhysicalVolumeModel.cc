// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicalVolumeModel.cc,v 1.12 2000-10-18 13:57:45 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "g4std/strstream"

G4PhysicalVolumeModel::G4PhysicalVolumeModel
(G4VPhysicalVolume*          pVPV,
 G4int                       soughtDepth,
 const G4Transform3D& modelTransformation,
 const G4ModelingParameters* pMP,
 G4bool useFullExtent):
  G4VModel       (modelTransformation, pMP),
  fpTopPV        (pVPV),
  fTopPVName     (pVPV -> GetName ()),
  fTopPVCopyNo   (pVPV -> GetCopyNo ()),
  fSoughtDepth   (soughtDepth),
  fCurrentDepth  (0),
  fpCurrentPV    (0),
  fpCurrentLV    (0),
  fpCurrentDepth (0),
  fppCurrentPV   (0),
  fppCurrentLV   (0)
{
  const int len = 8; char a [len];
  G4std::ostrstream o (a, len); o.seekp (G4std::ios::beg);
  o << fpTopPV -> GetCopyNo () << G4std::ends;
  fGlobalTag = fpTopPV -> GetName () + "." + a;
  fGlobalDescription = "G4PhysicalVolumeModel " + fGlobalTag;

  if (useFullExtent) {
    fExtent = fpTopPV -> GetLogicalVolume () -> GetSolid () -> GetExtent ();
  }
  else {
    G4BoundingSphereScene bsScene;
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
    fExtent = bsScene.GetBoundingSphereExtent ();
    fpMP = tempMP;
  }
}

G4PhysicalVolumeModel::~G4PhysicalVolumeModel () {}

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
				fSoughtDepth,
				startingTransformation,
				sceneHandler);

    // Clear current data and working space (via pointers to working space).
    fCurrentDepth = 0;
    fpCurrentPV   = 0;
    fpCurrentLV   = 0;
    if (fpCurrentDepth) *fpCurrentDepth = fCurrentDepth;
    if (fppCurrentPV)   *fppCurrentPV   = fpCurrentPV;
    if (fppCurrentLV)   *fppCurrentLV   = fpCurrentLV;

    sceneHandler.DecommissionSpecials (*this);

    // Clear pointers to working space.
    fpCurrentDepth = 0;
    fppCurrentPV   = 0;
    fppCurrentLV   = 0;
  }
}

G4String G4PhysicalVolumeModel::GetCurrentTag () const {
  const int len = 8; char a [len];
  G4std::ostrstream o (a, len); o.seekp (G4std::ios::beg);
  if (fpCurrentPV) {
    o << fpCurrentPV -> GetCopyNo () << G4std::ends;
    return fpCurrentPV -> GetName () + "." + a;
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
 G4LogicalVolume**   ppCurrentLV) {
  fpCurrentDepth = pCurrentDepth;
  fppCurrentPV   = ppCurrentPV;
  fppCurrentLV   = ppCurrentLV;
}

void G4PhysicalVolumeModel::VisitGeometryAndGetVisReps
(G4VPhysicalVolume* pVPV,
 G4int soughtDepth,
 const G4Transform3D& theAT,
 G4VGraphicsScene& sceneHandler) {

  // Visits geometry structure to a given depth (soughtDepth), starting
  //   at given physical volume with given starting transformation and
  //   describes volumes to the scene handler.
  // soughtDepth < 0 (default) implies full visit.
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
    DescribeAndDescend (pVPV, soughtDepth, pLV, pSol, pMaterial,
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
	// pVPV -> SetCopyNo (n);  // Uncertain of effect of this.
	DescribeAndDescend (pVPV, soughtDepth, pLV, pSol, pMaterial,
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
	  G4cerr <<
	    "G4PhysicalVolumeModel::VisitGeometryAndGetVisReps: WARNING:"
	    "\n  built-in replicated volumes replicated in radius are not yet"
	    "\n  properly visualizable."
	       << G4endl;
	  break;
	case kPhi:
	  rotation.rotateZ (-(offset+n*width));
	  // Minus Sign because for the physical volume we need the
	  // coordinate system rotation.
	  pRotation = &rotation;
	  break;
	} 
	pVPV -> SetTranslation (translation);
	pVPV -> SetRotation    (pRotation);
	// pVPV -> SetCopyNo (n); // Has no effect and might even be
	// dangerous.
	pSol = pLV -> GetSolid ();
	pMaterial = pLV -> GetMaterial ();
	DescribeAndDescend (pVPV, soughtDepth, pLV, pSol, pMaterial,
			    theAT, sceneHandler);
      }
    }
  }

  return;
}

void G4PhysicalVolumeModel::DescribeAndDescend
(G4VPhysicalVolume* pVPV,
 G4int soughtDepth,
 G4LogicalVolume* pLV,
 G4VSolid* pSol,
 const G4Material* pMaterial,
 const G4Transform3D& theAT,
 G4VGraphicsScene& sceneHandler) {

  // Maintain data members and store in working space (via pointers to
  // working space).
  if (fpCurrentDepth) *fpCurrentDepth = fCurrentDepth;
  fpCurrentPV = pVPV;
  fpCurrentLV = pLV;
  if (fppCurrentPV) *fppCurrentPV = fpCurrentPV;
  if (fppCurrentLV) *fppCurrentLV = fpCurrentLV;

  const HepRotation* pObjectRotation = pVPV -> GetObjectRotation ();
  const Hep3Vector&  translation     = pVPV -> GetTranslation ();
  G4Transform3D theLT (G4Transform3D (*pObjectRotation, translation));
  G4Transform3D theNewAT (theAT * theLT);

  /********************************************************
  G4cout << "G4PhysicalVolumeModel::DescribeAndDescend: "
	 << pVPV -> GetName () << "." << pVPV -> GetCopyNo ();
  G4cout << "\n  theAT: ";
  G4cout << "\n    Rotation: ";
  HepRotation rotation = theAT.getRotation ();
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

  // First check if mother covers...

  // This is only effective in surface drawing style, and then only if
  // the volumes are visible and opaque, and then only if no sections
  // or cutways are in operation.
  G4bool cullDaughter = thisToBeDrawn && IsDaughterCulled (pLV);
  if (!cullDaughter) {
    // OK, now let's check for daughters...
    if (soughtDepth != 0) {
      int nDaughters = pLV -> GetNoDaughters ();
      if (nDaughters) {
	for (int iDaughter = 0; iDaughter < nDaughters; iDaughter++) {
	  G4VPhysicalVolume* pVPV = pLV -> GetDaughter (iDaughter);
	  // Descend the geometry structure recursively...
	  fCurrentDepth++;
	  VisitGeometryAndGetVisReps
	    (pVPV, soughtDepth - 1, theNewAT, sceneHandler);
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
  sceneHandler.PreAddThis (theAT, *pVisAttribs);
  pSol -> DescribeYourselfTo (sceneHandler);
  sceneHandler.PostAddThis ();
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
	// &&
	// true // ...and opaque (transparency parameter not yet implemented).
	)
       )
      ;
  }
  else {
    return false;
  }
}

G4bool G4PhysicalVolumeModel::Validate () {
  G4VPhysicalVolume* world =
    G4TransportationManager::GetTransportationManager ()
    -> GetNavigatorForTracking () -> GetWorldVolume ();
  // The idea now is to seek a PV with the same name and copy no
  // in the hope it's the same one!!
  G4cout << "G4PhysicalVolumeModel::Validate() called." << G4endl;
  G4PhysicalVolumeSearchScene searchScene (fTopPVName, fTopPVCopyNo);
  G4PhysicalVolumeModel searchModel (world);
  searchModel.DescribeYourselfTo (searchScene);
  G4VPhysicalVolume* foundVolume = searchScene.GetFoundVolume ();
  if (foundVolume) {
    G4cout << "  Volume of the same name and copy number (\""
	   << fTopPVName << "\", copy " << fTopPVCopyNo
	   << ") still exists and is being used."
      "\n  Be warned that this does not necessarily guarantee it's the same"
      "\n  volume you originally specified in /vis/scene/add/."
	   << G4endl;
    fpTopPV = foundVolume;
    return true;
  }
  else {
    G4cout << "  A volume of the same name and copy number (\""
	   << fTopPVName << "\", copy " << fTopPVCopyNo
	   << ") no longer exists."
	   << G4endl;
    return false;
  }
}

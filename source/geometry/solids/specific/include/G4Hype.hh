// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Hype.hh,v 1.2 2000-04-11 16:03:39 johna Exp $
// $Original: G4Hype.hh,v 1.0 1998/06/09 16:57:50 safai Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4Hype
//
// Class description:
//
//   This class implements in G4 the volume equivalent to the
//   HYPE volume in Geant 3.21, i.e. a tube with hyperbolic profile.
//
//   For further informations, please read G4Hype.history and G4Hype.doc.
//
//   An hyperbolic volume  with curved sides parallel to the z-axis.
//   The Hype has a specified half-length along the z axis, about which
//   it is centred, and a given minimum and maximum radius.
//   A minimum radius of 0 signifies a filled Hype (with hyperbolical
//   inner surface). To have a filled Hype the user must specify 
//   inner radius = 0 AND inner stereo angle = 0.
// 
//   The inner and outer hyperbolical surfaces can have different
//   stereo angles. A stereo angle of 0 gives a cylindrical surface.
//
// Member functions:
//
//   As inherited from G4VSolid,
//
//   G4Hype(const G4String      &pName
//          const G4double      innerRadius
//          const G4double      outerRadius
//          const G4double      innerStereo
//          const G4double      outerStereo
//          const G4double      halfLenZ )
//
//   Construct an hype with the given name and dimensions.
//   The provided angles are in radians.
//
//
//   Protected:
//
//   G4ThreeVectorList*
//   CreateRotatedVertices(const G4AffineTransform& pTransform) const
//
//   Create the List of transformed vertices in the format required
//   for G4VSolid:: ClipCrossSection and ClipBetweenSections.

// Authors: 
//      Ernesto Lamanna (Ernesto.Lamanna@roma1.infn.it) &
//      Francesco Safai Tehrani (Francesco.SafaiTehrani@roma1.infn.it)
//      Rome, INFN & University of Rome "La Sapienza",  9 June 1998.
//
// History:
//      Updated Feb 2000 D.C. Williams
//
// --------------------------------------------------------------------

#ifndef G4HYPE_HH
#define G4HYPE_HH

#include "G4CSGSolid.hh"
#include "G4ThreeVector.hh"

class G4SolidExtentList;
class G4ClippablePolygon;

class G4Hype : public G4CSGSolid {
public:

  G4Hype(const G4String &pName,
	 const G4double newInnerRadius,
	 const G4double newOuterRadius,
	 const G4double newInnerStereo,
	 const G4double newOuterStereo,
	 const G4double newHalfLenZ);

  virtual ~G4Hype();
    
  void ComputeDimensions(G4VPVParameterisation* p,
			 const G4int n,
			 const G4VPhysicalVolume* pRep);

  G4bool CalculateExtent(const EAxis pAxis,
			 const G4VoxelLimits& pVoxelLimit,
			 const G4AffineTransform& pTransform,
			 G4double& pmin, G4double& pmax) const;

  G4double    GetInnerRadius () const { return innerRadius; }
  G4double    GetOuterRadius () const { return outerRadius; }
  G4double    GetZHalfLength () const { return halfLenZ; }
  G4double    GetInnerStereo () const { return innerStereo; }
  G4double    GetOuterStereo () const { return outerStereo; }

  void        SetInnerRadius (G4double newIRad) 
  { 
    innerRadius= newIRad;
    innerRadius2= newIRad*newIRad;
    endInnerRadius2=HypeInnerRadius2(halfLenZ);
    endInnerRadius=sqrt(endInnerRadius2);
  }

  void        SetOuterRadius (G4double newORad) 
  { 
    outerRadius= newORad;
    outerRadius2=newORad*newORad;
    endOuterRadius2=HypeOuterRadius2(halfLenZ);
    endOuterRadius=sqrt(endOuterRadius2);
  } 

  void        SetZHalfLength (G4double newHLZ)  { halfLenZ   = newHLZ ; }

  void        SetInnerStereo (G4double newISte) 
  { 
    innerStereo= fabs(newISte);
    tanInnerStereo=tan(innerStereo);
    tanInnerStereo2=tanInnerStereo*tanInnerStereo;
    endInnerRadius2=HypeInnerRadius2(halfLenZ);
    endInnerRadius=sqrt(endInnerRadius2);
  }

  void        SetOuterStereo (G4double newOSte)
  { 
    outerStereo= fabs(newOSte);
    tanOuterStereo=tan(outerStereo);
    tanOuterStereo2=tanOuterStereo*tanOuterStereo;
    endOuterRadius2=HypeOuterRadius2(halfLenZ);
    endOuterRadius=sqrt(endOuterRadius2);
  }

  EInside Inside(const G4ThreeVector& p) const;

  G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;

  G4double DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const;
  G4double DistanceToIn(const G4ThreeVector& p) const;
  G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			 const G4bool calcNorm=G4bool(false),
			 G4bool *validNorm=0,G4ThreeVector *n=0) const;
  G4double DistanceToOut(const G4ThreeVector& p) const;

  virtual G4GeometryType  GetEntityType() const { return G4String("G4Hype"); }

  void                DescribeYourselfTo (G4VGraphicsScene& scene) const;
  G4Polyhedron*       CreatePolyhedron   () const;
  G4NURBS*            CreateNURBS        () const;

protected:
  
  // whether we have an inner surface or not
  G4bool InnerSurfaceExists() const { return (innerRadius > DBL_MIN) || (innerStereo != 0); }

  // approximate isotropic distance to hyperbolic surface 
  static G4double ApproxDistOutside( const G4double pr, const G4double pz,
  			  	     const G4double r0, const G4double tanPhi );
  static G4double ApproxDistInside( const G4double pr, const G4double pz,
  			  	    const G4double r0, const G4double tan2Phi );

  // values of hype radius at a given Z
  G4double HypeInnerRadius2(const G4double zVal) const { return (tanInnerStereo2*zVal*zVal+innerRadius2); }
  G4double HypeOuterRadius2(const G4double zVal) const { return (tanOuterStereo2*zVal*zVal+outerRadius2); }

  // intersection with hyperbolic surface
  static G4int IntersectHype( const G4ThreeVector &p, const G4ThreeVector &v, 
                              const G4double r2, const G4double tan2Phi, G4double s[2] );

  static void AddPolyToExtent( const G4ThreeVector &v0,
  			       const G4ThreeVector &v1,
			       const G4ThreeVector &w1,
			       const G4ThreeVector &w0,
			       const G4VoxelLimits &voxelLimit,
			       const EAxis axis,
			       G4SolidExtentList &extentList );
	
 
  G4double innerRadius; // variable names are quite self explanative
  G4double outerRadius;
  G4double halfLenZ;
  G4double innerStereo;
  G4double outerStereo;

  // precalculated parameters, squared quantities
  G4double tanInnerStereo;
  G4double tanOuterStereo;
  G4double tanInnerStereo2; // squared tan of Inner Stereo angle
  G4double tanOuterStereo2; // squared tan of Outer Stereo angle
  G4double innerRadius2;    // squared Inner Radius
  G4double outerRadius2;    // squared Outer Radius
  G4double endInnerRadius2; // squared endcap Inner Radius
  G4double endOuterRadius2; // squared endcap Outer Radius
  G4double endInnerRadius; // endcap Inner Radius
  G4double endOuterRadius; // endcap Outer Radius
  
  // Used by distanceToOut
  enum ESide {outerFace,innerFace,leftCap, rightCap};
};
   	
#endif

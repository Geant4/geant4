// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Torus.hh,v 1.1 1999-01-07 16:07:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Torus
//
// A torus or torus segment with curved sides parallel to
// the z-axis. The torus has a specified swept radius
//  about which it is centred, and a given
// minimum and maximum radius. A minimum radius of 0
// signifies a filled torus . The torus segment is
// specified by starting and delta
// angles for phi, with 0 being the +x axis, PI/2
// the +y axis. A delta angle of 2PI signifies a
// complete, unsegmented torus/cylindr.
//
// Member functions:
//
// As inherited from G4CSGSolid+
//
// G4Torus(const G4String      &pName
//          G4double      pRmin
//          G4double      pRmax
//          G4double      pRtor
//          G4double      pSPhi
//          G4double      pDPhi )
//
//   Construct a torus with the given name and dimensions.
//   The angles are provided is radians. pRtor >= pRmax
//
//
// Protected:
//
// G4ThreeVectorList*
// CreateRotatedVertices(const G4AffineTransform& pTransform) const
//
//   Create the List of transformed vertices in the format required
//   for G4VSolid:: ClipCrossSection and ClipBetweenSections.
//   
// Member Data:
//
//	fRmin	Inside radius
//	fRmax	Outside radius
//	fRtor	swept radius of torus
//
//	fSPhi	The starting phi angle in radians,
//              adjusted such the fSPhi+fDPhi<=2PI,
//              fSPhi>-2PI
//
//	fDPhi	Delta angle of the segment in radians
//
//    You could find very often in G4Torus:: functions the values like pt or
//    it . These are the distances from p or i G4ThreeVector points in the
//    plane (Z axis points p or i) to fRtor point in XY plane. This value is
//    similar to rho for G4Tubs and is used for definiton of the point
//    relative to fRmin and fRmax, i.e. for solution of inside/outside
//    problems
//
// 
// History:
// 30.10.96 V.Grichine     First version of G4Torus
// 21.04.98 J.Apostolakis  Added SetAllParameters function

#ifndef G4Torus_HH
#define G4Torus_HH

#include "G4CSGSolid.hh"

class G4Torus : public G4CSGSolid {
public:
G4Torus(const G4String &pName,
	  G4double pRmin,
	  G4double pRmax,
	  G4double pRtor,
	  G4double pSPhi,
	  G4double pDPhi);

    virtual ~G4Torus();
    
    void SetAllParameters(G4double pRmin, G4double pRmax, G4double pRtor,
	       G4double pSPhi, G4double pDPhi);
 
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);
			   
    G4int TorusRoots(G4double Ri,
		     const G4ThreeVector& p,
		     const G4ThreeVector& v) const ;

    G4bool CalculateExtent(const EAxis pAxis,
			   const G4VoxelLimits& pVoxelLimit,
			   const G4AffineTransform& pTransform,
			   G4double& pmin, G4double& pmax) const;

    G4double    GetRmin() const { return fRmin ; }
    G4double    GetRmax() const { return fRmax ; } 
    G4double    GetRtor  () const { return fRtor   ; }
    G4double    GetSPhi() const { return fSPhi ; }
    G4double    GetDPhi() const { return fDPhi ; }

    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			   const G4bool calcNorm=G4bool(false),
			   G4bool *validNorm=0,G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    // Naming method (pseudo-RTTI : run-time type identification)
    virtual G4GeometryType  GetEntityType() const { return G4String("G4Torus"); }

    // Visualisation functions
    void                DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent         GetExtent          () const;
    G4Polyhedron*       CreatePolyhedron   () const;
    G4NURBS*            CreateNURBS        () const;

protected:

    G4int SolveBiQuadratic(double c[], double s[]  ) const ;
    G4int SolveCubic(double c[], double s[]  ) const ;
    G4int SolveQuadratic(double c[], double s[]  ) const ;
    
    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform,
                          G4int& noPolygonVertices) const;

    G4double fRmin,fRmax,fRtor,fSPhi,fDPhi;

  // Used by distanceToOut
  enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi};
  // used by normal
  enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi};

};
   	
#endif

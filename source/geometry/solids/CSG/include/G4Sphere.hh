// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Sphere.hh,v 1.3 2000-04-07 12:55:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Sphere
//
// Class description:
//
//   A G4Sphere is, in the general case, section of a spherical shell,
//   between specified phii and theta angles
//
//   The phi and theta segments are described by a starting angle,
//   and the +ve delta angle for the shape.
//   If the delta angle is >=2*M_PI, or >=M_PI the shape is treated as
//   continuous in phi or theta respectively.
//
//   Theta must lie between 0-PI (incl).
//
//   Member Data:
//
//	fRmin	inner radius
//	fRmax	outer radius
//
//	fSPhi	starting angle of the segment in radians
//	fDPhi	delta angle of the segment in radians
//
//	fSTheta	starting angle of the segment in radians
//	fDTheta	delta angle of the segment in radians
//
//     
//   Note:
//      Internally fSPhi & fDPhi are adjusted so that fDPhi<=2PI,
//      and fDPhi+fSPhi<=2PI. This enables simpler comparisons to be
//      made with (say) Phi of a point.

// History:
// 28.3.94 P.Kent Old C++ code converted to tolerant geometry
// 17.9.96 V.Grichine Final modifications to commit
// --------------------------------------------------------------------

#ifndef G4Sphere_HH
#define G4Sphere_HH

#include "G4CSGSolid.hh"

class G4Sphere : public G4CSGSolid {
public:
          G4Sphere(const G4String& pName,
	           G4double pRmin, G4double pRmax,
	           G4double pSPhi, G4double pDPhi,
	           G4double pSTheta, G4double pDTheta);
		   
	  virtual ~G4Sphere() ;
	  
	  // Access functions
		   
	  G4double    GetInsideRadius   () const { return fRmin;   }
	  G4double    GetOuterRadius  () const { return fRmax;   }
	  G4double    GetStartPhiAngle  () const { return fSPhi;   }
	  G4double    GetDeltaPhiAngle  () const { return fDPhi;   }
	  G4double    GetStartThetaAngle() const { return fSTheta; }
	  G4double    GetDeltaThetaAngle() const { return fDTheta; }

	  void   SetInsideRadius   (G4double newRmin)   { fRmin= newRmin; }
	  void   SetOuterRadius  (G4double newRmax)   { fRmax= newRmax; }
	  void   SetStartPhiAngle  (G4double newSphi)   { fSPhi= newSphi; }
	  void   SetDeltaPhiAngle  (G4double newDphi)   { fDPhi= newDphi; }
	  void   SetStartThetaAngle(G4double newSTheta) { fSTheta=newSTheta; }
	  void   SetDeltaThetaAngle(G4double newDTheta) { fDTheta=newDTheta; }

          void ComputeDimensions(G4VPVParameterisation* p,
                                 const G4int n,
                                 const G4VPhysicalVolume* pRep);

          G4bool CalculateExtent(const EAxis pAxis,
			         const G4VoxelLimits& pVoxelLimit,
			         const G4AffineTransform& pTransform,
			         G4double& pmin, G4double& pmax) const;
				 
          EInside Inside(const G4ThreeVector& p) const;

          G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

          G4double DistanceToIn(const G4ThreeVector& p,
	                        const G4ThreeVector& v) const;
	  
          G4double DistanceToIn(const G4ThreeVector& p) const;
	  
          G4double DistanceToOut(const G4ThreeVector& p,
	                         const G4ThreeVector& v,
			         const G4bool calcNorm=G4bool(false),
			         G4bool *validNorm=0,
				 G4ThreeVector *n=0) const;
				 
          G4double DistanceToOut(const G4ThreeVector& p) const;

             // Naming method (pseudo-RTTI : run-time type identification)

          virtual G4GeometryType  GetEntityType() const { return G4String("G4Sphere"); }

              // Visualisation functions
  
          void   DescribeYourselfTo(G4VGraphicsScene& scene) const;
       
          G4VisExtent GetExtent() const;
       
          G4Polyhedron* CreatePolyhedron() const;
       
          G4NURBS*      CreateNURBS() const;
       
          // Old access functions
          G4double  GetRmin()   const { return GetInsideRadius   (); }
          G4double  GetRmax()   const { return GetOuterRadius    (); }
          G4double  GetSPhi()   const { return GetStartPhiAngle  (); }
          G4double  GetDPhi()   const { return GetDeltaPhiAngle  (); }
	  G4double  GetSTheta() const { return GetStartThetaAngle(); }
          G4double  GetDTheta() const { return GetDeltaThetaAngle(); }

protected:
 
          G4ThreeVectorList*
          CreateRotatedVertices(const G4AffineTransform& pTransform,
	                        G4int& noPolygonVertices) const;
	
  // Used by distanceToOut
  
          enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kSTheta,kETheta};
  
  // used by normal
  
          enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNSTheta,kNETheta};

private:

    G4double fRmin,fRmax,
             fSPhi,fDPhi,
	     fSTheta,fDTheta;
	     
};
   	
#endif

// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Cons.hh,v 1.3 2000-04-07 12:55:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Cons
//
// Class description:
//
//   A G4Cons is, in the general case, a Phi segment of a cone, with
//   half-length fDz, inner and outer radii specified at -fDz and +fDz.
//   The Phi segment is described by a starting fSPhi angle, and the
//   +fDPhi delta angle for the shape.
//   If the delta angle is >=2*M_PI, the shape is treated as continuous
//   in Phi
//
//   Member Data:
//
//	fRmin1	inside radius at  -fDz
//	fRmin2	inside radius at  +fDz
//	fRmax1	outside radius at -fDz
//	fRmax2	outside radius at +fDz
//	fDz	half length in z
//
//	fSPhi	starting angle of the segment in radians
//	fDPhi	delta angle of the segment in radians
//
//   Note:
//      Internally fSPhi & fDPhi are adjusted so that fDPhi<=2PI,
//      and fDPhi+fSPhi<=2PI. This enables simpler comparisons to be
//      made with (say) Phi of a point.

// History:
// 19.3.94 P.Kent Old C++ code converted to tolerant geometry
// 13.9.96 V.Grichine Final modifications to commit
// --------------------------------------------------------------------

#ifndef G4Cons_HH
#define G4Cons_HH

#include "G4CSGSolid.hh"

class G4Cons : public G4CSGSolid
{
public:
        G4Cons(const G4String& pName,
	       G4double pRmin1, G4double pRmax1,
               G4double pRmin2, G4double pRmax2,
               G4double pDz,
	       G4double pSPhi, G4double pDPhi);
	       
	virtual ~G4Cons() ;

// Access functions	
        G4double    GetInnerRadiusMinusZ() const { return fRmin1 ; }
        G4double    GetOuterRadiusMinusZ() const { return fRmax1 ; } 
        G4double    GetInnerRadiusPlusZ()  const { return fRmin2 ; }
        G4double    GetOuterRadiusPlusZ()  const { return fRmax2 ; }
	
        G4double    GetZHalfLength()       const { return fDz    ; }
	
        G4double    GetStartPhiAngle () const { return fSPhi; }
        G4double    GetDeltaPhiAngle () const { return fDPhi; }
	
  // Modifier functions
        void    SetInnerRadiusMinusZ( G4double Rmin1 )  { fRmin1= Rmin1 ; }
        void    SetOuterRadiusMinusZ( G4double Rmax1 )  { fRmax1= Rmax1 ; } 
        void    SetInnerRadiusPlusZ ( G4double Rmin2 )  { fRmin2= Rmin2 ; }
        void    SetOuterRadiusPlusZ ( G4double Rmax2 )  { fRmax2= Rmax2 ; }
	       
        void    SetZHalfLength      ( G4double newDz )  { fDz=    newDz ; }
        void    SetStartPhiAngle    ( G4double newSPhi) { fSPhi=  newSPhi; }
        void    SetDeltaPhiAngle    ( G4double newDPhi) { fDPhi=  newDPhi; }

  // Other methods
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
			       
        G4double DistanceToOut(const G4ThreeVector& p) const ;
			       
        virtual G4GeometryType  GetEntityType() const { return G4String("G4Cons"); }
			   
        void   DescribeYourselfTo(G4VGraphicsScene& scene) const;
       
        G4VisExtent GetExtent() const;
       
        G4Polyhedron* CreatePolyhedron() const;
       
        G4NURBS*      CreateNURBS() const;
       
        //  Old access functions
        G4double    GetRmin1() const { return GetInnerRadiusMinusZ(); }
        G4double    GetRmax1() const { return GetOuterRadiusMinusZ(); }
        G4double    GetRmin2() const { return GetInnerRadiusPlusZ();  }
        G4double    GetRmax2() const { return GetOuterRadiusPlusZ();  }
	
        G4double    GetDz()    const { return GetZHalfLength() ; } // fDz  
	
        G4double    GetSPhi() const { return GetStartPhiAngle(); }  // fSPhi 
        G4double    GetDPhi() const { return GetDeltaPhiAngle(); }  // fDPhi

protected:
 
        G4ThreeVectorList*
        CreateRotatedVertices(const G4AffineTransform& pTransform) const;
	
  // Used by distanceToOut
  
        enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};
  
  // used by normal
  
        enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ};

private:

        G4double fRmin1,fRmin2,
	         fRmax1,fRmax2,
		 fDz,
		 fSPhi,fDPhi;
    
};
   	
#endif

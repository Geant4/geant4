// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Trap.hh,v 1.1 1999-01-07 16:07:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Trap
//
// A G4Trap is a general trapezoid: The faces perpendicular to the z planes
// are tapezia, and their centres are not necessarily on a line parallel to
// the z axis.
//
// Note that of the 11 parameters desribed below, only 9 are really
// independent - a check for planarity is made in the calculation of the
// equation for each plane. If the planes are not parallel, a call to
// G4Exception is made.
//
//      pDz     Half-length along the z-axis
//      pTheta  Polar angle of the line joining the centres of the faces
//              at -/+pDz
//      pPhi    Azimuthal angle of the line joing the centre of the face at
//              -pDz to the centre of the face at +pDz
//      pDy1     Half-length along y of the face at -pDz
//      pDx1    Half-length along x of the side at y=-pDy1 of the face at -pDz
//      pDx2    Half-length along x of the side at y=+pDy1 of the face at -pDz
//      pAlp1   Angle with respect to the y axis from the centre of the side
//              at y=-pDy1 to the centre at y=+pDy1 of the face at -pDz
//
//      pDy2     Half-length along y of the face at +pDz
//      pDx3    Half-length along x of the side at y=-pDy2 of the face at +pDz
//      pDx4    Half-length along x of the side at y=+pDy2 of the face at +pDz
//      pAlp2   Angle with respect to the y axis from the centre of the side
//              at y=-pDy2 to the centre at y=+pDy2 of the face at +pDz
//
//
// Member Data:
//
//      fDz     Half-length along the z axis
//      fTthetaCphi = tan(pTheta)*cos(pPhi)   These combinations are suitable for
//      fTthetaSphi = tan(pTheta)*sin(pPhi)   creation of the trapezoid corners
//
//      fDy1     Half-length along y of the face at -fDz
//      fDx1    Half-length along x of the side at y=-fDy1 of the face at -fDz
//      fDx2    Half-length along x of the side at y=+fDy1 of the face at -fDz
//      fTalpha1   Tan of Angle with respect to the y axis from the centre of
//                 the side at y=-fDy1 to the centre at y=+fDy1 of the face at -fDz
//
//      fDy2     Half-length along y of the face at +fDz
//      fDx3    Half-length along x of the side at y=-fDy2 of the face at +fDz
//      fDx4    Half-length along x of the side at y=+fDy2 of the face at +fDz
//      fTalpha2   Tan of Angle with respect to the y axis from the centre of
//                 the side at y=-fDy1 to the centre at y=+fDy1 of the face at -fDz
//
//      TrapSidePlane fPlanes[4] Plane equations of the faces not at +/-fDz
//                             *** order is important !!! : 
//  
//
// Member functions:
//
// As inherited from G4CSGSolid +  Constructors
//
//   G4Trap( const G4String& pName,
//            G4double pDz,
//	      G4double pTheta, G4double pPhi,
//	      G4double pDy1, G4double pDx1, G4double pDx2,
//	      G4double pAlp1,
//	      G4double pDy2, G4double pDx3, G4double pDx4,
//	      G4double pAlp2)
//   which prepare plane equations and corner coordinates from parameters
// &
//   G4Trap( const G4string& pName,
//           const G4ThreeVector pt[8]) 
//   which prepare plane equations and parameters from corner coordinates
// &
//     G4Trap( const G4String& pName,
//              G4double pZ,
//	        G4double pY,
//	        G4double pX, G4double pLTX);
//   for Right Angular Wedge from STEP
//
// & Constructor for G4Trd	     
//	     
//     G4Trap( const G4String& pName,
//              G4double pDx1,  G4double pDx2,
//	        G4double pDy1,  G4double pDy2,
//              G4double pDz);
//
// & Constructor for G4Para
//	     
//     G4Trap(const G4String& pName,
//	     G4double pDx, G4double pDy, G4double pDz,
//	     G4double pAlpha, G4double pTheta, G4double pPhi);
//
// + Access functions that return the respective parameter
//
// G4double GetZHalfLength()  const
// G4ThreeVector GetSymAxis() const Returns coordinates of unit vector along straight
//                                  line joining centers of -/+fDz planes   
// G4double GetYHalfLength1() const
// G4double GetXHalfLength1() const
// G4double GetXHalfLength2() const
// G4double GetTanAlpha1()    const
//
// G4double GetYHalfLength2() const
// G4double GetXHalfLength3() const
// G4double GetXHalfLength4() const
// G4double GetTanAlpha2()    const
//
// TrapSidePlane GetSidePlane(G4int n ) const  ;  n = 0,1,2,3
//
// Protected:
//    G4bool MakePlane( G4ThreeVector& p1,
//                      G4ThreeVector& p2,
//		        G4ThreeVector& p3, 
//		        G4ThreeVector& p4,
//		        TrapSidePlane& plane ) 
//
// 
//   G4ThreeVectorList*
//   CreateRotatedVertices(const G4Transform& pTransform) const
// 
//   Create the List of transformed vertices in the format required
//   for G4CSGSolid:: ClipCrossSection and ClipBetweenSections.
//
//
// History:
//
// 23.3.94 P.Kent: Old C++ code converted to tolerant geometry
// 9.9.96  V.Grichine: Final modifications before to commit 
// 1.11.96 V.Grichine Costructors for Right Angular Wedge from STEP & G4Trd/Para
// 8.12.97 J.Allison Added "nominal" contructor and method SetAllParameters.


#ifndef G4Trap_HH
#define G4Trap_HH

#include "G4CSGSolid.hh"

struct TrapSidePlane
{
    G4double a,b,c,d;		// Normal unit vector (a,b,c)  and offset (d)
				// => Ax+By+Cz+D=0  
};

class G4Trap : public G4CSGSolid {

public:

// The most general constructor for G4Trap     

    G4Trap( const G4String& pName,
             G4double pDz,
	     G4double pTheta, G4double pPhi,
	     G4double pDy1, G4double pDx1, G4double pDx2,
	     G4double pAlp1,
	     G4double pDy2, G4double pDx3, G4double pDx4,
	     G4double pAlp2);
	                              
     G4Trap( const G4String& pName,
	     const G4ThreeVector pt[8]) ;

// Constructor for Right Angular Wedge from STEP

     G4Trap( const G4String& pName,
              G4double pZ,
	      G4double pY,
	      G4double pX, G4double pLTX);

// Constructor for G4Trd	     
	     
     G4Trap( const G4String& pName,
               G4double pDx1,  G4double pDx2,
	       G4double pDy1,  G4double pDy2,
               G4double pDz);

// Constructor for G4Para
	     
     G4Trap(const G4String& pName,
	     G4double pDx, G4double pDy, G4double pDz,
	     G4double pAlpha, G4double pTheta, G4double pPhi);
	     
	     virtual ~G4Trap() ;
    
  // Constructor for "nominal" G4Trap whose parameters are to be set
  // by a G4VPramaterisation later.

     G4Trap(const G4String& pName);

                                      // Access functions

    G4double GetZHalfLength() const
    {
	return fDz ;
    }
    
    G4ThreeVector GetSymAxis() const
    {
     G4double cosTheta = 1.0/sqrt(1+fTthetaCphi*fTthetaCphi+fTthetaSphi*fTthetaSphi) ;
     
     return G4ThreeVector(fTthetaCphi*cosTheta,fTthetaSphi*cosTheta,cosTheta) ;
    }

    G4double GetYHalfLength1() const
    {
	return fDy1 ;
    }
    

    G4double GetXHalfLength1() const
    {
	return fDx1 ;
    }

    G4double GetXHalfLength2() const
    {
	return fDx2 ;
    }

    G4double GetTanAlpha1()    const
    {
        return fTalpha1 ; 
    }
    
    G4double GetYHalfLength2() const
    {
	return fDy2 ;
    }

    G4double GetXHalfLength3() const
    {
	return fDx3 ;
    }

    G4double GetXHalfLength4() const
    {
	return fDx4 ;
    }

    G4double GetTanAlpha2()    const
    {
        return fTalpha2 ; 
    }
    
    TrapSidePlane GetSidePlane(G4int n ) const
    {
        return fPlanes[n] ;
    }

                                                      // Set Methods
    void SetAllParameters ( G4double pDz,
			    G4double pTheta,
			    G4double pPhi,
			    G4double pDy1,
			    G4double pDx1,
			    G4double pDx2,
			    G4double pAlp1,
			    G4double pDy2,
			    G4double pDx3,
			    G4double pDx4,
			    G4double pAlp2);
	                              
                                                      // Methods
    
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const;    
        
    EInside Inside(const G4ThreeVector& p) const;
    
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const;
    
    G4double DistanceToIn(const G4ThreeVector& p) const;
    
    G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			   const G4bool calcNorm=false,
			   G4bool *validNorm=0,G4ThreeVector *n=0) const;
			   
    G4double DistanceToOut(const G4ThreeVector& p) const;

             // Naming method (pseudo-RTTI : run-time type identification

    virtual G4GeometryType  GetEntityType() const { return G4String("G4Trap"); }

                        // Visualisation functions

    void                DescribeYourselfTo (G4VGraphicsScene& scene) const;
    
    G4VisExtent         GetExtent          () const;
    
    G4Polyhedron* CreatePolyhedron   () const;
    
    G4NURBS*      CreateNURBS        () const;

protected:

    G4bool MakePlanes();
    G4bool MakePlane( const G4ThreeVector& p1,
                      const G4ThreeVector& p2,
		      const G4ThreeVector& p3, 
		      const G4ThreeVector& p4,
		      TrapSidePlane& plane ) ;

    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;

private:

    G4double fDz,fTthetaCphi,fTthetaSphi;
    G4double fDy1,fDx1,fDx2,fTalpha1;
    G4double fDy2,fDx3,fDx4,fTalpha2;
    TrapSidePlane fPlanes[4];

};

#endif


//  **************************** End of G4Trap.hh *****************************************

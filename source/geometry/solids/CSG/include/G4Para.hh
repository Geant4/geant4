// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Para.hh,v 1.1 1999-01-07 16:07:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Para
//
// A G4Parallepiped, essentially a box with half lengths dx,dy,dz `skewed'
// so that there are angles theta & phi of the polar line joining the faces at
// +-dz in z, and alpha formed by the y axis and the plane joinng the
// centre of the faces G4Parallel to the z-x plane at -dy and +dy.
//
// A G4Para is defined by:
//	dx,dy,dz	Half-length in x,y,z
//	alpha		Angle formed by the y axis and by the plane joining
//			the centre of the faces G4Parallel to the z-x plane
//			at -dy and +dy
//	theta		Polar angle of the line joining the centres of the
//			faces at -dz and +dz in z
//	phi		Azimuthal angle of the line joining the centres of the
//			faces at -dz and +dz in z
// Member data:
//
// Note that the angles parameters are not stored - precomputed trig is
// stored instead.
//
//      fDx   Half-length in x
//      fDy   Half-length in y
//      fDz   Half-length in z
//
//      fTalpha       Tan of alpha
//      fTthetaCphi   Tan theta * Cos phi
//      fTthetaSphi   Tan theta * Sin phi
//
// History:
// 21.3.94 P.Kent Old C++ code converted to tolerant geometry
// 31.10.96 V.Grichine Modifications according G4Box/Tubs before to commit

#ifndef G4Para_HH
#define G4Para_HH

#include "G4CSGSolid.hh"

class G4Para : public G4CSGSolid {
public:
          G4Para(const G4String& pName,
	         G4double pDx, G4double pDy, G4double pDz,
	         G4double pAlpha, G4double pTheta, G4double pPhi);
		 
          G4Para( const G4String& pName,
	          const G4ThreeVector pt[8]) ;
		 
	  virtual ~G4Para() ;
    
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

    G4double GetYHalfLength() const
    {
	return fDy ;
    }
    

    G4double GetXHalfLength() const
    {
	return fDx ;
    }


    G4double GetTanAlpha()    const
    {
        return fTalpha ; 
    }
    
                                             // Set  functions

    void SetXHalfLength(G4double val)
    {
	fDx= val;
    }
    void SetYHalfLength(G4double val) 
    {
	fDy= val;
    }
    void SetZHalfLength(G4double val) 
    {
	fDz= val;
    }
    void SetAlpha(double alpha)    
    {
        fTalpha= tan(alpha); 
    }
    void SetTanAlpha(double val)    
    {
        fTalpha= val; 
    }
    void SetThetaAndPhi(double pTheta, double pPhi)    
    {
	fTthetaCphi=tan(pTheta)*cos(pPhi);
	fTthetaSphi=tan(pTheta)*sin(pPhi);
    }
    void SetAllParameters(G4double pDx, G4double pDy, G4double pDz, 
            G4double pAlpha, G4double pTheta, G4double pPhi);
    
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

    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    
    G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			   const G4bool calcNorm=G4bool(false),
			   G4bool *validNorm=0,G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

             // Naming method (pseudo-RTTI : run-time type identification

    virtual G4GeometryType  GetEntityType() const { return G4String("G4Para"); }

                        // Visualisation functions

    void                DescribeYourselfTo (G4VGraphicsScene& scene) const;
    
    G4VisExtent         GetExtent          () const;
    
    G4Polyhedron* CreatePolyhedron   () const;
    
    G4NURBS*      CreateNURBS        () const;

protected:

    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;

private:
    G4double fDx,fDy,fDz;
    G4double fTalpha,fTthetaCphi,fTthetaSphi;
};
   	
#endif

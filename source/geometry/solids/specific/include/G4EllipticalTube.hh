// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EllipticalTube.hh,v 1.8 2000-11-02 16:54:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4EllipticalTube
//
// Class description:
//
//   Declaration of a CSG volume representing a tube with elliptical
//   cross section (geant3 solid 'ELTU'):
//   
//   G4EllipticalTube( const G4String& name, 
//                           G4double  Dx,
//                           G4double  Dy,
//                           G4double  Dz )
//
//   The equation of the surface in x/y is 1.0 = (x/dx)**2 + (y/dy)**2

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4EllipticalTube_hh
#define G4EllipticalTube_hh

#include "G4VSolid.hh"

class G4EllipticalTube : public G4VSolid
{
  public:

	G4EllipticalTube( const G4String &name, 
			  G4double theDx, G4double theDy, G4double theDz );
	virtual ~G4EllipticalTube();
			
	//
	// Standard CSG methods
	//
	virtual G4bool CalculateExtent(	const EAxis pAxis,
					const G4VoxelLimits& pVoxelLimit,
					const G4AffineTransform& pTransform,
					G4double& pmin, G4double& pmax) const;
	
	virtual EInside Inside( const G4ThreeVector& p) const;

	virtual G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

	virtual G4double DistanceToIn( const G4ThreeVector& p,const G4ThreeVector& v ) const;
	virtual G4double DistanceToIn( const G4ThreeVector& p ) const;
	virtual G4double DistanceToOut( const G4ThreeVector& p,const G4ThreeVector& v,
					const G4bool calcNorm=false,
					G4bool *validNorm=0,G4ThreeVector *n=0 ) const;
	virtual G4double DistanceToOut( const G4ThreeVector& p ) const;

	G4GeometryType GetEntityType() const { return G4String("G4EllipticalTube"); }

	//
	// Visualisation methods
	//
        virtual G4Polyhedron* CreatePolyhedron() const;
	virtual void DescribeYourselfTo( G4VGraphicsScene& scene ) const;
	virtual G4VisExtent GetExtent() const;


	//
	// Parameter access
	//
	G4double GetDx() const { return dx; }
	G4double GetDy() const { return dy; }
	G4double GetDz() const { return dz; }
	
	void SetDx( const G4double newDx ) { dx = newDx; }
	void SetDy( const G4double newDy ) { dy = newDy; }
	void SetDz( const G4double newDz ) { dz = newDz; }
	
  protected:

	G4double dx, dy, dz;
	
	//
	// Utility
	//
	inline G4double CheckXY( const G4double x, const G4double y, const G4double toler ) const {
		G4double rx = x/(dx+toler), ry = y/(dy+toler);
		return rx*rx + ry*ry;
	}
	inline G4double CheckXY( const G4double x, const G4double y ) const {
		G4double rx = x/dx, ry = y/dy;
		return rx*rx + ry*ry;
	}

	G4int IntersectXY( const G4ThreeVector &p,
			   const G4ThreeVector &v, G4double s[2] ) const;
};

#endif

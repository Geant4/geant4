// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polycone.hh,v 1.2 2000-06-26 15:35:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4Polycone
//
// Class description:
//
//   Class implementing a CSG type "PCON" Geant 3.21 volume,
//   inherited from  class G4VCSGSolid.

// --------------------------------------------------------------------

#ifndef G4Polycone_hh
#define G4Polycone_hh

#include "G4VCSGfaceted.hh"
#include "G4PolyconeSide.hh"


class G4EnclosingCylinder;
class G4ReduciblePolygon;
class G4VCSGface;

class G4Polycone : public G4VCSGfaceted 
{
	public:
	G4Polycone( G4String name, 
                    const G4double phiStart,		// initial phi starting angle
                    const G4double phiTotal,		// total phi angle
                    const G4int numZPlanes,		// number of z planes
                    const G4double zPlane[],		// position of z planes
                    const G4double rInner[],		// tangent distance to inner surface
                    const G4double rOuter[]  );		// tangent distance to outer surface

	G4Polycone( G4String name, 
		    const G4double phiStart,		// initial phi starting angle
                    const G4double phiTotal,		// total phi angle
		    const G4int    numRZ,		// number corners in r,z space
		    const G4double r[], 		// r coordinate of these corners
		    const G4double z[]       ); 	// z coordinate of these corners

	virtual ~G4Polycone();
	
	G4Polycone( const G4Polycone &source );
	const G4Polycone &operator=( const G4Polycone &source );

	//
	// A couple overrides to speed things up
	//
	EInside Inside( const G4ThreeVector &p ) const;
	G4double DistanceToIn( const G4ThreeVector &p, const G4ThreeVector &v ) const;
	G4double DistanceToIn( const G4ThreeVector &p ) const { return G4VCSGfaceted::DistanceToIn(p); }
	
	//
	// The usual G4VCSGface stuff
	//
	void ComputeDimensions(	G4VPVParameterisation* p,
				const G4int n,
				const G4VPhysicalVolume* pRep);

	virtual G4GeometryType  GetEntityType() const { return G4String("G4Polycone"); }

        G4Polyhedron* CreatePolyhedron() const;
	G4NURBS*      CreateNURBS() const;
	
	//
	// Access routines
	//
	inline G4double	GetStartPhi()		const { return startPhi; }
	inline G4double GetEndPhi()		const { return endPhi; }
	inline G4bool	IsOpen()		const { return phiIsOpen; }
	inline G4int	GetNumRZCorner()	const { return numCorner;}
	inline G4PolyconeSideRZ GetCorner( const G4int index ) const { return corners[index]; }	
	
	protected:
	//
	// Here are our parameters
	//
	G4double startPhi;		// Starting phi value (0 < phiStart < 2pi)
	G4double endPhi;		// end phi value (0 < endPhi-phiStart < 2pi)
	G4bool	 phiIsOpen;		// true if there is a phi segment
	G4int	 numCorner;		// number RZ points
	G4PolyconeSideRZ *corners;	// corner r,z points
	
	//
	// The following is temporary until graphics_reps is brought up to this design
	//
	class G4PolyconeHistorical;
	friend class G4PolyconeHistorical;
	
	struct G4PolyconeHistorical {
		G4PolyconeHistorical() {;}
		~G4PolyconeHistorical();
		G4PolyconeHistorical( const G4PolyconeHistorical &source );
	
		G4double Start_angle;
		G4double Opening_angle;
		G4int	 Num_z_planes;
		G4double *Z_values;
		G4double *Rmin;
		G4double *Rmax;
	};
	
	G4PolyconeHistorical	*original_parameters;

	//
	// Our quick test
	//
	G4EnclosingCylinder *enclosingCylinder;

	//
	// Generic initializer, called by all constructors
	//
	void Create( const G4double phiStart,	    // initial phi starting angle
            	     const G4double phiTotal,	    // total phi angle
		     G4ReduciblePolygon *rz    );   // r/z coordinate of these corners

	void CopyStuff( const G4Polycone &source );
};

#endif
